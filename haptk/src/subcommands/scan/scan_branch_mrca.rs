use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread;
use std::{collections::BTreeMap, path::PathBuf};

use color_eyre::eyre::ensure;
use color_eyre::Result;
use rayon::prelude::*;

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_recombination_file};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::Coord;
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::subcommands::scan::{read_tree_file, Limits};
use crate::utils::parse_coords;

use super::Hst;

type Row = [String; 14];

const HEADER: [&str; 14] = [
    "contig",
    "pos",
    "ref",
    "alt",
    "marker_id",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "mrca",
    "n_hom",
    "n_het",
    "centromere",
    "samples",
];

#[doc(hidden)]

pub fn run(
    args: ConciseArgs,
    limits: Limits,
    rec_rates: PathBuf,
    ad_hoc: bool,
    coords: Option<String>,
    step_size: usize,
    limit: f64,
) -> Result<()> {
    let rec_rates = read_recombination_file(rec_rates)?;

    let args = StandardArgs {
        file: args.file.clone(),
        output: args.output,
        coords: coords.clone().unwrap_or(String::from("")),
        selection: args.selection.clone(),
        prefix: args.prefix,
        samples: args.samples.clone(),
        no_alt: true,
    };

    let (tx, rx): (SyncSender<Row>, _) = sync_channel(2024);

    let cargs = args.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut output = cargs.output.clone();
        push_to_output(&cargs, &mut output, "branch_mrca_scan", "csv");
        let mut writer = open_csv_writer(output)?;
        writer.write_record(HEADER)?;

        while let Ok(row) = rx.recv() {
            writer.write_record(row)?;
        }
        Ok(())
    });

    match ad_hoc {
        true => {
            ensure!( coords.is_some(), "If run ad hoc, please give the contig (with --coords) and the wanted selection (--alleles all/longest-haplotype)");

            let (contig, start, stop) = parse_coords(&args.coords)?;
            let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None, true)?;

            tracing::info!("Starting the branch MRCA scan..");

            vcf.coords()
                .iter()
                .enumerate()
                .par_bridge()
                .filter(|(n, _)| *n % step_size == 0)
                .map(|(_, coord)| {
                    let start_idxs = match args.selection {
                        Selection::OnlyLongest => {
                            let res = vcf.only_longest_indexes_no_shard(coord);
                            if res.is_err() {
                                panic!("failed to find the longest-haplotype for coord {coord:?}");
                            }
                            res.ok()
                        }
                        _ => None,
                    };

                    (
                        coord,
                        construct_bhst_no_mut(&vcf, coord, limits.0, start_idxs).unwrap(),
                    )
                })
                .try_for_each(|(coord, hst)| {
                    find_optimized_value_and_create_csv_row(
                        hst,
                        coord.clone(),
                        limits,
                        &rec_rates,
                        vcf.samples(),
                        *vcf.ploidy,
                        limit,
                        tx.clone(),
                    )
                })?;
        }
        false => {
            let hsts = read_tree_file(args.file)?;

            let samples = hsts.metadata.samples.clone();
            let ploidy = *hsts.metadata.ploidy;

            tracing::info!("Starting the branch MRCA scan..");
            hsts.hsts.into_par_iter().try_for_each(|(coord, hst)| {
                find_optimized_value_and_create_csv_row(
                    hst,
                    coord,
                    limits,
                    &rec_rates,
                    &samples,
                    ploidy,
                    limit,
                    tx.clone(),
                )
            })?;
        }
    };

    tracing::info!("Finished the HST scan.");

    drop(tx);
    let _ = writer_handle.join();

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn find_optimized_value_and_create_csv_row(
    hst: Hst,
    coord: Coord,
    limits: Limits,
    rec_rates: &BTreeMap<u64, f32>,
    samples: &[String],
    ploidy: usize,
    limit: f64,
    tx: SyncSender<Row>,
) -> Result<()> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;

    let mut iterator = hst.node_indices();

    // Skip first node because it contains no haplotypes
    iterator.next();

    for idx in iterator {
        let node = hst.node_weight(idx).unwrap();

        if node.indexes.len() >= nmin_samples
            && node.indexes.len() <= nmax_samples
            && node.haplotype.len() >= nmin_variants
            && node.haplotype.len() <= nmax_variants
        {
            // OPTIMIZER CODE

            // Lowest MRCA
            let start_node = idx;
            let nodes = find_majority_nodes(&hst, start_node);

            // Iterate nodes from the leaf towards the root
            // NOTE: Sometimes genotyping data runs out and leaf nodes have more than 1 sample
            // NOTE: This is not really accounted for in the Gamma method without chromosome length corrections
            let mut nodes_iter = nodes.iter().map(|(_, idx)| idx).rev();
            let mut prev_idx = nodes_iter.next().unwrap();

            let mut used_indexes = vec![];
            let mut cm_sum = 0.;
            let mut n = 0;

            for idx in nodes_iter {
                let node = hst.node_weight(*idx).unwrap();
                let prev_node = hst.node_weight(*prev_idx).unwrap();

                for idx in &node.indexes {
                    if !used_indexes.contains(&idx) {
                        let (start, stop) = (prev_node.start.pos, prev_node.stop.pos);

                        let (_, start_cm) = rec_rates
                            .range(..=start)
                            .next_back()
                            .unwrap_or_else(|| rec_rates.range(start..).nth(1).unwrap());

                        let (_, stop_cm) = rec_rates
                            .range(stop..)
                            .next()
                            .unwrap_or_else(|| rec_rates.range(..stop).next_back().unwrap());

                        let distance = ((stop_cm - start_cm) / 100.).max(0.);

                        cm_sum += distance;
                        n += 1;

                        used_indexes.push(idx);
                    }
                }
                prev_idx = idx;
            }

            // Calculate MRCA estimate
            let n = n as f32;
            let b_c = (2.0 * n - 1.0) / (2.0 * n);
            let i_tau_hat = (b_c * 2.0 * n) / cm_sum;
            let value = i_tau_hat as f64;

            // Define the clause
            let clause = value < limit;

            // OPTIMIZER CODE END
            if clause {
                // Create a CSV row
                let (nhet, nhom) = node.zygosity(samples, ploidy);
                let is_centromere = node.check_for_centromere_hg38();
                let names = node.sample_name_list(samples, ploidy);

                tx.send([
                    coord.contig.clone(),
                    coord.pos.to_string(),
                    coord.reference.clone(),
                    coord.alt.clone(),
                    node.identifier(),
                    (node.stop.pos.saturating_sub(node.start.pos)).to_string(),
                    node.haplotype.len().to_string(),
                    node.start.pos.to_string(),
                    node.stop.pos.to_string(),
                    format!("{value:.4}"),
                    nhom.to_string(),
                    nhet.to_string(),
                    is_centromere.to_string(),
                    names,
                ])?;
            }
        }
    }
    Ok(())
}
