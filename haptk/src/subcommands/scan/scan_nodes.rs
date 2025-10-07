use std::sync::mpsc::sync_channel;
use std::sync::mpsc::SyncSender;
use std::thread;

use color_eyre::eyre::ensure;
use color_eyre::Result;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;

use crate::args::ConciseArgs;
use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::Coord;
use crate::io::open_csv_writer;
use crate::io::push_to_output;
use crate::read_vcf::read_vcf_to_matrix;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::subcommands::scan::read_tree_file;
use crate::subcommands::scan::Limits;
use crate::traits::OnlyLongest;
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
    "opt_value",
    "nhom",
    "nhet",
    "centromere",
    "samples",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    construct_hsts_ad_hoc: bool,
    coords: Option<String>,
    step_size: usize,
) -> Result<()> {
    let args = StandardArgs {
        file: args.file.clone(),
        output: args.output,
        coords: coords.clone().unwrap(),
        selection: args.selection.clone(),
        prefix: args.prefix,
        samples: args.samples.clone(),
        no_alt: true,
        list: None,
        include_indels: false,
    };

    let (tx, rx): (SyncSender<Row>, _) = sync_channel(2024);

    let cargs = args.clone();

    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut output = cargs.output.clone();
        push_to_output(&cargs, &mut output, "leaf_avg", "csv");
        let mut writer = open_csv_writer(output)?;
        writer.write_record(HEADER)?;

        while let Ok(row) = rx.recv() {
            writer.write_record(row)?;
        }

        Ok(())
    });

    match construct_hsts_ad_hoc {
        true => {
            ensure!( coords.is_some(), "If run ad hoc, please give the contig (with --coords) and the wanted selection (--alleles all/longest-haplotype)");

            let (contig, start, stop) = parse_coords(&args.coords)?;
            let vcf = read_vcf_to_matrix(&args, &contig, 0, Some((start, stop)), None, None, true)?;

            let samples = vcf.samples().clone();
            let ploidy = *vcf.ploidy;

            tracing::info!("Starting the leaf node sum scan..");
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
                .filter_map(|(coord, hst)| {
                    find_max_value_and_create_csv_row(hst, coord.clone(), limits, &samples, ploidy)
                })
                .try_for_each(|row| tx.send(row))?;
        }
        false => {
            let hsts = read_tree_file(args.file)?;

            let samples = hsts.metadata.samples.clone();
            let ploidy = *hsts.metadata.ploidy;

            tracing::info!("Starting the max length scan..");
            hsts.hsts
                .into_par_iter()
                .filter_map(|(coord, hst)| {
                    find_max_value_and_create_csv_row(hst, coord, limits, &samples, ploidy)
                })
                .try_for_each(|row| tx.send(row))?;
        }
    };

    tracing::info!("Finished the HST scan.");

    drop(tx);
    let _ = writer_handle.join();

    Ok(())
}

fn find_max_value_and_create_csv_row(
    hst: Hst,
    coord: Coord,
    limits: Limits,
    samples: &[String],
    ploidy: usize,
) -> Option<Row> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
    let mut top_node_idx = NodeIndex::new(0);

    // START VALUE
    let start_value = f64::MIN;

    let mut optimized_value = start_value;

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
            let node = hst.node_weight(idx).unwrap();

            // Longest haplotype in bp
            let start = node.start.pos;
            let stop = node.stop.pos;
            let value = stop
                .checked_sub(start)
                .expect("Overflowed subtraction, report to HAPTK github")
                as f64;

            // Longest haplotype in markers
            // let value = node.haplotype.len() as f64;

            // Clause
            let clause = value > optimized_value;

            // OPTIMIZER CODE END

            if clause {
                optimized_value = value;
                top_node_idx = idx;
            }
        }
    }

    if optimized_value == start_value {
        tracing::warn!(
            "No optimized value for the HST in contig {} at pos {}",
            coord.contig,
            coord.pos,
        );
        return None;
    }

    // Create CSV row

    let top_node = hst.node_weight(top_node_idx).unwrap();
    let is_centromere = top_node.check_for_centromere_hg38();

    let (nhet, nhom) = top_node.zygosity(samples, ploidy);
    let names = top_node.sample_name_list(samples, ploidy);

    Some([
        coord.contig,
        coord.pos.to_string(),
        coord.reference,
        coord.alt,
        top_node.identifier(),
        (top_node.stop.pos.saturating_sub(top_node.start.pos)).to_string(),
        top_node.haplotype.len().to_string(),
        top_node.start.pos.to_string(),
        top_node.stop.pos.to_string(),
        optimized_value.to_string(),
        nhom.to_string(),
        nhet.to_string(),
        is_centromere.to_string(),
        names,
    ])
}
