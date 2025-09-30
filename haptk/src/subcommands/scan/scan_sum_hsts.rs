use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread;
use std::{collections::BTreeMap, path::PathBuf};

use color_eyre::eyre::ensure;
use color_eyre::Result;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use rayon::prelude::*;

use super::Hst;
use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::io::{
    open_csv_writer, push_to_output, read_multiple_sample_ids, read_recombination_file,
};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::Coord;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::subcommands::scan::read_tree_file;
use crate::utils::{centromeres_hg38, parse_coords};

type Row = [String; 6];
const HEADER: [&str; 6] = ["contig", "pos", "ref", "alt", "avg", "centromere"];

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn run(
    args: ConciseArgs,
    rec_rates: PathBuf,
    centromere_cut_off: f32,
    construct_hsts_ad_hoc: bool,
    coords: Option<String>,
    step_size: usize,
    min_sample_size: usize,
    length_in_bp: bool,
    seg_samples: Option<Vec<PathBuf>>,
) -> Result<()> {
    let rec_rates = match length_in_bp {
        false => read_recombination_file(rec_rates)?,
        true => BTreeMap::new(),
    };

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
            let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None, true)?;

            let seg_samples = read_multiple_sample_ids(&seg_samples)?;
            let seg_samples = seg_samples.unwrap_or(vcf.samples().clone());

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
                        construct_bhst_no_mut(&vcf, coord, min_sample_size, start_idxs).unwrap(),
                    )
                })
                .map(|(coord, hst)| {
                    sum_and_create_csv_row(
                        vcf.samples(),
                        *vcf.ploidy,
                        hst,
                        coord.clone(),
                        &rec_rates,
                        &Some(&seg_samples),
                        centromere_cut_off,
                        length_in_bp,
                    )
                })
                .try_for_each(|row| tx.send(row))?;
        }
        false => {
            let samples = read_multiple_sample_ids(&args.samples)?;
            let hsts = read_tree_file(args.file)?;

            let vcf_samples = &hsts.metadata.samples;
            let ploidy = *hsts.metadata.ploidy;

            tracing::info!("Starting the leaf node sum scan..");
            hsts.hsts
                .into_par_iter()
                .map(|(coord, hst)| {
                    sum_and_create_csv_row(
                        vcf_samples,
                        ploidy,
                        hst,
                        coord,
                        &rec_rates,
                        &samples.as_ref(),
                        centromere_cut_off,
                        length_in_bp,
                    )
                })
                .try_for_each(|row| tx.send(row))?;
        }
    };

    tracing::info!("Finished the HST scan.");

    drop(tx);
    let _ = writer_handle.join();

    Ok(())
}

fn overlaps_centromeres(lengths: Vec<(u64, u64)>, contig: &str, centromere_cut_off: f32) -> bool {
    let (c_start, c_stop) = centromeres_hg38(contig);

    let lengths_len = lengths.len();

    let mut overlapping_centromere = 0;

    for (start, stop) in lengths.into_iter() {
        let c1 = start > c_start && start < c_stop;
        let c2 = stop > c_start && stop < c_stop;
        let c3 = start < c_start && stop > c_stop;

        if c1 || c2 || c3 {
            overlapping_centromere += 1;
        }
    }

    // If over x% of lengths overlap a centromere, flag as true
    overlapping_centromere as f32 / lengths_len as f32 > centromere_cut_off
}

#[allow(clippy::too_many_arguments)]
fn sum_and_create_csv_row(
    vcf_samples: &[String],
    ploidy: usize,
    hst: Hst,
    coord: Coord,
    rec_rates: &BTreeMap<u64, f32>,
    samples: &Option<&Vec<String>>,
    centromere_cut_off: f32,
    length_in_bp: bool,
) -> Row {
    let mut sum = 0.;
    let mut positions: Vec<(u64, u64)> = vec![];

    if length_in_bp {
        for node_idx in hst.node_indices() {
            let node = hst.node_weight(node_idx).unwrap();
            if hst
                .neighbors_directed(node_idx, Direction::Outgoing)
                .count()
                == 0
            {
                for idx in &node.indexes {
                    if let Some(samples) = samples {
                        let sample = &vcf_samples[idx / ploidy];

                        if samples.contains(sample) {
                            sum += (node.stop.pos as f32 - node.start.pos as f32).max(0.);
                            positions.push((node.start.pos, node.stop.pos));
                        }
                    } else {
                        sum += (node.stop.pos as f32 - node.start.pos as f32).max(0.);
                        positions.push((node.start.pos, node.stop.pos));
                    }
                }
            }
        }
    } else {
        for node_idx in hst.node_indices() {
            let node = hst.node_weight(node_idx).unwrap();
            if hst
                .neighbors_directed(node_idx, Direction::Outgoing)
                .count()
                == 0
            {
                let (_, start_cm) = rec_rates
                    .range(..=node.start.pos)
                    .next_back()
                    .unwrap_or_else(|| rec_rates.range(node.start.pos..).nth(1).unwrap());

                let (_, stop_cm) = rec_rates
                    .range(node.stop.pos..)
                    .next()
                    .unwrap_or_else(|| rec_rates.range(..node.stop.pos).next_back().unwrap());

                for idx in &node.indexes {
                    if let Some(samples) = samples {
                        let sample = &vcf_samples[idx / ploidy];

                        if samples.contains(sample) {
                            sum += (stop_cm - start_cm).max(0.);
                            positions.push((node.start.pos, node.stop.pos));
                        }
                    } else {
                        sum += (stop_cm - start_cm).max(0.);
                        positions.push((node.start.pos, node.stop.pos));
                    }
                }
            }
        }
    }

    let root = hst.node_weight(NodeIndex::new(0)).unwrap();

    let avg = if let Some(samples) = samples {
        let mut nhaplos = 0;
        for idx in &root.indexes {
            let sample = &vcf_samples[idx / ploidy];
            if samples.contains(sample) {
                nhaplos += 1;
            }
        }
        sum / nhaplos as f32
    } else {
        sum / root.indexes.len() as f32
    };

    let centromere = overlaps_centromeres(positions, &coord.contig, centromere_cut_off);

    [
        coord.contig,
        coord.pos.to_string(),
        coord.reference,
        coord.alt,
        avg.to_string(),
        centromere.to_string(),
    ]
}
