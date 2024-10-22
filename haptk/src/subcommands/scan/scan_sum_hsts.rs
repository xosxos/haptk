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

const HEADER: &[&str] = &["contig", "pos", "ref", "alt", "avg", "centromere"];

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn run(
    args: ConciseArgs,
    rec_rates: PathBuf,
    samples: Option<Vec<PathBuf>>,
    centromere_cut_off: f32,
    construct_hsts_ad_hoc: bool,
    coords: Option<String>,
    selection: Selection,
    step_size: usize,
    min_sample_size: usize,
    no_alt: bool,
) -> Result<()> {
    let rec_rates = read_recombination_file(rec_rates)?;

    // There was suprisingly very little in common with these if elses, this still needs to be refactored
    if construct_hsts_ad_hoc {
        ensure!( coords.is_some(), "If run ad hoc, please give the contig (with --coords) and the wanted selection (--alleles all/longest-haplotype)");

        let args = StandardArgs {
            file: args.file.clone(),
            output: args.output,
            info_limit: None,
            coords: coords.unwrap(),
            selection,
            prefix: args.prefix,
            samples: samples.clone(),
            no_alt,
        };

        let mut output = args.output.clone();
        push_to_output(&args, &mut output, "leaf_avg", "csv");
        let mut writer = open_csv_writer(output)?;

        let (contig, start, stop) = parse_coords(&args.coords)?;
        let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None, true)?;

        // Collect Vec::from_iter because BTreeSet requires `par_bridge` from Rayon which does not conserve order
        let results: Vec<(&Coord, f32, bool)> = Vec::from_iter(vcf.coords())
            .par_iter()
            .enumerate()
            .filter(|(n, _)| *n % step_size == 0)
            .map(|(_, &coord)| {
                let start_idxs = match args.selection == Selection::OnlyLongest {
                    true => {
                        let res = vcf.only_longest_indexes_no_shard(coord);
                        if res.is_err() {
                            panic!("failed to find the longest-haplotype for coord {coord:?}");
                        }
                        res.ok()
                    }
                    false => None,
                };

                (
                    coord,
                    construct_bhst_no_mut(&vcf, coord, min_sample_size, start_idxs).unwrap(),
                )
            })
            .map(|(coord, hst)| {
                sum(
                    vcf.samples(),
                    *vcf.ploidy,
                    &hst,
                    coord,
                    &rec_rates,
                    &Some(vcf.samples()),
                    centromere_cut_off,
                )
            })
            .collect();

        // Duplicate even the writer code because VCF and its Coord's, which we have a reference to, go out of scope at the end of this if
        writer.write_record(HEADER)?;

        for (coord, avg, centromere) in results {
            writer.write_record(vec![
                coord.contig.to_string(),
                coord.pos.to_string(),
                coord.reference.to_string(),
                coord.alt.to_string(),
                avg.to_string(),
                centromere.to_string(),
            ])?;
        }
    } else {
        let samples = read_multiple_sample_ids(&samples)?;
        let hsts = read_tree_file(args.file)?;

        let args = StandardArgs {
            file: PathBuf::new(),
            output: args.output,
            info_limit: None,
            coords: hsts.metadata.contig.clone(),
            selection: hsts.metadata.selection.clone(),
            prefix: args.prefix,
            samples: None,
            no_alt: false,
        };

        let mut output = args.output.clone();
        push_to_output(&args, &mut output, "leaf_avg", "csv");
        let mut writer = open_csv_writer(output)?;

        let vcf_samples = &hsts.metadata.samples;
        let ploidy = *hsts.metadata.ploidy;

        let results: Vec<(&Coord, f32, bool)> = hsts
            .hsts
            .par_iter()
            .map(|(coord, hst)| {
                sum(
                    vcf_samples,
                    ploidy,
                    hst,
                    coord,
                    &rec_rates,
                    &samples.as_ref(),
                    centromere_cut_off,
                )
            })
            .collect();

        writer.write_record(HEADER)?;

        for (coord, avg, centromere) in results {
            writer.write_record(vec![
                coord.contig.to_string(),
                coord.pos.to_string(),
                coord.reference.to_string(),
                coord.alt.to_string(),
                avg.to_string(),
                centromere.to_string(),
            ])?;
        }
    };

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

fn sum<'a>(
    vcf_samples: &[String],
    ploidy: usize,
    hst: &Hst,
    coord: &'a Coord,
    rec_rates: &BTreeMap<u64, f32>,
    samples: &Option<&Vec<String>>,
    centromere_cut_off: f32,
) -> (&'a Coord, f32, bool) {
    let mut sum = 0.;
    let mut positions: Vec<(u64, u64)> = vec![];

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

    (
        coord,
        avg,
        overlaps_centromeres(positions, &coord.contig, centromere_cut_off),
    )
}
