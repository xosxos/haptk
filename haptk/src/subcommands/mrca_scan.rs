use std::collections::BTreeMap;
use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use rayon::prelude::*;

use super::bhst::Node;
use super::mrca::mrca_independent;
use crate::args::{Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_recombination_file};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::{Coord, PhasedMatrix};
use crate::utils::{centromeres_hg38, parse_coords};

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn run(
    args: StandardArgs,
    rec_rates: PathBuf,
    step_size: usize,
    no_csv: bool,
    centromere_cut_off: f32,
) -> Result<()> {
    match args.selection {
        Selection::OnlyAlts | Selection::OnlyRefs | Selection::Unphased => {
            return Err(eyre!(
                "Running with selection {:?} is not supported.",
                args.selection
            ))
        }
        _ => (),
    };

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None, true)?;

    let rates = read_recombination_file(rec_rates)?;

    tracing::info!("Starting the MRCA scan..");

    let ages = find_ages(&vcf, &args, rates, step_size, centromere_cut_off)?;

    tracing::info!("Finished the MRCA scan.");

    if !no_csv {
        write_ages_to_csv(&args, &ages, args.output.clone(), &vcf)?;
    }

    Ok(())
}

fn find_ages(
    vcf: &PhasedMatrix,
    args: &StandardArgs,
    rates: BTreeMap<u64, f32>,
    step_size: usize,
    centromere_cut_off: f32,
) -> Result<Vec<(Coord, f64, bool)>> {
    match args.selection == Selection::OnlyLongest {
        true => Vec::from_iter(vcf.coords().clone())
            .par_iter()
            .enumerate()
            .filter(|(n, _)| *n % step_size == 0)
            .map(|(_, coord)| {
                let only_longest_lengths = vcf.only_longest_lengths_no_shard(coord)?;

                let check = check_for_centromeres(vcf, &only_longest_lengths, centromere_cut_off);

                let (i_tau_hat, _, _) = mrca_independent(only_longest_lengths, coord.pos, &rates);
                Ok((coord.clone(), i_tau_hat, check))
            })
            .collect(),
        false => Vec::from_iter(vcf.coords().clone())
            .par_iter()
            .enumerate()
            .filter(|(n, _)| *n % step_size == 0)
            .map(|(_, coord)| {
                let shared_lengths = vcf.get_lengths_from_uhst_no_mut(coord)?;
                let check = check_for_centromeres(vcf, &shared_lengths, centromere_cut_off);
                let (i_tau_hat, _, _) = mrca_independent(shared_lengths, coord.pos, &rates);
                Ok((coord.clone(), i_tau_hat, check))
            })
            .collect(),
    }
}

fn check_for_centromeres(
    vcf: &PhasedMatrix,
    lengths: &Vec<(Node, Node)>,
    centromere_cut_off: f32,
) -> bool {
    let (c_start, c_stop) = centromeres_hg38(vcf.get_contig());

    let mut overlapping_centromere = 0;

    for (left_node, right_node) in lengths {
        let start = left_node.start.pos;
        let stop = right_node.stop.pos;

        let c1 = start > c_start && start < c_stop;
        let c2 = stop > c_start && stop < c_stop;
        let c3 = start < c_start && stop > c_stop;

        if c1 || c2 || c3 {
            overlapping_centromere += 1;
        }
    }

    // If over x% of lengths overlap a centromere, flag as true
    overlapping_centromere as f32 / lengths.len() as f32 > centromere_cut_off
}

fn write_ages_to_csv(
    args: &StandardArgs,
    ages: &Vec<(Coord, f64, bool)>,
    mut output: PathBuf,
    vcf: &PhasedMatrix,
) -> Result<()> {
    push_to_output(args, &mut output, "mrca_scan", "csv");

    let mut writer = open_csv_writer(output)?;
    writer.write_record(vec!["contig", "pos", "ref", "alt", "mrca", "centromere"])?;
    for (coord, age, check) in ages {
        writer.write_record(vec![
            vcf.get_contig().to_string(),
            format!("{}", coord.pos),
            format!("{}", coord.reference),
            format!("{}", coord.alt),
            format!("{age:.5}"),
            format!("{check}"),
        ])?;
    }
    Ok(())
}
