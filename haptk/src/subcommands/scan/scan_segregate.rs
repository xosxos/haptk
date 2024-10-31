use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread;

use color_eyre::Result;
use fishers_exact::fishers_exact;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_coords_file, read_sample_ids};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::Coord;
use crate::subcommands::bhst::Node;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::subcommands::scan::{read_tree_file, Hst, Limits};

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
    "n",
    "nhet_cases",
    "nhet_ctrls",
    "nhom_cases",
    "nhom_ctrls",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    seg_samples: PathBuf,
    coord_list_path: Option<PathBuf>,
) -> Result<()> {
    let seg_samples = read_sample_ids(&Some(seg_samples))?.unwrap();

    let args = StandardArgs {
        file: args.file.clone(),
        output: args.output,
        info_limit: None,
        coords: String::from(""),
        selection: args.selection.clone(),
        prefix: args.prefix,
        samples: args.samples.clone(),
        no_alt: true,
    };

    let (tx, rx): (SyncSender<Row>, _) = sync_channel(2024);

    let cargs = args.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut output = cargs.output.clone();
        push_to_output(&cargs, &mut output, "segregate_scan", "csv");
        let mut writer = open_csv_writer(output)?;
        writer.write_record(HEADER)?;

        while let Ok(row) = rx.recv() {
            writer.write_record(row)?;
        }
        Ok(())
    });

    let coord_list: Option<Vec<Coord>> = coord_list_path.map(|v| {
        read_coords_file(&v).unwrap_or_else(|_| panic!("Could not read coords file {v:?}"))
    });

    match coord_list {
        Some(coord_list) => {
            tracing::info!("Reading {} coords to HSTs.", coord_list.len());
            let mut map = HashMap::new();

            for coord in coord_list {
                map.entry(coord.contig.clone())
                    .or_insert(Vec::new())
                    .push(coord);
            }

            let mut run_once = true;
            let mut seg_sample_idxs: Vec<usize> = vec![];

            for (contig, coords) in map {
                let args = StandardArgs {
                    file: args.file.clone(),
                    output: args.output.clone(),
                    info_limit: None,
                    coords: format!("{}", coords[0]),
                    selection: args.selection.clone(),
                    prefix: args.prefix.clone(),
                    samples: args.samples.clone(),
                    no_alt: true,
                };

                let vcf = read_vcf_to_matrix(&args, &contig, 0, None, None, None, true)?;

                if run_once {
                    // Filter out samples not present in the HST SCAN
                    let seg_samples: Vec<String> = seg_samples
                        .iter()
                        .filter(|&s| {
                            if vcf.samples().contains(s) {
                                true
                            } else {
                                tracing::warn!(
                                    "Sample {s:?} is on the ID list, but not in the HSTs."
                                );
                                false
                            }
                        })
                        .cloned()
                        .collect();

                    seg_sample_idxs = vcf.get_idxs_for_samples(&seg_samples).unwrap();
                    run_once = false;
                }
                let ncases = seg_sample_idxs.len();
                let ncontrols = vcf.nhaplotypes() - ncases;

                tracing::info!("Starting the case-ctrl scan for {contig}.");

                coords
                    .into_par_iter()
                    .map(|coord| {
                        let start_idxs = match args.selection {
                            Selection::OnlyLongest => {
                                let res = vcf.only_longest_indexes_no_shard(&coord);
                                if res.is_err() {
                                    panic!(
                                        "failed to find the longest-haplotype for coord {coord:?}"
                                    );
                                }
                                res.ok()
                            }
                            _ => None,
                        };
                        (
                            construct_bhst_no_mut(&vcf, &coord, limits.0, start_idxs).unwrap(),
                            coord,
                        )
                    })
                    .filter_map(|(hst, coord)| {
                        find_optimized_value_and_create_csv_row(
                            hst,
                            coord,
                            limits,
                            &seg_sample_idxs,
                            ncases,
                            ncontrols,
                        )
                    })
                    .try_for_each(|row| tx.send(row))?;
            }
        }
        None => {
            let hsts = read_tree_file(args.file)?;

            // Filter out samples not present in the HST SCAN
            let seg_samples: Vec<String> = seg_samples
                .into_iter()
                .filter(|s| {
                    if hsts.metadata.samples.contains(s) {
                        true
                    } else {
                        tracing::warn!("Sample {s:?} is on the ID list, but not in the HSTs.");
                        false
                    }
                })
                .collect();

            let seg_sample_idxs = hsts.get_sample_idxs(&seg_samples).unwrap();

            let ncases = seg_sample_idxs.len();
            let ncontrols = hsts.nhaplotypes() - ncases;

            tracing::info!("Starting the case-ctrl scan.");

            hsts.hsts
                .into_par_iter()
                .filter_map(|(coord, hst)| {
                    find_optimized_value_and_create_csv_row(
                        hst,
                        coord,
                        limits,
                        &seg_sample_idxs,
                        ncases,
                        ncontrols,
                    )
                })
                .try_for_each(|row| tx.send(row))?;
        }
    }

    tracing::info!("Finished the case-ctrl scan.");

    let _ = writer_handle.join();

    Ok(())
}

fn find_optimized_value_and_create_csv_row(
    hst: Hst,
    coord: Coord,
    limits: Limits,
    case_list: &[usize],
    ncases: usize,
    ncontrols: usize,
) -> Option<Row> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
    let mut top_node_idx = NodeIndex::new(0);

    // START VALUE
    let start_value = f64::MAX;

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

            // if let Some(value) = segregation_optimizer(case_list, names, optimized_value) {
            if let Some(value) =
                fisher_optimizer(ncases, ncontrols, case_list, &node.indexes, optimized_value)
            {
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

    let top_node = hst.node_weight(top_node_idx).unwrap();

    let (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls) =
        case_ctrl_zygosity_from_node(top_node, case_list);

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
        nhet_cases.to_string(),
        nhom_cases.to_string(),
        nhet_ctrls.to_string(),
        nhom_ctrls.to_string(),
    ])
}

#[allow(dead_code)]
fn segregation_optimizer(
    case_list: &[String],
    names: Vec<String>,
    optimized_value: f64,
) -> Option<f64> {
    // Check that all samples are only from the wanted list and return their amount
    let nsamples = match names.iter().all(|name| case_list.contains(name)) {
        true => names.len() as f64,
        false => 0.0,
    };

    match nsamples > optimized_value {
        true => Some(nsamples),
        false => None,
    }
}

fn fisher_optimizer(
    ncases: usize,
    ncontrols: usize,
    case_list: &[usize],
    indexes: &[usize],
    optimized_value: f64,
) -> Option<f64> {
    let cases_true = indexes
        .iter()
        .filter(|name| case_list.contains(name))
        .count();

    let controls_true = indexes.len() - cases_true;

    let cases_false = ncases - cases_true;
    let controls_false = ncontrols - controls_true;

    // println!(
    //     "{:?} {:?} {:?} {:?} {:?} {:?} {:?} {:?}",
    //     hsts.metadata.samples.len(),
    //     hsts.nhaplotypes(),
    //     ncases,
    //     ncontrols,
    //     cases_true,
    //     cases_false,
    //     controls_true,
    //     controls_false
    // );

    let pvalue = fishers_exact(&[
        cases_true as u32,
        controls_true as u32,
        cases_false as u32,
        controls_false as u32,
    ])
    .unwrap()
    .greater_pvalue;

    match pvalue < optimized_value {
        true => Some(pvalue),
        false => None,
    }
}

pub fn case_ctrl_zygosity_from_node(
    node: &Node,
    case_list: &[usize],
) -> (usize, usize, usize, usize) {
    let controls = node
        .indexes
        .iter()
        .filter(|i| !case_list.contains(i))
        .copied()
        .collect::<Vec<usize>>();

    let cases = node
        .indexes
        .iter()
        .filter(|i| case_list.contains(i))
        .copied()
        .collect::<Vec<usize>>();

    let (nhet_ctrls, nhom_ctrls) = calculate_zygote_n(controls);
    let (nhet_cases, nhom_cases) = calculate_zygote_n(cases);

    (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls)
}

fn calculate_zygote_n(mut indexes: Vec<usize>) -> (usize, usize) {
    let prior_len = indexes.len();
    indexes.sort();
    indexes.dedup();
    let post_len = indexes.len();

    let nhomozygotes = prior_len - post_len;
    let nheterozygotes = post_len - nhomozygotes;

    (nheterozygotes, nhomozygotes)
}
