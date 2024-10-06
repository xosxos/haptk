use std::{collections::BTreeMap, path::PathBuf, sync::Arc};

use color_eyre::Result;
use fishers_exact::fishers_exact;
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_sample_ids};
use crate::structs::Coord;
use crate::subcommands::scan::{
    read_tree_file, return_assoc, write_assoc, AssocRow, Hst, HstScan, Limits,
};
use crate::utils::parse_coords;

use super::hst_scan::{case_ctrl_zygosity_from_node, top_node_from_hsts};

const HEADER: &[&str] = &[
    "contig",
    "pos",
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
pub fn run(args: ConciseArgs, limits: Limits, samples: PathBuf) -> Result<()> {
    let seg_samples = read_sample_ids(&Some(samples))?.unwrap();

    let hsts = read_tree_file(args.file)?;
    let hsts = Arc::new(hsts);

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

    let args = StandardArgs {
        file: PathBuf::new(),
        output: args.output,
        info_limit: None,
        coords: hsts.metadata.contig.clone(),
        selection: hsts.metadata.selection.clone(),
        prefix: args.prefix,
        samples: None,
    };

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "segregate_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    // Define the minimizer/maximiser function for the HST
    let optimizer = |hsts: Arc<HstScan>,
                     hst: &Hst,
                     coord: &Coord,
                     limits: Limits|
     -> Option<(Coord, NodeIndex, f64)> {
        optimizer_inner(hsts, hst, coord, limits, &seg_samples)
    };

    // Define what information is wanted for each CSV row
    let rower = |hsts: Arc<HstScan>, coord: &Coord, top_node_idx, optimized_value| -> AssocRow {
        rower_inner(hsts, coord, top_node_idx, optimized_value, &seg_samples)
    };

    let write_bam = false;
    let assoc = return_assoc(&hsts, &args, limits, write_bam, optimizer, rower);

    write_assoc(HEADER, assoc, writer)?;

    Ok(())
}

fn rower_inner(
    hsts: Arc<HstScan>,
    coord: &Coord,
    top_node_idx: NodeIndex,
    optimized_value: f64,
    case_list: &[String],
) -> AssocRow {
    let node = top_node_from_hsts(&hsts.hsts, coord, top_node_idx);

    let case_list_idx = hsts.get_sample_idxs(case_list).unwrap();

    let (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls) =
        case_ctrl_zygosity_from_node(node, &case_list_idx);

    let mut assoc_hm = BTreeMap::new();
    assoc_hm.insert("nhet_cases".to_string(), nhet_cases.to_string());
    assoc_hm.insert("nhom_cases".to_string(), nhom_cases.to_string());
    assoc_hm.insert("nhet_ctrls".to_string(), nhet_ctrls.to_string());
    assoc_hm.insert("nhom_ctrls".to_string(), nhom_ctrls.to_string());

    AssocRow::new(hsts.clone(), coord, top_node_idx, optimized_value, assoc_hm)
}

fn optimizer_inner(
    hsts: Arc<HstScan>,
    hst: &Hst,
    coord: &Coord,
    limits: Limits,
    case_list: &[String],
) -> Option<(Coord, NodeIndex, f64)> {
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
            let names = hsts.get_sample_names(&node.indexes);

            // if let Some(value) = segregation_optimizer(case_list, names, optimized_value) {
            if let Some(value) = fisher_optimizer(hsts.clone(), case_list, names, optimized_value) {
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

    Some((coord.clone(), top_node_idx, optimized_value))
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
    hsts: Arc<HstScan>,
    case_list: &[String],
    names: Vec<String>,
    optimized_value: f64,
) -> Option<f64> {
    let ncases = case_list.len() * *hsts.metadata.ploidy;
    let ncontrols = hsts.nhaplotypes() - ncases;

    let cases_true = names.iter().filter(|name| case_list.contains(name)).count();

    let controls_true = names
        .iter()
        .filter(|name| !case_list.contains(name))
        .count();

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
