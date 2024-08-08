use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use color_eyre::{eyre::ensure, Result};
use fishers_exact::fishers_exact;
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::{open_csv_writer, push_to_output};
use crate::subcommands::hst_scan::{read_tree_file, HstScanRow};
use crate::utils::parse_coords;

use super::hst_scan::{return_assoc, write_assoc, AssocRow, Limits};
use super::HstScan;

const HEADER: &[&str] = &[
    "CHR",
    "POS",
    "MarkerID",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "opt_value",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    hst_scan_path: PathBuf,
    controls: Option<Vec<PathBuf>>,
) -> Result<()> {
    // ensure!(
    // controls.is_some(),
    // "Controls need to be provided by specifying control ids using --controls parameter."
    // );

    let hsts = read_tree_file(hst_scan_path)?;
    let hsts = Arc::new(hsts);

    let args = StandardArgs {
        file: PathBuf::new(),
        output: args.output,
        info_limit: hsts.metadata.info_limit,
        coords: hsts.metadata.contig.clone(),
        selection: hsts.metadata.selection.clone(),
        prefix: args.prefix,
        samples: None,
    };

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "hst_gwas", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    // let mut ctrls = read_multiple_sample_ids(&controls).unwrap().unwrap();
    // let mut samples = hsts.metadata.samples.clone();
    // samples.append(&mut ctrls);

    let write_bam = false;
    let assoc = return_assoc(hsts.clone(), &args, limits, write_bam, optimizer, rower);

    write_assoc(HEADER, assoc, writer)?;
    Ok(())
}

fn rower(
    hsts: Arc<HstScan>,
    tree_idx: usize,
    top_node_idx: NodeIndex,
    optimized_value: f64,
) -> AssocRow {
    // Insert to a HashMap for HST GWAS
    // cases_true,
    // ctrls_true,
    // n_case_hom,
    // n_ctrl_hom,
    // n_case_het,
    // n_ctrl_het,
    // af_case,
    // af_ctrl,

    AssocRow::new(
        hsts.clone(),
        tree_idx,
        top_node_idx,
        optimized_value,
        HashMap::new(),
    )
}

fn optimizer(
    hsts: Arc<HstScan>,
    hst_idx: usize,
    limits: Limits,
) -> Option<(usize, NodeIndex, f64)> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
    let mut top_node_idx = NodeIndex::new(0);

    // START VALUE
    let start_value = f64::MIN;

    let mut optimized_value = start_value;

    let hst: &HstScanRow = &hsts.hsts[hst_idx];
    let hst = &hst.hst;
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

            // HST GWAS
            // let (cases_true, controls_true, cases_false, controls_false) = find_from_node(node)
            // let p = fishers_exact(&[
            //     cases_true as u32,
            //     controls_true as u32,
            //     cases_false as u32,
            //     controls_false as u32,
            // ])
            // .unwrap();
            // let value = p.greater_pvalue;
            // let clause = value < optimized_value;

            // Longest haplotype
            let node = hst.node_weight(idx).unwrap();
            let value = node.haplotype.len() as f64;
            let clause = value > optimized_value;

            // OPTIMIZER CODE END

            if clause {
                optimized_value = value;
                top_node_idx = idx;
            }
        }
    }

    if optimized_value == start_value {
        println!("No optimized value for tree_idx {hst_idx:?}");
        return None;
    }

    Some((hst_idx, top_node_idx, optimized_value))
}
