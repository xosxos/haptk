use std::collections::HashMap;
use std::path::PathBuf;
use std::{collections::BTreeMap, sync::Arc};

use color_eyre::Result;
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::read_recombination_file;
use crate::io::{open_csv_writer, push_to_output};
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::hst_scan::{
    read_tree_file, return_assoc, write_assoc, AssocRow, HstScan, HstScanRow,
};
use crate::utils::parse_coords;

use super::hst_scan::Limits;

const HEADER: &[&str] = &[
    "CHR",
    "POS",
    "MarkerID",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "mrca",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    hst_scan_path: PathBuf,
    rec_rates: PathBuf,
) -> Result<()> {
    let rec_rates = read_recombination_file(rec_rates)?;
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
    push_to_output(&args, &mut output, "hst_mrca_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    let rower = |hsts: Arc<HstScan>, tree_idx, top_node_idx, optimized_value| -> AssocRow {
        AssocRow::new(
            hsts.clone(),
            tree_idx,
            top_node_idx,
            optimized_value,
            HashMap::new(),
        )
    };

    let optimizer =
        |hsts: Arc<HstScan>, hst_idx: usize, limits: Limits| -> Option<(usize, NodeIndex, f64)> {
            let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
            let mut top_node_idx = NodeIndex::new(0);

            // START VALUE
            let start_value = f64::MAX;

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

                    // Lowest MRCA
                    let start_node = idx;
                    let nodes = find_majority_nodes(&hst, start_node);

                    // Iterate nodes from the leaf towards the root
                    // NOTE: Sometimes genotyping data runs out and leaf nodes have more than 1 sample
                    // NOTE: This is not really accounted for in the Gamma method without chromosome length corrections
                    let mut nodes_iter = nodes.iter().map(|(_, idx)| idx).rev();
                    let mut prev_idx = nodes_iter.next().unwrap();

                    let (mut breakpoints, mut used_indexes) = (vec![], vec![]);
                    for idx in nodes_iter {
                        let node = hst.node_weight(*idx).unwrap();
                        let prev_node = hst.node_weight(*prev_idx).unwrap();

                        for idx in &node.indexes {
                            if !used_indexes.contains(&idx) {
                                let tuple = Breakpoint {
                                    start: hsts.get_pos(prev_node.start_idx),
                                    stop: hsts.get_pos(prev_node.stop_idx),
                                };
                                breakpoints.push(tuple);
                                used_indexes.push(idx);
                            }
                        }
                        prev_idx = idx;
                    }

                    // Calculate MRCA estimate
                    let value = bhst_mrca_independent(breakpoints, &rec_rates).unwrap() as f64;
                    let clause = value < optimized_value;

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
        };

    let write_bam = false;
    let assoc = return_assoc(hsts.clone(), &args, limits, write_bam, optimizer, rower);

    write_assoc(HEADER, assoc, writer)?;
    Ok(())
}

pub struct Breakpoint {
    pub start: u64,
    pub stop: u64,
}

///
/// The original R algorithm by Gandolfo et al translated to Rust without confidence intervals
/// <https://github.com/bahlolab/DatingRareMutations>
///
pub fn bhst_mrca_independent(
    breakpoints: Vec<Breakpoint>,
    rec_rates: &BTreeMap<u64, f32>,
) -> Result<f32> {
    let lr_sum = breakpoints
        .iter()
        .map(|point| {
            let stop_cm = match rec_rates.range(point.stop..).next() {
                Some(cm) => cm,
                None => rec_rates.range(..point.stop).next_back().unwrap(),
            };

            let start_cm = match rec_rates.range(..point.start).next_back() {
                Some(cm) => cm,
                None => rec_rates.range(point.start..).next().unwrap(),
            };

            (stop_cm.1 - start_cm.1) / 100.
        })
        .sum::<f32>();

    let n = breakpoints.len() as f32;
    let b_c = (2.0 * n - 1.0) / (2.0 * n);
    let i_tau_hat = (b_c * 2.0 * n) / lr_sum;

    Ok(i_tau_hat)
}
