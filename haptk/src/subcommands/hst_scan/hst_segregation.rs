use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use color_eyre::Result;
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::{open_csv_writer, push_to_output};
use crate::utils::parse_coords;

use super::hst_scan::Limits;
use super::{read_tree_file, return_assoc, write_assoc, AssocRow, HstScan, HstScanRow};

const HEADER: &[&str] = &[
    "CHR",
    "POS",
    "MarkerID",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "n",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    hst_scan_path: PathBuf,
    samples: PathBuf,
) -> Result<()> {
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
    push_to_output(&args, &mut output, "hst_segregation_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    let rower = |hsts: Arc<HstScan>, tree_idx, top_node_idx, optimized_value| -> AssocRow {
        // need these into HashMap
        // n_case_hom,
        // n_case_het,
        // let (n_case_hom, n_case_het) = find_homozygosity_no_ctrl(top_node, vcf.clone(), selection);
        AssocRow::new(
            hsts.clone(),
            tree_idx,
            top_node_idx,
            optimized_value,
            HashMap::new(),
        )
    };

    let start_value = f64::MIN;

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

                    let node = hst.node_weight(idx).unwrap();
                    let names = hsts.get_sample_names(&node.indexes);

                    // Check that all samples are only from the wanted list and return their amount
                    let nsamples = match names.iter().all(|name| wanted_ids.contains(name)) {
                        true => Some(names.len() as f64),
                        false => None,
                    };
                    if let Some(value) = nsamples {
                        if value > optimized_value {
                            optimized_value = value;
                            top_node_idx = idx;
                        }
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
