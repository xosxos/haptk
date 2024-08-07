use std::path::PathBuf;
use std::sync::Arc;
use std::{borrow::Borrow, collections::HashMap};

use color_eyre::{eyre::ensure, Result};
use fishers_exact::fishers_exact;
use ndarray::s;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output};
use crate::structs::{PhasedMatrix, Ploidy};
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::bhst::Node;
use crate::subcommands::hst_scan::{read_tree_file, HstScanRow};
use crate::utils::parse_coords;

use super::hst_scan::{return_assoc, write_assoc, AssocRow};

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
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    hst_scan_path: PathBuf,
    controls: Option<Vec<PathBuf>>,
) -> Result<()> {
    // ensure!(
    // controls.is_some(),
    // "Controls need to be provided by specifying control ids using --controls parameter."
    // );

    let hsts = read_tree_file(hst_scan_path)?;

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

    let hsts = Arc::new(hsts);

    let rower = |tree_idx, top_node_idx, optimized_value| -> AssocRow {
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
    };

    let optimizer = |hst_idx| -> Option<(usize, NodeIndex, f64)> {
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
    };

    let assoc = return_assoc(hsts.clone(), &args, write_bam, optimizer, rower);

    write_assoc(HEADER, assoc, writer)?;
    Ok(())
}

pub fn find_homozygosity<T: Borrow<PhasedMatrix>>(
    top_node: &Node,
    vcf: T,
    ctrl_vcf: &PhasedMatrix,
    ctrls_true: usize,
) -> (usize, usize, usize, usize) {
    let vcf = vcf.borrow();
    let (n_case_hom, n_ctrl_hom, n_case_het, n_ctrl_het) = match vcf.ploidy {
        Ploidy::Diploid => {
            let mut n_case_hom = 0;
            for i in &top_node.indexes {
                if top_node.indexes.contains(&(i + vcf.nsamples() / 2)) {
                    n_case_hom += 1;
                }
            }

            let ht = vcf.matrix.slice(s![
                top_node.indexes[0],
                top_node.start_idx..top_node.stop_idx + 1
            ]);

            let mut n_ctrl_hom = 0;
            for j in (0..ctrl_vcf.matrix.nrows()).step_by(2) {
                let ctrl_slice = ctrl_vcf
                    .matrix
                    .slice(s![j, top_node.start_idx..=top_node.stop_idx]);
                if ht == ctrl_slice {
                    let other_allele = ctrl_vcf.matrix.slice(s![
                        j + 1,
                        top_node.start_idx..=top_node.stop_idx
                    ]);
                    if ht == other_allele {
                        n_ctrl_hom += 1;
                    }
                }
            }

            let n_case_het = top_node.indexes.len() - n_case_hom * 2;
            let n_ctrl_het = ctrls_true - n_ctrl_hom * 2;
            (n_case_hom, n_ctrl_hom, n_case_het, n_ctrl_het)
        }
        // haploid so only homozygotes
        Ploidy::Haploid => (top_node.indexes.len(), ctrls_true, 0, 0),
        Ploidy::Mixed => panic!("Calculating homozygozity after a selection of a variant with --select only-refs or only-alts is not supported"),
    };
    (n_case_hom, n_ctrl_hom, n_case_het, n_ctrl_het)
}

pub fn find_homozygosity_no_ctrl<T: Borrow<PhasedMatrix>>(
    top_node: &Node,
    vcf: T,
    selection: &Selection,
) -> (usize, usize) {
    let (n_case_hom, n_case_het) = if selection == &Selection::All {
        let mut n_case_hom = 0;
        for i in &top_node.indexes {
            if top_node
                .indexes
                .contains(&(i + vcf.borrow().nsamples() / 2))
            {
                n_case_hom += 1;
            }
        }
        let n_case_het = top_node.indexes.len() - n_case_hom * 2;
        (n_case_hom, n_case_het)
    } else {
        // haploid so only homozygotes
        (top_node.indexes.len(), 0)
    };
    (n_case_hom, n_case_het)
}

pub fn return_mbah_lengths(hsts: &Vec<HstScanRow>) -> Vec<(usize, Node)> {
    hsts.par_iter().map(find_mbah_node).collect()
}

pub fn find_mbah_node(hst_scan_row: &HstScanRow) -> (usize, Node) {
    let nodes = find_majority_nodes(&hst_scan_row.hst);
    let mut iter = nodes.into_iter().rev();
    iter.next().unwrap();
    let (last_node, _last_node_idx) = iter.next().unwrap();
    (hst_scan_row.idx, last_node.clone())
}
