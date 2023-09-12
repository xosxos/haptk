use std::borrow::Borrow;
use std::path::PathBuf;
use std::{collections::BTreeMap, sync::Arc};

use color_eyre::{
    eyre::{eyre, WrapErr},
    Result,
};
use petgraph::graph::NodeIndex;
use petgraph::{Direction, Graph};

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::core::{open_csv_writer, parse_coords};
use crate::io::read_recombination_file;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::hst_scan::Trees;
use crate::utils::push_to_output;

use super::hst_gwas::{self, get_marker_id, get_sender, return_assoc, write_assoc, Assoc};
use super::{bhst::Node, hst_gwas::find_homozygosity_no_ctrl};

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    trees: PathBuf,
    rec_rates: PathBuf,
) -> Result<()> {
    let Trees { metadata, trees } = read_tree_file(trees)?;
    let mut args = StandardArgs {
        file: args.file,
        output: args.output,
        info_limit: metadata.info_limit,
        coords: metadata.coords,
        selection: metadata.selection,
        prefix: args.prefix,
        samples: None,
    };

    let rec_rates = read_recombination_file(rec_rates)?;

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "hst_mrca_scan", "csv");
    let writer = open_csv_writer(output)?;

    let samples = metadata.samples;

    let vcf = read_vcfs(&mut args, samples)?;

    let write_bam = true;
    let vcf = Arc::new(vcf);
    let tx = get_sender(write_bam, &args, vcf.clone());

    let vcfr = vcf.clone();
    let rower = |top_node, minimized_value, top_node_idx, tree_idx| -> Box<dyn Assoc + Send> {
        Box::new(MrcaAssocRow::new(
            top_node,
            vcfr.clone(),
            &args.selection,
            minimized_value,
            top_node_idx,
            tree_idx,
        ))
    };

    let vcfr = vcf.clone();
    let minimizer = |node_idx, min_value, tree| -> Option<f64> {
        minimize_mrca(node_idx, min_value, tree, &rec_rates, vcfr.clone())
    };

    let assoc = return_assoc(
        &trees,
        (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
        vcf,
        minimizer,
        tx,
        rower,
    );

    write_assoc(MrcaAssocRow::header(), assoc, writer)?;
    Ok(())
}

pub fn read_vcfs(args: &mut StandardArgs, samples: Vec<String>) -> Result<PhasedMatrix> {
    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(args, contig, 0, Some((start, stop)), Some(samples))?
        }
        _ => read_vcf_to_matrix(args, contig, 0, None, Some(samples))?,
    };

    Ok(vcf)
}

pub fn read_tree_file(path: PathBuf) -> Result<Trees> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;
    let rows = serde_json::from_reader(reader)?;

    Ok(rows)
}

#[derive(Clone, Debug)]
pub struct MrcaAssocRow {
    pub contig: String,
    pub pos: u64,
    pub marker_id: String,
    pub bp_len: f32,
    pub marker_len: usize,
    pub start: u64,
    pub start_idx: usize,
    pub stop: u64,
    pub stop_idx: usize,
    pub cases_true: usize,
    pub n_case_hom: usize,
    pub n_case_het: usize,
    pub opt_value: f64,
    pub first_sample_idx: usize,
    pub node_idx: NodeIndex,
    // pub pos_idx: usize,
}

hst_gwas::assoc!(MrcaAssocRow);
impl MrcaAssocRow {
    pub fn new(
        top_node: &Node,
        vcf: Arc<PhasedMatrix>,
        selection: &Selection,
        minimized_value: f64,
        top_node_idx: NodeIndex,
        tree_idx: usize,
    ) -> Self {
        let (n_case_hom, n_case_het) = find_homozygosity_no_ctrl(top_node, vcf.clone(), selection);

        Self {
            contig: vcf.get_contig().to_string(),
            marker_id: get_marker_id(top_node, vcf.clone()),
            pos: vcf.get_pos(tree_idx),
            bp_len: top_node.block_len,
            marker_len: top_node.stop_idx - top_node.start_idx + 1,
            start: vcf.get_pos(top_node.start_idx),
            stop: vcf.get_pos(top_node.stop_idx),
            start_idx: top_node.start_idx,
            stop_idx: top_node.stop_idx,
            cases_true: top_node.indexes.len(),
            n_case_hom,
            n_case_het,
            opt_value: minimized_value,
            first_sample_idx: top_node.indexes[0],
            node_idx: top_node_idx,
            // pos_idx: tree_row.idx,
        }
    }

    pub fn header<'a>() -> Vec<&'a str> {
        vec![
            "CHR",
            "POS",
            "MarkerID",
            "bp_len",
            "marker_len",
            "AF_case",
            "N_case",
            "N_case_hom",
            "N_case_het",
            "mrca",
            "SampleID",
        ]
    }
}

impl From<&MrcaAssocRow> for Vec<String> {
    fn from(row: &MrcaAssocRow) -> Vec<String> {
        vec![
            row.contig.to_string(),
            format!("{}", row.pos),
            row.marker_id.to_string(),
            format!("{}", row.bp_len),
            format!("{}", row.marker_len),
            format!("{}", row.cases_true),
            format!("{}", row.n_case_hom),
            format!("{}", row.n_case_het),
            format!("{}", row.opt_value),
        ]
    }
}

impl std::fmt::Display for MrcaAssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "pos: {}\nlen: {}\nstart: {}\nstop: {}\ncases true: {}\nmrca: {}\nnmarkers: {}",
            self.pos,
            self.bp_len,
            self.start,
            self.stop,
            self.cases_true,
            self.opt_value,
            self.marker_len
        );
        write!(f, "{line}")
    }
}

pub struct Breakpoint {
    pub sample: String,
    pub start: u64,
    pub stop: u64,
}

pub fn minimize_mrca<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    minimized_value: f64,
    tree: &Graph<Node, u8>,
    rec_rates: &BTreeMap<u64, f32>,
    vcf: T,
) -> Option<f64> {
    let mrca = calculate_mrca(idx, tree, rec_rates, vcf);
    match mrca < minimized_value {
        true => Some(mrca),
        false => None,
    }
}

pub fn calculate_mrca<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    tree: &Graph<Node, u8>,
    rec_rates: &BTreeMap<u64, f32>,
    vcf: T,
) -> f64 {
    let vcf = vcf.borrow();
    let mut idxs = vec![idx];

    loop {
        let children = tree.neighbors_directed(*idxs.last().unwrap(), Direction::Outgoing);

        if children.clone().count() < 2 {
            break;
        }

        let new_idx = children
            .max_by_key(|i| tree.node_weight(*i).unwrap().indexes.len())
            .unwrap();

        idxs.push(new_idx);
    }

    let mut breakpoints = vec![];
    let mut used_indexes = vec![];

    let mut nodes_iter = idxs.iter().rev();
    let mut prev_idx = nodes_iter.next().unwrap();

    for idx in nodes_iter {
        let node = tree.node_weight(*idx).unwrap();
        let prev_node = tree.node_weight(*prev_idx).unwrap();

        for idx in &node.indexes {
            if !used_indexes.contains(&idx) {
                let tuple = Breakpoint {
                    sample: vcf.get_sample_name(*idx),
                    start: vcf.get_pos(prev_node.start_idx),
                    stop: vcf.get_pos(prev_node.stop_idx),
                };
                breakpoints.push(tuple);
                used_indexes.push(idx);
            }
        }
        prev_idx = idx;
    }

    bhst_mrca_independent(breakpoints, rec_rates).unwrap() as f64
}

///
/// The original R algorithm by Gandolfo et al translated to Rust.
/// <https://github.com/bahlolab/DatingRareMutations>
///
// fn bhst_mrca_independent(
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
