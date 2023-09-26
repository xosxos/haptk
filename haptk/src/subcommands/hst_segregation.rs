use std::borrow::Borrow;
use std::path::PathBuf;
use std::sync::Arc;

use color_eyre::Result;
use petgraph::graph::NodeIndex;
use petgraph::Graph;

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::core::{open_csv_writer, parse_coords};
use crate::io::read_sample_ids;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::hst_gwas::{
    self, find_homozygosity_no_ctrl, get_marker_id, get_sender, return_assoc, write_assoc, Assoc,
};
use crate::subcommands::hst_scan::{read_tree_file, Node, Trees};
use crate::utils::push_to_output;

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    trees: PathBuf,
    samples: PathBuf,
) -> Result<()> {
    let Trees { metadata, trees } = read_tree_file(trees)?;

    let args = StandardArgs {
        file: args.file,
        output: args.output,
        info_limit: metadata.info_limit,
        coords: metadata.coords,
        selection: metadata.selection,
        prefix: args.prefix,
        samples: None,
    };
    let seg_samples = read_sample_ids(&Some(samples))?.unwrap();

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "hst_segregation_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => read_vcf_to_matrix(
            &args,
            contig,
            0,
            Some((start, stop)),
            Some(metadata.samples),
        )?,
        _ => read_vcf_to_matrix(&args, contig, 0, None, Some(metadata.samples))?,
    };

    let write_bam = false;
    let vcf = Arc::new(vcf);
    let tx = get_sender(write_bam, &args, vcf.clone());

    let vcfr = vcf.clone();
    let rower = |top_node, maximized_value, top_node_idx, tree_idx| -> Box<dyn Assoc + Send> {
        Box::new(SegAssocRow::new(
            top_node,
            vcfr.clone(),
            &args.selection,
            maximized_value,
            top_node_idx,
            tree_idx,
        ))
    };

    let vcfr = vcf.clone();
    let maximizer = |node_idx, max_value, tree| -> Option<f64> {
        maximize_wanted_samples(node_idx, max_value, tree, &seg_samples, vcfr.clone())
    };

    let start_value = f64::MIN;
    let assoc = return_assoc(
        &trees,
        (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
        vcf,
        maximizer,
        tx,
        rower,
        start_value,
    );

    write_assoc(SegAssocRow::header(), assoc, writer)?;
    Ok(())
}

#[derive(Clone, Debug)]
pub struct SegAssocRow {
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

hst_gwas::assoc!(SegAssocRow);
impl SegAssocRow {
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
            "N_case",
            "N_case_hom",
            "N_case_het",
            "n",
        ]
    }
}

impl From<&SegAssocRow> for Vec<String> {
    fn from(row: &SegAssocRow) -> Vec<String> {
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

impl std::fmt::Display for SegAssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "pos: {}\nlen: {}\nstart: {}\nstop: {}\ncases true: {}\nn: {}\nnmarkers: {}",
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

pub fn maximize_wanted_samples<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    optimized_value: f64,
    tree: &Graph<Node, u8>,
    wanted_ids: &Vec<String>,
    vcf: T,
) -> Option<f64> {
    let vcf = vcf.borrow();

    let node = tree.node_weight(idx).unwrap();
    let names = vcf.get_sample_names(&node.indexes);

    // Check that all samples are only from the wanted list and return their amount
    let nsamples = match names.iter().all(|name| wanted_ids.contains(name)) {
        true => Some(names.len() as f64),
        false => None,
    };

    // Check if n is larger than the biggest value up until now, and return the value if it is
    nsamples.map_or(None, |n| match n > optimized_value {
        true => Some(n),
        false => None,
    })
}
