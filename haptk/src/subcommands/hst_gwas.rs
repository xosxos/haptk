use std::borrow::Borrow;
use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;
use std::path::PathBuf;
use std::sync::mpsc::Sender;
use std::sync::Arc;

use color_eyre::{eyre::ensure, Result};
use fishers_exact::fishers_exact;
use ndarray::s;
use petgraph::graph::NodeIndex;
use petgraph::Graph;
use rayon::prelude::*;
use rust_htslib::bam::{
    header::HeaderRecord,
    record::{Aux, AuxArray, Cigar, CigarString},
    Format, Header, Record, Writer,
};

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::core::{open_csv_writer, parse_coords};
use crate::read_vcf::{
    read_vcf_to_ctrl_matrix, read_vcf_to_haploid_ctrl_matrix, read_vcf_to_matrix,
};
use crate::structs::{HapVariant, PhasedMatrix, Ploidy};
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::hst_scan::{read_tree_file, Node, TreeRow, Trees};
use crate::utils::push_to_output;

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    trees: PathBuf,
    controls: Option<PathBuf>,
    ctrl_file: Option<PathBuf>,
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

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "hst_gwas", "csv");
    let writer = open_csv_writer(output)?;

    let samples = metadata.samples;

    let (vcf, ctrl_vcf) = read_vcfs(&mut args, controls, ctrl_file, samples)?;

    let write_bam = true;

    let vcf = Arc::new(vcf);
    let tx = get_sender(write_bam, &args, vcf.clone());
    let assoc = return_binary_assoc(
        &trees,
        vcf,
        (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
        &ctrl_vcf,
        tx,
    );
    write_assoc(BinaryAssocRow::header(), assoc, writer)?;
    Ok(())
}

pub fn return_binary_assoc(
    trees: &[TreeRow],
    vcf: Arc<PhasedMatrix>,
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    ctrl_vcf: &PhasedMatrix,
    tx: Option<std::sync::mpsc::Sender<(Node, f64)>>,
) -> Vec<Box<dyn Assoc + Send>> {
    let vcfr = vcf.clone();
    let rower = |top_node, lowest_pvalue, top_node_idx, tree_idx| -> Box<dyn Assoc + Send> {
        Box::new(BinaryAssocRow::new(
            top_node,
            vcfr.clone(),
            ctrl_vcf,
            lowest_pvalue,
            top_node_idx,
            tree_idx,
        ))
    };

    let vcfr = vcf.clone();
    let minimizer = |node_idx, min_value, tree| -> Option<f64> {
        minimize_groups(node_idx, min_value, tree, vcfr.clone(), ctrl_vcf)
    };

    let start_value = f64::MAX;
    return_assoc(
        trees,
        (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
        vcf,
        minimizer,
        tx,
        rower,
        start_value,
    )
}

pub fn get_sender(
    write_bam: bool,
    args: &StandardArgs,
    vcf: Arc<PhasedMatrix>,
) -> Option<Sender<(Node, f64)>> {
    match write_bam {
        true => {
            let mut bam_output = args.output.clone();
            push_to_output(args, &mut bam_output, "hst_scan", "bam");
            let header = return_header(vcf.get_contig());
            let mut out = Writer::from_path(bam_output, &header, Format::Bam).unwrap();
            let (tx, rx) = std::sync::mpsc::channel();

            std::thread::spawn(move || {
                while let Ok((node, pvalue)) = rx.recv() {
                    let record = create_bam_record(&node, vcf.clone(), pvalue);
                    out.write(&record).unwrap();
                }
            });
            Some(tx)
        }
        false => None,
    }
}

pub fn read_vcfs(
    args: &mut StandardArgs,
    controls: Option<PathBuf>,
    ctrl_file: Option<PathBuf>,
    samples: Vec<String>,
) -> Result<(PhasedMatrix, PhasedMatrix)> {
    ensure!(controls.is_some() || ctrl_file.is_some(), "Controls need to be provided either by specifying control ids using --controls parameter or by giving a control vcf using --control-vcf parameter.");

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(args, contig, 0, Some((start, stop)), Some(samples))?
        }
        _ => read_vcf_to_matrix(args, contig, 0, None, Some(samples))?,
    };

    args.samples = controls;

    if let Some(ctrl_file) = ctrl_file {
        args.file = ctrl_file;
    }

    tracing::debug!("Reading in control data to a sample-variant matrix.");
    let ctrl_vcf = match (start, stop) {
        (Some(start), Some(stop)) => match args.selection {
            Selection::All => {
                read_vcf_to_ctrl_matrix(args, contig, 0, Some((start, stop)), vcf.coords())?
            }
            Selection::Haploid => {
                read_vcf_to_haploid_ctrl_matrix(args, contig, 0, Some((start, stop)), vcf.coords())?
            }
            _ => unreachable!(),
        },
        _ => match args.selection {
            Selection::All => read_vcf_to_ctrl_matrix(args, contig, 0, None, vcf.coords())?,
            Selection::Haploid => {
                read_vcf_to_haploid_ctrl_matrix(args, contig, 0, None, vcf.coords())?
            }
            _ => unreachable!(),
        },
    };
    Ok((vcf, ctrl_vcf))
}

pub fn return_assoc<'a, F, U>(
    trees: &'a [TreeRow],
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    vcf: Arc<PhasedMatrix>,
    optimizer: F,
    tx: Option<std::sync::mpsc::Sender<(Node, f64)>>,
    rower: U,
    start_value: f64,
) -> Vec<Box<dyn Assoc + Send>>
where
    F: Fn(NodeIndex, f64, &'a Graph<Node, u8>) -> Option<f64> + Sync + Send + Clone,
    U: Fn(&'a Node, f64, NodeIndex, usize) -> Box<dyn Assoc + Send> + Sync + Send,
{
    trees
        .par_iter()
        .filter_map(|tree_row| {
            optimize_tree(
                tree_row,
                (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
                vcf.clone(),
                optimizer.clone(),
                start_value,
            )
        })
        .map_with(tx, |tx, tuple| {
            if let Some(tx) = tx {
                let node = tuple.0.clone();
                tx.send((node, tuple.1)).unwrap();
            }
            tuple
        })
        .map(|(top_node, lowest_pvalue, top_node_idx, tree_idx)| {
            rower(top_node, lowest_pvalue, top_node_idx, tree_idx)
        })
        .collect()
}

pub fn optimize_tree<'a, F: Fn(NodeIndex, f64, &'a Graph<Node, u8>) -> Option<f64>>(
    tree_row: &'a TreeRow,
    (nmin_samples, nmax_samples, nmin_variants, nmax_variants): (usize, usize, usize, usize),
    vcf: Arc<PhasedMatrix>,
    optimizer: F,
    start_value: f64,
) -> Option<(&'a Node, f64, NodeIndex, usize)> {
    let mut top_node = tree_row.tree.node_weight(NodeIndex::new(0)).unwrap();
    let mut top_node_idx = NodeIndex::new(0);

    let mut optimized_value = start_value;

    let mut iterator = tree_row.tree.node_indices();

    // Skip first node because it contains no haplotypes
    iterator.next();

    for idx in iterator {
        let node = tree_row.tree.node_weight(idx).unwrap();

        if node.indexes.len() >= nmin_samples && node.indexes.len() <= nmax_samples {
            let index = node.indexes[0];
            let ht = vcf
                .matrix
                .slice(s![index, node.start_idx..node.stop_idx + 1]);

            if ht.len() >= nmin_variants && ht.len() <= nmax_variants {
                if let Some(value) = optimizer(idx, optimized_value, &tree_row.tree) {
                    optimized_value = value;
                    top_node_idx = idx;
                    top_node = node;
                }
            }
        }
    }
    if optimized_value == start_value {
        println!("found none at {}", vcf.get_pos(vcf.variant_idx()));
        return None;
    }
    Some((top_node, optimized_value, top_node_idx, tree_row.idx))
}

pub fn minimize_groups<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    minimized_value: f64,
    tree: &Graph<Node, u8>,
    vcf: Arc<PhasedMatrix>,
    ctrl_vcf: T,
) -> Option<f64> {
    let value = compare_groups(idx, tree, vcf, ctrl_vcf);

    match value < minimized_value {
        true => Some(value),
        false => None,
    }
}

pub fn compare_groups<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    tree: &Graph<Node, u8>,
    vcf: Arc<PhasedMatrix>,
    ctrl_vcf: T,
) -> f64 {
    let ctrl_vcf = ctrl_vcf.borrow();
    let node = tree.node_weight(idx).unwrap();

    let start_idx = node.start_idx;
    let stop_idx = node.stop_idx;
    let index = node.indexes[0];

    let ht = vcf.matrix.slice(s![index, start_idx..stop_idx + 1]);

    let mut hits = 0;
    for j in 0..ctrl_vcf.matrix.nrows() {
        let ctrl_slice = ctrl_vcf.matrix.slice(s![j, start_idx..stop_idx + 1]);
        if ht == ctrl_slice {
            hits += 1;
        }
    }

    let cases_true = node.indexes.len();
    let cases_false = vcf.matrix.nrows() - cases_true;
    let controls_false = ctrl_vcf.matrix.nrows() - hits;

    let p = fishers_exact(&[
        cases_true as u32,
        hits as u32,
        cases_false as u32,
        controls_false as u32,
    ])
    .unwrap();
    p.greater_pvalue
}

pub fn return_mbah_lengths(trees: &Vec<TreeRow>) -> Vec<(usize, Node)> {
    trees.par_iter().map(find_mbah_node).collect()
}

pub trait Assoc {
    fn pos(&self) -> u64;
    fn bp_len(&self) -> f32;
    fn marker_len(&self) -> usize;
    fn node_idx(&self) -> NodeIndex;
    fn start_idx(&self) -> usize;
    fn stop_idx(&self) -> usize;
    fn first_sample_idx(&self) -> usize;
    fn opt_value(&self) -> f64;
    fn show(&self) -> String;
    fn to_csv_row(&self) -> Vec<String>;
}

macro_rules! assoc {
    ($impl:ident) => {
        impl crate::subcommands::hst_gwas::Assoc for $impl {
            fn pos(&self) -> u64 {
                self.pos
            }
            fn bp_len(&self) -> f32 {
                self.bp_len
            }
            fn marker_len(&self) -> usize {
                self.marker_len
            }
            fn node_idx(&self) -> NodeIndex {
                self.node_idx
            }
            fn start_idx(&self) -> usize {
                self.start_idx
            }
            fn stop_idx(&self) -> usize {
                self.stop_idx
            }
            fn first_sample_idx(&self) -> usize {
                self.first_sample_idx
            }
            fn opt_value(&self) -> f64 {
                self.opt_value
            }
            fn show(&self) -> String {
                self.to_string()
            }
            fn to_csv_row(&self) -> Vec<String> {
                self.into()
            }
        }
    };
}
pub(crate) use assoc;

#[derive(Clone, Debug)]
pub struct BinaryAssocRow {
    pub contig: String,
    pub pos: u64,
    pub marker_id: String,
    pub bp_len: f32,
    pub marker_len: usize,
    pub start: u64,
    pub start_idx: usize,
    pub stop: u64,
    pub stop_idx: usize,
    pub af_case: f64,
    pub af_ctrl: f64,
    pub cases_true: usize,
    pub ctrls_true: usize,
    pub n_case_hom: usize,
    pub n_case_het: usize,
    pub n_ctrl_hom: usize,
    pub n_ctrl_het: usize,
    pub opt_value: f64,
    pub node_idx: NodeIndex,
    pub first_sample_idx: usize,
    // pub pos_idx: usize,
}

assoc!(BinaryAssocRow);
impl BinaryAssocRow {
    pub fn new(
        top_node: &Node,
        vcf: Arc<PhasedMatrix>,
        ctrl_vcf: &PhasedMatrix,
        opt_value: f64,
        top_node_idx: NodeIndex,
        tree_idx: usize,
    ) -> BinaryAssocRow {
        let ht = vcf.matrix.slice(s![
            top_node.indexes[0],
            top_node.start_idx..top_node.stop_idx + 1
        ]);

        let mut ctrls_true = 0;
        for j in 0..ctrl_vcf.matrix.nrows() {
            let ctrl_slice = ctrl_vcf
                .matrix
                .slice(s![j, top_node.start_idx..top_node.stop_idx + 1]);
            if ht == ctrl_slice {
                ctrls_true += 1;
            }
        }
        let (n_case_hom, n_ctrl_hom, n_case_het, n_ctrl_het) =
            find_homozygosity(top_node, vcf.clone(), ctrl_vcf, ctrls_true);

        let cases_true = top_node.indexes.len();

        let start_idx = top_node.start_idx;
        let stop_idx = top_node.stop_idx;
        let marker_len = stop_idx - start_idx + 1;

        let start = vcf.get_pos(start_idx);
        let stop = vcf.get_pos(stop_idx);

        let af_case = top_node.indexes.len() as f64 / vcf.matrix.nrows() as f64;
        let af_ctrl = ctrls_true as f64 / ctrl_vcf.matrix.nrows() as f64;

        BinaryAssocRow {
            contig: vcf.get_contig().to_string(),
            pos: vcf.get_pos(tree_idx),
            // pos_idx: tree_idx,
            marker_id: get_marker_id(top_node, vcf),
            bp_len: top_node.block_len,
            start_idx,
            stop_idx,
            start,
            stop,
            marker_len,
            af_case,
            af_ctrl,
            cases_true,
            ctrls_true,
            n_case_hom,
            n_case_het,
            n_ctrl_hom,
            n_ctrl_het,
            opt_value,
            node_idx: top_node_idx,
            first_sample_idx: top_node.indexes[0],
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
            "AF_ctrl",
            "N_case",
            "N_ctrl",
            "N_case_hom",
            "N_case_het",
            "N_ctrl_hom",
            "N_ctrl_het",
            "p.value",
        ]
    }
}

impl From<&BinaryAssocRow> for Vec<String> {
    fn from(row: &BinaryAssocRow) -> Vec<String> {
        vec![
            row.contig.to_string(),
            row.pos.to_string(),
            row.marker_id.to_string(),
            row.bp_len.to_string(),
            row.marker_len.to_string(),
            row.af_case.to_string(),
            row.af_ctrl.to_string(),
            row.cases_true.to_string(),
            row.ctrls_true.to_string(),
            row.n_case_hom.to_string(),
            row.n_case_het.to_string(),
            row.n_ctrl_hom.to_string(),
            row.n_ctrl_het.to_string(),
            row.opt_value.to_string(),
        ]
    }
}

impl std::fmt::Display for BinaryAssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!("pos: {}\nlen: {}\nstart: {}\nstop: {}\ncases true: {}\nctrls true: {}\npvalue: {}\nnmarkers: {}", self.pos, self.bp_len, self.start, self.stop, self.cases_true, self.ctrls_true, self.opt_value, self.marker_len);
        write!(f, "{line}")
    }
}

pub fn get_marker_id<T: Borrow<PhasedMatrix>>(top_node: &Node, vcf: T) -> String {
    let vcf = vcf.borrow();
    let start = vcf.get_pos(top_node.start_idx);
    let stop = vcf.get_pos(top_node.stop_idx);

    let ht = vcf
        .matrix
        .slice(s![
            top_node.indexes[0],
            top_node.start_idx..top_node.stop_idx + 1
        ])
        .to_vec();

    let ht = ht.iter().fold(
        format!("{}_{}_{}_", vcf.get_contig(), start, stop),
        |acc, e| format!("{acc}{e}"),
    );

    let mut hasher = DefaultHasher::new();
    hasher.write(ht.as_bytes());

    format!("{:02x}", hasher.finish())
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

pub fn find_mbah_node(tree_row: &TreeRow) -> (usize, Node) {
    let nodes = find_majority_nodes(&tree_row.tree);
    let mut iter = nodes.into_iter().rev();
    iter.next().unwrap();
    let (last_node, _last_node_idx) = iter.next().unwrap();
    (tree_row.idx, last_node.clone())
}

pub fn write_assoc(
    header: Vec<&str>,
    rows: Vec<Box<dyn Assoc + Send>>,
    mut writer: csv::Writer<Box<dyn std::io::Write>>,
) -> Result<()> {
    writer.write_record(header)?;

    rows.into_iter().try_for_each(|row| -> Result<()> {
        let record = row.to_csv_row();
        writer.write_record(record)?;
        Ok(())
    })
}

pub fn return_header(contig: &str) -> Header {
    let mut header = Header::new();
    let mut header_record = HeaderRecord::new(b"SQ");
    header_record.push_tag(b"SN", &contig);
    header_record.push_tag(b"LN", &0);
    header.push_record(&header_record);
    header
}

pub fn create_bam_record(top_node: &Node, vcf: Arc<PhasedMatrix>, lowest_pvalue: f64) -> Record {
    let ht = vcf.find_haplotype_for_sample(
        vcf.get_contig(),
        top_node.start_idx..=top_node.stop_idx,
        top_node.indexes[0],
    );

    let mut record = record_from_ht(ht);
    add_aux_data(&mut record, lowest_pvalue, top_node);
    record
}

fn add_aux_data(record: &mut Record, lowest_pvalue: f64, top_node: &Node) {
    // Add an integer field
    let aux_integer_field = Aux::Float(lowest_pvalue as f32);
    record.push_aux(b"PV", aux_integer_field).unwrap();

    let aux_array: &Vec<u32> = &top_node.indexes.iter().map(|v| *v as u32).collect();
    let aux_array: AuxArray<u32> = aux_array.into();
    let aux_array_field = Aux::ArrayU32(aux_array);
    record.push_aux(b"ID", aux_array_field).unwrap();
}

fn record_from_ht(ht: Vec<HapVariant>) -> Record {
    let start = ht.first().unwrap().pos;
    let stop = ht.last().unwrap().pos;

    let hash_str = ht
        .iter()
        .map(|ht| ht.genotype())
        .fold(format!("{}_{}_{}_", "chr9", start, stop), |acc, e| {
            format!("{acc}{e}")
        });

    let mut hasher = DefaultHasher::new();
    hasher.write(hash_str.as_bytes());

    let qname = format!("{:02x}", hasher.finish());

    let mut cigar = vec![];
    let mut seq = vec![];

    for (i, item) in ht.iter().enumerate() {
        let curr_pos = ht[i].pos;
        let last_pos = ht.get(i - 1).map_or(0, |last| last.pos);

        if last_pos != 0 {
            let distance = curr_pos as u32 - last_pos as u32 - 1;
            cigar.push(Cigar::RefSkip(distance));
        }

        seq.push(item.genotype().as_bytes()[0]);

        match item.gt {
            0 => cigar.push(Cigar::Equal(1)),
            1 => cigar.push(Cigar::Diff(1)),
            _ => unreachable!(),
        }
    }

    let mut record = Record::new();
    record.set_pos(start as i64);
    record.set(
        qname.as_bytes(),
        Some(&CigarString(cigar)),
        &seq,
        &vec![255; seq.len()],
    );
    record
}
