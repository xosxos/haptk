use std::borrow::Borrow;
use std::path::PathBuf;
use std::{collections::HashMap, sync::Arc};

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use petgraph::{graph::NodeIndex, Graph};

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::core::{open_csv_writer, parse_coords};
use crate::io::read_variable_data_file;
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
    var_data: PathBuf,
    var_name: String,
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

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "hst_qt_scan", "csv");
    let writer = open_csv_writer(output)?;

    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let mut vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None)?
        }
        _ => read_vcf_to_matrix(&args, contig, 0, None, None)?,
    };

    // Set clinical data
    vcf.set_variable_data(read_variable_data_file(var_data)?)?;

    let write_bam = true;

    let vcf = Arc::new(vcf);
    let tx = get_sender(write_bam, &args, vcf.clone());
    let hm = construct_variable_hm(vcf.clone(), var_name)?;

    let vcfr = vcf.clone();
    let rower = |top_node, minimized_value, top_node_idx, tree_idx| -> Box<dyn Assoc + Send> {
        Box::new(QtAssocRow::new(
            top_node,
            vcfr.clone(),
            &args.selection,
            minimized_value,
            top_node_idx,
            tree_idx,
            &hm,
        ))
    };

    let vcfr = vcf.clone();
    let minimizer = |node_idx, min_value, tree| -> Option<f64> {
        minimize_traits(node_idx, min_value, tree, &hm, vcfr.clone())
    };

    let start_value = f64::MAX;
    let assoc = return_assoc(
        &trees,
        (nmin_samples, nmax_samples, nmin_variants, nmax_variants),
        vcf,
        minimizer,
        tx,
        rower,
        start_value,
    );

    write_assoc(QtAssocRow::header(), assoc, writer)?;

    Ok(())
}

pub fn minimize_traits<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    minimized_value: f64,
    tree: &Graph<Node, u8>,
    hm: &HashMap<usize, f64>,
    vcf: T,
) -> Option<f64> {
    let value = compare_traits(idx, tree, hm, vcf);
    match value < minimized_value {
        true => Some(value),
        false => None,
    }
}

pub fn compare_traits<T: Borrow<PhasedMatrix>>(
    idx: NodeIndex,
    tree: &Graph<Node, u8>,
    hm: &HashMap<usize, f64>,
    vcf: T,
) -> f64 {
    let node = tree.node_weight(idx).unwrap();
    let other_indexes: Vec<usize> = (0..vcf.borrow().nsamples())
        .filter(|i| !node.indexes.contains(i))
        .collect();

    let mut v1 = get_variable_data(&node.indexes, hm);
    let v2 = get_variable_data(&other_indexes, hm);

    let mut genotypes = vec![1.0; node.indexes.len()];
    genotypes.extend(vec![0.0; other_indexes.len()]);

    v1.extend(v2.iter());

    // crate::stats::linear_regression(v1, genotypes)
    crate::stats::two_tail_welch_t_test(&v1, &v2)
}

#[derive(Clone, Debug)]
pub struct QtAssocRow {
    pub contig: String,
    pub pos: u64,
    // pub pos_idx: usize,
    pub marker_id: String,
    pub bp_len: f32,
    pub marker_len: usize,
    pub start: u64,
    pub start_idx: usize,
    pub stop: u64,
    pub stop_idx: usize,
    pub af_case: f64,
    pub cases_true: usize,
    pub n_case_hom: usize,
    pub n_case_het: usize,
    pub opt_value: f64,
    pub cases_mean: f64,
    pub ctrls_mean: f64,
    pub first_sample_idx: usize,
    pub node_idx: NodeIndex,
    // pub pos_idx: usize,
}

hst_gwas::assoc!(QtAssocRow);
impl QtAssocRow {
    pub fn new(
        top_node: &Node,
        vcf: Arc<PhasedMatrix>,
        selection: &Selection,
        minimized_value: f64,
        top_node_idx: NodeIndex,
        tree_idx: usize,
        hm: &HashMap<usize, f64>,
    ) -> Self {
        let other_indexes: Vec<usize> = (0..vcf.nsamples())
            .filter(|i| !top_node.indexes.contains(i))
            .collect();

        let v1 = get_variable_data(&top_node.indexes, hm);
        let v2 = get_variable_data(&other_indexes, hm);

        let cases_mean = v1.iter().sum::<f64>() / v1.len() as f64;
        let ctrls_mean = v2.iter().sum::<f64>() / v2.len() as f64;

        let cases_true = top_node.indexes.len();
        let start_idx = top_node.start_idx;
        let stop_idx = top_node.stop_idx;
        let marker_len = stop_idx - start_idx + 1;
        let af_case = cases_true as f64 / vcf.matrix.nrows() as f64;
        let (n_case_hom, n_case_het) = find_homozygosity_no_ctrl(top_node, vcf.clone(), selection);

        Self {
            contig: vcf.get_contig().to_string(),
            pos: vcf.get_pos(tree_idx),
            marker_id: get_marker_id(top_node, vcf.clone()),
            bp_len: top_node.block_len,
            marker_len,
            start: vcf.get_pos(start_idx),
            stop: vcf.get_pos(stop_idx),
            start_idx,
            stop_idx,
            af_case,
            cases_true,
            n_case_hom,
            n_case_het,
            opt_value: minimized_value,
            cases_mean,
            ctrls_mean,
            node_idx: top_node_idx,
            first_sample_idx: top_node.indexes[0],
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
            "cases_mean",
            "ctrls_mean",
            "p.value",
        ]
    }
}

impl From<&QtAssocRow> for Vec<String> {
    fn from(row: &QtAssocRow) -> Vec<String> {
        vec![
            row.contig.to_string(),
            format!("{}", row.pos),
            row.marker_id.to_string(),
            format!("{}", row.bp_len),
            format!("{}", row.marker_len),
            format!("{}", row.af_case),
            format!("{}", row.cases_true),
            format!("{}", row.n_case_hom),
            format!("{}", row.n_case_het),
            format!("{}", row.cases_mean),
            format!("{}", row.ctrls_mean),
            format!("{}", row.opt_value),
        ]
    }
}

impl std::fmt::Display for QtAssocRow {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "pos: {}\nlen: {}\nstart: {}\nstop: {}\ncases true: {}\npvalue: {}\nnmarkers: {}",
            self.pos,
            self.bp_len,
            self.start,
            self.stop,
            self.cases_true,
            self.opt_value,
            self.marker_len,
        );
        write!(f, "{line}")
    }
}

pub fn construct_variable_hm(
    vcf: Arc<PhasedMatrix>,
    col_name: String,
) -> Result<HashMap<usize, f64>> {
    let mut hm = HashMap::new();

    let df = vcf.clinical_data().as_ref().unwrap();

    let data: Result<Vec<f64>> = if let Ok(vec) = df[col_name.as_str()].f64() {
        Ok(vec.into_iter().flatten().collect())
    } else if let Ok(vec) = df[col_name.as_str()].i64() {
        Ok(vec.into_iter().flatten().map(|v| v as f64).collect())
    } else {
        Err(eyre!(
                "Error processing the data column {col_name:?}. The data is not floats or integers format"
            ))
    };

    let data = data?;
    let ids: Vec<&str> = df["parsed_id"].utf8()?.into_iter().flatten().collect();

    for (id, value) in ids.iter().zip(data.iter()) {
        hm.insert(vcf.get_sample_idx(id)?, *value);
    }
    Ok(hm)
}

fn get_variable_data(indexes: &[usize], hm: &HashMap<usize, f64>) -> Vec<f64> {
    indexes.iter().filter_map(|v| hm.get(v)).copied().collect()
}
