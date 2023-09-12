use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use haptk::core::open_csv_writer;
use haptk::structs::PhasedMatrix;
use haptk::subcommands::bhst::write_haplotype;
use haptk::subcommands::bhst::Node;
use haptk::subcommands::hst_gwas::Assoc;

pub fn assoc_vec_to_bm(vec: Vec<Box<dyn Assoc + Send>>) -> BTreeMap<u64, Box<dyn Assoc + Send>> {
    let mut hm = BTreeMap::new();
    for i in vec {
        hm.insert(i.pos(), i);
    }
    hm
}

pub fn write_indexes_to_file(
    node: &Node,
    vcf: Arc<PhasedMatrix>,
    tree_idx: usize,
    ctrls_true: usize,
) {
    if let Err(err) = std::fs::create_dir("results_gwas_haplotypes") {
        tracing::warn!("Error creating dir {err}");
    }
    let path = PathBuf::from(format!(
        "results_gwas_haplotypes/{}_{}_{}_{}_{}.ids",
        vcf.get_pos(tree_idx),
        tree_idx,
        node.stop_idx - node.start_idx + 1,
        node.indexes.len(),
        ctrls_true
    ));

    let mut writer = open_csv_writer(path).unwrap();
    for i in &node.indexes {
        let name = vcf.get_sample_name(*i);
        writer.write_record(vec![name]).unwrap();
    }
}

pub fn write_tp_ht_to_file(
    node: &Node,
    vcf: Arc<PhasedMatrix>,
    tree_idx: usize,
    ctrls_true: usize,
) {
    let sample_idx = node.indexes.get(0).unwrap();
    let ht = vcf.find_haplotype_for_sample(
        vcf.get_contig(),
        node.start_idx..node.stop_idx + 1,
        *sample_idx,
    );
    if let Err(err) = std::fs::create_dir("results_gwas_haplotypes") {
        tracing::warn!("Error creating dir {err}");
    }
    let path = PathBuf::from(format!(
        "results_gwas_haplotypes/{}_{}_{}_{}_{}.csv",
        vcf.get_pos(tree_idx),
        tree_idx,
        node.stop_idx - node.start_idx + 1,
        node.indexes.len(),
        ctrls_true
    ));

    let writer = open_csv_writer(path).unwrap();
    write_haplotype(ht, writer).unwrap();
}
