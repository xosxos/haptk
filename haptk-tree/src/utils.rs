use std::path::PathBuf;
use std::rc::Rc;

use haptk::{
    io::{get_output, get_vcf_writer, open_csv_writer},
    structs::{PhasedMatrix, Ploidy},
    subcommands::bhst::{write_haplotype, Node},
    subcommands::haplotype_to_vcf::{header, VcfRow},
};

pub fn write_indexes_to_file<S>(name: S, node: &Node, vcf: Rc<PhasedMatrix>)
where
    S: AsRef<str>,
{
    let path = PathBuf::from(format!(
        "results_tree_haplotypes/{}_{}_{}.ids",
        name.as_ref(),
        vcf.variant_idx_pos(),
        node.indexes.len(),
    ));

    let mut writer = open_csv_writer(path).unwrap();
    for i in &node.indexes {
        let name = vcf.get_sample_name(*i);
        writer.write_record(vec![name]).unwrap();
    }
}

pub fn write_tp_ht_to_csv<S>(name: S, node: &Node, vcf: Rc<PhasedMatrix>)
where
    S: AsRef<str>,
{
    let (start_idx, stop_idx) = match node.start_idx > node.stop_idx {
        true => (node.stop_idx, node.start_idx),
        false => (node.start_idx, node.stop_idx),
    };

    let sample_idx = node.indexes.first().unwrap();
    let ht = vcf.find_haplotype_for_sample(vcf.get_contig(), start_idx..=stop_idx, *sample_idx);
    let path = PathBuf::from(format!(
        "results_tree_haplotypes/{}_{}_{}.csv",
        name.as_ref(),
        vcf.variant_idx_pos(),
        node.indexes.len(),
    ));

    let writer = open_csv_writer(path).unwrap();
    write_haplotype(ht, writer).unwrap();
}

pub fn _write_tp_ht_to_vcf<S, T>(name: S, node: &Node, vcf: Rc<PhasedMatrix>)
where
    S: AsRef<str>,
{
    let path = PathBuf::from(format!(
        "results_tree_haplotypes/{}_{}_{}.vcf",
        name.as_ref(),
        vcf.variant_idx_pos(),
        node.indexes.len(),
    ));

    let (samples, gts) = match vcf.ploidy {
        Ploidy::Haploid => {
            let samples = vcf.get_sample_names(&node.indexes);
            let gts = vec![(1, 1); samples.len()];
            (samples, gts)
        }

        Ploidy::Diploid | Ploidy::Mixed => {
            let mut samples = vcf.get_sample_names(&node.indexes);
            samples.sort_unstable();
            let (u, d) = samples.partition_dedup();

            let gts = u
                .iter()
                .map(|u| match d.contains(u) {
                    true => (1, 1),
                    false => (1, 0),
                })
                .collect::<Vec<(u8, u8)>>();

            (u.to_vec(), gts)
        }
    };

    let pos = vcf.variant_idx_pos();
    let seqid = vcf.get_contig().to_string();
    let reference = "ref".to_string();
    let alt = "ht".to_string();
    let row = VcfRow::new(seqid, pos, reference, alt, gts, vcf.ploidy.clone());

    let header = header(samples, vcf.get_contig());

    let output = get_output(Some(path)).unwrap();
    let mut wrtr = get_vcf_writer(output);

    for line in header {
        wrtr.write_record(vec![line]).unwrap();
    }

    wrtr.write_record(Into::<Vec<String>>::into(row)).unwrap();
}
