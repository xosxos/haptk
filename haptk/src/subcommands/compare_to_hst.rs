use std::{collections::HashMap, path::PathBuf};

use color_eyre::{eyre::eyre, Result};
use ndarray::parallel::prelude::*;
use petgraph::Graph;
use petgraph::{prelude::NodeIndex, Direction};

use crate::{
    args::{Selection, StandardArgs},
    io::push_to_output,
    read_vcf::read_vcf_to_matrix,
    structs::PhasedMatrix,
    subcommands::bhst::{read_hst_file, write_hst_file, Hst, Node},
    utils::parse_snp_coord,
};

#[doc(hidden)]
pub fn run(args: StandardArgs, hst_path: PathBuf) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    // File reads
    let hst_import = read_hst_file(hst_path)?;

    // VCF read
    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;
    let (start, end) = (
        hst_import.metadata.coords.first().unwrap(),
        hst_import.metadata.coords.last().unwrap(),
    );
    tracing::info!("Reading coord range: {} - {}", start, end);

    let vcf = match args.selection {
        Selection::All | Selection::Haploid | Selection::OnlyAlts | Selection::OnlyRefs => {
            read_vcf_to_matrix(
                &args,
                contig,
                variant_pos,
                Some((Some(start.pos), Some(end.pos))),
                None,
                None,
            )?
        }
        Selection::OnlyLongest => {
            let mut vcf = read_vcf_to_matrix(&args, contig, variant_pos, None, None, None)?;
            vcf.select_only_longest_no_shard()?;

            vcf
        }
        Selection::Unphased => unreachable!(),
    };

    if hst_import.metadata.coords.len() > vcf.ncoords() {
        tracing::warn!(
            "The HST has more variants than the given VCF {} vs {}",
            hst_import.metadata.coords.len(),
            vcf.ncoords()
        );
    }

    if hst_import.metadata.coords.len() < vcf.ncoords() {
        tracing::warn!(
            "The HST has less variants than the given VCF {} vs {}",
            hst_import.metadata.coords.len(),
            vcf.ncoords()
        );
    }

    let count = hst_import
        .metadata
        .coords
        .iter()
        .filter(|r| {
            vcf.coords()
                .iter()
                .any(|s| s.reference == r.reference && s.alt == r.alt && s.pos == r.pos)
        })
        .count();

    tracing::info!(
        "The HST and the given VCF have {} variants in common. The HST has {} variants in total.",
        count,
        hst_import.metadata.coords.len(),
    );

    let hst_type = hst_import.metadata.hst_type.clone();

    let match_hst = populate_imported_hst(&vcf, hst_import);

    // Write .hst to file
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "match_hst", "hst.gz");

    write_hst_file(match_hst, &vcf, hst_output, false, args, hst_type)?;

    Ok(())
}

pub fn populate_imported_hst(vcf: &PhasedMatrix, mut original_hst: Hst) -> Graph<Node, ()> {
    let node_idx = NodeIndex::new(0);

    let nodes_without_contradictory_ht = (0..vcf.nhaplotypes())
        .into_par_iter()
        .map(|sample_idx| {
            // Return a list of node_idxs where the haplotype has no contradictory genotypes for the sample
            let mut nodes_without_contradictory_ht = vec![];

            recursive_haplotype_comparison(
                vcf,
                node_idx,
                &original_hst,
                sample_idx,
                &mut nodes_without_contradictory_ht,
            );

            (sample_idx, nodes_without_contradictory_ht)
        })
        .collect::<Vec<(usize, Vec<NodeIndex>)>>();

    let mut map = HashMap::new();

    for (sample_idx, nodes) in nodes_without_contradictory_ht {
        for node_idx in nodes {
            map.entry(node_idx).or_insert(Vec::new()).push(sample_idx);
        }
    }

    for node in original_hst.hst.node_weights_mut() {
        node.indexes = vec![];
    }

    for (key, value) in map.into_iter() {
        let node = original_hst.hst.node_weight_mut(key).unwrap();
        node.indexes = value;
    }

    original_hst.hst
}

// Travel all nodes and check if haplotypes are concordant with sample haplotypes
pub fn recursive_haplotype_comparison(
    vcf: &PhasedMatrix,
    node_idx: NodeIndex,
    hst: &Hst,
    sample_idx: usize,
    nodes_without_contradictory_ht: &mut Vec<NodeIndex>,
) {
    hst.hst
        // Find outgoing nodes from node idx
        .neighbors_directed(node_idx, Direction::Outgoing)
        .for_each(|child_idx| {
            let child_node = hst.hst.node_weight(child_idx).unwrap();

            if no_contradictory_genotypes_in_ht(vcf, hst, child_node, sample_idx) {
                nodes_without_contradictory_ht.push(child_idx);
            }

            recursive_haplotype_comparison(
                vcf,
                child_idx,
                hst,
                sample_idx,
                nodes_without_contradictory_ht,
            );
        })
}

// Check if the node haplotype has contradictory genotypes to the sample haplotype
pub fn no_contradictory_genotypes_in_ht(
    vcf: &PhasedMatrix,
    hst: &Hst,
    node: &Node,
    sample_idx: usize,
) -> bool {
    let (start, stop) = (&node.start, &node.stop);

    let sample_haplotype = vcf.find_haplotype_for_sample(start..=stop, sample_idx);

    let ref_haplotype = hst.get_haplotype(node);

    let has_contradictory_markers = ref_haplotype.iter().any(|r| {
        sample_haplotype
            .iter()
            .any(|s| s.reference == r.reference && s.alt == r.alt && s.pos == r.pos && s.gt != r.gt)
    });

    !has_contradictory_markers
}
