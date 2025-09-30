use std::{collections::HashMap, path::PathBuf};

use color_eyre::{eyre::eyre, Result};
use ndarray::parallel::prelude::*;
use petgraph::Graph;
use petgraph::{prelude::NodeIndex, Direction};

use crate::io::get_htslib_contig_len;
use crate::read_vcf::read_vcf_to_matrix_by_indexes;
use crate::{
    args::{Selection, StandardArgs},
    io::push_to_output,
    read_vcf::read_vcf_to_matrix,
    structs::PhasedMatrix,
    subcommands::bhst::{read_hst_file, write_hst_file, Hst, Node},
    utils::parse_snp_coord,
};

#[doc(hidden)]
pub fn run(args: StandardArgs, hst_path: PathBuf, only_longest_leafs: bool) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    if args.selection == Selection::List {
        return Err(eyre!(
            "Running with list selection is not supported for now."
        ));
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
                false,
            )?
        }
        Selection::OnlyLongest => {
            let (only_longest_lookups, vcf) = if get_htslib_contig_len(&args.file, contig).is_ok() {
                let mut vcf = read_vcf_to_matrix(
                    &args,
                    contig,
                    variant_pos,
                    None,
                    None,
                    Some(5_000_000),
                    false,
                )?;
                let lookups = vcf.get_only_longest_lookups()?;
                (lookups, vcf)
            } else {
                let vcf = read_vcf_to_matrix(&args, contig, variant_pos, None, None, None, false)?;
                let lookups = vcf.get_only_longest_lookups_no_shard()?;
                (lookups, vcf)
            };

            read_vcf_to_matrix_by_indexes(
                &args.file,
                variant_pos,
                contig,
                Some((Some(start.pos), Some(end.pos))),
                vcf.samples().clone(),
                vcf.metadata.indexes,
                only_longest_lookups,
                None,
                args.no_alt,
                &Selection::Haploid,
                false,
                args.only_snv,
            )?
        }
        Selection::Unphased | Selection::List => unreachable!(),
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
        .filter_map(|r| vcf.coords().get(r))
        .count();

    tracing::info!(
        "The HST and the given VCF have {} variants in common. The HST has {} variants in total.",
        count,
        hst_import.metadata.coords.len(),
    );

    let hst_type = hst_import.metadata.hst_type.clone();

    let match_hst = populate_imported_hst(&vcf, hst_import, only_longest_leafs);

    // Write .hst to file
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "match_hst", "hst.gz");

    write_hst_file(match_hst, &vcf, hst_output, false, args, hst_type)?;

    Ok(())
}

pub fn populate_imported_hst(
    vcf: &PhasedMatrix,
    mut original_hst: Hst,
    only_longest: bool,
) -> Graph<Node, ()> {
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

    for (sample_idx, nodes) in nodes_without_contradictory_ht.into_iter() {
        let nodes = match only_longest {
            true => remove_non_max_leaf_nodes(&original_hst, nodes),
            false => nodes,
        };

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

pub fn remove_non_max_leaf_nodes(hst: &Hst, nodes: Vec<NodeIndex>) -> Vec<NodeIndex> {
    // Helper function
    let is_leaf_node = |node_idx| {
        hst.hst
            .neighbors_directed(node_idx, Direction::Outgoing)
            .count()
            == 0
    };

    let mut max_bp_len = 0;

    // Find max bp_len for leaf nodes
    for node_idx in &nodes {
        let node = hst.hst.node_weight(*node_idx).unwrap();
        let bp_len = node.stop.pos.saturating_sub(node.start.pos);

        if is_leaf_node(*node_idx) && bp_len > max_bp_len {
            max_bp_len = bp_len;
        }
    }

    // Many leaf nodes can have max bp_len, so re-iterate to get max_nodes
    let max_nodes: Vec<NodeIndex> = nodes
        .iter()
        .filter(|&node_idx| {
            let node = hst.hst.node_weight(*node_idx).unwrap();
            let bp_len = node.stop.pos.saturating_sub(node.start.pos);

            is_leaf_node(*node_idx) && bp_len == max_bp_len
        })
        .copied()
        .collect();

    // Return nodes
    nodes
        .into_iter()
        // Keep, if node == max_node, or if the node is not a leaf node
        .filter(|idx| max_nodes.contains(idx) || !is_leaf_node(*idx))
        .collect()
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

    let has_contradictory_markers = hst
        .metadata
        .coords
        // This gets the hst haplotype
        .range(start..=stop)
        .enumerate()
        .map(|(hap_index, coord)| (coord, node.haplotype[hap_index]))
        // And compares all the markers in it
        .any(|(r, other_gt)| {
            sample_haplotype.iter().any(|s| {
                s.reference == r.reference && s.alt == r.alt && s.pos == r.pos && s.gt != other_gt
            })
        });

    !has_contradictory_markers
}
