#![allow(clippy::comparison_chain)]
use std::{collections::HashMap, path::PathBuf};

use color_eyre::{
    eyre::{eyre, WrapErr},
    Result,
};
use ndarray::parallel::prelude::*;
use petgraph::Graph;
use petgraph::{prelude::NodeIndex, Direction};

use crate::args::Selection;
use crate::core::parse_snp_coord;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::bhst::read_hst_file;
use crate::utils::push_to_output;
use crate::{args::StandardArgs, subcommands::bhst::write_hst_file};

use super::bhst::{Hst, Node};

#[doc(hidden)]
pub fn run(args: StandardArgs, hst_path: PathBuf) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    // File reads
    let hst_import = read_hst_file(hst_path)?;

    // IMG output
    let mut img_output = args.output.clone();
    push_to_output(&args, &mut img_output, "hst_comparison", "svg");

    svg::save(&img_output, &svg::Document::new())
        .wrap_err(eyre!("Failed writing to {:?}", img_output))?;

    // VCF read
    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;
    let (start, end) = (
        hst_import.coords.first().unwrap(),
        hst_import.coords.last().unwrap(),
    );
    tracing::info!("Reading coord range: {} - {}", start, end);

    let vcf = match args.selection {
        Selection::All | Selection::Haploid => {
            let vcf =
                read_vcf_to_matrix(&args, contig, variant_pos, Some((start.pos, end.pos)), None)?;
            vcf
        }
        Selection::OnlyAlts | Selection::OnlyRefs => {
            let mut vcf =
                read_vcf_to_matrix(&args, contig, variant_pos, Some((start.pos, end.pos)), None)?;

            // Select carriers before switching to the given haplotype as the reference
            vcf.select_carriers(variant_pos, &args.selection)?;
            vcf
        }
        Selection::OnlyLongest => {
            let mut vcf = read_vcf_to_matrix(&args, contig, variant_pos, None, None)?;
            vcf.select_only_longest();

            // Use start and end from the haplotype to celect columns from the matrix by range
            let start = vcf.get_first_idx_on_left_by_pos(start.pos);
            let mut stop = vcf.get_first_idx_on_right_by_pos(end.pos);
            if stop != vcf.ncoords() {
                stop += 1;
            }
            vcf.select_columns_by_range(
                // inclusive range not possible due to ndarray range generics being so difficult
                start..stop,
            );
            vcf
        }
        Selection::Unphased => unreachable!(),
    };

    if hst_import.coords.len() > vcf.ncoords() {
        tracing::warn!(
            "The HST has more variants than the given VCF {} vs {}",
            hst_import.coords.len(),
            vcf.ncoords()
        );
    } else if hst_import.coords.len() < vcf.ncoords() {
        tracing::warn!(
            "The HST has less variants than the given VCF {} vs {}",
            hst_import.coords.len(),
            vcf.ncoords()
        );
    }

    let count = hst_import
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
        hst_import.coords.len(),
    );

    let match_hst = create_match_hst(&vcf, hst_import);

    // Write .hst to file
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "match_hst", "hst.gz");
    write_hst_file(match_hst, &vcf, hst_output, false, args)?;

    Ok(())
}

pub fn init_match_hst(vcf: &PhasedMatrix, hst: &Hst) -> Graph<Node, u8> {
    let mut match_hst: Graph<Node, u8> = Graph::new();
    let node_idx = NodeIndex::new(0);
    let node = hst.hst.node_weight(node_idx).unwrap();

    let new_node = Node {
        start_idx: node.start_idx,
        stop_idx: node.stop_idx,
        indexes: vec![],
        haplotype: vec![],
    };
    let copy_idx = match_hst.add_node(new_node);

    recursive_copy(copy_idx, node_idx, hst, &mut match_hst, vcf);

    fn recursive_copy(
        copy_idx: NodeIndex,
        node_idx: NodeIndex,
        hst: &Hst,
        match_hst: &mut Graph<Node, u8>,
        vcf: &PhasedMatrix,
    ) {
        let children = hst.hst.neighbors_directed(node_idx, Direction::Outgoing);
        for child_idx in children {
            let child = hst.hst.node_weight(child_idx).unwrap();

            // let ref_samples = vec!["s".to_string(); child.indexes.len()];
            let new_node = Node {
                start_idx: child.start_idx,
                stop_idx: child.stop_idx,
                indexes: vec![],
                haplotype: vec![],
            };
            let copy_child_idx = match_hst.add_node(new_node);

            match_hst.add_edge(copy_idx, copy_child_idx, 0);
            recursive_copy(copy_child_idx, child_idx, hst, match_hst, vcf)
        }
    }
    match_hst
}

pub fn create_match_hst(vcf: &PhasedMatrix, hst: Hst) -> Graph<Node, u8> {
    let node_idx = NodeIndex::new(0);

    let node_idx_vecs = (0..vcf.nrows())
        .into_par_iter()
        .map(|sample_idx| {
            // Return a list of node_idxs where the haplotype has no contradictory genotypes for the sample
            let mut node_idxs = vec![];
            recursive_haplotype_comparison(vcf, node_idx, &hst, sample_idx, &mut node_idxs);
            node_idxs
        })
        .collect::<Vec<Vec<NodeIndex>>>();

    let mut match_hst = init_match_hst(vcf, &hst);

    let mut map = HashMap::new();
    for (sample_idx, vec) in node_idx_vecs.into_iter().enumerate() {
        for node_idx in vec {
            map.entry(node_idx).or_insert(Vec::new()).push(sample_idx);
        }
    }

    for (key, value) in map.into_iter() {
        let node = match_hst.node_weight_mut(key).unwrap();
        node.indexes.extend(value);
    }
    // println!("{match_hst:?}");

    match_hst
}

// Travel all nodes and check if haplotypes are concordant with sample haplotypes
pub fn recursive_haplotype_comparison(
    vcf: &PhasedMatrix,
    node_idx: NodeIndex,
    hst: &Hst,
    sample_idx: usize,
    found_node_idxs: &mut Vec<NodeIndex>,
) {
    let children = hst.hst.neighbors_directed(node_idx, Direction::Outgoing);
    for child_idx in children {
        if !haplotype_has_contradictory_genotypes(vcf, hst, child_idx, sample_idx) {
            found_node_idxs.push(child_idx);
        }

        recursive_haplotype_comparison(vcf, child_idx, hst, sample_idx, found_node_idxs);
    }
}

// Check if the node haplotype has any contradictory genotypes to the sample haplotype
pub fn haplotype_has_contradictory_genotypes(
    vcf: &PhasedMatrix,
    hst: &Hst,
    node_idx: NodeIndex,
    sample_idx: usize,
) -> bool {
    let child = hst.hst.node_weight(node_idx).unwrap();

    // How to handle direction and missing variants here?
    let start_idx = child.start_idx;
    let stop_idx = child.stop_idx;
    if start_idx == stop_idx {
        return false;
    }

    let start_pos = hst.get_pos(child.start_idx);
    let stop_pos = hst.get_pos(child.stop_idx);
    let start = vcf.get_first_idx_on_left_by_pos(start_pos);
    let mut stop = vcf.get_first_idx_on_right_by_pos(stop_pos);
    if stop != vcf.ncoords() {
        stop += 1;
    }

    let sample_haplotype = vcf.find_haplotype_for_sample(vcf.get_contig(), start..stop, sample_idx);

    let ref_haplotype = hst.get_haplotype(node_idx);
    // for ht in &ref_haplotype {
    // println!("{ht}");
    // }

    let has_contradictory_markers = ref_haplotype.iter().any(|r| {
        sample_haplotype
            .iter()
            .any(|s| s.reference == r.reference && s.alt == r.alt && s.pos == r.pos && s.gt != r.gt)
    });

    // let contradictory_marker_pos = ref_haplotype.iter().position(|r| {
    // sample_haplotype
    // .iter()
    // .any(|s| s.reference == r.reference && s.alt == r.alt && s.pos == r.pos && s.gt != r.gt)
    // });
    // let samples_gt: Vec<usize> = ref_haplotype
    // .iter()
    // .filter_map(|r| {
    // sample_haplotype.iter().position(|s| {
    // s.reference == r.reference && s.alt == r.alt && s.pos == r.pos && s.gt != r.gt
    // })
    // })
    // .collect();
    // if let Some(contradictory_marker_pos) = contradictory_marker_pos {
    // println!("here {}", ref_haplotype[contradictory_marker_pos]);
    // println!("heree {}", sample_haplotype[samples_gt[0]]);
    // }
    // println!(
    // "comparing {} vs {}",
    // sample_haplotype.len(),
    // ref_haplotype.len()
    // );
    has_contradictory_markers
}
