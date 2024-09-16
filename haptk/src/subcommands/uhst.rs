use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use ndarray::s;
use petgraph::graph::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;

use crate::{
    args::{Selection, StandardArgs},
    io::{open_csv_writer, push_to_output, write_haplotype},
    structs::{HapVariant, PhasedMatrix},
    subcommands::bhst::{self, write_hst_file, HstType, Node},
};

#[doc(hidden)]
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone)]
pub enum LocDirection {
    Left,
    Right,
}

impl std::fmt::Display for LocDirection {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            LocDirection::Right => write!(f, "right"),
            LocDirection::Left => write!(f, "left"),
        }
    }
}

#[doc(hidden)]
pub fn run(args: StandardArgs, min_size: usize, publish: bool) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    let vcf = bhst::read_vcf_with_selections(&args)?;

    ensure!(
        vcf.nrows() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nrows()
    );

    let vec = vec![LocDirection::Left, LocDirection::Right];
    let first_and_mbah_nodes = vec
        .par_iter()
        .map(|direction| -> Result<(Node, Node)> {
            let uhst = construct_uhst(&vcf, direction, vcf.variant_idx(), min_size, false);

            // Find first, second last and last nodes on the majority branch
            // for downstream analyses
            let start_idx = NodeIndex::new(0);
            let nodes = bhst::find_majority_nodes(&uhst, start_idx);

            ensure!(
                nodes.len() > 2,
                "The majority branch has only less than 3 nodes."
            );

            let first_maj_node = nodes.get(1).unwrap().0.clone();
            let second_last_node = nodes.get(nodes.len() - 2).unwrap().0;
            let last_node = nodes.last().unwrap().0.clone();

            tracing::info!(
                "{}-side majority-based ancestral samples:{}",
                direction,
                second_last_node
                    .indexes
                    .iter()
                    .fold(String::new(), |cur, acc| {
                        let acc = vcf.get_sample_name(*acc);
                        format!("{cur} {acc}")
                    }),
            );

            // Write to .hst
            let mut hst_output = args.output.clone();
            push_to_output(
                &args,
                &mut hst_output,
                &format!("uhst_{direction}"),
                "hst.gz",
            );

            let hst_type = match direction {
                LocDirection::Left => HstType::UhstLeft,
                LocDirection::Right => HstType::UhstRight,
            };

            write_hst_file(uhst, &vcf, hst_output, publish, args.clone(), hst_type)?;

            Ok((first_maj_node, last_node))
        })
        .collect::<Result<Vec<(Node, Node)>>>()?;

    tracing::debug!("Finished constructing unilateral HSTs");

    tracing::debug!("Finished majority branches");

    // Shared core haplotype
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "uhst_shared_core_haplotype", "csv");
    let writer = open_csv_writer(sh_output)?;
    let core_haplotype = combine_node_haplotypes(
        &[
            first_and_mbah_nodes[0].0.clone(),
            first_and_mbah_nodes[1].0.clone(),
        ],
        &vcf,
    );
    write_haplotype(core_haplotype, writer)?;
    tracing::debug!("Finished writing core haplotype");

    // Majority based haplotype
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "uhst_mbah", "csv");
    let writer = open_csv_writer(sh_output)?;
    let mbah = combine_node_haplotypes(
        &[
            first_and_mbah_nodes[0].1.clone(),
            first_and_mbah_nodes[1].1.clone(),
        ],
        &vcf,
    );
    write_haplotype(mbah, writer)?;
    tracing::debug!("Finished writing majority based ancestral haplotype");

    Ok(())
}

#[doc(hidden)]
pub fn construct_uhst(
    vcf: &PhasedMatrix,
    direction: &LocDirection,
    idx: usize,
    min_size: usize,
    only_majority: bool,
) -> Graph<Node, u8> {
    let mut uhst = Graph::<_, _>::new();

    uhst.add_node(Node {
        start_idx: idx,
        stop_idx: idx,
        indexes: (0..vcf.matrix.nrows()).collect(),
        haplotype: vec![],
    });

    let mut min_size_blacklist = vec![];
    loop {
        // Filter indices
        let indices = uhst
            .node_indices()
            .filter(|node_idx| !min_size_blacklist.contains(node_idx))
            .filter(|node_idx| {
                let node = uhst.node_weight(*node_idx).unwrap();
                let count = uhst
                    .neighbors_directed(*node_idx, Direction::Outgoing)
                    .count();
                count == 0 && node.indexes.len() > 1
            })
            .collect::<Vec<NodeIndex>>();

        // Multithread horizontally all childless nodes
        let nodes = indices
            .into_par_iter()
            .filter_map(|node_idx| {
                find_contradictory_gt(vcf, &uhst, idx, node_idx, direction)
                    .unwrap()
                    .map(|(node1, node2)| (node_idx, node1, node2))
            })
            .collect::<Vec<(NodeIndex, Node, Node)>>();

        // Terminate if no new nodes will be added to the tree
        if nodes.is_empty() {
            break;
        }

        if only_majority {
            // This is only used for VCF struct only-longest selection
            for (idx, node1, node2) in nodes {
                match node1.indexes.len().cmp(&node2.indexes.len()) {
                    std::cmp::Ordering::Greater | std::cmp::Ordering::Equal => {
                        let node_idx = uhst.add_node(node1);
                        uhst.add_edge(idx, node_idx, 0);
                    }
                    std::cmp::Ordering::Less => {
                        let node_idx = uhst.add_node(node2);
                        uhst.add_edge(idx, node_idx, 0);
                    }
                }
            }
        } else {
            // Add nodes to the tree
            for (idx, node1, node2) in nodes {
                let uhst_node_count_before = uhst.node_count();

                if node1.indexes.len() >= min_size {
                    let node_idx1 = uhst.add_node(node1);
                    uhst.add_edge(idx, node_idx1, 0);
                }

                if node2.indexes.len() >= min_size {
                    let node_idx2 = uhst.add_node(node2);
                    uhst.add_edge(idx, node_idx2, 0);
                }

                if uhst.node_count() == uhst_node_count_before {
                    min_size_blacklist.push(idx);
                }
            }
        }
    }
    uhst
}

#[doc(hidden)]
fn find_contradictory_gt(
    vcf: &PhasedMatrix,
    uhst: &Graph<Node, u8>,
    idx: usize,
    node_idx: NodeIndex,
    direction: &LocDirection,
) -> Result<Option<(Node, Node)>> {
    let node = uhst.node_weight(node_idx).unwrap();

    let mut var_idx = node.start_idx;

    if node_idx == NodeIndex::new(0) {
        var_idx = match direction {
            LocDirection::Left => idx + 1,
            LocDirection::Right => idx.saturating_sub(1),
        };
    }

    let (mut z, mut o);
    loop {
        (z, o) = (vec![], vec![]);
        match direction {
            LocDirection::Left => {
                var_idx = var_idx.saturating_sub(1);
            }
            LocDirection::Right => {
                if var_idx < vcf.matrix.ncols() - 1 {
                    var_idx += 1;
                }
            }
        }

        for i in node.indexes.iter() {
            let gt = vcf.matrix.slice(s![*i, var_idx]);

            match gt.into_scalar() {
                0 => z.push(*i),
                1 => o.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        if (!z.is_empty() && !o.is_empty()) || var_idx == vcf.matrix.ncols() - 1 || var_idx == 0 {
            break;
        }
    }

    if !z.is_empty() && !o.is_empty() {
        let (start_idx, stop_idx) = match direction {
            LocDirection::Left => (var_idx, idx),
            LocDirection::Right => (idx, var_idx),
        };

        let n1 = Node {
            start_idx,
            stop_idx,
            haplotype: vcf.find_u8_haplotype_for_sample(start_idx..stop_idx + 1, z[0]),
            indexes: z,
        };

        let n2 = Node {
            start_idx,
            stop_idx,
            haplotype: vcf.find_u8_haplotype_for_sample(start_idx..stop_idx + 1, o[0]),
            indexes: o,
        };

        return Ok(Some((n1, n2)));
    }

    Ok(None)
}

pub fn combine_node_haplotypes(nodes: &[Node], vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let left_pos = nodes[0].start_idx;
    let left_sample = nodes[0].indexes[0];

    let left_range = left_pos + 1..vcf.variant_idx;

    let mut left_ht = vcf.find_haplotype_for_sample(vcf.get_contig(), left_range, left_sample);

    let right_pos = nodes[1].stop_idx;
    let right_sample = nodes[1].indexes[0];

    let right_range = vcf.variant_idx..right_pos;

    let right_ht = vcf.find_haplotype_for_sample(vcf.get_contig(), right_range, right_sample);

    left_ht.extend(right_ht);
    left_ht
}
