use std::{cmp::Ordering, collections::HashMap, sync::mpsc::sync_channel};

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use petgraph::Graph;
use petgraph::{graph::NodeIndex, Direction};
use rayon::prelude::*;

use crate::{
    args::{Selection, StandardArgs},
    error::HaptkError,
    io::{open_csv_writer, push_to_output, write_haplotype},
    libs::structs::{CoordDataSlot, PhasedMatrix},
    structs::{Coord, HapVariant},
    subcommands::{
        bhst::{self, HstType, Node},
        compare_to_haplotype::{
            find_shared_haplotype_ranges, range_length_avg, range_length_median,
            transform_gt_matrix_to_match_matrix, write_ranges_to_csv,
        },
    },
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
pub fn run(args: StandardArgs, min_size: usize, publish: bool, window: u64) -> Result<()> {
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    let mut vcf = bhst::read_vcf_with_selections(&args, Some(window))?;

    ensure!(
        vcf.nhaplotypes() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nhaplotypes()
    );

    let start = vcf.start_coord().clone();
    let vec = [LocDirection::Left, LocDirection::Right];
    let first_and_mbah_nodes = vec
        // .par_iter()
        .iter()
        .map(|direction| -> Result<(Node, Node)> {
            let uhst = construct_uhst(&mut vcf, direction, &start, min_size, false)?;

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
            push_to_output(&args, &mut hst_output, &format!("{direction}"), "hst.gz");

            let hst_type = match direction {
                LocDirection::Left => HstType::HstLeft,
                LocDirection::Right => HstType::HstRight,
            };

            bhst::write_hst_file(uhst, &vcf, hst_output, publish, args.clone(), hst_type)?;

            Ok((first_maj_node, last_node))
        })
        .collect::<Result<Vec<(Node, Node)>>>()?;

    tracing::debug!("Finished constructing unilateral HSTs");

    tracing::debug!("Finished majority branches");

    // Shared core haplotype
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "core_haplotype", "csv");
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
    push_to_output(&args, &mut sh_output, "ancestral_haplotype", "csv");
    let writer = open_csv_writer(sh_output)?;
    let ancestral_haplotype = combine_node_haplotypes(
        &[
            first_and_mbah_nodes[0].1.clone(),
            first_and_mbah_nodes[1].1.clone(),
        ],
        &vcf,
    );

    write_haplotype(ancestral_haplotype.clone(), writer)?;
    tracing::debug!("Finished writing majority based ancestral haplotype");

    // Vec into Map for speed up
    let ht: HashMap<Coord, HapVariant> = ancestral_haplotype
        .into_iter()
        .map(|v| (v.clone().into(), v))
        .collect();

    let variant_pos = vcf.variant_idx_pos();
    let vcf = transform_gt_matrix_to_match_matrix(vcf, &ht, variant_pos)?;

    let shared_ranges = find_shared_haplotype_ranges(&vcf);

    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "ancestral_segments", "csv");
    let mut writer = open_csv_writer(sh_output)?;

    write_ranges_to_csv(&vcf, &shared_ranges, None, None, &mut writer)?;

    tracing::info!(
        "Sample haplotypes: {}, average length: {}, median length: {}",
        shared_ranges.len(),
        range_length_avg(&shared_ranges),
        range_length_median(&shared_ranges)
    );
    tracing::debug!("Finished writing ancestral segment lengths");

    Ok(())
}

#[doc(hidden)]
pub fn construct_uhst(
    vcf: &mut PhasedMatrix,
    direction: &LocDirection,
    start_coord: &Coord,
    min_size: usize,
    only_majority: bool,
) -> Result<Graph<Node, ()>> {
    let mut hst = bhst::initiate_hst(vcf, start_coord, None);

    let mut blacklist_nodes = vec![];

    loop {
        let blacklist = insert_nodes_to_uhst(
            vcf,
            &mut hst,
            &blacklist_nodes,
            min_size,
            direction,
            only_majority,
            start_coord,
        );

        match blacklist {
            Err(e) => {
                vcf.read_more(e)?;
                continue;
            }
            Ok(blacklist) => {
                if blacklist.is_empty() {
                    break;
                }
                if only_majority && blacklist[0].is_empty() {
                    break;
                }
                blacklist_nodes.extend(blacklist.iter().flatten());
            }
        }
    }
    Ok(hst)
}

pub fn insert_nodes_to_uhst(
    vcf: &PhasedMatrix,
    hst: &mut Graph<Node, ()>,
    blacklist_nodes: &[NodeIndex],
    min_size: usize,
    direction: &LocDirection,
    only_majority: bool,
    start_coord: &Coord,
) -> std::result::Result<Vec<Vec<NodeIndex>>, HaptkError> {
    let (tx, rx) = sync_channel(2024);

    // Loop all nodes of the HST
    let blacklist = hst
        .node_indices()
        .par_bridge()
        // Filter out nodes that were blacklisted (genotyping data ran out or node haplotypes min_size was reached )
        .filter(|n| !blacklist_nodes.contains(n))
        // Filter in leaf nodes with 2 or more haplotypes
        .filter_map(|node_idx| {
            let node = hst.node_weight(node_idx).unwrap();

            let count = hst
                .neighbors_directed(node_idx, Direction::Outgoing)
                .count();

            match count == 0 && node.indexes.len() >= min_size.max(2) {
                true => Some((node_idx, node.clone())),
                false => None,
            }
        })
        // Search for contradictory genotypes by extending the node haplotype
        .map(|(parent_idx, node)| {
            match find_contradictory_gt_uhst(vcf, start_coord, &node, direction) {
                Err(e) => Err(e),
                Ok(Some((node1, node2))) => {
                    let black_list_nodes: Vec<NodeIndex> = if only_majority {
                        match node1.indexes.len().cmp(&node2.indexes.len()) {
                            Ordering::Greater | Ordering::Equal => {
                                tx.send((parent_idx, node1)).unwrap()
                            }
                            Ordering::Less => tx.send((parent_idx, node2)).unwrap(),
                        }
                        vec![parent_idx]
                    } else {
                        [node1, node2]
                            .into_iter()
                            .filter_map(|insert_node| {
                                if insert_node.indexes.len() >= min_size {
                                    tx.send((parent_idx, insert_node)).unwrap();
                                    None
                                } else {
                                    Some(parent_idx)
                                }
                            })
                            .collect()
                    };

                    Ok(black_list_nodes)
                }
                _ => Ok(vec![parent_idx]),
            }
        })
        .collect::<std::result::Result<Vec<Vec<NodeIndex>>, HaptkError>>();

    drop(tx);

    // Receive nodes for the iterator and add them to the tree
    // We use a channel because if genotyping data needes to be dynamically extended, we dont want to lose the already computed nodes
    while let Ok((parent_idx, insert_node)) = rx.recv() {
        let idx = hst.add_node(insert_node.clone());
        hst.add_edge(parent_idx, idx, ());
    }

    blacklist
}

#[doc(hidden)]
pub fn find_contradictory_gt_uhst(
    vcf: &PhasedMatrix,
    start_coord: &Coord,
    node: &Node,
    direction: &LocDirection,
) -> std::result::Result<Option<(Node, Node)>, HaptkError> {
    let next_contradictory_idx = match direction {
        LocDirection::Left => vcf.prev_contradictory(&node.start, &node.indexes),
        LocDirection::Right => vcf.next_contradictory(&node.start, &node.indexes),
    };

    let next_contradictory_idx = match (next_contradictory_idx, direction) {
        (Err(_), LocDirection::Left) => return Err(HaptkError::HstLeftEndError),
        (Err(_), LocDirection::Right) => return Err(HaptkError::HstRightEndError),
        (Ok(idx), _) => idx,
    };

    if let Some((next_coord, genotypes)) = next_contradictory_idx {
        let (mut ones, mut zeroes) = (vec![], vec![]);

        for i in &node.indexes {
            match genotypes[*i] == 1 {
                true => ones.push(*i),
                false => zeroes.push(*i),
            }
        }

        let (start, stop) = match direction {
            LocDirection::Left => (next_coord, start_coord),
            LocDirection::Right => (start_coord, next_coord),
        };

        let n1 = Node {
            haplotype: vcf.find_u8_haplotype_for_sample(start..=stop, zeroes[0]),
            // start: Cow::Borrowed(start),
            // stop: Cow::Borrowed(stop),
            start: start.clone(),
            stop: stop.clone(),
            indexes: zeroes,
        };

        let n2 = Node {
            haplotype: vcf.find_u8_haplotype_for_sample(start..=stop, ones[0]),
            start: start.clone(),
            stop: stop.clone(),
            indexes: ones,
        };

        return Ok(Some((n1, n2)));
    }

    if !vcf.is_genome_wide() {
        tracing::warn!(
            "Genotyping data ran out on {} side for samples {:?} starting from position {}",
            direction,
            vcf.get_sample_names(&node.indexes),
            node.start.pos
        );
    }

    // If data ran out and no contradictory genotypes were found, return none
    Ok(None)
}

pub fn combine_node_haplotypes(nodes: &[Node], vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let left = &nodes[0].start;
    let left_plus = vcf.coords().range(&nodes[0].start..).nth(1).unwrap();

    let left = match left_plus > vcf.start_coord() {
        true => left,
        false => left_plus,
    };

    let left_sample = nodes[0].indexes[0];

    let left_range = left..vcf.start_coord();

    let mut left_ht = vcf.find_haplotype_for_sample(left_range, left_sample);

    let right = &nodes[1].stop;
    let right_sample = nodes[1].indexes[0];

    let right_range = vcf.start_coord()..right;

    let right_ht = vcf.find_haplotype_for_sample(right_range, right_sample);

    left_ht.extend(right_ht);
    left_ht
}
