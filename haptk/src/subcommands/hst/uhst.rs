use std::cmp::Ordering;
use std::collections::HashMap;
use std::sync::mpsc::sync_channel;

use color_eyre::eyre::ensure;
use color_eyre::eyre::eyre;
use color_eyre::Result;
use petgraph::graph::NodeIndex;
use petgraph::Direction;
use petgraph::Graph;
use rayon::prelude::*;

use crate::subcommands::hst::find_majority_nodes;
use crate::subcommands::hst::pair_wise::calculate_pair_wise;
use crate::subcommands::hst::pair_wise::PairWiseMatrix;
use crate::subcommands::hst::read_vcf_with_selections;
use crate::subcommands::hst::CoordDataSlot;
use crate::subcommands::hst::HstType;
use crate::subcommands::hst::Node;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::Coord;
use crate::core::HapVariant;
use crate::core::PhasedMatrix;
use crate::error::Error;
use crate::io::open_csv_writer;
use crate::io::push_to_output;
use crate::io::write_haplotype;

use crate::subcommands::compare_to_haplotype::find_shared_haplotype_ranges;
use crate::subcommands::compare_to_haplotype::range_length_avg;
use crate::subcommands::compare_to_haplotype::range_length_median;
use crate::subcommands::compare_to_haplotype::transform_gt_matrix_to_match_matrix;
use crate::subcommands::compare_to_haplotype::write_ranges_to_csv;

use super::Hst;

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

    let mut vcf = read_vcf_with_selections(&args, Some(window))?;

    // Less haplotypes available than the minimum requested node_size
    ensure!(
        vcf.nhaplotypes() >= min_size,
        Error::MinNodes {
            n_haplotypes: vcf.nhaplotypes(),
            min_size
        }
    );

    let start = vcf.start_coord().clone();
    let vec = [LocDirection::Left, LocDirection::Right];
    let first_and_mbah_nodes = vec
        // .par_iter()
        .iter()
        .map(|direction| -> Result<(Node, Node, PairWiseMatrix)> {
            let hst_type = match direction {
                LocDirection::Left => HstType::HstLeft,
                LocDirection::Right => HstType::HstRight,
            };

            let mut uhst = Hst::new(&vcf, &args, hst_type);

            populate_uhst(&mut vcf, &mut uhst.hst, direction, &start, min_size, false)?;

            // Calulate pair-wise sharing matrix
            let pair_wise: PairWiseMatrix = calculate_pair_wise(&vcf, &uhst);

            // Find first, second last and last nodes on the majority branch
            // for downstream analyses
            let start_idx = NodeIndex::new(0);
            let nodes = find_majority_nodes(&uhst.hst, start_idx);

            // Check that theres at least 3 samples
            ensure!(nodes.len() > 2, Error::HstTooSmall);

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

            uhst.write_to_file(hst_output, publish)?;

            Ok((first_maj_node, last_node, pair_wise))
        })
        .collect::<Result<Vec<(Node, Node, PairWiseMatrix)>>>()?;

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

    // Pair-wise matrix
    let pairwise_left = first_and_mbah_nodes[0].2.clone();
    let pairwise_right = first_and_mbah_nodes[1].2.clone();

    let mut header = vec!["id".to_string()];

    header.extend(
        pairwise_left
            .iter()
            .map(|(id, _v)| vcf.get_sample_name(*id)),
    );

    let mut rows = vec![header];
    for ((id, left), (_, right)) in pairwise_left.iter().zip(pairwise_right.iter()) {
        let name = vcf.get_sample_name(*id);
        let mut sum = vec![name];

        for ((_, (start, _, v1)), (_, (_, stop, v2))) in left.iter().zip(right.iter()) {
            let value = if start == stop && *v1 && *v2 {
                1
            } else if start == stop {
                0
            } else {
                stop.saturating_sub(*start) + 1
            };

            sum.push(value.to_string());
        }

        rows.push(sum);
    }
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "pair-wise", "csv");
    let mut writer = open_csv_writer(sh_output)?;

    for row in rows {
        writer.write_record(row)?;
    }

    // Shared haplotype ranges
    //
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
pub fn populate_uhst(
    vcf: &mut PhasedMatrix,
    hst: &mut Graph<Node, ()>,
    direction: &LocDirection,
    start_coord: &Coord,
    min_size: usize,
    only_majority: bool,
) -> Result<()> {
    // let mut hst = initiate_hst(vcf, start_coord, None);

    let mut blacklist_nodes = vec![];

    loop {
        let blacklist = insert_nodes_to_uhst(
            vcf,
            hst,
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
    Ok(())
}

pub fn insert_nodes_to_uhst(
    vcf: &PhasedMatrix,
    hst: &mut Graph<Node, ()>,
    blacklist_nodes: &[NodeIndex],
    min_size: usize,
    direction: &LocDirection,
    only_majority: bool,
    start_coord: &Coord,
) -> std::result::Result<Vec<Vec<NodeIndex>>, Error> {
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

            let n_children = hst
                .neighbors_directed(node_idx, Direction::Outgoing)
                .count();

            match n_children == 0 && node.indexes.len() >= min_size.max(2) {
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
        .collect::<std::result::Result<Vec<Vec<NodeIndex>>, Error>>();

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
) -> std::result::Result<Option<(Node, Node)>, Error> {
    let next_contradictory_idx = match direction {
        LocDirection::Left => vcf.prev_contradictory(&node.start, &node.indexes),
        LocDirection::Right => vcf.next_contradictory(&node.start, &node.indexes),
    };

    let next_contradictory_idx = match (next_contradictory_idx, direction) {
        (Err(_), LocDirection::Left) => return Err(Error::HstLeftEnd),
        (Err(_), LocDirection::Right) => return Err(Error::HstRightEnd),
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
