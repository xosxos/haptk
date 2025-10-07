use std::sync::mpsc::sync_channel;

use color_eyre::eyre::ensure;
use color_eyre::Result;
use petgraph::prelude::NodeIndex;
use petgraph::Direction;
use petgraph::Graph;
use rayon::prelude::*;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::Coord;
use crate::core::HapVariant;
use crate::core::PhasedMatrix;
use crate::error::Error;
use crate::io::open_csv_writer;
use crate::io::push_to_output;
use crate::io::write_haplotype;

use crate::subcommands::hst::find_majority_nodes;
use crate::subcommands::hst::find_shared_haplotype;
use crate::subcommands::hst::read_vcf_with_selections;
use crate::subcommands::hst::CoordDataSlot;
use crate::subcommands::hst::HstType;
use crate::subcommands::hst::Node;

use super::Hst;

#[doc(hidden)]
pub fn run(args: StandardArgs, min_size: usize, publish: bool, window: Option<u64>) -> Result<()> {
    ensure!(
        args.selection != Selection::Unphased,
        "Running with unphased data is not supported."
    );

    let mut vcf = read_vcf_with_selections(&args, window)?;

    ensure!(
        vcf.nhaplotypes() >= min_size,
        "VCF has less haplotypes than the required minimum node size ({} < {min_size})",
        vcf.nhaplotypes()
    );

    // Construct the bidirectional HST
    let mut bhst = construct_bhst(&mut vcf, &args, min_size)?;
    tracing::info!("Finished HST construction.");

    // Find the majority based haplotype
    let mbah = find_ancestral_haplotype(&bhst, &vcf)?;

    // Write it to file
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_ancestral_haplotype", "csv");
    write_haplotype(mbah, open_csv_writer(sh_output)?)?;

    // Find the shared core haplotype
    let ht = find_shared_haplotype(&bhst.hst, &vcf);

    // Write it to file
    let mut sh_output = args.output.clone();
    push_to_output(&args, &mut sh_output, "bhst_core_haplotype", "csv");
    write_haplotype(ht, open_csv_writer(sh_output)?)?;

    // Write the HST to file
    let mut hst_output = args.output.clone();
    push_to_output(&args, &mut hst_output, "bhst", "hst.gz");
    bhst.write_to_file(hst_output, publish)?;

    Ok(())
}

#[doc(hidden)]
pub fn construct_bhst(vcf: &mut PhasedMatrix, args: &StandardArgs, min_size: usize) -> Result<Hst> {
    // Initate the HST by inserting the root node, and two child nodes if there is a contradictory genotype right at the starting variant

    let mut hst = Hst::new(vcf, args, HstType::Bhst);
    // let mut hst = initiate_hst(vcf, start_coord, None);

    let mut blacklist_nodes: Vec<NodeIndex> = vec![];

    loop {
        // Find contradictory genotypes for samples in leaf nodes
        // Insert new nodes into a the HST using a parallel iterator and  channels
        // Return a blacklist of nodes to not exclude from iteration each round
        let blacklist = insert_nodes_to_bhst(vcf, &mut hst.hst, &blacklist_nodes, min_size);

        match blacklist {
            Err(e) => {
                // If data runs, and the file has not ended, read more from disk
                vcf.read_more(e)?;
                continue;
            }
            Ok(blacklist) => {
                if blacklist.is_empty() {
                    break;
                }
                blacklist_nodes.extend(blacklist.iter().flatten());
            }
        }
    }

    Ok(hst)
}

pub fn insert_nodes_to_bhst(
    vcf: &PhasedMatrix,
    hst: &mut Graph<Node, ()>,
    blacklist_nodes: &[NodeIndex],
    min_size: usize,
) -> std::result::Result<Vec<Vec<NodeIndex>>, Error> {
    let (tx, rx) = sync_channel(2024);

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
        .map(
            |(parent_idx, node)| match find_contradictory_gt_bhst(vcf, &node) {
                Err(e) => Err(e),
                Ok(Some(nodes)) => {
                    let black_list_nodes: Vec<NodeIndex> = nodes
                        .into_iter()
                        .filter_map(|insert_node| {
                            if insert_node.indexes.len() >= min_size {
                                tx.send((parent_idx, insert_node)).unwrap();
                                None
                            } else {
                                Some(parent_idx)
                            }
                        })
                        .collect();
                    Ok(black_list_nodes)
                }
                _ => Ok(vec![parent_idx]),
            },
        )
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
pub fn find_contradictory_gt_bhst(
    vcf: &PhasedMatrix,
    node: &Node,
) -> std::result::Result<Option<Vec<Node>>, Error> {
    let (left_coord, right_coord) = (&node.start, &node.stop);

    // NOTE: Remains to be seen what the overhead here would be
    // let (prev, next) = thread::scope(|s| {
    //     let prev_handle = s.spawn(|| vcf.prev_contradictory(left_coord, &node.indexes));
    //     let next = vcf.next_contradictory(right_coord, &node.indexes);
    //     let prev = prev_handle.join().unwrap();
    //     (prev, next)
    // });

    let prev = vcf.prev_contradictory(left_coord, &node.indexes);
    let next = vcf.next_contradictory(right_coord, &node.indexes);

    let (prev, next) = match (prev, next) {
        (Err(_), Err(_)) => return Err(Error::HstBothEnd),
        (_, Err(_)) => return Err(Error::HstRightEnd),
        (Err(_), _) => return Err(Error::HstLeftEnd),
        (Ok(prev), Ok(next)) => (prev, next),
    };

    // NOTE: Benchmarking required
    // Allocate Vecs with capacity so no reallocation is required
    let (mut zo, mut zz, mut oz, mut oo) = (
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
        Vec::with_capacity(node.indexes.len()),
    );

    // Fill the genotype buckets and then create HST nodes from the buckets
    match (prev, next) {
        (Some((left, left_vec)), Some((right, right_vec))) => {
            for i in node.indexes.iter() {
                let left_bit = left_vec[*i] == 1;
                let right_bit = right_vec[*i] == 1;

                match (left_bit, right_bit) {
                    (false, false) => zz.push(*i),
                    (false, true) => zo.push(*i),
                    (true, false) => oz.push(*i),
                    (true, true) => oo.push(*i),
                }
            }
            // Create a new node for each bucket if it is not empty
            let nodes = create_nodes_from_buckets(vcf, left, right, oo, oz, zo, zz);
            Ok(Some(nodes))
        }
        // Handle if data ran out on the right side
        (Some((left, left_vec)), None) => {
            if !vcf.is_genome_wide() {
                tracing::warn!(
                    "Genotyping data ran out on the right with samples {:?}",
                    vcf.get_sample_names(&node.indexes)
                );
            }

            for i in node.indexes.iter() {
                match left_vec[*i] == 1 {
                    true => oz.push(*i),
                    false => zz.push(*i),
                }
            }
            let nodes =
                create_nodes_from_buckets(vcf, left, vcf.coords().last().unwrap(), oo, oz, zo, zz);
            Ok(Some(nodes))
        }
        // Handle if data ran out on the left side
        (None, Some((right, right_vec))) => {
            if !vcf.is_genome_wide() {
                tracing::warn!(
                    "Genotyping data ran out on the left with samples {:?}",
                    vcf.get_sample_names(&node.indexes)
                );
            }

            for i in node.indexes.iter() {
                match right_vec[*i] == 1 {
                    true => zo.push(*i),
                    false => zz.push(*i),
                }
            }
            let nodes = create_nodes_from_buckets(
                vcf,
                vcf.coords().first().unwrap(),
                right,
                oo,
                oz,
                zo,
                zz,
            );
            Ok(Some(nodes))
        }
        // Handle if data ran out on both sides
        (None, None) => Ok(None),
    }
}

fn create_nodes_from_buckets(
    vcf: &PhasedMatrix,
    left: &Coord,
    right: &Coord,
    oo: Vec<usize>,
    oz: Vec<usize>,
    zo: Vec<usize>,
    zz: Vec<usize>,
) -> Vec<Node> {
    let mut nodes = Vec::with_capacity(4);
    if !zz.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zz[0]),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, zo[0]),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oz[0]),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            start: left.clone(),
            stop: right.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(left..=right, oo[0]),
            indexes: oo,
        };
        nodes.push(node);
    }
    nodes
}

pub fn find_ancestral_haplotype(hst: &Hst, vcf: &PhasedMatrix) -> Result<Vec<HapVariant>> {
    let start_idx = NodeIndex::new(0);
    let nodes = find_majority_nodes(&hst.hst, start_idx);

    // Check that theres at least 3 samples
    ensure!(nodes.len() > 2, Error::HstTooSmall);

    let mut iter = nodes.iter().rev();
    let (last_node, _last_node_idx) = iter.next().unwrap();
    let (second_last_node, _) = iter.next().unwrap();

    tracing::info!(
        "Majority-based ancestral samples:{}",
        second_last_node
            .indexes
            .iter()
            .fold(String::new(), |cur, acc| {
                let acc = vcf.get_sample_name(*acc);
                format!("{cur} {acc}")
            }),
    );

    // Start and stop idxs are the first deviating genotype and it cannot be a part of the haplotype
    let start = vcf.coords().range(&last_node.start..).nth(1).unwrap();

    let start = match start > &last_node.stop {
        true => &last_node.start,
        false => start,
    };

    Ok(vcf.find_haplotype_for_sample(start..&last_node.stop, last_node.indexes[0]))
}
