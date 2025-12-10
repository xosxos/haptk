use std::path::PathBuf;

use petgraph::graph::NodeIndex;

use crate::error::Error;
use crate::subcommands::hst::{Hst, Node};

#[cfg_attr(feature = "clap", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct Args {
    /// A file path to a HST file
    hst: PathBuf,

    /// Algorithm for finding subtrees
    #[cfg_attr(feature = "clap", arg(short = 'a', long, value_enum, default_value_t = Algorithm::MajorityBranch))]
    algorithm: Algorithm,

    /// Max amount of indexes in subtree
    #[cfg_attr(feature = "clap", arg(long, default_value_t = usize::MAX))]
    max_node_size: usize,

    /// Min amount of indexes in subtree
    #[cfg_attr(feature = "clap", arg(long, default_value_t = 10))]
    min_node_size: usize,

    /// Min haplotype length
    #[cfg_attr(feature = "clap", arg(long, default_value_t = 0))]
    min_ht_len: u64,

    /// Max haplotype length
    #[cfg_attr(feature = "clap", arg(long, default_value_t = u64::MAX))]
    max_ht_len: u64,

    /// Set the root node
    #[cfg_attr(feature = "clap", arg(long, default_value_t = 0))]
    root_node: usize,
}

#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Default)]
pub enum Algorithm {
    #[default]
    MajorityBranch,
    LeafNodes,
    Size,
    All,
}

#[doc(hidden)]
pub fn run(args: Args) -> Result<(), Error> {
    let mut hst = Hst::from_file(&args.hst)?;

    if args.root_node != 0 {
        hst.subtree(NodeIndex::new(args.root_node));
    }

    match args.algorithm {
        Algorithm::MajorityBranch => majority_branch(&hst, &args),
        Algorithm::LeafNodes => leaf_nodes(&hst),
        Algorithm::Size => size(&hst, &args),
        Algorithm::All => all(&hst, &args),
    }

    Ok(())
}

fn majority_branch(hst: &Hst, args: &Args) {
    println!("node_idx,bp_size,maj_branch_size,min_branch_size");

    let mut parent_node_indexes = hst.nhaplotypes();

    for idx in hst.majority_branch() {
        let node_data = hst.node_weight(idx).unwrap();
        let indexes_len = node_data.indexes.len();

        let branch_size = parent_node_indexes.saturating_sub(indexes_len);

        if branch_size >= args.min_node_size && branch_size <= args.max_node_size {
            let siblings = hst.get_siblings(idx);

            if !siblings.is_empty() {
                let sibling = siblings[0];
                let sibling = hst.node_weight(sibling).unwrap();

                let ht_length = node_data.stop.pos.saturating_sub(node_data.start.pos);
                let names = sibling.sample_names(
                    &hst.metadata.samples,
                    *hst.metadata.ploidy,
                    &hst.metadata.selected_haplotypes,
                );

                println!(
                    "{:?},{:?},{:?},{:?},{}",
                    idx.index(),
                    ht_length,
                    node_data.indexes.len(),
                    sibling.indexes.len(),
                    names,
                );
            }
        }

        parent_node_indexes = indexes_len;
    }
}

fn leaf_nodes(hst: &Hst) {
    println!("node_idx,bp_size,n_haplotypes,samples");

    for (idx, node) in hst.nodes() {
        if node.indexes.len() == 2 {
            let ht_length = node.stop.pos.saturating_sub(node.start.pos);

            let names = node.sample_names(
                &hst.metadata.samples,
                *hst.metadata.ploidy,
                &hst.metadata.selected_haplotypes,
            );

            println!(
                "{:?},{:?},{:?},{}",
                idx.index(),
                ht_length,
                node.indexes.len(),
                names
            );
        }
    }
}

fn all(hst: &Hst, args: &Args) {
    println!("node_idx,bp_size,n_haplotypes,samples");

    let nodes: Vec<_> = hst
        .nodes()
        // Filter in nodes
        .filter(|(_idx, node)| {
            let ht_len = node.stop.pos.saturating_sub(node.start.pos);

            node.indexes.len() >= args.min_node_size
                && node.indexes.len() <= args.max_node_size
                && ht_len >= args.min_ht_len
                && ht_len <= args.max_ht_len
        })
        .collect();

    for (idx, node) in nodes {
        let ht_length = node.stop.pos.saturating_sub(node.start.pos);

        let names = node.sample_names(
            &hst.metadata.samples,
            *hst.metadata.ploidy,
            &hst.metadata.selected_haplotypes,
        );

        println!(
            "{:?},{:?},{:?},{}",
            idx.index(),
            ht_length,
            node.indexes.len(),
            names
        );
    }
}

// First draft of this
//
// We use Vec for now as HSTs are small
//
// but remember:
//
// https://stackoverflow.com/questions/64226562/check-if-vec-contains-all-elements-from-another-vec
//
// Speed in nanoseconds for Vec vs HashSet contains
// size      Vec          HashSet
// 10        14           386
// 100       1754         3187
// 1000      112_306      31233
// 10000     2_821_867    254_801
//

fn size(hst: &Hst, args: &Args) {
    println!("node_idx,bp_size,length,pair_wise_avg,samples");

    let pair_wise = hst.calculate_pair_wise();

    let mut nodes: Vec<_> = hst
        .nodes()
        // Filter in nodes
        .filter(|(_idx, node)| {
            let ht_len = node.stop.pos.saturating_sub(node.start.pos);

            node.indexes.len() >= args.min_node_size
                && node.indexes.len() <= args.max_node_size
                && ht_len >= args.min_ht_len
                && ht_len <= args.max_ht_len
        })
        .collect();

    // Sort nodes by size
    nodes.sort_by_key(|(_idx, data)| data.indexes.len());

    // Deduplicated nodes
    let mut dedup: Vec<(NodeIndex, &Node)> = vec![];

    // Iterate sorted nodes in reverse to get size by ascending order
    'outer: for item in nodes.into_iter().rev() {
        for bigs in &dedup {
            // If any of the samples are shared between nodes
            // there must be a child-parent relationship
            //
            // Because the list is sorted so that the largest come first
            // We can just filter out all nodes like so
            if item
                .1
                .indexes
                .iter()
                // Could be HashSet
                .any(|idx| bigs.1.indexes.contains(idx))
            {
                // Jump to next node
                continue 'outer;
            }
        }

        // If there was no sample overlap, push to deduplicated nodes
        dedup.push(item);
    }

    for (idx, node) in dedup {
        let ht_length = node.stop.pos.saturating_sub(node.start.pos);

        let names = node.sample_names(
            &hst.metadata.samples,
            *hst.metadata.ploidy,
            &hst.metadata.selected_haplotypes,
        );

        let mean = pair_wise.mean(&node.indexes);

        println!(
            "{},{},{},{:.0},{}",
            idx.index(),
            node.indexes.len(),
            ht_length,
            mean,
            names
        );
    }
}
