use std::path::PathBuf;

use crate::error::Error;
use crate::subcommands::hst::Hst;

#[cfg_attr(feature = "clap", derive(clap::Args))]
#[derive(Debug, Clone)]
pub struct Args {
    /// A file path to a HST file
    hst: PathBuf,

    /// Algorithm for finding subtrees
    #[cfg_attr(feature = "clap", arg(long, value_enum, default_value_t = Algorithm::MajorityBranch))]
    algorithm: Algorithm,

    /// Max amount of indexes in subtree
    #[cfg_attr(feature = "clap", arg(long, default_value_t = usize::MAX))]
    max_node_size: usize,

    /// Min amount of indexes in subtree
    #[cfg_attr(feature = "clap", arg(long, default_value_t = 10))]
    min_node_size: usize,
}

#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
#[derive(Debug, Clone, Default)]
pub enum Algorithm {
    #[default]
    MajorityBranch,
}

#[doc(hidden)]
pub fn run(args: Args) -> Result<(), Error> {
    let hst = Hst::from_file(&args.hst)?;

    match args.algorithm {
        Algorithm::MajorityBranch => majority_branch(&hst, &args),
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

        if branch_size > args.min_node_size && branch_size < args.max_node_size {
            let siblings = hst.get_siblings(idx);

            if !siblings.is_empty() {
                let sibling = siblings[0];
                let sibling = hst.node_weight(sibling).unwrap();

                let ht_length = node_data.stop.pos.saturating_sub(node_data.start.pos);

                println!(
                    "{:?},{:?},{:?},{:?}",
                    idx.index(),
                    ht_length,
                    node_data.indexes.len(),
                    sibling.indexes.len()
                );
            }
        }

        parent_node_indexes = indexes_len;
    }
}
