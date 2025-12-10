use std::slice::Iter;

use indexmap::IndexMap;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;

use super::Hst;

#[derive(Clone)]
pub struct PairWiseMatrix(Vec<(usize, IndexMap<usize, (u64, u64, bool)>)>);

impl PairWiseMatrix {
    pub fn iter(&self) -> Iter<'_, (usize, IndexMap<usize, (u64, u64, bool)>)> {
        self.0.iter()
    }

    pub fn mean(&self, idxs: &[usize]) -> f32 {
        let total = idxs.len();

        let mut count = 0;

        let sum: f32 = self
            .iter()
            .filter(|(idx, _map)| idxs.contains(idx))
            .inspect(|_| count += 1)
            .map(|(_idx, map)| {
                let sum: u64 = idxs
                    .iter()
                    .flat_map(|idx| map.get(idx))
                    .map(|(start, stop, is_snp)| match (start == stop, is_snp) {
                        (true, true) => 1,
                        (true, false) => 0,
                        (false, false) => stop.saturating_sub(*start),
                        // (true, true) |
                        // (true, false) |
                        (false, true) => unreachable!(),
                    })
                    .sum();
                sum as f32 / total as f32
            })
            .sum();

        assert_eq!(count, total);

        sum / total as f32
    }
}

impl Hst {
    pub fn calculate_pair_wise(&self) -> PairWiseMatrix {
        calculate_pair_wise_inner(&self)
    }
}

// Recurse up branches from all leaf nodes to find how many bp is shared between the leaf node
// sample and the rest
//
pub fn calculate_pair_wise_inner(hst: &Hst) -> PairWiseMatrix {
    // Iterate all nodes
    let mut rows: Vec<(usize, _)> = hst
        .node_indices()
        // Filter in only leaf nodes
        .filter(|n| hst.n_children(*n) == 0)
        // Flatten the sample indexes of each node
        .flat_map(|n| {
            hst.node_weight(n)
                .unwrap()
                .indexes
                .iter()
                .map(move |v| (n, v))
        })
        .par_bridge()
        .map(|(node_idx, sample_idx)| {
            // Store in all shared lengths for this sample
            let mut shared_lengths = IndexMap::new();

            let data = hst.node_weight(node_idx).unwrap();

            // Get parent node
            let parent_node = hst.get_parent(node_idx).unwrap();
            let parent_data = hst.node_weight(parent_node).unwrap();

            // Iterate over parent indexes
            for parent_sample_idx in &parent_data.indexes {
                //
                // Separate logic if the parent is a root node
                // and the haplotype is just 1 nucleotide long
                if data.stop == data.start {
                    // Samples share the same nucleotide
                    let is_same = data.indexes.contains(parent_sample_idx);

                    shared_lengths.entry(*parent_sample_idx).or_insert((
                        data.start.pos,
                        data.stop.pos,
                        is_same,
                    ));
                } else {
                    let stop = hst.coords().range(..&data.stop).next_back().unwrap();
                    let start = hst
                        .coords()
                        .range(&data.start..)
                        .nth(1)
                        .unwrap_or(hst.coords().range(&data.start..).nth(0).unwrap());

                    shared_lengths
                        .entry(*parent_sample_idx)
                        .or_insert((start.pos, stop.pos, false));
                }
            }

            // Start recursion up the branch to get shared lengths between
            // all samples and this sample
            recurse_branch(hst, parent_node, &mut shared_lengths);

            // Sort by id
            shared_lengths.sort_by_key(|id, _length| *id);

            (*sample_idx, shared_lengths)
        })
        .collect();

    // Sort by id
    rows.sort_by_key(|(id, _row)| *id);

    PairWiseMatrix(rows)
}

fn recurse_branch(
    hst: &Hst,
    node: NodeIndex,
    shared_lengths: &mut IndexMap<usize, (u64, u64, bool)>,
) {
    let data = hst.node_weight(node).unwrap();

    let stop = hst.coords().range(..&data.stop).next_back().unwrap();
    let start = hst
        .coords()
        .range(&data.start..)
        .nth(1)
        .unwrap_or(hst.coords().range(&data.start..).nth(0).unwrap());

    // If there is no parent node i.e. we are a root node
    // check lengths from the node itself and not the parent
    let indexes = match hst.get_parent(node) {
        Some(parent_node) => &hst.node_weight(parent_node).unwrap().indexes,
        None => &data.indexes,
    };

    for parent_sample_idx in indexes {
        if data.stop == data.start {
            let is_same = data.indexes.contains(parent_sample_idx);

            shared_lengths.entry(*parent_sample_idx).or_insert((
                data.start.pos,
                data.stop.pos,
                is_same,
            ));
        } else {
            shared_lengths
                .entry(*parent_sample_idx)
                .or_insert((start.pos, stop.pos, false));
        }
    }

    if let Some(parent_node) = hst.get_parent(node) {
        recurse_branch(hst, parent_node, shared_lengths);
    }
}
