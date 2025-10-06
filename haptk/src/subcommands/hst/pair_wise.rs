use std::sync::mpsc::channel;

use indexmap::IndexMap;
use petgraph::graph::NodeIndex;

use crate::core::PhasedMatrix;
use rayon::prelude::*;

use super::Hst;

pub type PairWiseMatrix = Vec<(usize, IndexMap<usize, (u64, u64, bool)>)>;

// Length is 1 SNP too much to both directions
pub fn calculate_pair_wise(vcf: &PhasedMatrix, hst: &Hst) -> PairWiseMatrix {
    //
    let (tx, rx) = channel();

    // Iterate all leaf nodes
    hst.node_indices()
        .par_bridge()
        .filter(|n| hst.n_children(*n) == 0)
        .for_each(|leaf| {
            let data = hst.node_weight(leaf).unwrap();

            // Get parent node
            let parent_node = hst.get_parent(leaf).unwrap();

            let parent_data = hst.node_weight(parent_node).unwrap();

            // Collect all leaf node idxs
            for leaf_sample_idx in &data.indexes {
                let mut row = IndexMap::new();

                // Iterate leaf node indexes
                for other_sample_idx in &parent_data.indexes {
                    // Insert the length if the pair does not already have a length

                    if data.stop == data.start {
                        let ht_other = vcf.find_u8_haplotype_for_sample(
                            &data.start..=&data.stop,
                            *other_sample_idx,
                        );
                        let ht_this = vcf.find_u8_haplotype_for_sample(
                            &data.start..=&data.stop,
                            *leaf_sample_idx,
                        );
                        row.entry(*other_sample_idx).or_insert((
                            data.start.pos,
                            data.stop.pos,
                            ht_this == ht_other,
                        ));
                    }

                    let stop = vcf.coords().range(..&data.stop).next_back().unwrap();
                    let start = vcf
                        .coords()
                        .range(&data.start..)
                        .nth(1)
                        .unwrap_or(vcf.coords().range(&data.start..).nth(0).unwrap());

                    row.entry(*other_sample_idx)
                        .or_insert((start.pos, stop.pos, false));
                }

                // Start recursion up the branch
                recurse_branch(vcf, hst, parent_node, *leaf_sample_idx, &mut row);

                row.sort_by_key(|id, _length| *id);

                let _ = tx.send((*leaf_sample_idx, row));
            }
        });

    drop(tx);

    let mut rows: PairWiseMatrix = vec![];

    while let Ok(row) = rx.recv() {
        rows.push(row);
    }

    rows.sort_by_key(|(id, _row)| *id);

    rows
}

fn recurse_branch(
    vcf: &PhasedMatrix,
    hst: &Hst,
    node: NodeIndex,
    leaf_idx: usize,
    row: &mut IndexMap<usize, (u64, u64, bool)>,
) {
    let data = hst.node_weight(node).unwrap();

    let stop = vcf.coords().range(..&data.stop).next_back().unwrap();
    let start = vcf
        .coords()
        .range(&data.start..)
        .nth(1)
        .unwrap_or(vcf.coords().range(&data.start..).nth(0).unwrap());

    if let Some(parent_node) = hst.get_parent(node) {
        let parent_data = hst.node_weight(parent_node).unwrap();

        for other_sample_idx in &parent_data.indexes {
            if data.stop == data.start {
                let ht_other =
                    vcf.find_u8_haplotype_for_sample(&data.start..=&data.stop, *other_sample_idx);
                let ht_this = vcf.find_u8_haplotype_for_sample(&data.start..=&data.stop, leaf_idx);

                row.entry(*other_sample_idx).or_insert((
                    data.start.pos,
                    data.stop.pos,
                    ht_other == ht_this,
                ));
            } else {
                row.entry(*other_sample_idx)
                    .or_insert((start.pos, stop.pos, false));
            }
        }

        recurse_branch(vcf, hst, parent_node, leaf_idx, row);
    } else {
        for other_sample_idx in &data.indexes {
            let ht_other =
                vcf.find_u8_haplotype_for_sample(&data.start..=&data.stop, *other_sample_idx);
            let ht_this = vcf.find_u8_haplotype_for_sample(&data.start..=&data.stop, leaf_idx);

            if data.stop == data.start {
                row.entry(*other_sample_idx).or_insert((
                    data.start.pos,
                    data.stop.pos,
                    ht_other == ht_this,
                ));
            } else {
                row.entry(*other_sample_idx)
                    .or_insert((start.pos, stop.pos, false));
            }
        }
    }
}
