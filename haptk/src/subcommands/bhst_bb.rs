use std::fmt;

use color_eyre::eyre::ensure;
use ndarray::s;
use petgraph::prelude::NodeIndex;
use petgraph::{Direction, Graph};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::{Selection, StandardArgs};
use crate::core::parse_snp_coord;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::{Coord, HapVariant, PhasedMatrix};

#[derive(Serialize, Deserialize, Default, Clone, Debug, PartialEq, PartialOrd)]
pub struct Node {
    pub indexes: Vec<usize>,
    pub haplotype: Vec<HapVariant>,
}

impl std::fmt::Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "n: {}\nnmarkers: {}\nblock length (bp): {}\n",
            self.indexes.len(),
            self.haplotype.len(),
            (self.haplotype.last().unwrap().pos - self.haplotype.first().unwrap().pos) as f32
        );
        write!(f, "{line}")
    }
}

pub trait Indexes {
    fn indexes(&self) -> &Vec<usize>;
    fn block_len(&self) -> f32;
    fn show(&self) -> String;
}

impl Indexes for Node {
    fn indexes(&self) -> &Vec<usize> {
        &self.indexes
    }
    fn block_len(&self) -> f32 {
        (self.haplotype.last().unwrap().pos - self.haplotype.first().unwrap().pos) as f32
    }

    fn show(&self) -> String {
        self.to_string()
    }
}

type Result<T> = std::result::Result<T, HstError>;

enum HstError {
    PrematureEndOfSequencingData,
    EndOfSequencingData,
    EndOfLeafNodes,
}

impl fmt::Display for HstError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::PrematureEndOfSequencingData => write!(f, "Premature end of sequencing data"),
            Self::EndOfSequencingData => write!(f, "End of sequencing data"),
            Self::EndOfLeafNodes => write!(f, "End of leaf nodes"),
        }
    }
}

#[doc(hidden)]
pub fn run(args: StandardArgs) -> color_eyre::Result<()> {
    ensure!(
        args.selection == Selection::All,
        "For the biobank version of hst, only selecting all is supported for now."
    );

    let (contig, pos) = parse_snp_coord(&args.coords)?;

    let mut vcf = read_vcf_to_matrix(&args, contig, pos, None, None)?;
    let min_size = 1;

    let mut bhst = Graph::<_, _>::new();

    bhst.add_node(Node {
        indexes: (0..vcf.matrix.nrows()).collect(),
        haplotype: vec![],
    });

    loop {
        // Filter out nodes with children and nodes with less indexes than min_size
        let indices = bhst
            .node_indices()
            .filter(|node_idx| {
                let node = bhst.node_weight(*node_idx).unwrap();
                let count = bhst
                    .neighbors_directed(*node_idx, Direction::Outgoing)
                    .count();

                count == 0 && node.indexes.len() > min_size
            })
            .collect::<Vec<NodeIndex>>();

        // Multithread horizontally all childless nodes
        let nodes = indices
            .into_par_iter()
            .map(|node_idx| {
                find_contradictory_gt(&vcf, &bhst, node_idx).map(|nodes| (node_idx, nodes))
            })
            .collect::<Result<Vec<(NodeIndex, Vec<Node>)>>>();

        let nodes: Vec<(NodeIndex, Vec<Node>)> = match nodes {
            Ok(nodes) => nodes,
            Err(HstError::PrematureEndOfSequencingData) => {
                vcf = read_more_vcf();
                continue;
            }
            Err(HstError::EndOfLeafNodes) => break,
            Err(HstError::EndOfSequencingData) => unreachable!(),
        };

        // Add nodes to the tree
        for (node_idx, new_nodes) in nodes {
            for new_node in new_nodes {
                let new_node_idx = bhst.add_node(new_node.clone());

                bhst.add_edge(node_idx, new_node_idx, 0);
            }
        }
    }

    tracing::debug!("Finished HST construction: {bhst:?}");

    Ok(())
}

fn find_previous_idx_from_hapvariant(vcf: &PhasedMatrix, hap: &HapVariant) -> Result<usize> {
    todo!()
}

fn find_next_idx_from_hapvariant(vcf: &PhasedMatrix, hap: &HapVariant) -> Result<usize> {
    todo!()
}

fn find_coord(vcf: &PhasedMatrix, idx: usize) -> Coord {
    todo!()
}

fn coord_to_hapvariant(coord: Coord, gt: u8) -> HapVariant {
    todo!()
}

fn find_haplotype(haplotype: &Vec<HapVariant>, left: u8, right: u8) -> Vec<HapVariant> {
    todo!()
}

fn read_more_vcf() -> PhasedMatrix {
    todo!()
}

fn find_contradictory_gt(
    vcf: &PhasedMatrix,
    bhst: &Graph<Node, u8>,
    node_idx: NodeIndex,
) -> Result<Vec<Node>> {
    let node = bhst.node_weight(node_idx).unwrap();

    let (mut left, mut right) = if node_idx == NodeIndex::new(0) {
        let left = find_coord(vcf, vcf.variant_idx());
        // Minus 1 to account for the starting variant itself as well
        let right = find_coord(vcf, vcf.variant_idx() - 1);
        (coord_to_hapvariant(left, 1), coord_to_hapvariant(right, 1))
    } else {
        (
            node.haplotype.first().unwrap().clone(),
            node.haplotype.last().unwrap().clone(),
        )
    };

    let (mut l_haplotype, mut r_haplotype) = (vec![], vec![]);

    let (mut rz, mut ro);
    // Expand to the right until a contradictory genotype is found or the end of sequencing data
    loop {
        (rz, ro) = (vec![], vec![]);
        let idx = find_next_idx_from_hapvariant(vcf, &right);

        let idx = match idx {
            Ok(idx) => idx,
            Err(e) => match e {
                HstError::EndOfSequencingData => break,
                HstError::PrematureEndOfSequencingData => {
                    return Err(HstError::PrematureEndOfSequencingData)
                }
                HstError::EndOfLeafNodes => unreachable!(),
            },
        };

        for i in node.indexes.iter() {
            let right_gt = vcf.matrix.slice(s![*i, idx]);

            match right_gt.into_scalar() {
                0 => rz.push(*i),
                1 => ro.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        if !rz.is_empty() && !ro.is_empty() {
            break;
        }

        let coord = find_coord(vcf, idx);
        let hapvariant = if rz.is_empty() {
            coord_to_hapvariant(coord, 1)
        } else {
            coord_to_hapvariant(coord, 0)
        };
        r_haplotype.push(hapvariant.clone());
        right = hapvariant;
    }

    let (mut lz, mut lo);
    // Expand to the left until a contradictory genotype is found or the end of sequencing data
    loop {
        (lz, lo) = (vec![], vec![]);
        let idx = find_previous_idx_from_hapvariant(vcf, &left);

        let idx = match idx {
            Ok(idx) => idx,
            Err(e) => match e {
                HstError::EndOfSequencingData => break,
                HstError::PrematureEndOfSequencingData => {
                    return Err(HstError::PrematureEndOfSequencingData)
                }
                HstError::EndOfLeafNodes => unreachable!(),
            },
        };

        for i in node.indexes.iter() {
            let left_gt = vcf.matrix.slice(s![*i, idx]);

            match left_gt.into_scalar() {
                0 => lz.push(*i),
                1 => lo.push(*i),
                _ => unreachable!("Other genotypes than 0 and 1 are present in the matrix"),
            }
        }

        if !lz.is_empty() && !lo.is_empty() {
            break;
        }

        let coord = find_coord(vcf, idx);
        let hapvariant = if lz.is_empty() {
            coord_to_hapvariant(coord, 1)
        } else {
            coord_to_hapvariant(coord, 0)
        };

        l_haplotype.push(hapvariant.clone());
        left = hapvariant;
    }

    l_haplotype.extend(r_haplotype);

    // Divide left and right side contradictory genotypes to buckets:
    // [0,1] [0,0] [1, 0] [0, 0]
    let (mut zo, mut zz, mut oz, mut oo) = (vec![], vec![], vec![], vec![]);

    for i in &lz {
        if rz.contains(i) {
            zz.push(*i);
        } else if ro.contains(i) {
            zo.push(*i)
        }
    }

    for i in &lo {
        if rz.contains(i) {
            oz.push(*i);
        } else if ro.contains(i) {
            oo.push(*i)
        }
    }

    // Create a new node for each bucket if it is not empty
    let mut nodes = Vec::with_capacity(4);

    if !zz.is_empty() {
        let node = Node {
            haplotype: find_haplotype(&l_haplotype, 0, 0),
            indexes: zz,
        };
        nodes.push(node);
    }

    if !zo.is_empty() {
        let node = Node {
            haplotype: find_haplotype(&l_haplotype, 0, 1),
            indexes: zo,
        };
        nodes.push(node);
    }

    if !oz.is_empty() {
        let node = Node {
            haplotype: find_haplotype(&l_haplotype, 1, 0),
            indexes: oz,
        };
        nodes.push(node);
    }

    if !oo.is_empty() {
        let node = Node {
            haplotype: find_haplotype(&l_haplotype, 1, 1),
            indexes: oo,
        };
        nodes.push(node);
    }

    // Return buckets if more than one bucket exists
    match nodes.len() > 1 {
        true => Ok(nodes),
        false => Err(HstError::EndOfLeafNodes),
    }
}
