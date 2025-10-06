use color_eyre::eyre::ensure;
use color_eyre::Result;
use petgraph::prelude::NodeIndex;
use petgraph::Direction;
use petgraph::Graph;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::read_vcf::get_sample_names;
use crate::core::read_vcf::read_vcf_to_matrix;
use crate::core::structs::CoordDataSlot;
use crate::core::PhasedMatrix;
use crate::error::Error;
use crate::structs::Coord;
use crate::structs::HapVariant;
use crate::utils::parse_snp_coord;

pub mod bhst;
pub mod immutable_hst;
pub mod pair_wise;
pub mod structs;
pub mod uhst;

pub use structs::Hst;
pub use structs::HstType;
pub use structs::Metadata;
pub use structs::Node;

pub fn read_vcf_with_selections(args: &StandardArgs, window: Option<u64>) -> Result<PhasedMatrix> {
    let (contig, pos) = parse_snp_coord(&args.coords)?;
    let (indexes, _) = get_sample_names(args, contig, None)?;

    // Check that theres at least 3 samples
    ensure!(indexes.len() > 2, Error::HstTooSmall);

    let mut vcf = read_vcf_to_matrix(args, contig, pos, None, None, window, false)?;

    if args.selection == Selection::OnlyLongest {
        vcf.select_only_longest()?
    }

    Ok(vcf)
}

pub fn initiate_hst(
    vcf: &PhasedMatrix,
    start_coord: &Coord,
    indexes: Option<Vec<usize>>,
) -> Graph<Node, ()> {
    let mut hst = Graph::new();

    let indexes = indexes.unwrap_or_else(|| (0..vcf.nhaplotypes()).collect());

    let root = Node {
        start: start_coord.clone(),
        stop: start_coord.clone(),
        indexes: indexes.clone(),
        haplotype: vec![],
    };

    let _ = hst.add_node(root.clone());

    // If the first variant is already contradictory, split into two nodes
    if vcf.is_contradictory(start_coord, &indexes) {
        let genotypes = vcf.matrix_column(start_coord);
        let (mut ones, mut zeroes) = (vec![], vec![]);

        for i in indexes.iter() {
            match genotypes[*i] == 1 {
                true => ones.push(*i),
                false => zeroes.push(*i),
            }
        }

        let node1 = Node {
            start: start_coord.clone(),
            stop: start_coord.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, ones[0]),
            indexes: ones,
        };

        let node2 = Node {
            start: start_coord.clone(),
            stop: start_coord.clone(),
            haplotype: vcf.find_u8_haplotype_for_sample(start_coord..=start_coord, zeroes[0]),
            indexes: zeroes,
        };

        let node_idx1 = hst.add_node(node1.clone());
        let node_idx2 = hst.add_node(node2.clone());

        hst.add_edge(NodeIndex::new(0), node_idx1, ());
        hst.add_edge(NodeIndex::new(0), node_idx2, ());
    };

    hst
}

pub fn find_shared_haplotype(g: &Graph<Node, ()>, vcf: &PhasedMatrix) -> Vec<HapVariant> {
    let start_idx = NodeIndex::new(0);
    let nodes = find_majority_nodes(g, start_idx);
    let (second_node, _second_node_idx) = nodes[1];

    let start = vcf.coords().range(&second_node.start..).nth(1).unwrap();

    let start = match start > &second_node.stop {
        true => &second_node.start,
        false => start,
    };

    vcf.find_haplotype_for_sample(start..&second_node.stop, second_node.indexes[0])
}

pub fn find_majority_nodes(g: &Graph<Node, ()>, start_idx: NodeIndex) -> Vec<(&Node, NodeIndex)> {
    let mut idx = start_idx;
    let start_node = g.node_weight(idx).unwrap();
    let mut majority_nodes = vec![(start_node, idx)];

    // Max by returns None if the iterator is empty, thus when a leaf node is reached
    while let Some(max) = g
        .neighbors_directed(idx, Direction::Outgoing)
        .map(|idx| (g.node_weight(idx).unwrap(), idx))
        .max_by(|x, y| x.0.indexes.len().cmp(&y.0.indexes.len()))
    {
        majority_nodes.push(max);
        idx = max.1;
    }

    majority_nodes
}
