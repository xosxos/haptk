use color_eyre::eyre::ensure;
use color_eyre::Result;
use petgraph::prelude::NodeIndex;
use petgraph::Direction;
use petgraph::Graph;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::read_vcf::read_vcf_to_matrix;
use crate::core::Coord;
use crate::core::HapVariant;
use crate::core::PhasedMatrix;
use crate::error::Error;
use crate::io::get_indexes_and_sample_ids_from_vcf;
use crate::traits::OnlyLongest;
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
    let (indexes, _) = get_indexes_and_sample_ids_from_vcf(&args.file, &args.samples, None)?;

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

use ndarray::ArrayView1;

use crate::core::Matrix;
use crate::read_vcf::read_shard_of_vcf;
use crate::subcommands::uhst::LocDirection;

pub trait CoordDataSlot {
    fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool;

    fn prev_contradictory(
        &'_ self,
        coord: &Coord,
        sample_idxs: &[usize],
    ) -> std::result::Result<Option<(&'_ Coord, ArrayView1<'_, u8>)>, Error>;

    fn next_contradictory(
        &self,
        coord: &Coord,
        sample_idxs: &[usize],
    ) -> std::result::Result<Option<(&'_ Coord, ArrayView1<'_, u8>)>, Error>;

    fn is_file_end(&self, side: LocDirection) -> bool;
    fn read_more(&mut self, error: Error) -> Result<()>;
}

pub fn is_contradictory_by_idx(matrix: &Matrix, positions: &[usize], variant_idx: usize) -> bool {
    let genotypes = matrix.genotypes(variant_idx);
    let first = positions[0];
    genotypes.len() > 1 && positions.iter().any(|x| genotypes[*x] != genotypes[first])
}

impl CoordDataSlot for PhasedMatrix {
    // A column is contractory if it contains both 0 and 1.
    fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool {
        let slot = self.matrix_column(coord);
        let first = positions[0];
        slot.len() > 1 && positions.iter().any(|x| slot[*x] != slot[first])
    }

    fn prev_contradictory(
        &'_ self,
        coord: &Coord,
        positions: &[usize],
    ) -> std::result::Result<Option<(&'_ Coord, ArrayView1<'_, u8>)>, Error> {
        if positions.len() < 2 {
            return Ok(None);
        }

        let (matrix_key, variant_idx) = &self.indexer[coord];

        let (mut travel, mut variant_idx, mut past_first) = (0, *variant_idx as isize, false);

        for (_key, prev_matrix) in self.matrix.range(..=matrix_key.clone()).rev() {
            match past_first {
                true => variant_idx = prev_matrix.ncols() as isize - 1,
                false => variant_idx -= 1,
            }

            while variant_idx >= 0 {
                match is_contradictory_by_idx(prev_matrix, positions, variant_idx as usize) {
                    true => {
                        let new_coord = self.coords().range(..coord).rev().nth(travel).unwrap();

                        let genotypes = prev_matrix.genotypes(variant_idx as usize);

                        return Ok(Some((new_coord, genotypes)));
                    }
                    false => travel += 1,
                }
                variant_idx -= 1;
            }

            past_first = true;
        }

        match self.is_file_end(LocDirection::Left) {
            true => Ok(None),
            false => Err(Error::HstEnd),
        }
    }

    fn next_contradictory(
        &'_ self,
        coord: &Coord,
        positions: &[usize],
    ) -> std::result::Result<Option<(&'_ Coord, ArrayView1<'_, u8>)>, Error> {
        if positions.len() < 2 {
            return Ok(None);
        }
        let (matrix_key, variant_idx) = &self.indexer[coord];

        let (mut travel, mut variant_idx, mut past_first) = (0, *variant_idx, false);

        for (_key, next_matrix) in self.matrix.range(matrix_key.clone()..) {
            match past_first {
                true => variant_idx = 0,
                false => variant_idx += 1,
            }

            while variant_idx < next_matrix.ncols() {
                match is_contradictory_by_idx(next_matrix, positions, variant_idx) {
                    true => {
                        let new_coord = self.coords().range(coord..).skip(1).nth(travel).unwrap();

                        let genotypes = next_matrix.genotypes(variant_idx);

                        return Ok(Some((new_coord, genotypes)));
                    }
                    false => travel += 1,
                };
                variant_idx += 1;
            }

            past_first = true;
        }

        match self.is_file_end(LocDirection::Right) {
            true => Ok(None),
            false => Err(Error::HstEnd),
        }
    }

    fn is_file_end(&self, side: LocDirection) -> bool {
        match side {
            LocDirection::Left => match !self.metadata.sharded {
                true => true,
                false => self.metadata.fetch_range.0 == 0 || !self.metadata.sharded,
            },
            LocDirection::Right => match !self.metadata.sharded {
                true => true,
                false => {
                    self.metadata.fetch_range.1 >= self.metadata.contig_len.unwrap()
                        || !self.metadata.sharded
                }
            },
        }
    }

    fn read_more(&mut self, error: Error) -> Result<()> {
        let start = self.coords().first().unwrap().clone();
        let stop = self.coords().last().unwrap().clone();

        match error {
            Error::HstLeftEnd => {
                read_shard_of_vcf(self, start.pos.saturating_sub(1_000_000), start.pos - 1)?
            }
            Error::HstRightEnd => read_shard_of_vcf(self, stop.pos + 1, stop.pos + 1_000_000)?,
            Error::HstBothEnd => {
                read_shard_of_vcf(self, start.pos.saturating_sub(1_000_000), start.pos - 1)?;
                read_shard_of_vcf(self, stop.pos + 1, stop.pos + 1_000_000)?
            }
            _ => unreachable!(),
        }

        Ok(())
    }
}
