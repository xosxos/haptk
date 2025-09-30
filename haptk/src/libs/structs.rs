use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::HashMap;
use std::ops::Bound;
use std::ops::RangeBounds;
use std::path::PathBuf;
use std::sync::Arc;

use color_eyre::{eyre::ensure, Result};
use ndarray::iter::AxisIter;
use ndarray::ArrayView2;
use ndarray::{s, Array2, ArrayView1, Axis};
use petgraph::graph::NodeIndex;
use petgraph::Graph;
use serde::{Deserialize, Serialize};

use crate::args::Selection;
use crate::error::Error;
use crate::read_vcf::read_shard_of_vcf;
use crate::subcommands::bhst::{find_majority_nodes, Node};
use crate::subcommands::immutable_hst;
use crate::subcommands::uhst::{self, LocDirection};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Coord {
    pub contig: String,
    pub pos: u64,
    #[serde(rename(serialize = "ref"), alias = "ref")]
    pub reference: String,
    pub alt: String,
}

impl std::hash::Hash for Coord {
    fn hash<H>(&self, state: &mut H)
    where
        H: std::hash::Hasher,
    {
        self.pos.hash(state);
        self.alt.hash(state);
    }
}

impl Ord for Coord {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.pos.cmp(&other.pos) {
            Ordering::Less => Ordering::Less,
            Ordering::Equal => {
                let cmp = self.reference.cmp(&other.reference);
                if cmp == Ordering::Equal {
                    if self.alt == "-" || other.alt == "-" {
                        Ordering::Equal
                    } else {
                        self.alt.cmp(&other.alt)
                    }
                } else {
                    cmp
                }
            }
            Ordering::Greater => Ordering::Greater,
        }
    }
}

impl PartialOrd for Coord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for Coord {
    fn eq(&self, other: &Self) -> bool {
        (&self.contig, self.pos, &self.reference, &self.alt)
            == (&other.contig, other.pos, &other.reference, &other.alt)
    }
}

impl Eq for Coord {}

impl From<HapVariant> for Coord {
    fn from(hap: HapVariant) -> Self {
        Self {
            contig: hap.contig,
            pos: hap.pos,
            reference: hap.reference,
            alt: hap.alt,
        }
    }
}

#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct HapVariant {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub gt: u8,
}

impl HapVariant {
    pub fn genotype(&self) -> &String {
        match self.gt {
            0 => &self.reference,
            1 => &self.alt,
            _ => panic!("Please normalize the vcf using bcftools norm"),
        }
    }
}

impl std::fmt::Display for HapVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "{}_{}_{}_{}_{}",
            self.contig, self.pos, self.reference, self.alt, self.gt
        );
        write!(f, "{line}")
    }
}

impl PartialEq<HapVariant> for Coord {
    fn eq(&self, other: &HapVariant) -> bool {
        if other.alt == "-" && other.gt == 0 {
            self.contig == other.contig
                && self.pos == other.pos
                && self.reference == other.reference
        } else {
            self.contig == other.contig
                && self.pos == other.pos
                && self.reference == other.reference
                && self.alt == other.alt
        }
    }
}

impl PartialEq<Coord> for HapVariant {
    fn eq(&self, other: &Coord) -> bool {
        if self.alt == "-" && self.gt == 0 {
            self.contig == other.contig
                && self.pos == other.pos
                && self.reference == other.reference
        } else {
            self.contig == other.contig
                && self.pos == other.pos
                && self.reference == other.reference
                && self.alt == other.alt
        }
    }
}

impl std::fmt::Display for Coord {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "{}_{}_{}_{}",
            self.contig, self.pos, self.reference, self.alt
        );
        write!(f, "{line}")
    }
}

#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CigarVariant {
    pub pos: u64,
    pub alt: char,
}

impl Ord for CigarVariant {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.pos < other.pos {
            return Ordering::Less;
        }

        if self.pos > other.pos {
            return Ordering::Greater;
        }

        return self.alt.cmp(&other.alt);
    }
}

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub enum Ploidy {
    Mixed,
    Haploid,
    #[default]
    Diploid,
}

impl From<Ploidy> for usize {
    fn from(ploidy: Ploidy) -> Self {
        match ploidy {
            Ploidy::Haploid => 1,
            Ploidy::Mixed => 1,
            Ploidy::Diploid => 2,
        }
    }
}

impl From<&Ploidy> for usize {
    fn from(ploidy: &Ploidy) -> Self {
        match ploidy {
            Ploidy::Haploid => 1,
            Ploidy::Mixed => 1,
            Ploidy::Diploid => 2,
        }
    }
}

impl std::ops::Deref for Ploidy {
    type Target = usize;
    fn deref(&self) -> &Self::Target {
        match self {
            Ploidy::Haploid => &1,
            Ploidy::Mixed => &1,
            Ploidy::Diploid => &2,
        }
    }
}

impl From<&Selection> for Ploidy {
    fn from(selection: &Selection) -> Self {
        match selection {
            Selection::Haploid | Selection::OnlyRefs | Selection::OnlyAlts | Selection::List => {
                Ploidy::Haploid
            }
            _ => Ploidy::Diploid,
        }
    }
}

impl From<Selection> for Ploidy {
    fn from(selection: Selection) -> Self {
        match selection {
            Selection::Haploid | Selection::OnlyRefs | Selection::OnlyAlts | Selection::List => {
                Ploidy::Haploid
            }
            _ => Ploidy::Diploid,
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReadMetadata {
    pub file_path: PathBuf,
    pub fetch_range: (u64, u64),
    pub lookups: Vec<[bool; 2]>,
    pub indexes: Vec<usize>,
    pub contig_len: Option<u64>,
    pub sharded: bool,
    pub remove_no_alt: bool,
    pub include_indels: bool,
    pub is_genome_wide: bool,
}

#[derive(Debug, Default, Clone)]
pub struct PhasedMatrix {
    pub start_coord: Coord,
    pub variant_idx: usize,
    pub matrix: BTreeMap<Arc<Coord>, Array2<u8>>,
    samples: Vec<String>,
    coords: BTreeSet<Coord>,
    pub indexer: HashMap<Coord, (Arc<Coord>, usize)>,
    pub ploidy: Ploidy,
    pub metadata: ReadMetadata,
}

pub type Offset = (Coord, Coord);

impl PhasedMatrix {
    pub fn new<T: AsRef<Selection>>(
        variant_idx: usize,
        start_coord: Coord,
        matrix: Array2<u8>,
        samples: Vec<String>,
        coords: BTreeSet<Coord>,
        selection: T,
        metadata: ReadMetadata,
    ) -> Self {
        let start = Arc::new(coords.first().unwrap().clone());

        let indexer: HashMap<Coord, (Arc<Coord>, usize)> = coords
            .iter()
            .enumerate()
            .map(|(i, v)| (v.clone(), (start.clone(), i)))
            .collect();

        let mut map = BTreeMap::new();
        map.insert(start.clone(), matrix);

        Self {
            variant_idx,
            start_coord,
            matrix: map,
            samples,
            coords,
            ploidy: selection.as_ref().into(),
            indexer,
            metadata,
        }
    }

    pub fn nsamples(&self) -> usize {
        self.samples.len()
    }

    pub fn nhaplotypes(&self) -> usize {
        self.samples.len() * *self.ploidy
    }

    pub fn matrix_nrows(&self) -> usize {
        self.matrix.values().map(|m| m.nrows()).sum()
    }

    pub fn matrix_ncols(&self) -> usize {
        self.matrix.values().map(|m| m.ncols()).sum()
    }

    pub fn matrix_point_coord(&self, x: usize, coord: &Coord) -> u8 {
        match self.indexer.get(coord) {
            Some((matrix_key, index)) => {
                let matrix = self.matrix.get(matrix_key).unwrap();
                matrix[[x, *index]]
            }
            None => {
                let nearest_coord = self.get_nearest_coord_by_pos(coord.pos);
                tracing::warn!(
                    "Could not find coord {}, using {} instead",
                    coord,
                    nearest_coord
                );
                let (matrix_key, index) = self.indexer.get(nearest_coord).unwrap();
                let matrix = self.matrix.get(matrix_key).unwrap();
                matrix[[x, *index]]
            }
        }
    }

    pub fn matrix_column(&self, coord: &'_ Coord) -> ArrayView1<'_, u8> {
        let (matrix, index) = self
            .indexer
            .get(coord)
            .unwrap_or_else(|| panic!("Could not find coord {coord} from the indexer"));

        let matrix = self.matrix.get(matrix).unwrap();

        matrix.index_axis(Axis(1), *index)
    }

    pub fn insert_matrix(&mut self, start_coord: Coord, matrix: Array2<u8>) {
        self.matrix.insert(Arc::new(start_coord), matrix);
    }

    // Sort inside select_rows to avoid bugs down the line
    pub fn select_rows(&mut self, mut to_keep: Vec<usize>) {
        to_keep.sort();
        self.matrix_select(0, &to_keep);

        self.samples = to_keep
            .iter()
            .map(|index| self.get_sample_name(*index))
            .collect();
    }

    pub fn matrix_select(&mut self, axis: usize, to_keep: &[usize]) {
        if axis == 0 {
            for (_, matrix) in self.matrix.iter_mut() {
                *matrix = matrix.select(Axis(0), to_keep)
            }
        } else {
            unreachable!("matrix select on columns is not yet implemented")
        }
    }

    pub fn ncoords(&self) -> usize {
        self.coords.len()
    }

    pub fn samples(&self) -> &Vec<String> {
        &self.samples
    }

    pub fn is_genome_wide(&self) -> bool {
        self.metadata.is_genome_wide
    }

    pub fn variant_idx(&self) -> usize {
        self.variant_idx
    }

    pub fn set_variant_idx(&mut self, variant_idx: usize) {
        let new_start_coord = self.get_coord(variant_idx).clone();
        self.set_start_coord(new_start_coord);

        self.variant_idx = variant_idx;
    }

    pub fn variant_idx_pos(&self) -> u64 {
        if !self.coords.is_empty() {
            self.start_coord.pos
        } else {
            0
        }
    }

    pub fn get_pos(&self, idx: usize) -> u64 {
        self.coords.iter().nth(idx).unwrap().pos
    }

    pub fn get_coord(&self, idx: usize) -> &Coord {
        self.coords.iter().nth(idx).unwrap()
    }

    pub fn coords(&self) -> &BTreeSet<Coord> {
        &self.coords
    }

    pub fn coords_mut(&mut self) -> &mut BTreeSet<Coord> {
        &mut self.coords
    }

    pub fn set_coords(&mut self, coords: BTreeSet<Coord>) {
        self.coords = coords;
    }

    pub fn start_coord(&self) -> &Coord {
        &self.start_coord
    }

    pub fn get_coord_idx(&self, coord: &Coord) -> usize {
        self.coords.range(..=coord).count() - 1
    }

    pub fn set_start_coord(&mut self, coord: Coord) {
        self.start_coord = coord;
    }

    pub fn get_contig(&self) -> &String {
        &self.coords.first().unwrap().contig
    }

    pub fn get_nearest_coord_by_pos(&self, pos: u64) -> &Coord {
        if let Some(idx) = self.coords.iter().position(|c| c.pos >= pos) {
            let coord = self.coords.iter().nth(idx).unwrap();

            if idx == 0 || coord.pos == pos {
                return coord;
            }

            let before = self.coords.iter().nth(idx - 1).unwrap();

            match pos - before.pos >= coord.pos - pos {
                true => coord,
                false => before,
            }
        } else {
            // return last element if no element is larger
            self.coords.last().unwrap()
        }
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.samples.get(index / *self.ploidy).unwrap().clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_idxs_for_samples(&self, samples: &[String]) -> Result<Vec<usize>> {
        let idxs: Vec<_> = self
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| samples.contains(s))
            .flat_map(|(i, _)| {
                ((i * *self.ploidy)..(i * *self.ploidy) + *self.ploidy).collect::<Vec<usize>>()
            })
            .collect();

        ensure!(
            !idxs.is_empty(),
            "None of the control samples are found in the vcf."
        );
        Ok(idxs)
    }

    pub fn select_only_longest(&mut self) -> Result<()> {
        let longest_indexes = self.only_longest_indexes()?;

        // Update lookups
        let mut lookups = vec![];

        for idx in &longest_indexes {
            let idxs = self.get_idxs_for_samples(&[self.get_sample_name(*idx)])?;
            let pos = idxs.iter().position(|i| idx == i).unwrap();
            let lookup = match pos {
                0 => [true, false],
                1 => [false, true],
                _ => unreachable!("Only diploid genotypes are supported"),
            };
            lookups.push(lookup);
        }

        self.metadata.lookups = lookups;

        self.select_rows(longest_indexes);
        self.ploidy = Ploidy::Haploid;

        tracing::info!("Finished only-longest selection.");
        Ok(())
    }

    pub fn select_only_longest_no_shard(&mut self) -> Result<()> {
        let longest_indexes = self.only_longest_indexes_no_shard(self.start_coord())?;

        self.select_rows(longest_indexes);
        self.ploidy = Ploidy::Haploid;

        tracing::info!("Finished only-longest selection.");
        Ok(())
    }

    pub fn only_longest_indexes(&mut self) -> Result<Vec<usize>> {
        let start = self.start_coord().clone();
        let lengths = self.get_lengths_from_uhst(&start)?;

        Ok(self.lengths_to_indexes(lengths))
    }

    pub fn get_only_longest_lookups(&mut self) -> Result<Vec<[bool; 2]>> {
        let indexes = self.only_longest_indexes()?;

        let mut lookups = vec![];

        for idx in indexes {
            let idxs = self.get_idxs_for_samples(&[self.get_sample_name(idx)])?;
            let pos = idxs.iter().position(|i| idx == *i).unwrap();
            let lookup = match pos {
                0 => [true, false],
                1 => [false, true],
                _ => unreachable!("Only diploid genotypes are supported"),
            };
            lookups.push(lookup);
        }

        Ok(lookups)
    }

    pub fn get_only_longest_lookups_no_shard(&self) -> Result<Vec<[bool; 2]>> {
        let indexes = self.only_longest_indexes_no_shard(self.start_coord())?;

        let mut lookups = vec![];

        for idx in indexes {
            let idxs = self.get_idxs_for_samples(&[self.get_sample_name(idx)])?;
            let pos = idxs.iter().position(|i| idx == *i).unwrap();
            let lookup = match pos {
                0 => [true, false],
                1 => [false, true],
                _ => unreachable!("Only diploid genotypes are supported"),
            };
            lookups.push(lookup);
        }

        Ok(lookups)
    }

    pub fn only_longest_indexes_no_shard(&self, start: &Coord) -> Result<Vec<usize>> {
        let lengths = self.get_lengths_from_uhst_no_mut(start)?;

        Ok(self.lengths_to_indexes(lengths))
    }

    fn lengths_to_indexes(&self, lengths: Vec<(Node, Node)>) -> Vec<usize> {
        let lengths = lengths
            .iter()
            .map(|(lnode, rnode)| {
                if lnode.start == rnode.stop {
                    0
                } else {
                    let stop = self.coords().range(..=&rnode.stop).next_back().unwrap();

                    let start = self.coords().range(&lnode.start..).nth(1).unwrap();
                    stop.pos - start.pos + 1
                }
            })
            .collect::<Vec<u64>>();

        (0..self.nsamples())
            .map(|i| {
                let lengths = lengths
                    .iter()
                    .enumerate()
                    .skip(i * *self.ploidy)
                    .take(*self.ploidy);

                let (max_idx, max_len) = lengths.clone().max_by_key(|(_, l)| *l).unwrap();

                if !self.is_genome_wide() && lengths.filter(|(_, l)| *l == max_len).count() > 1 {
                    tracing::warn!(
                        "Sample {} has two equally long haplotypes in only-longest selection.",
                        self.get_sample_name(i)
                    );
                }

                max_idx
            })
            .collect()
    }

    pub fn only_longest_lengths(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let lengths = self.get_lengths_from_uhst(start_coord)?;

        Ok(self.nodes_to_lengths(lengths))
    }

    pub fn only_longest_lengths_no_shard(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let lengths = self.get_lengths_from_uhst_no_mut(start_coord)?;

        Ok(self.nodes_to_lengths(lengths))
    }

    fn nodes_to_lengths(&self, lengths: Vec<(Node, Node)>) -> Vec<(Node, Node)> {
        let calculate_len = |(lnode, rnode): &(Node, Node)| {
            if lnode.start == rnode.stop {
                0
            } else {
                let stop = self.coords().range(..&rnode.stop).next_back().unwrap();
                let start = self.coords().range(&lnode.start..).nth(1).unwrap();
                stop.pos.saturating_sub(start.pos + 1)
            }
        };

        (0..self.nsamples())
            .map(|i| {
                let lengths = lengths
                    .iter()
                    .enumerate()
                    .skip(i * *self.ploidy)
                    .take(*self.ploidy);

                let (_, max_nodes) = lengths
                    .clone()
                    .max_by_key(|(_, l)| calculate_len(l))
                    .unwrap();

                let max_len = calculate_len(max_nodes);

                if !self.is_genome_wide()
                    && lengths.filter(|(_, l)| calculate_len(l) == max_len).count() > 1
                {
                    tracing::warn!(
                        "Sample {} has two equally long haplotypes.",
                        self.get_sample_name(i)
                    );
                }

                max_nodes.clone()
            })
            .collect()
    }

    pub fn get_lengths_from_uhst(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let uhst_left = uhst::construct_uhst(self, &LocDirection::Left, start_coord, 1, true)?;
        let uhst_right = uhst::construct_uhst(self, &LocDirection::Right, start_coord, 1, true)?;

        Ok(self.get_lengths(uhst_left, uhst_right))
    }

    pub fn get_lengths_from_uhst_no_mut(&self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let uhst_right =
            immutable_hst::construct_uhst_no_mut(self, &LocDirection::Right, start_coord, 1, true)?;
        let uhst_left =
            immutable_hst::construct_uhst_no_mut(self, &LocDirection::Left, start_coord, 1, true)?;

        Ok(self.get_lengths(uhst_left, uhst_right))
    }

    fn get_lengths(
        &self,
        uhst_left: Graph<Node, ()>,
        uhst_right: Graph<Node, ()>,
    ) -> Vec<(Node, Node)> {
        let start_idx = NodeIndex::new(0);
        let lmaj_branch = find_majority_nodes(&uhst_left, start_idx);
        let rmaj_branch = find_majority_nodes(&uhst_right, start_idx);

        (0..self.nhaplotypes())
            .map(|idx| {
                let mut lnode = lmaj_branch.len()
                    - lmaj_branch
                        .iter()
                        .rev()
                        .position(|(node, _)| node.indexes.contains(&idx))
                        .unwrap();
                let mut rnode = rmaj_branch.len()
                    - rmaj_branch
                        .iter()
                        .rev()
                        .position(|(node, _)| node.indexes.contains(&idx))
                        .unwrap();
                if rnode == rmaj_branch.len() {
                    rnode -= 1;
                }

                if lnode == lmaj_branch.len() {
                    lnode -= 1;
                }

                (
                    idx,
                    lmaj_branch.get(lnode).unwrap(),
                    rmaj_branch.get(rnode).unwrap(),
                )
            })
            .map(|(_idx, (lnode, _), (rnode, _))| ((*lnode).clone(), (*rnode).clone()))
            .collect()
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_haplotype_for_sample<R: RangeBounds<Coord> + std::fmt::Debug>(
        &self,
        range: R,
        sample: usize,
    ) -> Vec<HapVariant> {
        self.coords()
            .range(range)
            .map(|coord| {
                let gt = self.matrix_point_coord(sample, coord);

                HapVariant {
                    contig: coord.contig.to_string(),
                    pos: coord.pos,
                    alt: coord.alt.clone(),
                    reference: coord.reference.clone(),
                    gt,
                }
            })
            .collect()
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_u8_haplotype_for_sample<R: RangeBounds<Coord>>(
        &self,
        col_range: R,
        sample: usize,
    ) -> Vec<u8> {
        let (col_first_idx, col_last_idx) = match (col_range.start_bound(), col_range.end_bound()) {
            (Bound::Included(first), Bound::Included(last)) => (first, last),
            _ => unreachable!(
                "Range for find_u8_haplotype_for_sample should always be inclusive x..=y"
            ),
        };
        let (matrix_key_first, index_first) = &self.indexer[col_first_idx];
        let (matrix_key_last, index_last) = &self.indexer[col_last_idx];

        if matrix_key_first == matrix_key_last {
            let matrix = self.matrix.get(matrix_key_first).unwrap();
            return matrix
                .slice(s![sample, *index_first..=*index_last])
                .to_vec();
        }

        self.matrix
            .range(matrix_key_first.clone()..matrix_key_last.clone())
            .flat_map(|(key, next_matrix)| {
                let index_stop = match key == matrix_key_last {
                    true => *index_last,
                    false => next_matrix.ncols() - 1,
                };

                let index_start = match key == matrix_key_first {
                    true => *index_first,
                    false => 0,
                };

                next_matrix
                    .slice(s![sample, index_start..=index_stop])
                    .to_vec()
            })
            .collect()
    }

    // NOTE: Legacy compare-to-haplotype methods, will be deprecated when index based access to matrices is completely removed
    pub fn matrix_axis_iter(&self, axis: usize) -> AxisIter<'_, u8, ndarray::Dim<[usize; 1]>> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix.axis_iter(Axis(axis))
    }

    pub fn set_matrix(&mut self, start_coord: Coord, matrix: Array2<u8>) {
        let mut map = BTreeMap::new();

        map.insert(Arc::new(start_coord), matrix);
        self.matrix = map;
    }

    pub fn slice_cols<R: RangeBounds<usize>>(&'_ self, col_range: R) -> ArrayView2<'_, u8> {
        let (col_first_idx, col_last_idx) = match (col_range.start_bound(), col_range.end_bound()) {
            (Bound::Included(first), Bound::Excluded(last)) => (*first, last.saturating_sub(1)),
            (Bound::Included(first), Bound::Included(last)) => (*first, *last),
            _ => unreachable!(
                "Range for find_u8_haplotype_for_sample should always be x..=y or x..y"
            ),
        };

        let (_, matrix) = self.matrix.iter().nth(0).unwrap();

        matrix.slice(s![.., col_first_idx..=col_last_idx])
    }

    pub fn matrix(&self) -> &Array2<u8> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix
    }
}

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

pub fn is_contradictory_by_idx(matrix: &Array2<u8>, positions: &[usize], index: usize) -> bool {
    let slot = matrix.index_axis(Axis(1), index);
    let first = positions[0];
    slot.len() > 1 && positions.iter().any(|x| slot[*x] != slot[first])
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

        let (matrix_key, index) = &self.indexer[coord];

        let (mut travel, mut index, mut past_first) = (0, *index as isize, false);

        for (_key, prev_matrix) in self.matrix.range(..=matrix_key.clone()).rev() {
            match past_first {
                true => index = prev_matrix.ncols() as isize - 1,
                false => index -= 1,
            }

            while index >= 0 {
                match is_contradictory_by_idx(prev_matrix, positions, index as usize) {
                    true => {
                        let new_coord = self.coords().range(..coord).rev().nth(travel).unwrap();
                        let view = prev_matrix.index_axis(Axis(1), index as usize);
                        return Ok(Some((new_coord, view)));
                    }
                    false => travel += 1,
                }
                index -= 1;
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
        let (matrix_key, index) = &self.indexer[coord];

        let (mut travel, mut index, mut past_first) = (0, *index, false);

        for (_key, next_matrix) in self.matrix.range(matrix_key.clone()..) {
            match past_first {
                true => index = 0,
                false => index += 1,
            }

            while index < next_matrix.ncols() {
                match is_contradictory_by_idx(next_matrix, positions, index) {
                    true => {
                        let new_coord = self.coords().range(coord..).skip(1).nth(travel).unwrap();
                        let view = next_matrix.index_axis(Axis(1), index);
                        return Ok(Some((new_coord, view)));
                    }
                    false => travel += 1,
                };
                index += 1;
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
        let start = self.coords.first().unwrap().clone();
        let stop = self.coords.last().unwrap().clone();

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coords_displays() {
        let coord = Coord {
            contig: "chr9".to_string(),
            pos: 25,
            reference: "G".to_string(),
            alt: "T".to_string(),
        };

        assert_eq!("chr9_25_G_T".to_string(), format!("{coord}"))
    }

    #[test]
    fn test_hapvariant_genotype_getter() {
        let mut hv = HapVariant {
            contig: "chr9".to_string(),
            pos: 25,
            reference: "G".to_string(),
            alt: "T".to_string(),
            gt: 1,
        };

        assert_eq!(&"T".to_string(), hv.genotype());

        hv.gt = 0;
        assert_eq!(&"G".to_string(), hv.genotype());
    }
}
