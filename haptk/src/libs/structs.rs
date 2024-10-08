use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::ops::Bound;
use std::ops::RangeBounds;
use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use ndarray::iter::AxisIter;
use ndarray::ArrayView2;
use ndarray::{s, Array2, ArrayView1, Axis};
use petgraph::graph::NodeIndex;
use polars::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::Selection;
use crate::error::HaptkError;
use crate::read_vcf::read_shard_of_vcf;
use crate::subcommands::bhst_shard::find_majority_nodes;
use crate::subcommands::bhst_shard::Node;
use crate::subcommands::uhst_shard;
use crate::subcommands::uhst_shard::LocDirection;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Coord {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
}

impl Ord for Coord {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.pos.cmp(&other.pos) {
            Ordering::Less => Ordering::Less,
            Ordering::Equal => {
                if self.alt == "-" || other.alt == "-" {
                    Ordering::Equal
                } else {
                    self.alt.cmp(&other.alt)
                }
            }
            Ordering::Greater => Ordering::Greater,
        }
        // self.pos.cmp(&other.pos)
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

pub type ClinicalData = (f64, f64);

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
            Selection::Haploid | Selection::OnlyRefs | Selection::OnlyAlts => Ploidy::Haploid,
            _ => Ploidy::Diploid,
        }
    }
}

impl From<Selection> for Ploidy {
    fn from(selection: Selection) -> Self {
        match selection {
            Selection::Haploid | Selection::OnlyRefs | Selection::OnlyAlts => Ploidy::Haploid,
            _ => Ploidy::Diploid,
        }
    }
}

pub enum MatrixSlice {
    All,
    Point(usize),
}

#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReadMetadata {
    pub file_path: PathBuf,
    pub fetch_range: (u64, u64),
    pub lookups: Vec<[bool; 2]>,
    pub indexes: Vec<usize>,
    pub contig_len: Option<u64>,
    pub sharded: bool,
}

//TODO: Most of these fields, if not all, should be private to ensure the correctness of
// for example common selects and slicings performed by the user
// Slice indexing matrices is also an anti-pattern in Polars
#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhasedMatrix {
    pub start_coord: Coord,
    pub variant_idx: usize,
    matrix: BTreeMap<Offset, Array2<u8>>,
    samples: Vec<String>,
    coords: BTreeSet<Coord>,
    clinical_data: Option<DataFrame>,
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
        let start = coords.first().unwrap().to_owned();
        let end = coords.last().unwrap().to_owned();

        let mut map = BTreeMap::new();
        map.insert((start, end), matrix);

        Self {
            variant_idx,
            start_coord,
            matrix: map,
            samples,
            coords,
            clinical_data: None,
            ploidy: selection.as_ref().into(),
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

    pub fn matrix_point(&self, x: usize, y: usize) -> u8 {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix[[x, y]]
    }

    pub fn matrix_point_coord(&self, x: usize, coord: &Coord) -> u8 {
        let (col_index, matrix) = self.find_matrix_and_col_idx(coord);
        matrix[[x, col_index]]
    }

    pub fn find_matrix_and_col_idx(&self, coord: &Coord) -> (usize, &Array2<u8>) {
        let matrix_idx = self
            .matrix
            .iter()
            .position(|((start, end), _)| coord >= start && end >= coord)
            .unwrap();

        let ((start, _), matrix) = self.matrix.iter().nth(matrix_idx).unwrap();

        let index = self.get_coord_idx(coord) - self.get_coord_idx(start);
        (index, matrix)
    }

    pub fn matrix_column(&self, coord: &Coord) -> ArrayView1<u8> {
        let (index, matrix) = self.find_matrix_and_col_idx(coord);

        matrix.index_axis(Axis(1), index)
    }

    // pub fn matrix_row(&self, coord: &Coord) -> ArrayView1<u8> {
    //     let index = self.get_coord_idx(coord);
    //     // Box::new(
    //     //     self.matrix
    //     //         .iter()
    //     //         .map(move |(_, m)| m.index_axis(Axis(1), index).into_iter().copied())
    //     //         .flatten(),
    //     // )
    //     let (_, matrix) = self.matrix.iter().nth(0).unwrap();
    // }

    pub fn matrix_slice(&self, row: MatrixSlice, col: MatrixSlice) -> ArrayView1<u8> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();

        let slice = match (row, col) {
            (MatrixSlice::Point(x), MatrixSlice::All) => matrix.slice(s![x, ..]),
            (MatrixSlice::All, MatrixSlice::Point(y)) => matrix.slice(s![.., y]),
            _ => panic!(),
        };
        slice
    }

    pub fn matrix_axis_iter(&self, axis: usize) -> AxisIter<'_, u8, ndarray::Dim<[usize; 1]>> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix.axis_iter(Axis(axis))
    }

    pub fn slice_cols<R: RangeBounds<usize>>(&self, col_range: R) -> ArrayView2<u8> {
        let (col_first_idx, col_last_idx) = match (col_range.start_bound(), col_range.end_bound()) {
            (Bound::Included(first), Bound::Excluded(last)) => (*first, last.saturating_sub(1)),
            (Bound::Included(first), Bound::Included(last)) => (*first, *last),
            _ => panic!("Range problem"),
        };

        let (_, matrix) = self.matrix.iter().nth(0).unwrap();

        matrix.slice(s![.., col_first_idx..=col_last_idx])
    }

    pub fn matrix(&self) -> &Array2<u8> {
        let (_, matrix) = self.matrix.iter().nth(0).unwrap();
        matrix
    }

    pub fn insert_matrix(&mut self, (start_coord, end_coord): (Coord, Coord), matrix: Array2<u8>) {
        self.matrix.insert((start_coord, end_coord), matrix);
    }

    pub fn set_matrix(&mut self, start_coord: Coord, end_coord: Coord, matrix: Array2<u8>) {
        let mut map = BTreeMap::new();

        map.insert((start_coord, end_coord), matrix);
        self.matrix = map;
    }

    pub fn matrix_select(&mut self, axis: usize, to_keep: &[usize]) {
        if axis == 0 {
            for (_, matrix) in self.matrix.iter_mut() {
                *matrix = matrix.select(Axis(0), to_keep)
            }
        } else {
            panic!("matrix select on columns not implemented")
        }
    }

    pub fn ncoords(&self) -> usize {
        self.coords.len()
    }

    pub fn samples(&self) -> &Vec<String> {
        &self.samples
    }

    pub fn has_samples(&self) -> bool {
        !self.samples.is_empty()
    }

    pub fn clinical_data(&self) -> &Option<DataFrame> {
        &self.clinical_data
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
        self.coords.iter().position(|c| c == coord).unwrap()
    }

    pub fn try_coord_by_hapvariant(&self, hv: &HapVariant) -> Option<usize> {
        self.coords.iter().position(|c| c == hv)
    }

    pub fn set_start_coord(&mut self, coord: Coord) {
        self.start_coord = coord;
    }

    pub fn get_contig(&self) -> &String {
        &self.coords.get(&self.start_coord).unwrap().contig
    }

    pub fn get_first_idx_on_right_by_pos(&self, pos: u64) -> usize {
        match self.coords.iter().position(|c| c.pos >= pos) {
            Some(idx) => idx,
            None => self.coords.len() - 1,
        }
    }

    pub fn get_first_idx_on_left_by_pos(&self, pos: u64) -> usize {
        match self.coords.iter().position(|c| c.pos > pos) {
            Some(idx) => idx.saturating_sub(1),
            None => 0,
        }
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

    pub fn idx_by_pos(&self, pos: u64) -> Option<usize> {
        self.coords.iter().position(|c| c.pos == pos)
    }

    pub fn idx_by_coord(&self, coord: &Coord) -> Option<usize> {
        self.coords.iter().position(|c| c == coord)
    }

    pub fn idx_by_hapvariant(&self, hap: &HapVariant) -> Option<usize> {
        self.coords.iter().position(|coord| coord == hap)
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.samples.get(index / *self.ploidy).unwrap().clone()
        // self.samples.get(index).unwrap().clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_sample_idxs(&self, samples: &[String]) -> Result<Vec<usize>> {
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

    // pub fn select_carriers(&mut self, variant_pos: u64, selection: &Selection) -> Result<()> {
    //     ensure!(
    //         self.variant_idx_pos() == variant_pos,
    //         "Cannot select variant carriers: no variant at position {variant_pos}"
    //     );

    //     let coord = self.get_nearest_coord_by_pos(variant_pos);
    //     let coord_idx = self.get_coord_idx(coord);

    //     let slice = self.matrix_slice(MatrixSlice::All, MatrixSlice::Point(coord_idx));

    //     let indexes = slice
    //         .iter()
    //         .enumerate()
    //         .filter(|(_, gt)| match selection {
    //             Selection::OnlyAlts => **gt == 1,
    //             Selection::OnlyRefs => **gt == 0,
    //             _ => panic!("Invalid selection method for alleles"),
    //         })
    //         .map(|(i, _)| i)
    //         .collect::<Vec<usize>>();

    //     match selection {
    //         Selection::OnlyAlts => tracing::info!("Selected {} ALT alleles", indexes.len()),
    //         Selection::OnlyRefs => tracing::info!("Selected {} REF alleles", indexes.len()),
    //         _ => panic!("Invalid selection method for alleles"),
    //     };

    //     self.select_rows(indexes);
    //     self.ploidy = Ploidy::Mixed;

    //     tracing::info!("Finished selecting by alleles.");
    //     Ok(())
    // }

    pub fn get_lengths_from_uhst(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let uhst_right = uhst_shard::construct_uhst(
            self,
            &uhst_shard::LocDirection::Right,
            start_coord,
            1,
            true,
        )?;
        let uhst_left = uhst_shard::construct_uhst(
            self,
            &uhst_shard::LocDirection::Left,
            start_coord,
            1,
            true,
        )?;

        let start_idx = NodeIndex::new(0);
        let lmaj_branch = find_majority_nodes(&uhst_left, start_idx);
        let rmaj_branch = find_majority_nodes(&uhst_right, start_idx);

        let lengths = (0..self.nhaplotypes())
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
            .collect();
        Ok(lengths)
    }

    pub fn select_only_longest(&mut self) -> Result<()> {
        let longest_indexes = self.only_longest_indexes()?;

        self.select_rows(longest_indexes);
        self.ploidy = Ploidy::Haploid;

        tracing::info!("Finished only-longest selection.");
        Ok(())
    }

    pub fn only_longest_indexes(&mut self) -> Result<Vec<usize>> {
        let start = self.start_coord().clone();
        let lengths = self.get_lengths_from_uhst(&start)?;
        let lengths = lengths
            .iter()
            .map(|(lnode, rnode)| {
                if lnode.start == rnode.stop {
                    0
                } else {
                    let stop_idx = self.get_coord_idx(&rnode.stop);
                    let stop = self.coords.iter().nth(stop_idx - 1).unwrap();

                    let start_idx = self.get_coord_idx(&lnode.start);
                    let start = self.coords.iter().nth(start_idx + 1).unwrap();
                    stop.pos - start.pos + 1
                }
            })
            .collect::<Vec<u64>>();

        let idxs = (0..self.nsamples())
            .map(|i| {
                let lengths = lengths
                    .iter()
                    .enumerate()
                    .skip(i * *self.ploidy)
                    .take(*self.ploidy);

                let (max_idx, max_len) = lengths.clone().max_by_key(|(_, l)| *l).unwrap();
                if lengths.filter(|(_, l)| *l == max_len).count() > 1 {
                    tracing::warn!(
                        "Sample {} has two equally long haplotypes in only-longest selection.",
                        self.get_sample_name(i)
                    );
                }
                max_idx
            })
            .collect::<Vec<usize>>();
        Ok(idxs)
    }

    pub fn only_longest_lengths(&mut self, start_coord: &Coord) -> Result<Vec<(Node, Node)>> {
        let lengths = self.get_lengths_from_uhst(start_coord)?;

        let calculate_len = |(lnode, rnode): &(Node, Node)| {
            if lnode.start == rnode.stop {
                0
            } else {
                let stop_idx = self.get_coord_idx(&rnode.stop);
                let stop = self.coords.iter().nth(stop_idx - 1).unwrap();

                let start_idx = self.get_coord_idx(&lnode.start);
                let start = self.coords.iter().nth(start_idx + 1).unwrap();
                stop.pos - start.pos + 1
            }
        };

        let lengths = (0..self.nsamples())
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
                if lengths.filter(|(_, l)| calculate_len(l) == max_len).count() > 1 {
                    tracing::warn!(
                        "Sample {} has two equally long haplotypes.",
                        self.get_sample_name(i)
                    );
                }
                max_nodes.clone()
            })
            .collect::<Vec<(Node, Node)>>();
        Ok(lengths)
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

    // pub fn select_columns_by_range<R: RangeBounds<Coord>>(&mut self, range: R) {
    //     let (first, last) = match (range.start_bound(), range.end_bound()) {
    //         (Bound::Included(first), Bound::Excluded(last)) => {
    //             let last_idx = self.get_coord_idx(last);
    //             (
    //                 first.clone(),
    //                 self.get_coord(last_idx.saturating_sub(1)).clone(),
    //             )
    //         }
    //         (Bound::Included(first), Bound::Included(last)) => (first.clone(), last.clone()),
    //         _ => panic!("Range problem"),
    //     };

    //     let first_idx = self.get_coord_idx(&first);
    //     let last_idx = self.get_coord_idx(&last);

    //     let last_idx = last_idx.min(self.ncoords() - 1);

    //     self.coords = self.coords.range(first..=last).cloned().collect();

    //     let offset_coord_start = self.coords.first().unwrap().clone();
    //     let offset_coord_end = self.coords.last().unwrap().clone();

    //     self.set_matrix(
    //         offset_coord_start,
    //         offset_coord_end,
    //         self.slice_cols(first_idx..last_idx).to_owned(),
    //     );

    //     let coord_idx = self.get_coord_idx(self.start_coord());
    //     self.start_coord = self.get_coord(coord_idx - first_idx).clone();
    // }

    pub fn select_columns_by_range_idx<R: RangeBounds<usize>>(&mut self, range: R) {
        let (first_idx, last_idx) = match (range.start_bound(), range.end_bound()) {
            (Bound::Included(first), Bound::Excluded(last)) => (*first, last.saturating_sub(1)),
            (Bound::Included(first), Bound::Included(last)) => (*first, *last),
            _ => panic!("Range problem"),
        };

        let first = self.get_coord(first_idx).clone();

        let last_idx = last_idx.min(self.ncoords() - 1);
        let last = self.get_coord(last_idx).clone();

        self.coords = self.coords.range(first..=last).cloned().collect();

        let offset_coord_start = self.coords.first().unwrap().clone();
        let offset_coord_end = self.coords.last().unwrap().clone();

        self.set_matrix(
            offset_coord_start,
            offset_coord_end,
            self.slice_cols(first_idx..last_idx).to_owned(),
        );

        let new_var_idx = self.variant_idx - first_idx;
        self.set_variant_idx(new_var_idx);
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
        range: R,
        sample: usize,
    ) -> Vec<u8> {
        self.coords()
            .range(range)
            .map(|coord| self.matrix_point_coord(sample, coord))
            .collect()
    }

    pub fn set_variable_data(&mut self, df: DataFrame) -> Result<()> {
        let clinical_samples: Vec<String> = if let Ok(raw_utf8) = df["id"].utf8() {
            raw_utf8
                .into_iter()
                .flatten()
                .map(|v| v.to_string())
                .collect()
        } else if let Ok(i64) = df["id"].i64() {
            i64.into_iter().flatten().map(|v| v.to_string()).collect()
        } else {
            return Err(eyre!("Make sure variable file has an id column"));
        };

        let mut counter = 0;
        for sample in &self.samples {
            if !clinical_samples.contains(sample) {
                tracing::warn!("Sample {sample} is not present in clinical data");
            } else {
                counter += 1;
            }
        }
        if counter == 0 {
            return Err(eyre!(
                "None of the wanted vcf samples match sample ids in the variable data file"
            ));
        }
        let ids = Series::new("parsed_id", clinical_samples);
        let df_lazy = df.lazy().with_columns([ids.lit()]).collect().unwrap();

        let names = lit(Series::from_iter(self.samples.clone()));

        let df = df_lazy
            // .clone()
            .lazy()
            .filter(col("parsed_id").is_in(names))
            .collect()?;

        df["parsed_id"].utf8().wrap_err(eyre!("Error processing the `id` column of the variables data file. Either the column does not exist or id's are not utf8 compliant"))?;

        self.clinical_data = Some(df);
        Ok(())
    }

    pub fn get_variable_data_mean(
        &self,
        indexes: &[usize],
        column_names: &[String],
    ) -> Result<Option<Vec<f64>>> {
        let samples = self.get_sample_names(indexes);
        let names = lit(Series::from_iter(samples.clone()));

        if let Some(df) = &self.clinical_data {
            let mut means = Vec::with_capacity(column_names.len());
            for column_name in column_names {
                tracing::debug!("Analyzing variable column {column_name:?}..");
                let check_that_not_empty = df
                    .clone()
                    .lazy()
                    .filter(col("parsed_id").is_in(names.clone()))
                    .select([col(column_name)])
                    .collect()?;

                // Reminder: This check might need to be evaluated in the future to speed things up
                // It is here, because the mean function returns 0.0 for an empty vector
                // This can be considered a bug in under circumstances
                if check_that_not_empty.shape().0 == 0 {
                    return Err(eyre!(
                        "Samples {samples:?} have no data for `{column_name}`"
                    ));
                }

                let mean = df
                    .clone()
                    .lazy()
                    .filter(col("parsed_id").is_in(names.clone()))
                    .select([col(column_name).mean()])
                    .collect()?;

                let mean = mean.column(column_name)?.sum::<f64>().ok_or(eyre!(
                    "A value {mean:?} cannot be parsed as an f64 in the column {column_name}."
                ))?;
                means.push(mean)
            }
            Ok(Some(means))
        } else {
            Ok(None)
        }
    }

    pub fn get_variable_data_vecs(
        &self,
        indexes: &[usize],
        column_names: &[String],
    ) -> Result<Option<Vec<Vec<f64>>>> {
        let names_vec = self.get_sample_names(indexes);
        let names = lit(Series::from_iter(names_vec));

        if let Some(df) = &self.clinical_data {
            let mut vecs = Vec::with_capacity(column_names.len());
            for column_name in column_names {
                let df = df
                    .clone()
                    .lazy()
                    .filter(col("parsed_id").is_in(names.clone()))
                    .select([col(column_name)])
                    .collect()?;

                if let Ok(vec) = df[column_name.as_str()].f64() {
                    let vec: Vec<f64> = vec.into_iter().flatten().collect();
                    vecs.push(vec)
                } else if let Ok(vec) = df[column_name.as_str()].i64() {
                    let vec: Vec<f64> = vec.into_iter().flatten().map(|v| v as f64).collect();
                    vecs.push(vec)
                }
            }
            Ok(Some(vecs))
        } else {
            Ok(None)
        }
    }
}

pub trait CoordDataSlot {
    fn read_more(&mut self, error: HaptkError) -> Result<()>;
    fn is_file_end(&self, side: LocDirection) -> bool;
    fn get_slot(&self, coord: &Coord) -> ArrayView1<u8>;
    fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool;

    fn prev_contradictory(
        &self,
        coord: &Coord,
        sample_idxs: &[usize],
    ) -> std::result::Result<Option<&Coord>, HaptkError>;

    fn next_contradictory(
        &self,
        coord: &Coord,
        sample_idxs: &[usize],
    ) -> std::result::Result<Option<&Coord>, HaptkError>;
}

impl CoordDataSlot for PhasedMatrix {
    fn get_slot(&self, coord: &Coord) -> ArrayView1<u8> {
        self.matrix_column(coord)
    }
    // Slot is contractory if it contains both 0 and 1.
    fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool {
        let slot = self.get_slot(coord);
        let mut iter = positions.iter();
        let first = iter.next().unwrap();
        slot.len() > 1 && iter.any(|x| slot[*x] != slot[*first])
    }

    fn prev_contradictory(
        &self,
        coord: &Coord,
        positions: &[usize],
    ) -> std::result::Result<Option<&Coord>, HaptkError> {
        if positions.len() < 2 {
            return Ok(None);
        }

        match self
            .coords
            .range(..coord)
            .rev()
            .find(|&coord| self.is_contradictory(coord, positions))
        {
            Some(coord) => Ok(Some(coord)),
            None => match self.is_file_end(LocDirection::Left) {
                true => Ok(None),
                false => Err(HaptkError::HstEndError),
            },
        }
    }

    fn next_contradictory(
        &self,
        coord: &Coord,
        positions: &[usize],
    ) -> std::result::Result<Option<&Coord>, HaptkError> {
        if positions.len() < 2 {
            return Ok(None);
        }

        match self
            .coords
            .range(coord..)
            .skip(1)
            .find(|&coord| self.is_contradictory(coord, positions))
        {
            Some(coord) => Ok(Some(coord)),
            None => match self.is_file_end(LocDirection::Right) {
                true => Ok(None),
                false => Err(HaptkError::HstEndError),
            },
        }
    }

    fn is_file_end(&self, side: LocDirection) -> bool {
        match side {
            LocDirection::Left => self.metadata.fetch_range.0 == 0 || !self.metadata.sharded,
            LocDirection::Right => {
                self.metadata.fetch_range.1 >= self.metadata.contig_len.unwrap()
                    || !self.metadata.sharded
            }
        }
    }

    fn read_more(&mut self, error: HaptkError) -> Result<()> {
        let start = self.coords.first().unwrap().clone();
        let stop = self.coords.last().unwrap().clone();

        match error {
            HaptkError::HstLeftEndError => {
                read_shard_of_vcf(self, start.pos.saturating_sub(1_000_000), start.pos)
            }
            HaptkError::HstRightEndError => read_shard_of_vcf(self, stop.pos, stop.pos + 1_000_000),
            HaptkError::HstBothEndError => {
                read_shard_of_vcf(self, start.pos.saturating_sub(1_000_000), start.pos)?;
                read_shard_of_vcf(self, stop.pos, stop.pos + 1_000_000)
            }
            _ => unreachable!(),
        }
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
