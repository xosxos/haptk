use std::cmp::Ordering;
use std::collections::btree_set::Range;
use std::collections::BTreeSet;
use std::ops::Bound;
// use std::ops::Range;
use std::ops::RangeBounds;
use std::ops::RangeInclusive;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use ndarray::{s, Array2, ArrayView1, Axis};
use petgraph::graph::NodeIndex;
use polars::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::Selection;
use crate::structs::Coord;
use crate::structs::HapVariant;
use crate::structs::Ploidy;

use crate::subcommands::bhst_shard::{find_majority_nodes, Node};
use crate::subcommands::uhst_shard;

// //TODO: Most of these fields, if not all, should be private to ensure the correctness of
// // for example common selects and slicings performed by the user
// // Slice indexing matrices is also an anti-pattern in Polars
// #[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
// pub struct PhasedMatrix {
//     pub start_coord: Coord,
//     pub variant_idx: usize,
//     pub matrix: Array2<u8>,
//     samples: Vec<String>,
//     coords: BTreeSet<Coord>,
//     clinical_data: Option<DataFrame>,
//     pub ploidy: Ploidy,
// }

// impl PhasedMatrix {
//     pub fn new<T: AsRef<Selection>>(
//         variant_idx: usize,
//         start_coord: Coord,
//         matrix: Array2<u8>,
//         samples: Vec<String>,
//         coords: BTreeSet<Coord>,
//         selection: T,
//     ) -> Self {
//         Self {
//             variant_idx,
//             start_coord,
//             matrix,
//             samples,
//             coords,
//             clinical_data: None,
//             ploidy: selection.as_ref().into(),
//         }
//     }

//     pub fn nsamples(&self) -> usize {
//         self.samples.len()
//     }

//     pub fn nhaplotypes(&self) -> usize {
//         self.samples.len() * *self.ploidy
//     }

//     pub fn nrows(&self) -> usize {
//         self.matrix.nrows()
//     }

//     pub fn ncoords(&self) -> usize {
//         self.coords.len()
//     }

//     pub fn samples(&self) -> &Vec<String> {
//         &self.samples
//     }

//     pub fn clinical_data(&self) -> &Option<DataFrame> {
//         &self.clinical_data
//     }

//     pub fn variant_idx(&self) -> usize {
//         self.variant_idx
//     }

//     pub fn set_variant_idx(&mut self, variant_idx: usize) {
//         self.variant_idx = variant_idx;
//     }

//     pub fn variant_idx_pos(&self) -> u64 {
//         if !self.coords.is_empty() {
//             self.start_coord.pos
//         } else {
//             0
//         }
//     }

//     pub fn coords(&self) -> &BTreeSet<Coord> {
//         &self.coords
//     }

//     pub fn coords_mut(&mut self) -> &mut BTreeSet<Coord> {
//         &mut self.coords
//     }

//     pub fn set_coords(&mut self, coords: BTreeSet<Coord>) {
//         self.coords = coords;
//     }

//     pub fn start_coord(&self) -> Coord {
//         // self.variant_idx
//         self.start_coord.clone()
//     }

//     pub fn get_coord_idx(&self, coord: &Coord) -> usize {
//         self.coords.iter().position(|c| c == coord).unwrap()
//     }

//     pub fn set_start_coord(&mut self, coord: Coord) {
//         self.start_coord = coord;
//     }

//     pub fn get_contig(&self) -> &String {
//         &self.coords.get(&self.start_coord).unwrap().contig
//     }

//     pub fn get_first_idx_on_right_by_pos(&self, pos: u64) -> usize {
//         match self.coords.iter().position(|c| c.pos >= pos) {
//             Some(idx) => idx,
//             None => self.coords.len(),
//         }
//     }

//     pub fn get_first_idx_on_left_by_pos(&self, pos: u64) -> usize {
//         match self.coords.iter().position(|c| c.pos > pos) {
//             Some(idx) => idx.saturating_sub(1),
//             None => 0,
//         }
//     }

//     pub fn get_nearest_coord_by_pos(&self, pos: u64) -> &Coord {
//         if let Some(idx) = self.coords.iter().position(|c| c.pos >= pos) {
//             let coord = self.coords.iter().nth(idx).unwrap();

//             if idx == 0 || coord.pos == pos {
//                 return coord;
//             }

//             let before = self.coords.iter().nth(idx - 1).unwrap();

//             match pos - before.pos >= coord.pos - pos {
//                 true => coord,
//                 false => before,
//             }
//         } else {
//             // return last element if no element is larger
//             self.coords.last().unwrap()
//         }
//     }

//     pub fn idx_by_pos(&self, pos: u64) -> Option<usize> {
//         self.coords.iter().position(|c| c.pos == pos)
//     }

//     pub fn idx_by_coord(&self, coord: &Coord) -> Option<usize> {
//         self.coords.iter().position(|c| c == coord)
//     }

//     pub fn idx_by_hapvariant(&self, hap: &HapVariant) -> Option<usize> {
//         self.coords.iter().position(|coord| coord == hap)
//     }

//     pub fn matrix(&self) -> &Array2<u8> {
//         &self.matrix
//     }

//     pub fn get_sample_name(&self, index: usize) -> String {
//         self.samples.get(index / *self.ploidy).unwrap().clone()
//         // self.samples.get(index).unwrap().clone()
//     }

//     pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
//         indexes.iter().map(|i| self.get_sample_name(*i)).collect()
//     }

//     pub fn get_sample_idxs(&self, samples: &[String]) -> Result<Vec<usize>> {
//         let idxs: Vec<_> = self
//             .samples
//             .iter()
//             .enumerate()
//             .filter(|(_, s)| samples.contains(s))
//             .flat_map(|(i, _)| {
//                 ((i * *self.ploidy)..(i * *self.ploidy) + *self.ploidy).collect::<Vec<usize>>()
//             })
//             .collect();

//         ensure!(
//             !idxs.is_empty(),
//             "None of the control samples are found in the vcf."
//         );
//         Ok(idxs)
//     }

//     pub fn select_carriers(&mut self, variant_pos: u64, selection: &Selection) -> Result<()> {
//         ensure!(
//             self.variant_idx_pos() == variant_pos,
//             "Cannot select variant carriers: no variant at position {variant_pos}"
//         );

//         let coord = self.get_nearest_coord_by_pos(variant_pos);
//         let coord_idx = self.get_coord_idx(coord);

//         let slice = self.matrix.slice(s![0..self.matrix.nrows(), coord_idx]);

//         let indexes = slice
//             .iter()
//             .enumerate()
//             .filter(|(_, gt)| match selection {
//                 Selection::OnlyAlts => **gt == 1,
//                 Selection::OnlyRefs => **gt == 0,
//                 _ => panic!("Invalid selection method for alleles"),
//             })
//             .map(|(i, _)| i)
//             .collect::<Vec<usize>>();

//         match selection {
//             Selection::OnlyAlts => tracing::info!("Selected {} ALT alleles", indexes.len()),
//             Selection::OnlyRefs => tracing::info!("Selected {} REF alleles", indexes.len()),
//             _ => panic!("Invalid selection method for alleles"),
//         };

//         self.select_rows(indexes);
//         self.ploidy = Ploidy::Mixed;

//         tracing::info!("Finished selecting by alleles.");
//         Ok(())
//     }

//     pub fn get_lengths_from_uhst(&self, start_coord: &Coord) -> Vec<(Node, Node)> {
//         let uhst_right = uhst_shard::construct_uhst(
//             self,
//             &uhst_shard::LocDirection::Right,
//             start_coord.clone(),
//             1,
//             true,
//         );
//         let uhst_left = uhst_shard::construct_uhst(
//             self,
//             &uhst_shard::LocDirection::Left,
//             start_coord.clone(),
//             1,
//             true,
//         );

//         let start_idx = NodeIndex::new(0);
//         let lmaj_branch = find_majority_nodes(&uhst_left, start_idx);
//         let rmaj_branch = find_majority_nodes(&uhst_right, start_idx);

//         (0..self.matrix.nrows())
//             .map(|idx| {
//                 let mut lnode = lmaj_branch.len()
//                     - lmaj_branch
//                         .iter()
//                         .rev()
//                         .position(|(node, _)| node.indexes.contains(&idx))
//                         .unwrap();
//                 let mut rnode = rmaj_branch.len()
//                     - rmaj_branch
//                         .iter()
//                         .rev()
//                         .position(|(node, _)| node.indexes.contains(&idx))
//                         .unwrap();
//                 if rnode == rmaj_branch.len() {
//                     rnode -= 1;
//                 }

//                 if lnode == lmaj_branch.len() {
//                     lnode -= 1;
//                 }

//                 (
//                     idx,
//                     lmaj_branch.get(lnode).unwrap(),
//                     rmaj_branch.get(rnode).unwrap(),
//                 )
//             })
//             .map(|(_idx, (lnode, _), (rnode, _))| ((*lnode).clone(), (*rnode).clone()))
//             .collect()
//     }

//     pub fn select_only_longest(&mut self) {
//         let longest_indexes = self.only_longest_indexes();

//         self.select_rows(longest_indexes);
//         self.ploidy = Ploidy::Haploid;

//         tracing::info!("Finished only-longest selection.");
//     }

//     pub fn only_longest_indexes(&self) -> Vec<usize> {
//         let lengths = self.get_lengths_from_uhst(&self.start_coord());
//         let lengths = lengths
//             .iter()
//             .map(|(lnode, rnode)| {
//                 if lnode.start == rnode.stop {
//                     0
//                 } else {
//                     let stop_idx = self.get_coord_idx(&rnode.stop);
//                     let stop = self.coords.iter().nth(stop_idx - 1).unwrap();

//                     let start_idx = self.get_coord_idx(&lnode.start);
//                     let start = self.coords.iter().nth(start_idx + 1).unwrap();
//                     stop.pos - start.pos + 1
//                 }
//             })
//             .collect::<Vec<u64>>();

//         (0..self.nsamples())
//             .map(|i| {
//                 let lengths = lengths
//                     .iter()
//                     .enumerate()
//                     .skip(i * *self.ploidy)
//                     .take(*self.ploidy);

//                 let (max_idx, max_len) = lengths.clone().max_by_key(|(_, l)| *l).unwrap();
//                 if lengths.filter(|(_, l)| *l == max_len).count() > 1 {
//                     tracing::warn!(
//                         "Sample {} has two equally long haplotypes in only-longest selection.",
//                         self.get_sample_name(i)
//                     );
//                 }
//                 max_idx
//             })
//             .collect::<Vec<usize>>()
//     }

//     pub fn only_longest_lengths(&self, start_coord: &Coord) -> Vec<(Node, Node)> {
//         let lengths = self.get_lengths_from_uhst(start_coord);

//         let calculate_len = |(lnode, rnode): &(Node, Node)| {
//             if lnode.start == rnode.stop {
//                 0
//             } else {
//                 let stop_idx = self.get_coord_idx(&rnode.stop);
//                 let stop = self.coords.iter().nth(stop_idx - 1).unwrap();

//                 let start_idx = self.get_coord_idx(&lnode.start);
//                 let start = self.coords.iter().nth(start_idx + 1).unwrap();
//                 stop.pos - start.pos + 1
//             }
//         };

//         (0..self.nsamples())
//             .map(|i| {
//                 let lengths = lengths
//                     .iter()
//                     .enumerate()
//                     .skip(i * *self.ploidy)
//                     .take(*self.ploidy);

//                 let (_, max_nodes) = lengths
//                     .clone()
//                     .max_by_key(|(_, l)| calculate_len(l))
//                     .unwrap();

//                 let max_len = calculate_len(max_nodes);
//                 if lengths.filter(|(_, l)| calculate_len(l) == max_len).count() > 1 {
//                     tracing::warn!(
//                         "Sample {} has two equally long haplotypes.",
//                         self.get_sample_name(i)
//                     );
//                 }
//                 max_nodes.clone()
//             })
//             .collect::<Vec<(Node, Node)>>()
//     }

//     // Sort inside select_rows to avoid bugs down the line
//     pub fn select_rows(&mut self, mut to_keep: Vec<usize>) {
//         to_keep.sort();
//         self.matrix = self.matrix.select(Axis(0), &to_keep);

//         self.samples = to_keep
//             .iter()
//             .map(|index| self.get_sample_name(*index))
//             .collect();
//     }

//     pub fn select_columns_by_idx(&mut self, to_keep: &mut [usize]) {
//         to_keep.sort();
//         self.matrix = self.matrix.select(Axis(1), to_keep);

//         self.coords = to_keep
//             .iter()
//             .map(|index| self.coords.iter().nth(*index).unwrap().clone())
//             .collect();
//     }

//     pub fn select_columns_by_range<R: RangeBounds<Coord>>(&mut self, range: R) {
//         // Take the first element of the range and subtract from variant_idx
//         if let (Bound::Included(first), Bound::Excluded(last)) =
//             (range.start_bound(), range.end_bound())
//         {
//             let first_idx = self.get_coord_idx(first);
//             let last_idx = self.get_coord_idx(last);

//             self.start_coord = first.clone();

//             self.matrix = self.matrix.slice(s![.., first_idx..last_idx]).to_owned();

//             self.coords = self.coords.range(range).cloned().collect();
//         } else {
//             panic!("range problem")
//         }
//     }

//     // Inclusive ranges not supported so remember to add + 1 to stop_idx
//     pub fn find_haplotype_for_sample<R: RangeBounds<Coord>>(
//         &self,
//         range: R,
//         sample: usize,
//     ) -> Vec<HapVariant> {
//         self.coords()
//             .range(range)
//             .map(|coord| {
//                 let index = self.get_coord_idx(coord);
//                 let gt = *self.matrix.slice(ndarray::s![sample, index]).into_scalar();

//                 HapVariant {
//                     contig: coord.contig.to_string(),
//                     pos: coord.pos,
//                     alt: coord.alt.clone(),
//                     reference: coord.reference.clone(),
//                     gt,
//                 }
//             })
//             .collect()
//     }

//     // Inclusive ranges not supported so remember to add + 1 to stop_idx
//     pub fn find_u8_haplotype_for_sample<R: RangeBounds<Coord>>(
//         &self,
//         range: R,
//         sample: usize,
//     ) -> Vec<u8> {
//         self.coords()
//             .range(range)
//             .map(|coord| {
//                 let index = self.get_coord_idx(coord);

//                 *self.matrix.slice(ndarray::s![sample, index]).into_scalar()
//             })
//             .collect()
//     }

//     pub fn set_variable_data(&mut self, df: DataFrame) -> Result<()> {
//         let clinical_samples: Vec<String> = if let Ok(raw_utf8) = df["id"].utf8() {
//             raw_utf8
//                 .into_iter()
//                 .flatten()
//                 .map(|v| v.to_string())
//                 .collect()
//         } else if let Ok(i64) = df["id"].i64() {
//             i64.into_iter().flatten().map(|v| v.to_string()).collect()
//         } else {
//             return Err(eyre!("Make sure variable file has an id column"));
//         };

//         let mut counter = 0;
//         for sample in &self.samples {
//             if !clinical_samples.contains(sample) {
//                 tracing::warn!("Sample {sample} is not present in clinical data");
//             } else {
//                 counter += 1;
//             }
//         }
//         if counter == 0 {
//             return Err(eyre!(
//                 "None of the wanted vcf samples match sample ids in the variable data file"
//             ));
//         }
//         let ids = Series::new("parsed_id", clinical_samples);
//         let df_lazy = df.lazy().with_columns([ids.lit()]).collect().unwrap();

//         let names = lit(Series::from_iter(self.samples.clone()));

//         let df = df_lazy
//             // .clone()
//             .lazy()
//             .filter(col("parsed_id").is_in(names))
//             .collect()?;

//         df["parsed_id"].utf8().wrap_err(eyre!("Error processing the `id` column of the variables data file. Either the column does not exist or id's are not utf8 compliant"))?;

//         self.clinical_data = Some(df);
//         Ok(())
//     }

//     pub fn get_variable_data_mean(
//         &self,
//         indexes: &[usize],
//         column_names: &[String],
//     ) -> Result<Option<Vec<f64>>> {
//         let samples = self.get_sample_names(indexes);
//         let names = lit(Series::from_iter(samples.clone()));

//         if let Some(df) = &self.clinical_data {
//             let mut means = Vec::with_capacity(column_names.len());
//             for column_name in column_names {
//                 tracing::debug!("Analyzing variable column {column_name:?}..");
//                 let check_that_not_empty = df
//                     .clone()
//                     .lazy()
//                     .filter(col("parsed_id").is_in(names.clone()))
//                     .select([col(column_name)])
//                     .collect()?;

//                 // Reminder: This check might need to be evaluated in the future to speed things up
//                 // It is here, because the mean function returns 0.0 for an empty vector
//                 // This can be considered a bug in under circumstances
//                 if check_that_not_empty.shape().0 == 0 {
//                     return Err(eyre!(
//                         "Samples {samples:?} have no data for `{column_name}`"
//                     ));
//                 }

//                 let mean = df
//                     .clone()
//                     .lazy()
//                     .filter(col("parsed_id").is_in(names.clone()))
//                     .select([col(column_name).mean()])
//                     .collect()?;

//                 let mean = mean.column(column_name)?.sum::<f64>().ok_or(eyre!(
//                     "A value {mean:?} cannot be parsed as an f64 in the column {column_name}."
//                 ))?;
//                 means.push(mean)
//             }
//             Ok(Some(means))
//         } else {
//             Ok(None)
//         }
//     }

//     pub fn get_variable_data_vecs(
//         &self,
//         indexes: &[usize],
//         column_names: &[String],
//     ) -> Result<Option<Vec<Vec<f64>>>> {
//         let names_vec = self.get_sample_names(indexes);
//         let names = lit(Series::from_iter(names_vec));

//         if let Some(df) = &self.clinical_data {
//             let mut vecs = Vec::with_capacity(column_names.len());
//             for column_name in column_names {
//                 let df = df
//                     .clone()
//                     .lazy()
//                     .filter(col("parsed_id").is_in(names.clone()))
//                     .select([col(column_name)])
//                     .collect()?;

//                 if let Ok(vec) = df[column_name.as_str()].f64() {
//                     let vec: Vec<f64> = vec.into_iter().flatten().collect();
//                     vecs.push(vec)
//                 } else if let Ok(vec) = df[column_name.as_str()].i64() {
//                     let vec: Vec<f64> = vec.into_iter().flatten().map(|v| v as f64).collect();
//                     vecs.push(vec)
//                 }
//             }
//             Ok(Some(vecs))
//         } else {
//             Ok(None)
//         }
//     }
// }

// pub trait CoordDataSlot {
//     fn get_slot(&self, coord: &Coord) -> ArrayView1<u8>;
//     fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool;
//     fn prev_contradictory(&self, coord: &Coord, positions: &[usize]) -> Option<&Coord>;
//     fn next_contradictory(&self, coord: &Coord, positions: &[usize]) -> Option<&Coord>;
// }

// impl CoordDataSlot for PhasedMatrix {
//     fn get_slot(&self, coord: &Coord) -> ArrayView1<u8> {
//         let index = self.get_coord_idx(coord);
//         self.matrix.index_axis(Axis(1), index)
//     }
//     // Slot is contractory if it contains both 0 and 1.
//     fn is_contradictory(&self, coord: &Coord, positions: &[usize]) -> bool {
//         let slot = self.get_slot(coord);
//         let mut iter = positions.iter();
//         let first = iter.next().unwrap();
//         slot.len() > 1 && iter.any(|x| slot[*x] != slot[*first])
//     }

//     fn prev_contradictory(&self, coord: &Coord, positions: &[usize]) -> Option<&Coord> {
//         let index = self.get_coord_idx(coord);

//         if index == 0 || positions.len() < 2 {
//             return None;
//         }

//         let mut idx = index as isize - 1;

//         while idx >= 0 {
//             let coord = self.coords.iter().nth(idx as usize).unwrap();
//             if self.is_contradictory(coord, positions) {
//                 return Some(coord);
//             }
//             idx -= 1;
//         }
//         None
//     }

//     fn next_contradictory(&self, coord: &Coord, positions: &[usize]) -> Option<&Coord> {
//         let index = self.get_coord_idx(coord);

//         if index == self.matrix.ncols() - 1 || positions.len() < 2 {
//             return None;
//         }

//         let mut idx = index + 1;

//         while idx < self.matrix.ncols() {
//             let coord = self.coords.iter().nth(idx).unwrap();
//             if self.is_contradictory(coord, positions) {
//                 return Some(coord);
//             }
//             idx += 1;
//         }
//         None
//     }
// }

// #[cfg(test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test_coords_displays() {
//         let coord = Coord {
//             contig: "chr9".to_string(),
//             pos: 25,
//             reference: "G".to_string(),
//             alt: "T".to_string(),
//         };

//         assert_eq!("chr9_25_G_T".to_string(), format!("{coord}"))
//     }

//     #[test]
//     fn test_hapvariant_genotype_getter() {
//         let mut hv = HapVariant {
//             contig: "chr9".to_string(),
//             pos: 25,
//             reference: "G".to_string(),
//             alt: "T".to_string(),
//             gt: 1,
//         };

//         assert_eq!(&"T".to_string(), hv.genotype());

//         hv.gt = 0;
//         assert_eq!(&"G".to_string(), hv.genotype());
//     }
// }
