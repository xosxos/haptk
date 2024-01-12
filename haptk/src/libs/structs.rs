use std::ops::Range;
use std::ops::RangeBounds;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use ndarray::{s, Array2, ArrayView1, Axis};
use polars::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::Selection;
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::bhst::Node;
use crate::subcommands::uhst;

#[derive(Debug, Default, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Coord {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
}

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

#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Eq, Serialize, Deserialize)]
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
        self.contig == other.contig
            && self.pos == other.pos
            && self.reference == other.reference
            && self.alt == other.alt
    }
}

impl PartialEq<Coord> for HapVariant {
    fn eq(&self, other: &Coord) -> bool {
        self.contig == other.contig
            && self.pos == other.pos
            && self.reference == other.reference
            && self.alt == other.alt
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
            Selection::Haploid => Ploidy::Haploid,
            _ => Ploidy::Diploid,
        }
    }
}

impl From<Selection> for Ploidy {
    fn from(selection: Selection) -> Self {
        match selection {
            Selection::Haploid => Ploidy::Haploid,
            _ => Ploidy::Diploid,
        }
    }
}

//TODO: Most of these fields, if not all, should be private to ensure the correctness of
// for example common selects and slicings performed by the user
// Slice indexing matrices is also an anti-pattern in Polars
#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhasedMatrix {
    pub variant_idx: usize,
    pub matrix: Array2<u8>,
    samples: Vec<String>,
    coords: Vec<Coord>,
    clinical_data: Option<DataFrame>,
    pub ploidy: Ploidy,
}

impl PhasedMatrix {
    pub fn new<T: AsRef<Selection>>(
        variant_idx: usize,
        matrix: Array2<u8>,
        samples: Vec<String>,
        coords: Vec<Coord>,
        selection: T,
    ) -> Self {
        Self {
            variant_idx,
            matrix,
            samples,
            coords,
            clinical_data: None,
            ploidy: selection.as_ref().into(),
        }
    }

    pub fn nsamples(&self) -> usize {
        self.samples.len()
    }

    pub fn nrows(&self) -> usize {
        self.matrix.nrows()
    }

    pub fn ncoords(&self) -> usize {
        self.coords.len()
    }

    pub fn samples(&self) -> &Vec<String> {
        &self.samples
    }

    pub fn clinical_data(&self) -> &Option<DataFrame> {
        &self.clinical_data
    }

    pub fn coords(&self) -> &Vec<Coord> {
        &self.coords
    }

    pub fn coords_mut(&mut self) -> &mut Vec<Coord> {
        &mut self.coords
    }

    pub fn set_coords(&mut self, coords: Vec<Coord>) {
        self.coords = coords;
    }

    pub fn variant_idx(&self) -> usize {
        self.variant_idx
    }

    pub fn set_variant_idx(&mut self, variant_idx: usize) {
        self.variant_idx = variant_idx;
    }

    pub fn variant_idx_pos(&self) -> u64 {
        if !self.coords.is_empty() {
            self.coords[self.variant_idx].pos
        } else {
            0
        }
    }

    pub fn get_coord(&self, idx: usize) -> &Coord {
        &self.coords[idx]
    }

    pub fn get_pos(&self, idx: usize) -> u64 {
        self.coords[idx].pos
    }

    pub fn get_contig(&self) -> &String {
        &self.coords[0].contig
    }

    pub fn get_first_idx_on_right_by_pos(&self, pos: u64) -> usize {
        match self.coords.iter().position(|c| c.pos >= pos) {
            Some(idx) => idx,
            None => self.coords.len(),
        }
    }

    pub fn get_first_idx_on_left_by_pos(&self, pos: u64) -> usize {
        match self.coords.iter().position(|c| c.pos > pos) {
            Some(idx) => idx.saturating_sub(1),
            None => 0,
        }
    }

    pub fn get_nearest_idx_by_pos(&self, pos: u64) -> usize {
        if let Some(idx) = self.coords.iter().position(|c| c.pos >= pos) {
            if idx == 0 {
                return idx;
            }
            if self.get_pos(idx) == pos {
                return idx;
            }

            let one_before_idx = idx - 1;
            let before = self.get_pos(one_before_idx);
            let after = self.get_pos(idx);
            if pos - before >= after - pos {
                idx
            } else {
                one_before_idx
            }
        } else {
            // return last element if no element is larger
            self.coords.len() - 1
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

    pub fn matrix(&self) -> &Array2<u8> {
        &self.matrix
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.samples.get(index / *self.ploidy).unwrap().clone()
        // self.samples.get(index).unwrap().clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_sample_idx(&self, sample: &str) -> Result<usize> {
        self.samples
            .iter()
            .position(|s| s == sample)
            .ok_or_else(|| eyre!("sample {sample} not found in the vcf"))
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

    pub fn select_carriers(&mut self, variant_pos: u64, selection: &Selection) -> Result<()> {
        ensure!(
            self.variant_idx_pos() == variant_pos,
            "Cannot select variant carriers: no variant at position {variant_pos}"
        );

        let slice = self
            .matrix
            .slice(s![0..self.matrix.nrows(), self.variant_idx]);

        let indexes = slice
            .iter()
            .enumerate()
            .filter(|(_, gt)| match selection {
                Selection::OnlyAlts => **gt == 1,
                Selection::OnlyRefs => **gt == 0,
                _ => panic!("Invalid selection method for alleles"),
            })
            .map(|(i, _)| i)
            .collect();

        self.select_rows(indexes);
        self.ploidy = Ploidy::Mixed;

        tracing::info!("Finished selecting by alleles.");
        Ok(())
    }

    pub fn get_lengths_from_uhst(&self) -> Vec<(Node, Node)> {
        let uhst_right = uhst::construct_uhst(
            self,
            &uhst::LocDirection::Right,
            self.variant_idx(),
            1,
            true,
        );
        let uhst_left =
            uhst::construct_uhst(self, &uhst::LocDirection::Left, self.variant_idx(), 1, true);

        let lmaj_branch = find_majority_nodes(&uhst_left);
        let rmaj_branch = find_majority_nodes(&uhst_right);

        (0..self.matrix.nrows())
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

    pub fn select_only_longest(&mut self) {
        let longest_indexes = self.only_longest_indexes();

        self.select_rows(longest_indexes);
        self.ploidy = Ploidy::Haploid;

        tracing::info!("Finished only-longest selection.");
    }

    pub fn only_longest_indexes(&self) -> Vec<usize> {
        let lengths = self.get_lengths_from_uhst();
        let lengths = lengths
            .iter()
            .map(|(lnode, rnode)| {
                if lnode.start_idx == rnode.stop_idx {
                    0
                } else {
                    self.get_pos(rnode.stop_idx - 1) - self.get_pos(lnode.start_idx + 1) + 1
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
                if lengths.filter(|(_, l)| *l == max_len).count() > 1 {
                    tracing::warn!(
                        "Sample {} has two equally long haplotypes in only-longest selection.",
                        self.get_sample_name(i)
                    );
                }
                max_idx
            })
            .collect::<Vec<usize>>()
    }

    // Sort inside select_rows to avoid bugs down the line
    pub fn select_rows(&mut self, mut to_keep: Vec<usize>) {
        to_keep.sort();
        self.matrix = self.matrix.select(Axis(0), &to_keep);

        self.samples = to_keep
            .iter()
            .map(|index| self.get_sample_name(*index))
            .collect();
    }

    pub fn select_columns_by_idx(&mut self, to_keep: &mut [usize]) {
        to_keep.sort();
        self.matrix = self.matrix.select(Axis(1), to_keep);

        self.coords = to_keep
            .iter()
            .map(|index| self.coords[*index].clone())
            .collect();
    }

    pub fn select_columns_by_range(&mut self, range: Range<usize>) {
        // Take the first element of the range and subtract from variant_idx
        self.variant_idx -= range.clone().next().unwrap();
        self.matrix = self.matrix.slice(s![.., range.clone()]).to_owned();
        self.coords = self.coords[range].to_vec();
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_haplotype_for_sample<R: RangeBounds<usize> + Iterator<Item = usize>>(
        &self,
        contig: &str,
        range: R,
        sample: usize,
    ) -> Vec<HapVariant> {
        range
            .map(|index| HapVariant {
                contig: contig.to_string(),
                pos: self.coords[index].pos,
                alt: self.coords[index].alt.clone(),
                reference: self.coords[index].reference.clone(),
                gt: *self.matrix.slice(ndarray::s![sample, index]).into_scalar(),
            })
            .collect()
    }

    // Inclusive ranges not supported so remember to add + 1 to stop_idx
    pub fn find_u8_haplotype_for_sample<R: RangeBounds<usize> + Iterator<Item = usize>>(
        &self,
        range: R,
        sample: usize,
    ) -> Vec<u8> {
        range
            .map(|index| *self.matrix.slice(ndarray::s![sample, index]).into_scalar())
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
    fn get_slot(&self, index: usize) -> ArrayView1<u8>;
    fn is_contradictory(&self, index: usize, positions: &Vec<usize>) -> bool;
    fn prev_contradictory(&self, index: usize, positions: &Vec<usize>) -> Option<usize>;
    fn next_contradictory(&self, index: usize, positions: &Vec<usize>) -> Option<usize>;
}

impl CoordDataSlot for PhasedMatrix {
    fn get_slot(&self, index: usize) -> ArrayView1<u8> {
        self.matrix.index_axis(Axis(1), index)
    }
    // Slot is contractory if it contains both 0 and 1.
    fn is_contradictory(&self, index: usize, positions: &Vec<usize>) -> bool {
        let slot = self.get_slot(index);
        let mut iter = positions.iter();
        let first = iter.next().unwrap();
        slot.len() > 1 && iter.any(|x| slot[*x] != slot[*first])
    }

    fn prev_contradictory(&self, index: usize, positions: &Vec<usize>) -> Option<usize> {
        if index == 0 || positions.len() < 2 {
            return None;
        }
        let mut idx = (index - 1) as isize;
        while idx >= 0 {
            if self.is_contradictory(idx as usize, positions) {
                return Some(idx as usize);
            }
            idx -= 1;
        }
        None
    }

    fn next_contradictory(&self, index: usize, positions: &Vec<usize>) -> Option<usize> {
        if index == self.matrix.ncols() - 1 || positions.len() < 2 {
            return None;
        }
        let mut idx = index + 1;
        while idx < self.matrix.ncols() {
            if self.is_contradictory(idx as usize, positions) {
                return Some(idx);
            }
            idx += 1;
        }
        None
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
