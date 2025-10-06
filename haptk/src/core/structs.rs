use std::cmp::Ordering;
use std::path::PathBuf;

use color_eyre::Result;
use ndarray::{Array2, ArrayView1, Axis};
use serde::{Deserialize, Serialize};

use crate::args::Selection;
use crate::error::Error;
use crate::read_vcf::read_shard_of_vcf;
use crate::subcommands::uhst::LocDirection;

use super::phased_matrix::PhasedMatrix;

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

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CigarVariant {
    pub pos: u64,
    pub alt: char,
}
impl PartialOrd for CigarVariant {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CigarVariant {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.pos < other.pos {
            return Ordering::Less;
        }

        if self.pos > other.pos {
            return Ordering::Greater;
        }

        self.alt.cmp(&other.alt)
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
