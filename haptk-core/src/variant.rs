use std::cmp::Ordering;

use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Coord {
    pub contig: String,
    pub pos: u64,
    #[serde(rename(serialize = "ref"), alias = "ref")]
    pub reference: String,
    pub alt: String,
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct HapVariant {
    pub contig: String,
    pub pos: u64,
    pub reference: String,
    pub alt: String,
    pub gt: u8,
}

#[derive(Debug, Default, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CigarVariant {
    pub pos: u64,
    pub alt: char,
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

impl PartialEq for Coord {
    fn eq(&self, other: &Self) -> bool {
        (&self.contig, self.pos, &self.reference, &self.alt)
            == (&other.contig, other.pos, &other.reference, &other.alt)
    }
}

impl PartialEq<HapVariant> for Coord {
    fn eq(&self, other: &HapVariant) -> bool {
        if other.alt == "-" && other.gt == 0 {
            return (&self.contig, self.pos, &self.reference)
                == (&other.contig, other.pos, &other.reference);
        }

        (&self.contig, self.pos, &self.reference, &self.alt)
            == (&other.contig, other.pos, &other.reference, &other.alt)
    }
}

impl PartialEq<Coord> for HapVariant {
    fn eq(&self, other: &Coord) -> bool {
        if self.alt == "-" && self.gt == 0 {
            return (&self.contig, self.pos, &self.reference)
                == (&other.contig, other.pos, &other.reference);
        }

        (&self.contig, self.pos, &self.reference, &self.alt)
            == (&other.contig, other.pos, &other.reference, &other.alt)
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

impl Ord for Coord {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.pos < other.pos {
            return Ordering::Less;
        }

        if self.pos > other.pos {
            return Ordering::Greater;
        }

        let cmp = self.reference.cmp(&other.reference);

        if cmp == Ordering::Equal {
            if self.alt == "-" || other.alt == "-" {
                return Ordering::Equal;
            } else {
                return self.alt.cmp(&other.alt);
            }
        }

        cmp
    }
}

impl Eq for Coord {}

impl std::hash::Hash for Coord {
    fn hash<H>(&self, state: &mut H)
    where
        H: std::hash::Hasher,
    {
        self.pos.hash(state);
        self.alt.hash(state);
    }
}

impl PartialOrd for Coord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialOrd for CigarVariant {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
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

impl std::fmt::Display for HapVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!(
            "{}_{}_{}_{}_{}",
            self.contig, self.pos, self.reference, self.alt, self.gt
        );
        write!(f, "{line}")
    }
}

impl std::fmt::Display for CigarVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let line = format!("{}_{}", self.pos, self.alt);
        write!(f, "{line}")
    }
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
