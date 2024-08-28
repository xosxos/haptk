use std::path::PathBuf;

use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, PartialEq)]
pub struct StandardArgs {
    pub file: PathBuf,
    pub coords: String,
    pub output: PathBuf,
    pub samples: Option<Vec<PathBuf>>,
    pub selection: Selection,
    pub info_limit: Option<f32>,
    pub prefix: Option<String>,
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct ConciseArgs {
    pub file: PathBuf,
    pub output: PathBuf,
    pub prefix: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GraphArgs {
    pub width: f32,
    pub height: f32,
    pub mark_locus: bool,
    pub font_size: f32,
    pub stroke_width: u32,
    pub color: String,
    pub background_color: String,
}

impl Default for GraphArgs {
    fn default() -> Self {
        Self {
            width: 2560.0,
            height: 1440.0,
            mark_locus: false,
            font_size: 20.0,
            stroke_width: 5,
            color: String::from("black"),
            background_color: String::from("white"),
        }
    }
}

impl AsRef<Selection> for Selection {
    fn as_ref(&self) -> &Selection {
        self
    }
}
#[derive(Clone, Debug)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum SortOption {
    Left,
    Right,
    Total,
}

#[derive(Serialize, Deserialize, Clone, Default, Debug, PartialEq, Eq)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Selection {
    #[default]
    /// Select all alleles from samples (as of now only diploid or haploid organisms are supported)
    All,
    /// Select only the allele per each sample sharing the most haplotype with the haplotypes of the other samples
    OnlyLongest,
    /// Select only the alleles containing the REF variant at given a coordinate
    OnlyRefs,
    /// Select only the alleles containing the ALT variant at given a coordinate
    OnlyAlts,
    /// Use for unphased data (currently not supported for most commands)
    Unphased,
    /// Use for haploid genotypes
    Haploid,
}
