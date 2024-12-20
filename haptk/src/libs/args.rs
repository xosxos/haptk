use std::path::PathBuf;

use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct StandardArgs {
    pub file: PathBuf,

    /// The starting coordinate, i.e. chr9:27573534
    #[arg(short = 'c', long)]
    pub coords: String,

    /// Output directory
    #[arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("./"))]
    pub output: PathBuf,

    /// List of samples for HST construction (one ID per row)
    #[arg(short = 'S', long, value_delimiter = ' ', num_args = 1.. )]
    pub samples: Option<Vec<PathBuf>>,

    #[arg(short = 'a', long = "alleles", value_enum, default_value_t = Selection::All)]
    pub selection: Selection,

    pub info_limit: Option<f32>,

    /// Output filename prefix
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,

    /// Do not include no ALTs
    #[arg(long)]
    pub no_alt: bool,
}

#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct ConciseArgs {
    /// Bilateral HST file obtained from bhst-scan
    pub file: PathBuf,

    /// Output directory
    #[arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("./"))]
    pub output: PathBuf,

    /// Output filename prefix
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,

    /// List of samples for HST construction (one ID per row)
    #[arg(short = 'S', long, value_delimiter = ' ', num_args = 1.. )]
    pub samples: Option<Vec<PathBuf>>,

    #[arg(short = 'a', long = "alleles", value_enum, default_value_t = Selection::All)]
    pub selection: Selection,
}

#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct GraphArgs {
    /// Graph width in px
    #[arg(long, default_value_t = 5000.0)]
    pub width: f32,

    /// Graph height in px
    #[arg(long, default_value_t = 7500.0)]
    pub height: f32,

    /// Mark the variant coordinate
    #[arg(long)]
    pub mark_locus: bool,

    // Font size
    #[arg(long, default_value_t = 75.0)]
    pub font_size: f32,

    // Line stroke width
    #[arg(long, default_value_t = 7)]
    pub stroke_width: u32,

    // Font color
    #[arg(long, default_value_t = String::from("black"))]
    pub color: String,

    // Background color
    #[arg(long, default_value_t = String::from("white"))]
    pub background_color: String,
}

impl Default for GraphArgs {
    fn default() -> Self {
        Self {
            width: 5000.0,
            height: 7500.0,
            mark_locus: false,
            font_size: 75.0,
            stroke_width: 7,
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
#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum SortOption {
    #[default]
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
    #[value(name = "longest-haplotype")]
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
