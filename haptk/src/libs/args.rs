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

#[derive(Serialize, Deserialize, Clone, Default, Debug, PartialEq, Eq)]
pub enum Selection {
    #[default]
    All,
    OnlyAlts,
    OnlyRefs,
    OnlyLongest,
    Unphased,
    Haploid,
}

impl AsRef<Selection> for Selection {
    fn as_ref(&self) -> &Selection {
        self
    }
}
