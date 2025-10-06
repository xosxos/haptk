use std::path::PathBuf;

use thiserror::Error as ThisError;

#[rustfmt::skip]
#[derive(ThisError, Debug)]
pub enum Error {
    #[error("Failed to parse coords: {coord}")]
    CoordParse { coord: String },

    #[error("Position: {value:?} is not an integer in coords {coord}")]
    PosParse { coord: String, value: String },

    #[error("Htslib error: {0}")]
    HtsLib(#[from] rust_htslib::errors::Error),

    #[error("Io error: {0} {1}")]
    Io(PathBuf, std::io::Error),

    #[error("Io error: {0}")]
    StdIo(#[from] std::io::Error),

    #[error("parse error: {0}")]
    Parse(#[from] noodles::sam::header::record::value::map::tag::ParseError),

    #[error("Unknown file extension: {0:?}")]
    UnknownExtension(PathBuf),
}
