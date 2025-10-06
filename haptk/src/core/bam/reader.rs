use std::ops::Range;
use std::path::Path;

use rust_htslib::bam::{IndexedReader, Read};

use color_eyre::eyre::WrapErr;
use color_eyre::Result;

use crate::error::Error;

use crate::core::bam::Record;

pub struct Reader {
    rdr: IndexedReader,
    contig: String,
    sample_id: String,
}

impl Reader {
    pub fn from_path(
        bam_path: &Path,
        ref_path: &Path,
        contig: &str,
        range: &Range<u32>,
        sample_id: &str,
    ) -> Result<Self> {
        // Fetch reads in a given contig from the bam
        let mut rdr = IndexedReader::from_path(bam_path).wrap_err(Error::Io {
            path: bam_path.to_path_buf(),
        })?;
        rdr.set_reference(ref_path)?;

        rdr.fetch((&contig, range.start, range.end))?;
        tracing::info!("Fetched {}..{} for {}", range.start, range.end, sample_id);

        Ok(Self {
            rdr,
            contig: contig.to_string(),
            sample_id: sample_id.to_string(),
        })
    }

    pub fn records(&mut self) -> impl Iterator<Item = Record> + '_ {
        self.rdr.records().filter_map(|v| match v {
            Ok(record) => Some(Record::new(record)),
            Err(_) => {
                tracing::warn!(
                    "Invalid record found for {} {}",
                    self.contig,
                    self.sample_id
                );
                None
            }
        })
    }
}
