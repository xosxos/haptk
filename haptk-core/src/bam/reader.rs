use std::ops::Range;
use std::path::Path;

use rust_htslib::bam::{IndexedReader, Read};
use std::result::Result;

use crate::bam;
use crate::error;
use crate::io::FileType;

pub struct Reader {
    rdr: IndexedReader,
    contig: String,
    sample_id: String,
}

impl Reader {
    pub fn from_path<P>(
        bam_path: &Path,
        ref_path: Option<P>,
        contig: &str,
        range: &Range<u64>,
        sample_id: &str,
    ) -> Result<Self, error::Error>
    where
        P: AsRef<Path>,
    {
        // Fetch reads in a given contig from the bam
        let mut rdr = IndexedReader::from_path(bam_path)?;

        let file_type = FileType::from_path(bam_path)?;

        if matches!(file_type, FileType::CRAM) {
            rdr.set_reference(
                ref_path
                    .ok_or_else(|| {
                        error::Error::New(format!("No reference file with CRAM file {bam_path:?}"))
                    })?
                    .as_ref(),
            )?;
        }

        rdr.fetch((&contig, range.start, range.end))?;
        tracing::info!("Fetched {}..{} for {}", range.start, range.end, sample_id);

        Ok(Self {
            rdr,
            contig: contig.to_string(),
            sample_id: sample_id.to_string(),
        })
    }

    pub fn records(&mut self) -> impl Iterator<Item = bam::Record> + '_ {
        self.rdr.records().filter_map(|v| match v {
            Ok(record) => Some(bam::Record::new(record)),
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
