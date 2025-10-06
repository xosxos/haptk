use std::path::Path;

use std::result::Result;

use rust_htslib::bam::CompressionLevel;
use rust_htslib::bam::Format;
use rust_htslib::bam::Header;
use rust_htslib::bam::HeaderView;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::Read;

use crate::error;

pub struct Writer {
    writer: rust_htslib::bam::Writer,
    header_view: HeaderView,
}

impl Writer {
    pub fn from_path(
        output: &Path,
        input_path: &Path,
        ref_path: &Path,
    ) -> Result<Self, error::Error> {
        let bam = IndexedReader::from_path(input_path)?;

        let header_view = bam.header();
        let header = Header::from_template(header_view);

        let mut cram_writer = rust_htslib::bam::Writer::from_path(output, &header, Format::Cram)?;
        cram_writer.set_reference(ref_path)?;
        cram_writer.set_compression_level(CompressionLevel::Fastest)?;

        Ok(Self {
            writer: cram_writer,
            header_view: header_view.clone(),
        })
    }

    pub fn write(&mut self, record_string: String) -> Result<(), error::Error> {
        let record =
            rust_htslib::bam::Record::from_sam(&self.header_view, record_string.as_bytes())?;

        Ok(self.writer.write(&record)?)
    }
}
