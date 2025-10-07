pub use htslib::Header;

mod htslib {
    use std::path::Path;

    use crate::Error;

    use rust_htslib::bcf;
    use rust_htslib::bcf::HeaderRecord;
    use rust_htslib::bcf::IndexedReader;
    use rust_htslib::bcf::Read;

    pub struct Header(bcf::header::HeaderView);

    impl Header {
        pub fn try_get(path: &Path) -> Result<Self, Error> {
            let bcf = IndexedReader::from_path(path)?;
            let header = bcf.header().clone();
            Ok(Self(header))
        }

        pub fn samples(&self) -> Result<Vec<String>, Error> {
            self.0
                .samples()
                .into_iter()
                .map(|sample| Ok(String::from_utf8_lossy(sample).to_string()))
                .collect()
        }

        pub fn reference_sequences(&self) -> Result<Vec<(String, u64)>, Error> {
            self.0
                .header_records()
                .iter()
                .filter(|r| matches!(r, HeaderRecord::Contig { .. }))
                .map(|r| match r {
                    HeaderRecord::Contig { values, .. } => {
                        let id = values.get("ID").ok_or_else(|| {
                            Error::New(format!(
                                "The input VCF has no ID for {values:?} in the header"
                            ))
                        })?;

                        Ok((
                            id.clone(),
                            values
                                .get("length")
                                .ok_or_else(|| Error::NoContigLength {
                                    contig: id.to_string(),
                                })?
                                .parse::<u64>()?,
                        ))
                    }
                    _ => {
                        unreachable!(
                            "
                    VCF header contig filtering failed for some reason.
                    Create an issue at Github: https://github.com/xosxos/haptk
                    "
                        )
                    }
                })
                .collect()
        }
    }
}

use crate::error::Error;
use crate::utils;
use std::collections::HashMap;

impl Header {
    #[allow(clippy::type_complexity)]
    pub fn filter_contigs(
        &self,
        contigs: Vec<String>,
    ) -> Result<
        (
            Vec<(String, Option<u64>, Option<u64>)>,
            HashMap<String, u64>,
        ),
        Error,
    > {
        // println!("{contigs:?}");

        let mut file_contigs_and_lengths: HashMap<String, u64> = HashMap::new();
        let mut file_contigs = vec![];

        for (seq_name, seq_len) in self.reference_sequences()? {
            file_contigs_and_lengths.insert(seq_name.clone(), seq_len);
            file_contigs.push(seq_name.to_string());
        }

        let mut contigs = contigs
            .into_iter()
            .map(utils::parse_coords)
            .inspect(|c| {
                if c.is_err() {
                    panic!("could not parse coord {c:?}")
                }
            })
            .map(|v| v.unwrap())
            .filter(|contig| file_contigs.contains(&contig.0))
            .collect::<Vec<(String, Option<u64>, Option<u64>)>>();

        if contigs.is_empty() {
            tracing::warn!(
                "None of the given contigs are in the bam file, using all contigs of the bam file"
            );
            contigs = file_contigs.into_iter().map(|c| (c, None, None)).collect();
        }

        Ok((contigs, file_contigs_and_lengths))
    }
}
