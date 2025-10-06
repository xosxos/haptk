use std::collections::HashMap;
use std::{fs::File, path::PathBuf};

use color_eyre::eyre::eyre;
use color_eyre::eyre::WrapErr;
use color_eyre::Result;

use crate::error::Error;
use crate::utils;

pub struct Header(noodles::sam::header::Header);

impl Header {
    pub fn try_get(path: &PathBuf) -> Result<Self> {
        let extension = path
            .extension()
            .unwrap_or_else(|| panic!("bam file {path:?} does not have a proper extension"))
            .to_str()
            .unwrap();

        let header: noodles::sam::Header = match extension {
            "cram" => File::open(path)
                .map(noodles::cram::io::Reader::new)
                .wrap_err(Error::Io { path: path.clone() })?
                .read_header()?,
            "bam" => File::open(path)
                .map(noodles::bam::io::Reader::new)
                .wrap_err(Error::Io { path: path.clone() })?
                .read_header()?,
            e => return Err(eyre!("unknown file extension {e} for {path:?}")),
        };

        Ok(Self(header))
    }

    pub fn read_group_samples(&self) -> Result<Vec<String>> {
        self.0
            .read_groups()
            .iter()
            .map(|(_, map)| {
                let sm =
                    noodles::sam::header::record::value::map::tag::Other::try_from([b'S', b'M'])?;

                let sample_name = map.other_fields().get(&sm).unwrap().to_string();
                Ok(sample_name)
            })
            .collect()
    }

    pub fn reference_sequences(&self) -> Vec<(String, u64)> {
        self.0
            .reference_sequences()
            .iter()
            .map(|(seqname, map)| (seqname.to_string(), map.length().get() as u64))
            .collect()
    }

    #[allow(clippy::type_complexity)]
    pub fn filter_contigs(
        &self,
        contigs: Vec<String>,
    ) -> Result<(
        Vec<(String, Option<u64>, Option<u64>)>,
        HashMap<String, u64>,
    )> {
        // println!("{contigs:?}");

        let mut bam_contigs_and_lengths: HashMap<String, u64> = HashMap::new();
        let mut bam_contigs = vec![];

        for (seq_name, seq_len) in self.reference_sequences() {
            bam_contigs_and_lengths.insert(seq_name.clone(), seq_len);
            bam_contigs.push(seq_name.to_string());
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
            .filter(|contig| bam_contigs.contains(&contig.0))
            .collect::<Vec<(String, Option<u64>, Option<u64>)>>();

        if contigs.is_empty() {
            tracing::warn!(
                "None of the given contigs are in the bam file, using all contigs of the bam file"
            );
            contigs = bam_contigs.into_iter().map(|c| (c, None, None)).collect();
        }

        Ok((contigs, bam_contigs_and_lengths))
    }
}
