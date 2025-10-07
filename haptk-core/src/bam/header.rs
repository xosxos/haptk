pub use htslib::Header;

mod htslib {
    use std::path::PathBuf;

    use crate::Error;

    use rust_htslib::bam::{self, IndexedReader, Read};

    pub struct Header(bam::HeaderView);

    impl Header {
        pub fn try_get(path: &PathBuf) -> Result<Self, Error> {
            let bam = IndexedReader::from_path(path)?;
            let header = bam.header().clone();
            Ok(Self(header))
        }

        pub fn read_group_samples(&self) -> Result<Vec<String>, Error> {
            let bam_header = bam::Header::from_template(&self.0);
            let bam_hm = bam_header.to_hashmap();

            if !bam_hm.contains_key("RG") {
                return Ok(vec![]);
            }

            let rg = bam_hm.get("RG").unwrap();

            let samples = rg
                .iter()
                .filter_map(|item| item.get("SM"))
                .cloned()
                .collect();

            Ok(samples)
        }

        pub fn reference_sequences(&self) -> Vec<(String, u64)> {
            let bam_header = bam::Header::from_template(&self.0);
            let bam_hm = bam_header.to_hashmap();

            let sq = bam_hm.get("SQ").unwrap();

            sq.iter()
                .map(|seq| {
                    let name = seq.get("SN").unwrap();
                    let len = seq.get("LN").unwrap();
                    (name.clone(), len.parse::<u64>().unwrap())
                })
                .collect()
        }

        pub fn to_vcf_header(
            &self,
            input_header: Option<String>,
            info_rows: Vec<String>,
            program: &str,
            contigs: &[(String, Option<u64>, Option<u64>)],
            samples_and_paths: &[(String, PathBuf)],
            dia_line: String,
        ) -> Result<Vec<String>, Error> {
            let bam_header = bam::Header::from_template(&self.0);
            let bam_hm = bam_header.to_hashmap();
            let mut header = vec![];

            let time = chrono::offset::Utc::now();

            header.push("##fileformat=VCFv4.2".into());
            header.push(format!("##fileDate={time}"));
            header.push(format!("##source={program}"));
            header.push("##FILTER=<ID=PASS,Description=\"All filters passed\">".into());
            header.push(
                "##INFO=<ID=TYPE,Number=.,Type=String,Description=\"Insertion types\">".into(),
            );

            if let Some(rows) = input_header {
                // Filter out empty lines
                for elem in rows.split("\n").filter(|&x| !x.is_empty()) {
                    header.push(elem.to_string());
                }
            }
            header.push("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".into());
            header.push("##FORMAT=<ID=N,Number=1,Type=Integer,Description=\"Number of reads containing the insertion\">".into());
            header.push("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">".into());

            for row in info_rows {
                header.push(row);
            }

            if let Some(sq) = bam_hm.get("SQ") {
                for i in sq {
                    let sn = i.get("SN").unwrap();
                    let ln = i.get("LN").unwrap();
                    if contigs.iter().any(|c| sn == &c.0) {
                        header.push(format!("##contig=<ID={sn},length={ln}>"));
                    }
                }
            }

            if let Some(pg) = bam_hm.get("PG") {
                for i in pg {
                    let cl = i.get("CL");

                    if let Some(cl) = cl {
                        header.push(format!("##PG=\"{}\"", cl))
                    }
                }
            }

            header.push(dia_line);

            let last_header_row = format!(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
                samples_and_paths
                    .iter()
                    .map(|s| s.0.clone())
                    .collect::<Vec<String>>()
                    .join("\t")
            );
            header.push(last_header_row);

            Ok(header)
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

#[cfg(feature = "noodles")]
#[allow(dead_code)]
mod noodles {
    use std::{fs::File, path::PathBuf};

    use crate::error;

    use noodles::sam::header;
    pub struct Header(header::Header);

    impl Header {
        pub fn try_get(path: &PathBuf) -> Result<Self, error::Error> {
            let extension = path
                .extension()
                .unwrap_or_else(|| panic!("bam file {path:?} does not have a proper extension"))
                .to_str()
                .unwrap();

            let header: header::Header = match extension {
                "cram" => File::open(path)
                    .map_err(|e| error::Error::Io(path.clone(), e))
                    .map(noodles::cram::io::Reader::new)?
                    .read_header()?,
                "bam" => File::open(path)
                    .map_err(|e| error::Error::Io(path.clone(), e))
                    .map(noodles::bam::io::Reader::new)?
                    .read_header()?,
                _ => return Err(error::Error::UnknownExtension(path.clone())),
            };

            Ok(Self(header))
        }

        pub fn read_group_samples(&self) -> Result<Vec<String>, error::Error> {
            self.0
                .read_groups()
                .iter()
                .map(|(_, map)| {
                    let sm = noodles::sam::header::record::value::map::tag::Other::try_from([
                        b'S', b'M',
                    ])?;

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
    }
}
