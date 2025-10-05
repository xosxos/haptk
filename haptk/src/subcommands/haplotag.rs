use std::path::PathBuf;

use crate::args::{Selection, StandardArgs};
use color_eyre::{eyre::ensure, Result};

#[doc(hidden)]
pub fn run(
    mut args: StandardArgs,
    bam_file: PathBuf,
    ref_file: PathBuf,
    contigs: Vec<String>,
    threads: usize,
) -> Result<()> {
    ensure!(
        args.selection == Selection::All,
        "Only running with phased data and all chromosomes is supported."
    );

    let output = args.output.clone();
    let conf = tagger::Configuration { threads };

    if args.include_indels {
        args.include_indels = false;
        tracing::info!("Setting `include_indels` to false for haplotagging");
    }

    tagger::run_haplotag(
        &args,
        vec![bam_file],
        &ref_file,
        Some(output),
        contigs,
        conf,
    )?;

    Ok(())
}

mod tagger {
    use std::collections::{BTreeMap, HashMap};
    use std::ops::Range;
    use std::path::{Path, PathBuf};
    use std::thread;

    use color_eyre::Result;
    use crossbeam_channel::{unbounded, Sender};
    use rayon::prelude::*;
    use rust_htslib::bam::{Read, Record};

    use crate::args::StandardArgs;
    use crate::read_vcf::read_vcf_to_matrix;
    use crate::structs::CigarVariant;
    use crate::subcommands::haplotag::cigar_iterator::{CigarIterType, CigarIterator};
    use crate::subcommands::haplotag::io::spawn_collector;
    use crate::subcommands::haplotag::utils::get_bam_header;

    use super::utils::{filter_contigs_by_bam, get_sample_name_and_bam_paths, RangeDivisions};

    pub struct Configuration {
        pub threads: usize,
    }

    // type ChannelObj = RecordBuf;
    pub type ChannelObj = Vec<String>;

    // static RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"(?i)gaa").unwrap());

    pub fn run_haplotag(
        args: &StandardArgs,
        paths: Vec<PathBuf>,
        ref_path: &Path,
        output: Option<PathBuf>,
        contigs: Vec<String>,
        conf: Configuration,
    ) -> Result<()> {
        let file_path = paths.first().unwrap().clone();

        let header = get_bam_header(&file_path).unwrap();

        // Read the first bam file and base the contigs and vcf header on this bam
        let (contigs, contigs_and_len) = filter_contigs_by_bam(contigs, &header)?;

        // Attach sample IDs to all bam file paths
        let sample_id_and_bam_paths = get_sample_name_and_bam_paths(paths)?;

        let (collector_tx, collector_rx) = unbounded();

        let collector_handle =
            spawn_collector(collector_rx, output, file_path, ref_path.to_path_buf());

        // Find insertions for all samples in each contig,
        thread::scope(|_s| -> Result<()> {
            // println!("thread {contigs:?}");
            contigs
                .into_par_iter()
                .filter(|(c, _, _)| c != "chrY")
                .filter(|(c, _, _)| c != "chrM")
                .try_for_each(|(contig, start, stop)| -> Result<()> {
                    let (start, stop) = match (start, stop) {
                        (Some(start), Some(stop)) => (start, stop),
                        (None, None) => (0, *contigs_and_len.get(&contig).unwrap()),
                        _ => unreachable!(),
                    };

                    let vcf = read_vcf_to_matrix(
                        &args,
                        &contig,
                        0,
                        Some((Some(start), Some(stop))),
                        Some(
                            sample_id_and_bam_paths
                                .iter()
                                .map(|(id, _)| id.clone())
                                .collect(),
                        ),
                        None,
                        false,
                    )?;

                    let ploidy = *vcf.ploidy;

                    let mut haplotypes = HashMap::new();

                    for (sample_id, _) in &sample_id_and_bam_paths {
                        let idxs = vcf.get_idxs_for_samples(&[sample_id.clone()]).unwrap();

                        let haps: BTreeMap<CigarVariant, Vec<u8>> = vcf
                            .coords()
                            .iter()
                            .flat_map(|coord| {
                                let variant = CigarVariant {
                                    pos: coord.pos,
                                    alt: coord.alt.chars().nth(0).unwrap(),
                                };

                                let gts: Vec<u8> = idxs
                                    .iter()
                                    .map(|sample_idx| vcf.matrix_point_coord(*sample_idx, coord))
                                    .collect();

                                // Allelles are homozygous, return
                                if gts.iter().all(|v| v == &gts[0]) {
                                    return None;
                                }

                                Some((variant, gts))
                            })
                            .collect();
                        haplotypes.insert(sample_id.clone(), haps);
                    }

                    let max_chunk_len = 5_000_000;
                    let min_chunk_len = 10_000_000;

                    let nsamples = sample_id_and_bam_paths.len();
                    let chunks = (conf.threads / nsamples).max(1);

                    let length = stop.saturating_sub(start);
                    let chunk_len = length / chunks as u64;

                    let should_multithread = || -> bool { chunk_len > max_chunk_len };

                    // Sometimes its best to single thread depending on I/O speeds
                    // the current heuristic is obv not a very sophisticated way to make this decision
                    match should_multithread() {
                        // Single thread (sometimes its best to single thread depending on I/O speeds)
                        true => {
                            tracing::info!("Chunk len > max_chunk_len, skipping par_iter");
                            let chunks = stop.saturating_sub(start) / max_chunk_len;
                            let coordinate_ranges: Vec<Range<u32>> = (start as u32..stop as u32)
                                .divide_evenly_into(chunks as usize)
                                .collect();

                            coordinate_ranges
                                .into_par_iter()
                                // .into_iter()
                                .try_for_each(|range| -> Result<()> {
                                    iterate_region(
                                        &sample_id_and_bam_paths,
                                        ref_path,
                                        &contig,
                                        range,
                                        &collector_tx,
                                        &haplotypes,
                                        ploidy,
                                    )
                                })?;
                        }
                        // Multi thread
                        false => {
                            let chunks =
                                stop.saturating_sub(start) / (chunk_len).max(min_chunk_len);
                            let coordinate_ranges: Vec<Range<u32>> = (start as u32..stop as u32)
                                .divide_evenly_into(chunks.max(1) as usize)
                                .collect();
                            tracing::info!("false chunks {}", chunks.max(1));

                            coordinate_ranges
                                .into_par_iter()
                                // .into_iter()
                                .try_for_each(|range| -> Result<()> {
                                    tracing::info!("checking range {:?}", range);
                                    iterate_region(
                                        &sample_id_and_bam_paths,
                                        ref_path,
                                        &contig,
                                        range,
                                        &collector_tx,
                                        &haplotypes,
                                        ploidy,
                                    )
                                })?;
                        }
                    }

                    Ok(())
                })?;
            Ok(())
        })?;

        drop(collector_tx);
        collector_handle.join().unwrap()?;

        Ok(())
    }

    type BamReader = rust_htslib::bam::IndexedReader;
    // type BamRecord = rust_htslib::errors::Result<rust_htslib::bam::Record>;

    pub fn get_bam_reader<'a>(
        bam_path: &'a PathBuf,
        ref_path: &Path,
        contig: &String,
        range: &Range<u32>,
        sample_id: &String,
    ) -> Result<BamReader> {
        // Fetch reads in a given contig from the bam
        let mut rdr = BamReader::from_path(bam_path)?;
        rdr.set_reference(ref_path)?;

        rdr.fetch((&contig, range.start, range.end))?;
        tracing::info!("Fetched {}..{} for {}", range.start, range.end, sample_id);

        Ok(rdr)
    }

    pub fn iterate_region<'a>(
        sample_id_and_bam_paths: &[(String, PathBuf)],
        ref_path: &Path,
        contig: &String,
        range: Range<u32>,
        collector_tx: &Sender<ChannelObj>,
        haplotypes: &HashMap<String, BTreeMap<CigarVariant, Vec<u8>>>,
        ploidy: usize,
    ) -> Result<()> {
        // Individual level merges and process
        sample_id_and_bam_paths
            .iter()
            .enumerate()
            // .par_bridge()
            .try_for_each(|(_sample_idx, (sample_id, bam_path))| -> Result<()> {
                tracing::info!("Searching insertions from {contig} for {sample_id}");

                let mut rdr = get_bam_reader(bam_path, ref_path, contig, &range, sample_id)?;

                let mut records = rdr.records();

                records.try_for_each(|record| match record {
                    Ok(record) => match_cigar_to_haplotype(
                        record,
                        contig,
                        sample_id,
                        collector_tx,
                        haplotypes,
                        ploidy,
                    ),
                    Err(_) => {
                        tracing::warn!("Invalid record found for {sample_id} {contig}");
                        Ok(())
                    }
                })?;

                tracing::info!("Unloading reads");
                Ok(())
            })?;

        Ok(())
    }

    fn match_cigar_to_haplotype<'a>(
        record: Record,
        contig: &'a str,
        sample_id: &'a str,
        collector_tx: &Sender<ChannelObj>,
        haplotypes: &HashMap<String, BTreeMap<CigarVariant, Vec<u8>>>,
        ploidy: usize,
    ) -> Result<()> {
        let seq_len = record.seq().len();

        if seq_len == 0 {
            return Ok(());
        }

        let haplotype_map = haplotypes.get(sample_id).unwrap();

        let mut haplotype_score = vec![0; ploidy];
        let mut positions_debug = vec![];

        CigarIterator::new(
            record.cigar().take().into(),
            record.pos() as u32,
            record.seq(),
        )
        .for_each(|cigar_item| match cigar_item {
            CigarIterType::Match((pos, seq)) => {
                // Position returned is the last pos before start of the match, so add 1
                // Or maybe this is 0 indexed positions from HTSLib
                let pos = pos + 1;

                // println!("match: {} {}", pos, seq);
                let start = CigarVariant {
                    pos: pos as u64,
                    alt: 'A',
                };

                let stop = CigarVariant {
                    pos: (pos as u64 + (seq.len() - 1) as u64),
                    alt: 'X',
                };

                for (variant, gts) in haplotype_map.range(start..=stop) {
                    // +1 cause indexing starts from 0
                    // println!("{}-{} vs {}", variant.pos, variant.alt, pos);

                    let query_pos = variant.pos as u32 - pos;
                    // println!("query pos {}", query_pos);

                    let allele = seq.chars().nth(query_pos as usize).unwrap();

                    for (i, gt) in gts.iter().enumerate() {
                        if gt == &1 && variant.alt == allele {
                            positions_debug.push(variant.pos);
                            haplotype_score[i] += 1;
                        }
                    }
                }
            }
            _ => (),
        });

        fn fold_func<S: ToString>(acc: String, cur: S) -> String {
            match acc.is_empty() {
                true => cur.to_string(),
                false => format!("{acc},{}", cur.to_string()),
            }
        }

        let aux: Vec<String> = record
            .aux_iter()
            .map(|v| v.unwrap())
            .map(|(tag, v)| {
                let tag = String::from_utf8(tag.to_vec()).unwrap();
                match v {
                    rust_htslib::bam::record::Aux::Char(i) => i.to_string(),
                    rust_htslib::bam::record::Aux::I8(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U8(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::I16(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U16(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::I32(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::U32(i) => format!("{tag}:i:{i}"),
                    rust_htslib::bam::record::Aux::Float(i) => format!("{tag}:f:{i}"),
                    rust_htslib::bam::record::Aux::Double(i) => format!("{tag}:f:{i}"),
                    rust_htslib::bam::record::Aux::String(i) => format!("{tag}:Z:{i}"),
                    rust_htslib::bam::record::Aux::HexByteArray(i) => format!("{tag}:H:{i}"),
                    rust_htslib::bam::record::Aux::ArrayI8(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU8(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayI16(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU16(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayI32(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayU32(i) => {
                        format!("{tag}:B:i,{}", i.iter().fold(String::new(), fold_func))
                    }
                    rust_htslib::bam::record::Aux::ArrayFloat(i) => {
                        format!("{tag}:B:f,{}", i.iter().fold(String::new(), fold_func))
                    }
                }
            })
            .collect();

        let seq = String::from_utf8(record.seq().as_bytes()).unwrap();

        let mut row: Vec<String> = vec![
            String::from_utf8(record.qname().to_vec()).unwrap(),
            record.flags().to_string(),
            contig.to_string(),
            // Convert back to 1-based positions
            (record.pos() + 1).to_string(),
            record.mapq().to_string(),
            record.cigar().to_string(),
            // record.mtid().to_string(),
            "*".to_string(),
            // record.mpos().to_string(),
            "0".to_string(),
            // seq.len().to_string(),
            "0".to_string(),
            seq,
            "*".to_string(),
            // String::from_utf8(record.qual().to_vec()).unwrap(),
        ];

        row.extend(aux);

        let find_best_match = |haplotype_score: &Vec<u32>| -> Option<&str> {
            let margin = 3;

            if haplotype_score[0] >= haplotype_score[1] + margin {
                return Some("0");
            }

            if haplotype_score[0] + margin <= haplotype_score[1] {
                return Some("1");
            }

            return None;
        };

        row.push(format!(
            "ht:Z:{}",
            find_best_match(&haplotype_score).unwrap_or(&"NA")
        ));

        row.push(format!(
            "sc:B:i,{},{}",
            haplotype_score[0], haplotype_score[1],
        ));

        collector_tx.send(row)?;

        Ok(())
    }
}

mod utils {
    use std::fs::File;
    use std::path::PathBuf;
    use std::{hash::BuildHasher, ops::Range};

    use color_eyre::{
        eyre::{eyre, OptionExt, WrapErr},
        Result,
    };
    use num_traits::{FromPrimitive, Num};
    use rayon::prelude::*;

    use gxhash::GxBuildHasher;
    use indexmap::IndexMap;

    use crate::error::Error;

    /// A `HashMap` using a (DOS-resistant) [`GxBuildHasher`].
    pub type HashMap<K, V> = IndexMap<K, V, GxBuildHasher>;

    /// A convenience trait that can be used together with the type aliases defined
    /// to get access to the `new()` and `with_capacity()` methods for the
    /// [`HashMap`] type alias.
    pub trait HashMapExt {
        /// Constructs a new HashMap.
        fn new() -> Self;
    }

    impl<K, V, S> HashMapExt for IndexMap<K, V, S>
    where
        S: BuildHasher + Default,
    {
        fn new() -> Self {
            IndexMap::with_hasher(S::default())
        }
    }

    #[rustfmt::skip]
    pub fn get_bam_header(path: &PathBuf) -> Result<noodles::sam::Header> {
        let extension = path
            .extension()
            .expect(&format!(
                "bam file {path:?} does not have a proper extension"
            ))
            .to_str()
            .unwrap();

        let header: noodles::sam::Header = match extension {
            "cram" => File::open(path).map(noodles::cram::io::Reader::new).wrap_err(Error::Io { path: path.clone()})?.read_header()?,
            "bam" => File::open(path).map(noodles::bam::io::Reader::new).wrap_err(Error::Io {path: path.clone()})?.read_header()?,
            e => return Err(eyre!("unknown file extension {e} for {path:?}")),
        };

        Ok(header)
    }

    pub fn get_sample_name_and_bam_paths(paths: Vec<PathBuf>) -> Result<Vec<(String, PathBuf)>> {
        let mut names_and_paths: Vec<(String, PathBuf)> = paths
        .into_par_iter()
        // .into_iter()
        .enumerate()
        .map(|(i, path)| {

            let header = get_bam_header(&path)?;
            let read_groups: &noodles::sam::header::ReadGroups = header.read_groups();

            let Some((_, map)) = read_groups.first() else  {
                tracing::warn!("No sample name in the .bam header for {path:?}. Naming the sample unnamed_{i}");
                let sample_name = format!("unnamed_{i}");

                return Ok((sample_name, path));
            };

            let sm = noodles::sam::header::record::value::map::tag::Other::try_from([b'S', b'M'])?;

            let sample_name = map.other_fields().get(&sm).unwrap().to_string();

            if read_groups.len() > 1 {
                tracing::warn!( "Multiple read groups, selecting the first one: {sample_name}", )
            }
            Ok((sample_name, path))
        })
        .collect::<Result<Vec<(String, PathBuf)>>>()?;

        names_and_paths.sort_by(|a, b| a.0.cmp(&b.0));

        let name_vec: Vec<String> = names_and_paths
            .iter()
            .map(|(name, _path)| name.clone())
            .collect();

        names_and_paths.iter_mut().enumerate().for_each(|(i, v)| {
            if name_vec[i + 1..].contains(&v.0) {
                v.0.push_str(&i.to_string());
            }
        });

        Ok(names_and_paths)
    }

    pub fn filter_contigs_by_bam(
        contigs: Vec<String>,
        header: &noodles::sam::Header,
    ) -> Result<(
        Vec<(String, Option<u64>, Option<u64>)>,
        HashMap<String, u64>,
    )> {
        // println!("{contigs:?}");
        let reference_sequences = header.reference_sequences();

        let mut bam_contigs_and_lengths: HashMap<String, u64> = HashMap::new();
        let mut bam_contigs = vec![];

        for (seqname, map) in reference_sequences {
            let seq_len = map.length();
            bam_contigs_and_lengths.insert(seqname.to_string(), seq_len.get() as u64);
            bam_contigs.push(seqname.to_string());
        }

        let mut contigs = contigs
            .into_iter()
            .map(parse_coords)
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

    // Coords are in the format [contig] OR [contig]:[start]-[stop]
    pub fn parse_coords(coords: String) -> Result<(String, Option<u64>, Option<u64>)> {
        let mut coord_split = coords.split(':');

        let contig = coord_split
            .next()
            .ok_or_eyre("cannot parse coord: {coords}")?;

        let positions = coord_split.next();

        if positions.is_none() {
            return Ok((contig.to_string(), None, None));
        }

        let mut pos_split = positions.unwrap().split('-');
        let (start, stop) = (pos_split.next(), pos_split.next());

        if let (Some(start), Some(stop)) = (start, stop) {
            let start = start
                .parse::<u64>()
                .wrap_err(eyre!("cannot start {start:?} of coord {coords:?}"))?;
            let stop = stop
                .parse::<u64>()
                .wrap_err(eyre!("cannot stop {stop:?} of coord {coords:?}"))?;

            Ok((contig.to_string(), Some(start), Some(stop)))
        } else {
            Err(eyre!("cannot parse coord: {coords}"))
        }
    }

    /// Split range into an iterator of smaller ranges
    pub trait RangeDivisions<T: Num + FromPrimitive + PartialOrd + Copy> {
        fn divide_evenly_into(self, divisions: usize) -> Even<T>;
    }

    impl<T: Num + FromPrimitive + PartialOrd + Copy> RangeDivisions<T> for Range<T> {
        fn divide_evenly_into(self, divisions: usize) -> Even<T> {
            Even {
                next_start: self.start,
                next_end: self.start,
                end: self.end,
                div_remaining: divisions,
            }
        }
    }

    pub struct Even<T: Num + FromPrimitive + PartialOrd + Copy> {
        next_start: T,
        next_end: T,
        end: T,
        div_remaining: usize,
    }

    impl<T: Num + FromPrimitive + PartialOrd + Copy> Iterator for Even<T> {
        type Item = Range<T>;

        fn next(&mut self) -> Option<Self::Item> {
            if self.div_remaining > 0 {
                self.next_start = self.next_end;
                if self.next_start < self.end {
                    self.next_end = self.next_start
                        + (self.end - self.next_start)
                            / T::from_usize(self.div_remaining)
                                .expect("divisions usize cannot be converted to range type");
                } else {
                    self.next_end = self.next_start
                        - (self.next_start - self.end)
                            / T::from_usize(self.div_remaining)
                                .expect("divisions usize cannot be converted to range type");
                }

                self.div_remaining -= 1;

                Some(self.next_start..self.next_end)
            } else {
                None
            }
        }
    }
}

mod cigar_iterator {
    use std::vec::IntoIter;

    use rust_htslib::bam::record::{Cigar, Seq};

    #[derive(Debug, Clone, Eq, PartialEq)]
    pub enum CigarIterType {
        Diff((u32, String)),
        Del((u32, u32)),
        Ins((u32, String)),
        LeadingSoftClip((u32, String)),
        TrailingSoftClip((u32, String)),
        HardClip(u32),
        Equal(u32),
        Match((u32, String)),
    }

    // impl CigarIterType {
    //     pub fn pos(&self) -> u32 {
    //         match self {
    //             CigarIterType::Diff((pos, _)) => *pos,
    //             CigarIterType::Del((pos, _)) => *pos,
    //             CigarIterType::Ins((pos, _)) => *pos,
    //             CigarIterType::Equal(pos) => *pos,
    //             CigarIterType::Match((pos, _)) => *pos,
    //             CigarIterType::HardClip(pos) => *pos,
    //             CigarIterType::LeadingSoftClip((pos, _)) => *pos,
    //             CigarIterType::TrailingSoftClip((pos, _)) => *pos,
    //         }
    //     }
    // }

    #[derive(Debug, Clone)]
    pub struct CigarIterator {
        read_start_pos: u32,
        current_pos: u32,
        cigar: IntoIter<Cigar>,
        seq: Vec<u8>,
        passed: usize,
    }

    impl CigarIterator {
        pub fn new(cigar: Vec<Cigar>, read_start_pos: u32, seq: Seq<'_>) -> Self {
            Self {
                read_start_pos,
                current_pos: 0,
                cigar: cigar.into_iter(),
                seq: seq.as_bytes(),
                passed: 0,
            }
        }
    }

    impl Iterator for CigarIterator {
        type Item = CigarIterType;

        fn next(&mut self) -> Option<Self::Item> {
            if let Some(cigar) = self.cigar.next() {
                let pos = self.read_start_pos + self.current_pos;

                let start = self.current_pos as usize;

                let variant = match cigar {
                    Cigar::RefSkip(_) => panic!("Refskip N present in Cigar strings"),
                    Cigar::Pad(_) => panic!("Padding P present in Cigar strings"),
                    Cigar::Match(bp_len) => {
                        self.current_pos += bp_len;
                        CigarIterType::Match((
                            pos,
                            to_string_seq(&self.seq[start..self.current_pos as usize]),
                        ))
                    }
                    // current_pos is not incremented in deletions
                    Cigar::Del(bp_len) => {
                        //self.current_pos += n;
                        self.read_start_pos += bp_len;
                        CigarIterType::Del((pos, bp_len))
                    }
                    Cigar::Equal(bp_len) => {
                        self.current_pos += bp_len;
                        CigarIterType::Equal(pos)
                    }
                    Cigar::Ins(bp_len) => {
                        self.read_start_pos = self.read_start_pos.saturating_sub(bp_len);
                        self.current_pos += bp_len;
                        CigarIterType::Ins((
                            pos,
                            to_string_seq(&self.seq[start..self.current_pos as usize]),
                        ))
                    }

                    Cigar::SoftClip(bp_len) => {
                        if self.current_pos == 0 {
                            // Htslib .pos() returns the read start position without leading softclips
                            self.read_start_pos = self.read_start_pos.saturating_sub(bp_len);

                            self.current_pos += bp_len;
                            CigarIterType::LeadingSoftClip((
                                self.read_start_pos,
                                to_string_seq(&self.seq[start..self.current_pos as usize]),
                            ))
                        } else {
                            self.current_pos += bp_len;
                            CigarIterType::TrailingSoftClip((
                                pos,
                                to_string_seq(&self.seq[start..self.current_pos as usize]),
                            ))
                        }
                    }
                    Cigar::HardClip(_bp_len) => {
                        CigarIterType::HardClip(pos)
                        // if self.current_pos == 0 {
                        // Htslib bam record .pos() returns a position that maybe disregards
                        // hardclip?
                        // TODO: check this out
                        // self.read_start_pos = self.read_start_pos.saturating_sub(bp_len);
                        // self.current_pos += bp_len;
                        // CigarIterType::HardClip(pos)
                        // } else {
                        // self.current_pos += bp_len;
                        // CigarIterType::HardClip(pos)
                        // }
                        //self.current_pos += n;
                    }

                    Cigar::Diff(bp_len) => {
                        self.current_pos += bp_len;
                        CigarIterType::Diff((
                            pos,
                            to_string_seq(&self.seq[start..self.current_pos as usize]),
                        ))
                    }
                };

                self.passed += 1;

                Some(variant)
            } else {
                if self.passed == 0 && self.current_pos as usize != self.seq.len() {
                    // unsafe {
                    //     let snapshot = String::from_utf8_unchecked(self.seq.clone())
                    //         .chars()
                    //         .take(70)
                    //         .collect::<String>();

                    //     let last_chars = String::from_utf8_unchecked(self.seq.clone())
                    //         .chars()
                    //         .rev()
                    //         .take(10)
                    //         .collect::<String>();

                    //     tracing::debug!(
                    //         "Unmapped read for {} pos: {} sequence: {}..{}..{}",
                    //         self.sample_id,
                    //         self.read_start_pos,
                    //         snapshot,
                    //         self.seq
                    //             .len()
                    //             .saturating_sub(snapshot.len())
                    //             .saturating_sub(last_chars.len()),
                    //         last_chars,
                    //     );
                    // }
                } else if self.passed > 0 {
                    assert_eq!(
                    self.current_pos as usize,
                    self.seq.len(),
                    "Incorrect cigar string: cigar: {:?}, iterated nucletotides {:?}, total sequence: {:?}",
                    self.cigar,
                    self.passed,
                    String::from_utf8(self.seq.clone()).unwrap(),
                    // to_string_seq(self.seq),
                );
                }
                None
            }
        }
    }

    pub fn to_string_seq(bytes: &[u8]) -> String {
        // String::from_utf8(bytes.to_vec()).unwrap()
        unsafe { String::from_utf8_unchecked(bytes.to_vec()) }
    }
}

mod io {
    use color_eyre::eyre::Context;
    use color_eyre::{eyre::eyre, Result};
    use crossbeam_channel::Receiver;
    use rust_htslib::bam::{CompressionLevel, Format, Header, IndexedReader, Read, Writer};
    use std::path::PathBuf;
    use std::thread::{self, JoinHandle};

    use crate::subcommands::haplotag::tagger::ChannelObj;

    // use crate::structs::Vcf;

    // use std::io;
    // pub fn get_output(filename: Option<PathBuf>) -> Result<Box<dyn io::Write>> {
    //     let output: Box<dyn io::Write> = match filename {
    //         Some(name) => match name.to_str() {
    //             Some("-") => Box::new(io::stdout()),
    //             Some(name) => Box::new(
    //                 match std::fs::File::options()
    //                     .create(true)
    //                     .write(true)
    //                     .truncate(true)
    //                     .open(name)
    //                 {
    //                     Ok(x) => x,
    //                     Err(err) => {
    //                         let msg = format!("failed to open \"{}\": {err}", name);
    //                         return Err(std::io::Error::new(std::io::ErrorKind::NotFound, msg))?;
    //                     }
    //                 },
    //             ),
    //             None => return Err(eyre!("unknown error 1")),
    //         },
    //         None => Box::new(io::stdout()),
    //     };
    //     Ok(output)
    // }

    pub fn spawn_collector(
        rx: Receiver<ChannelObj>,
        output: Option<PathBuf>,
        bam_path: PathBuf,
        ref_path: PathBuf,
    ) -> JoinHandle<Result<()>> {
        thread::spawn(move || -> Result<()> {
            let bam = IndexedReader::from_path(bam_path).wrap_err(eyre!("path:?"))?;
            let header_view = bam.header();
            // let output = get_output(output)?;

            let header = Header::from_template(header_view);

            let mut cram_writer = Writer::from_path(&output.unwrap(), &header, Format::Cram)?;
            cram_writer.set_reference(ref_path)?;
            cram_writer.set_compression_level(CompressionLevel::Fastest)?;

            // let mut writer = csv::WriterBuilder::new()
            //     .delimiter(b'\t')
            //     .has_headers(false)
            //     .flexible(true)
            //     .quote_style(csv::QuoteStyle::Never)
            //     .from_writer(output);

            // let stdout = io::stdout().lock();
            // let mut writer = noodles::sam::io::Writer::new(BufWriter::new(stdout));

            // if let Err(e) = out_header.iter().try_for_each(|s| writer.serialize(s)) {
            //     panic!("Collector thread panics due to writer error: {e:?}")
            // }

            while let Ok(record) = rx.recv() {
                let record = record
                    .into_iter()
                    .enumerate()
                    .fold(String::new(), |acc, (i, cur)| {
                        if i == 0 {
                            cur
                        } else {
                            format!("{acc}\t{cur}")
                        }
                    });

                // println!("{record}");

                let record = rust_htslib::bam::Record::from_sam(header_view, record.as_bytes())?;
                // writer.write_alignment_record(&header, &record)?;

                if let Err(e) = cram_writer.write(&record) {
                    panic!("Collector thread panics due to writer error: {e:?}")
                }
            }

            Ok(())
        })
    }
}
