use std::collections::BTreeMap;
use std::collections::HashMap;
use std::ops::Range;
use std::path::PathBuf;
use std::thread;
use std::thread::JoinHandle;

use color_eyre::{eyre::ensure, Result};
use crossbeam_channel::Receiver;
use crossbeam_channel::{unbounded, Sender};
use rayon::prelude::*;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::bam;
use crate::core::cigar_iterator::CigarIterType;
use crate::core::cigar_iterator::CigarIterator;
use crate::core::utils::RangeDivisions;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::CigarVariant;

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

    let conf = Configuration {
        threads,
        ref_path: ref_file,
        output: Some(args.output.clone()),
    };

    if args.include_indels {
        args.include_indels = false;
        tracing::info!("Setting `include_indels` to false for haplotagging");
    }

    run_haplotag(&args, vec![bam_file], contigs, conf)?;

    Ok(())
}

#[derive(Debug, Clone)]
pub struct Configuration {
    pub threads: usize,
    pub output: Option<PathBuf>,
    pub ref_path: PathBuf,
}

pub type ChannelObj = Vec<String>;

pub fn run_haplotag(
    args: &StandardArgs,
    paths: Vec<PathBuf>,
    contigs: Vec<String>,
    conf: Configuration,
) -> Result<()> {
    let file_path = paths.first().unwrap().clone();

    let header = bam::Header::try_get(&file_path)?;

    // Read the first bam file and base the contigs and vcf header on this bam
    let (contigs, contigs_and_len) = header.filter_contigs(contigs)?;

    // Attach sample IDs to all bam file paths
    let sample_id_and_bam_paths = get_sample_name_and_bam_paths(paths)?;

    let (collector_tx, collector_rx) = unbounded();

    let collector_handle = spawn_collector(collector_rx, file_path, conf.clone());

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
                    args,
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
                                // Take first character
                                alt: coord.alt.chars().next().unwrap(),
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
                        let coordinate_ranges: Vec<Range<u64>> =
                            (start..stop).divide_evenly_into(chunks as usize).collect();

                        coordinate_ranges
                            .into_par_iter()
                            // .into_iter()
                            .try_for_each(|range| -> Result<()> {
                                iterate_region(
                                    &sample_id_and_bam_paths,
                                    &conf,
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
                        let chunks = stop.saturating_sub(start) / (chunk_len).max(min_chunk_len);
                        let coordinate_ranges: Vec<Range<u64>> = (start..stop)
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
                                    &conf,
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

pub fn iterate_region(
    sample_id_and_bam_paths: &[(String, PathBuf)],
    conf: &Configuration,
    contig: &str,
    range: Range<u64>,
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

            let mut rdr =
                bam::Reader::from_path(bam_path, &conf.ref_path, contig, &range, sample_id)?;

            rdr.records().try_for_each(|record| {
                match_cigar_to_haplotype(
                    record,
                    contig,
                    sample_id,
                    collector_tx,
                    haplotypes,
                    ploidy,
                )
            })?;

            tracing::info!("Unloading reads");
            Ok(())
        })?;

    Ok(())
}

fn match_cigar_to_haplotype<'a>(
    record: bam::Record,
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
    // let mut positions_debug = vec![];

    let seq = record.seq();

    CigarIterator::new(
        record.cigar(),
        record.pos(),
        &seq,
        sample_id,
        false,
        u64::MAX,
    )
    .for_each(|cigar_item| {
        if let CigarIterType::Match(ref_pos, _, seq) = cigar_item {
            // 0-based positions from HTSLib, change to 1-based positions
            let ref_pos = ref_pos + 1;

            // println!("match: {} {}", pos, seq);
            let start = CigarVariant {
                pos: ref_pos,
                alt: 'A',
            };

            let stop = CigarVariant {
                pos: (ref_pos + (seq.len() - 1) as u64),
                alt: 'X',
            };

            for (variant, gts) in haplotype_map.range(start..=stop) {
                // Because the sequence in a Cigar::Match is identical to the reference sequence,
                // we can use ref_pos based indexes to index the Match sequence
                let query_pos = variant.pos - ref_pos;

                let allele = seq.chars().nth(query_pos as usize).unwrap();

                for (i, gt) in gts.iter().enumerate() {
                    if gt == &1 && variant.alt == allele {
                        haplotype_score[i] += 1;

                        // positions_debug.push(variant.pos);
                    }
                }
            }
        }
    });

    let seq = String::from_utf8(seq).unwrap();

    let mut row: Vec<String> = vec![
        record.qname(),
        record.flags().to_string(),
        contig.to_string(),
        // Convert back to 1-based positions
        (record.pos() + 1).to_string(),
        record.mapq().to_string(),
        record.cigar_string().to_string(),
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

    row.extend(record.aux());

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
        find_best_match(&haplotype_score).unwrap_or("NA")
    ));

    row.push(format!(
        "sc:B:i,{},{}",
        haplotype_score[0], haplotype_score[1],
    ));

    collector_tx.send(row)?;

    Ok(())
}

pub fn spawn_collector(
    rx: Receiver<ChannelObj>,
    input_path: PathBuf,
    conf: Configuration,
) -> JoinHandle<Result<()>> {
    thread::spawn(move || -> Result<()> {
        let mut cram_writer =
            bam::Writer::from_path(&conf.output.unwrap(), &input_path, &conf.ref_path)?;

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

            if let Err(e) = cram_writer.write(record) {
                panic!("Collector thread panics due to writer error: {e:?}")
            }
        }

        Ok(())
    })
}

pub fn get_sample_name_and_bam_paths(paths: Vec<PathBuf>) -> Result<Vec<(String, PathBuf)>> {
    let mut names_and_paths: Vec<(String, PathBuf)> = paths
        .into_par_iter()
        // .into_iter()
        .enumerate()
        .map(|(i, path)| {
            let header = bam::Header::try_get(&path)?;

            let samples = header.read_group_samples()?;

            if samples.is_empty() {
                tracing::warn!(
                    "No sample name in the .bam header for {path:?}. Naming the sample unnamed_{i}"
                );
                let sample_name = format!("unnamed_{i}");

                return Ok((sample_name, path));
            }

            if samples.len() > 1 {
                tracing::warn!(
                    "Multiple read groups, selecting the first one: {}",
                    &samples[0]
                )
            }

            Ok((samples[0].to_string(), path))
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
