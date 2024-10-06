use std::collections::BTreeSet;
#[doc(inline)]
use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use ndarray::{Array2, ShapeBuilder};
use rayon::prelude::*;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;
// use rust_htslib::bcf::Record;

use crate::{
    args::{Selection, StandardArgs},
    error::HatkError::{NormalizeError, PloidyError, SamplesNotFoundError},
    io::{get_htslib_contig_len, read_multiple_sample_ids},
    libs::structs::PhasedMatrix,
    structs::{Coord, Ploidy},
};

pub fn get_reader(
    path: &PathBuf,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
) -> Result<IndexedReader> {
    let mut reader = IndexedReader::from_path(path)?;
    let rid = reader.header().name2rid(contig.as_bytes())?;

    match range {
        // RUST-HTSLIB is 0-based so subtract 1
        Some((Some(start), Some(end))) => {
            reader.fetch(rid, start.saturating_sub(1), Some(end.saturating_sub(1)))?
        }
        Some((Some(start), None)) => reader.fetch(rid, start.saturating_sub(1), None)?,
        Some((None, Some(end))) => reader.fetch(rid, 0, Some(end.saturating_sub(1)))?,
        _ => reader.fetch(rid, 0, None)?,
    };

    Ok(reader)
}

pub fn filter_samples(samples: &[String], wanted: Option<Vec<String>>) -> Vec<usize> {
    if let Some(wanted) = wanted {
        for i in &wanted {
            if !samples.contains(i) {
                tracing::warn!("Wanted sample {i} is not in the VCF");
            }
        }

        samples
            .iter()
            .enumerate()
            .filter(|(_, s)| wanted.contains(s))
            .map(|(i, _)| i)
            .collect()
    } else {
        (0..samples.len()).collect()
    }
}

pub fn get_samples(header: &HeaderView) -> Result<Vec<String>> {
    header
        .samples()
        .into_iter()
        .map(|sample| Ok(String::from_utf8_lossy(sample).to_string()))
        .collect()
}

pub fn get_sample_names(
    args: &StandardArgs,
    contig: &str,
    wanted_samples: Option<Vec<String>>,
) -> Result<(Vec<usize>, Vec<String>)> {
    let reader = get_reader(&args.file, contig, None)?;

    let samples = reader
        .header()
        .samples()
        .into_iter()
        .map(|sample| Ok(String::from_utf8_lossy(sample).to_string()))
        .collect::<Result<Vec<String>>>()?;

    let wanted = match wanted_samples.clone() {
        Some(wanted) => Some(wanted),
        None => read_multiple_sample_ids(&args.samples)?,
    };

    let sample_indexes = filter_samples(&samples, wanted);

    ensure!(!sample_indexes.is_empty(), SamplesNotFoundError);
    let samples: Vec<String> = sample_indexes
        .iter()
        .map(|s| &samples[*s])
        .cloned()
        .collect();

    Ok((sample_indexes, samples))
}

pub fn read_vcf_to_matrix(
    args: &StandardArgs,
    contig: &str,
    variant_pos: u64,
    range: Option<(Option<u64>, Option<u64>)>,
    wanted_samples: Option<Vec<String>>,
) -> Result<PhasedMatrix> {
    let (indexes, samples) = get_sample_names(args, contig, wanted_samples)?;

    tracing::info!("Input VCF: {:?}", args.file);
    tracing::info!("Reading phased genotypes from {contig} with target position at {variant_pos}.");
    tracing::info!("Using {} samples from the VCF.", indexes.len());

    // If contig length not reported, read single thread
    if get_htslib_contig_len(&args.file, contig).is_err() {
        tracing::info!(
            "Contig {} length missing from the VCF. Reading genotypes single-threaded.",
            contig
        );
        let (markers, coords) = read_vcf_batch_to_matrix(args, contig, range, &indexes)?;

        tracing::info!("Finished reading all genotypes.");
        construct_phased_matrix(
            samples,
            markers,
            coords,
            &args.selection,
            variant_pos,
            contig,
        )
    } else {
        read_parallel(args, contig, variant_pos, range, samples, &indexes)
    }
}

fn read_parallel(
    args: &StandardArgs,
    contig: &str,
    variant_pos: u64,
    range: Option<(Option<u64>, Option<u64>)>,
    samples: Vec<String>,
    sample_indexes: &[usize],
) -> Result<PhasedMatrix> {
    let (first_pos, last_pos) = match range {
        Some((Some(start), Some(end))) => (start, end),
        _ => (0, get_htslib_contig_len(&args.file, contig)?),
    };
    let mut batches = vec![];
    let nthreads = rayon::current_num_threads();
    let batch_size = (last_pos - first_pos) / (nthreads * 5) as u64;
    let mut start = first_pos;

    loop {
        // plus one to evade rounding down
        let end = start + batch_size + 1;
        if end <= last_pos {
            batches.push((start, end));
            start = end + 1;
        } else {
            batches.push((start, last_pos));
            break;
        }
    }

    tracing::info!("Reading the VCF genotypes in {} batches.", batches.len());

    let matrices = batches
        .par_iter()
        .map(|(start, stop)| {
            read_vcf_batch_to_matrix(
                args,
                contig,
                Some((Some(*start), Some(*stop))),
                sample_indexes,
            )
        })
        .collect::<Result<Vec<(Vec<u8>, BTreeSet<Coord>)>>>()?;

    tracing::info!("Finished parallelized construction of genotype matrices.");

    let (markers, coords) = matrices
        .into_iter()
        .fold((vec![], BTreeSet::new()), |mut acc, cur| {
            acc.0.extend(cur.0);
            acc.1.extend(cur.1);
            acc
        });

    construct_phased_matrix(
        samples,
        markers,
        coords,
        &args.selection,
        variant_pos,
        contig,
    )
}

fn read_vcf_batch_to_matrix(
    args: &StandardArgs,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    sample_indexes: &[usize],
) -> Result<(Vec<u8>, BTreeSet<Coord>)> {
    let mut reader = get_reader(&args.file, contig, range)?;

    let mut coords = BTreeSet::new();
    let mut markers = vec![];

    let mut gt_buffer = rust_htslib::bcf::record::Buffer::new();

    for record in reader.records() {
        let record = record?;

        // HTSlib is 0-based so add 1
        let pos = (record.pos() + 1) as u64;

        tracing::trace!("Reading record at position {pos}");

        ensure!(record.alleles().len() == 2, NormalizeError(pos));

        let gts = record.genotypes_shared_buffer(&mut gt_buffer)?;

        let markers_len_before = markers.len();

        for i in sample_indexes {
            markers.extend(gts.get(*i).iter().filter_map(genotype_to_u8));
        }

        let diff = markers.len() - markers_len_before;
        let gts_per_sample = diff / sample_indexes.len();
        let gts_mod = diff % sample_indexes.len();

        ensure!(
            check_ploidy(gts_per_sample, gts_mod, args),
            PloidyError((pos, gts_per_sample))
        );

        coords.insert(construct_coord(&record, contig, pos));
    }

    Ok((markers, coords))
}

fn construct_phased_matrix(
    samples: Vec<String>,
    markers: Vec<u8>,
    coords: BTreeSet<Coord>,
    selection: &Selection,
    variant_pos: u64,
    contig: &str,
) -> Result<PhasedMatrix> {
    let ploidy: Ploidy = selection.as_ref().into();
    let matrix =
        Array2::from_shape_vec((samples.len() * ploidy as usize, coords.len()).f(), markers)?;

    let mut vcf = PhasedMatrix::new(0, Coord::default(), matrix, samples, coords, selection);

    vcf.variant_idx = vcf.get_first_idx_on_right_by_pos(variant_pos);
    vcf.start_coord = vcf.get_nearest_coord_by_pos(variant_pos).clone();

    if selection == &Selection::OnlyAlts || selection == &Selection::OnlyRefs {
        vcf.idx_by_pos(variant_pos).ok_or_else(|| {
            eyre!("Cannot select using a variant. Not found at position {contig}:{variant_pos}")
        })?;
    }

    tracing::info!(
        "Constructed a genotype matrix from {} records. Starting variant position: {}.",
        vcf.ncoords(),
        vcf.variant_idx_pos(),
    );
    Ok(vcf)
}

fn check_ploidy(length: usize, gts_mod: usize, args: &StandardArgs) -> bool {
    if gts_mod != 0 {
        return false;
    }
    match args.selection {
        Selection::Haploid => length == 1 || length == 2,
        _ => length == 2,
    }
}
fn genotype_to_u8(g: &GenotypeAllele) -> Option<u8> {
    match g {
        GenotypeAllele::Unphased(v) => Some(*v as u8),
        GenotypeAllele::Phased(v) => Some(*v as u8),
        GenotypeAllele::UnphasedMissing => None,
        GenotypeAllele::PhasedMissing => None,
    }
}

fn construct_coord(record: &rust_htslib::bcf::Record, contig: &str, pos: u64) -> Coord {
    let alleles = record.alleles();
    let reference = String::from_utf8_lossy(alleles.first().unwrap()).to_string();
    let alt = String::from_utf8_lossy(alleles.get(1).unwrap()).to_string();

    Coord {
        contig: contig.to_string(),
        reference,
        alt,
        pos,
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn sample_filtering() {
        let samples = vec!["foo".to_string(), "fii".to_string()];
        let wanted = vec!["foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![0,1]);

        let samples = vec!["foo".to_string(), "fii".to_string()];
        let wanted = vec!["fee".to_string(), "foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![0,1]);

        let samples = vec!["faa".to_string(), "foo".to_string(), "fii".to_string()];
        let wanted = vec!["fee".to_string(), "foo".to_string(), "fii".to_string()];
        let sample_indexes = filter_samples(&samples, Some(wanted));
        assert_eq!(sample_indexes, vec![1,2]);

        let sample_indexes = filter_samples(&samples, None);
        assert_eq!(sample_indexes, vec![0,1,2]);
    }
}
