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

use crate::{
    args::{Selection, StandardArgs},
    error::HaptkError::{NormalizeError, PloidyError, SamplesNotFoundError},
    io::{get_htslib_contig_len, read_multiple_sample_ids},
    structs::{Coord, PhasedMatrix, Ploidy, ReadMetadata},
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

fn genotype_to_u8(g: &GenotypeAllele) -> u8 {
    match g {
        GenotypeAllele::Unphased(v) => *v as u8,
        GenotypeAllele::Phased(v) => *v as u8,
        _ => panic!(),
    }
}

fn prune_by_gt(
    args: &StandardArgs,
    contig: &str,
    variant_pos: u64,
    indexes: &[usize],
    samples: Vec<String>,
    wanted_gt: u8,
) -> Result<(Vec<[bool; 2]>, Vec<String>)> {
    let mut reader = IndexedReader::from_path(&args.file)?;
    let rid = reader.header().name2rid(contig.as_bytes())?;

    // Positions are 0 based in HTSLib
    reader.fetch(
        rid,
        variant_pos.saturating_sub(1),
        Some(variant_pos.saturating_sub(1)),
    )?;

    let mut gt_buffer = rust_htslib::bcf::record::Buffer::new();

    let record = reader.records().next();
    let record = record.unwrap()?;
    let gts = record.genotypes_shared_buffer(&mut gt_buffer)?;

    // let pos = (record.pos() + 1) as u64;
    // let coord = construct_coord(&record, contig, pos);

    let lookups: Vec<[bool; 2]> = indexes
        .iter()
        .map(|i| {
            let sample_gts = gts.get(*i);

            let gt1 = genotype_to_u8(&sample_gts[0]);
            let gt2 = genotype_to_u8(&sample_gts[1]);
            [gt1 == wanted_gt, gt2 == wanted_gt]
        })
        .collect();

    let samples = samples
        .into_iter()
        .zip(lookups.iter())
        // cheap enough
        .flat_map(|(s, lookup)| match (lookup[0], lookup[1]) {
            (true, true) => vec![s.clone(), s],
            (true, false) | (false, true) => vec![s],
            (false, false) => vec![],
        })
        .collect();

    Ok((lookups, samples))
}

pub fn read_vcf_to_matrix(
    args: &StandardArgs,
    contig: &str,
    variant_pos: u64,
    mut range: Option<(Option<u64>, Option<u64>)>,
    wanted_samples: Option<Vec<String>>,
    sharded: bool,
) -> Result<PhasedMatrix> {
    let (indexes, samples) = get_sample_names(args, contig, wanted_samples)?;

    let (lookups, samples) = match args.selection {
        Selection::OnlyRefs => prune_by_gt(args, contig, variant_pos, &indexes, samples, 0)?,
        Selection::OnlyAlts => prune_by_gt(args, contig, variant_pos, &indexes, samples, 1)?,
        _ => (indexes.iter().map(|_| [true, true]).collect(), samples),
    };

    tracing::info!("Input VCF: {:?}", args.file);
    tracing::info!("Reading phased genotypes from {contig} with target position at {variant_pos}.");
    tracing::info!("Using {} samples from the VCF.", indexes.len());

    let window = 10_000_000;

    if sharded {
        ensure!(get_htslib_contig_len(&args.file, contig).is_ok(), "Cannot construct HST using partial reads of the VCF if the contig length is not defined in the header");

        range = Some((
            Some(variant_pos.saturating_sub(window)),
            Some(variant_pos + window),
        ));
        tracing::info!(
            "Reading first batch in range {}-{}",
            variant_pos.saturating_sub(window),
            variant_pos + window
        );
    }

    // let (markers, coords) =
    // read_vcf_batch_to_matrix(&args.file, contig, range, &indexes, &lookups, args.no_alt)?;
    let (markers, coords) = match get_htslib_contig_len(&args.file, contig).is_err() {
        true => {
            tracing::info!("No contig {contig} length in the VCF. Reading single-threaded.");
            read_vcf_batch_to_matrix(&args.file, contig, range, &indexes, &lookups, args.no_alt)?
        }
        false => read_parallel(&args.file, contig, range, &indexes, &lookups, args.no_alt)?,
    };

    let start = coords.first().unwrap();
    let end = coords.last().unwrap();
    let fetch_range = (start.pos, end.pos);

    let metadata = ReadMetadata {
        indexes,
        lookups,
        file_path: args.file.clone(),
        fetch_range,
        contig_len: get_htslib_contig_len(&args.file, contig).ok(),
        sharded,
        remove_no_alt: args.no_alt,
    };

    construct_phased_matrix(
        samples,
        markers,
        coords,
        &args.selection,
        variant_pos,
        contig,
        metadata,
    )
}

fn read_vcf_batch_to_matrix(
    file_path: &PathBuf,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    sample_indexes: &[usize],
    lookups: &[[bool; 2]],
    remove_no_alt: bool,
) -> Result<(Vec<u8>, BTreeSet<Coord>)> {
    let mut reader = get_reader(file_path, contig, range)?;

    let mut coords = BTreeSet::new();
    let mut markers = vec![];

    let mut gt_buffer = rust_htslib::bcf::record::Buffer::new();

    for record in reader.records() {
        let record = record?;

        // HTSlib is 0-based so add 1
        let pos = (record.pos() + 1) as u64;

        // There is a weird bug where the IndexedReader fetches records outside of the range when multithreading
        if let Some((start, stop)) = range {
            if let (Some(start), Some(stop)) = (start, stop) {
                if pos < start || pos > stop {
                    continue;
                }
            }
        }

        tracing::trace!("Reading record at position {pos}");

        ensure!(record.alleles().len() == 2, NormalizeError(pos));

        let gts = record.genotypes_shared_buffer(&mut gt_buffer)?;

        let markers_len_before = markers.len();

        for (i, sample_idx) in sample_indexes.iter().enumerate() {
            let sample_gts = gts.get(*sample_idx);

            match (lookups[i][0], lookups[i][1]) {
                (true, true) => markers.extend([
                    genotype_to_u8(&sample_gts[0]),
                    genotype_to_u8(&sample_gts[1]),
                ]),
                (true, false) => {
                    markers.push(genotype_to_u8(&sample_gts[0]));
                }
                (false, true) => {
                    markers.push(genotype_to_u8(&sample_gts[1]));
                }
                (false, false) => (),
            }
        }

        let diff = markers.len() - markers_len_before;
        // let gts_per_sample = diff / sample_indexes.len();
        // let gts_mod = diff % sample_indexes.len();

        // tracing::trace!(
        // "gts {gts_per_sample} mod {gts_mod} before_len {markers_len_before} len_now {}",
        // markers.len()
        // );
        // ensure!(
        //     check_ploidy(gts_per_sample, gts_mod, Selection::All),
        //     PloidyError((pos, gts_per_sample))
        // );

        if remove_no_alt {
            let slice = &markers[(markers.len() - diff)..markers.len()];
            if slice.iter().all(|gt| gt == &0u8) {
                let _ = markers.drain((markers.len() - diff)..);
                continue;
            }
        }

        if !coords.insert(construct_coord(&record, contig, pos)) {
            tracing::error!("Duplicate coord at {contig}:{pos}. Not adding the duplicate");
            let _ = markers.drain((markers.len() - diff)..);
        }
    }

    Ok((markers, coords))
}

fn read_parallel(
    file_path: &PathBuf,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    sample_indexes: &[usize],
    lookups: &[[bool; 2]],
    remove_no_alt: bool,
) -> Result<(Vec<u8>, BTreeSet<Coord>)> {
    let (first_pos, last_pos) = match range {
        Some((Some(start), Some(end))) => (start, end),
        _ => (0, get_htslib_contig_len(file_path, contig)?),
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
                file_path,
                contig,
                Some((Some(*start), Some(*stop))),
                sample_indexes,
                lookups,
                remove_no_alt,
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
    Ok((markers, coords))
}

#[allow(clippy::too_many_arguments)]
fn construct_phased_matrix(
    samples: Vec<String>,
    markers: Vec<u8>,
    coords: BTreeSet<Coord>,
    selection: &Selection,
    variant_pos: u64,
    contig: &str,
    metadata: ReadMetadata,
) -> Result<PhasedMatrix> {
    let ploidy: Ploidy = selection.as_ref().into();
    tracing::debug!(
        "Contructing matrix: {} x {} = {}",
        samples.len() * *ploidy,
        coords.len(),
        markers.len()
    );
    let matrix = Array2::from_shape_vec((samples.len() * *ploidy, coords.len()).f(), markers)?;

    let mut vcf = PhasedMatrix::new(
        0,
        Coord::default(),
        matrix,
        samples,
        coords,
        selection,
        metadata,
    );

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

fn check_ploidy(length: usize, gts_mod: usize, selection: Selection) -> bool {
    if gts_mod != 0 {
        return false;
    }

    match selection {
        Selection::Haploid => length == 1 || length == 2,
        Selection::OnlyAlts => length == 1,
        Selection::OnlyRefs => length == 1,
        _ => length == 2,
    }
}

// fn genotype_to_u8_missing(g: &GenotypeAllele) -> Option<u8> {
//     match g {
//         GenotypeAllele::Unphased(v) => Some(*v as u8),
//         GenotypeAllele::Phased(v) => Some(*v as u8),
//         GenotypeAllele::UnphasedMissing => None,
//         GenotypeAllele::PhasedMissing => None,
//     }
// }

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

pub fn read_shard_of_vcf(vcf: &mut PhasedMatrix, start: u64, stop: u64) -> Result<()> {
    let (indexes, lookups) = (&vcf.metadata.indexes, &vcf.metadata.lookups);

    let range = Some((Some(start), Some(stop)));

    tracing::info!(
        "Reading more from: {:?}, range: {:?}",
        vcf.metadata.file_path,
        (start, stop),
    );

    let (markers, coords) =
        match get_htslib_contig_len(&vcf.metadata.file_path, &vcf.start_coord.contig).is_err() {
            true => read_vcf_batch_to_matrix(
                &vcf.metadata.file_path,
                &vcf.start_coord.contig,
                range,
                indexes,
                lookups,
                vcf.metadata.remove_no_alt,
            )?,
            false => read_parallel(
                &vcf.metadata.file_path,
                &vcf.start_coord.contig,
                range,
                indexes,
                lookups,
                vcf.metadata.remove_no_alt,
            )?,
        };

    if start < vcf.metadata.fetch_range.0 {
        vcf.metadata.fetch_range.0 = start;
    }

    if stop > vcf.metadata.fetch_range.1 {
        vcf.metadata.fetch_range.1 = stop;
    }

    let matrix = Array2::from_shape_vec(
        (vcf.samples().len() * *vcf.ploidy, coords.len()).f(),
        markers,
    )?;

    let start = coords.first().unwrap();

    vcf.insert_matrix(start.clone(), matrix);

    vcf.coords_mut().extend(coords);
    Ok(())
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
