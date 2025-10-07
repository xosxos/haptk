use std::collections::BTreeSet;
use std::path::Path;

use ndarray::{Array2, ShapeBuilder};
use rayon::prelude::*;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;

use crate::error::Error;
use crate::phased_matrix::PhasedMatrix;
use crate::phased_matrix::ReadMetadata;
use crate::ploidy::Ploidy;
use crate::variant::Coord;
use crate::vcf::contig_len_from_vcf;
use crate::vcf::get_indexes_and_sample_ids_from_vcf;

pub fn get_vcf_reader(
    path: &Path,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
) -> Result<IndexedReader, Error> {
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

pub fn read_phased_matrix(
    path: &Path,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    samples: Option<Vec<String>>,
) -> Result<PhasedMatrix, Error> {
    let (indexes, samples) = get_indexes_and_sample_ids_from_vcf(path, &None, samples)?;

    let lookups = indexes.iter().map(|_| [true, true]).collect();

    let ploidy: Ploidy = Ploidy::Diploid;
    let variant_pos = 0;
    let window = None;
    let remove_no_alt = false;
    let include_indels = false;
    let is_genome_wide = false;

    read_vcf_to_matrix_by_indexes(
        path,
        variant_pos,
        contig,
        range,
        samples,
        indexes,
        lookups,
        window,
        remove_no_alt,
        ploidy,
        is_genome_wide,
        include_indels,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn read_vcf_to_matrix_by_indexes(
    file: &Path,
    variant_pos: u64,
    contig: &str,
    mut range: Option<(Option<u64>, Option<u64>)>,
    samples: Vec<String>,
    indexes: Vec<usize>,
    lookups: Vec<[bool; 2]>,
    window: Option<u64>,
    remove_no_alt: bool,
    ploidy: Ploidy,
    is_genome_wide: bool,
    include_indels: bool,
) -> Result<PhasedMatrix, Error> {
    tracing::info!("Input VCF: {:?}", file);
    tracing::info!("Reading phased genotypes from {contig} with target position at {variant_pos}.");
    tracing::info!("Using {} samples from the VCF.", indexes.len());

    if let Some(window) = window {
        assert!(variant_pos.checked_add(window).is_some(), "Error: Variant position: {variant_pos} added to the window: {window} is larger than the largest allowed 64-bit unsigned integer {}", u64::MAX);
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

    let (markers, coords) = match contig_len_from_vcf(file, contig).is_err() {
        true => {
            tracing::info!("No contig {contig} length in the VCF. Reading single-threaded.");
            read_vcf_batch_to_matrix(
                file,
                contig,
                range,
                &indexes,
                &lookups,
                remove_no_alt,
                include_indels,
            )?
        }
        false => read_parallel(
            file,
            contig,
            range,
            &indexes,
            &lookups,
            remove_no_alt,
            include_indels,
        )?,
    };

    let start = coords.first().unwrap();
    let end = coords.last().unwrap();
    let fetch_range = (start.pos, end.pos);

    let metadata = ReadMetadata {
        indexes,
        lookups,
        file_path: file.to_path_buf(),
        fetch_range,
        contig_len: contig_len_from_vcf(file, contig).ok(),
        sharded: window.is_some(),
        remove_no_alt,
        include_indels,
        is_genome_wide,
    };

    construct_phased_matrix(samples, markers, coords, ploidy, variant_pos, metadata)
}

fn read_vcf_batch_to_matrix(
    file_path: &Path,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    sample_indexes: &[usize],
    lookups: &[[bool; 2]],
    remove_no_alt: bool,
    include_indels: bool,
) -> Result<(Vec<u8>, BTreeSet<Coord>), Error> {
    let mut reader = get_vcf_reader(file_path, contig, range)?;

    let mut coords = BTreeSet::new();
    let mut markers = vec![];

    let mut gt_buffer = rust_htslib::bcf::record::Buffer::new();

    let mut prev_pos = 0;

    for record in reader.records() {
        let record = record?;

        // HTSlib is 0-based so add 1
        let pos = (record.pos() + 1) as u64;

        // Check that there are no multiallelics, else return error
        if record.alleles().len() != 2 {
            return Err(Error::Normalize { pos });
        }

        // Build the `Coord` struct
        let alleles = record.alleles();
        let reference = String::from_utf8_lossy(alleles.first().unwrap()).to_string();
        let alt = String::from_utf8_lossy(alleles.get(1).unwrap()).to_string();

        let coord = Coord {
            contig: contig.to_string(),
            reference,
            alt,
            pos,
        };

        // Check that the current is larger than or equals the previous position
        // Else return OrderError

        if pos < prev_pos {
            return Err(Error::Order {
                prev_pos,
                pos,
                coord: coord.to_string(),
            });
        }

        if !include_indels && (coord.alt.len() > 1 || coord.reference.len() > 1) {
            tracing::warn!("Only SNVs wanted, disregarding: {coord}");
            continue;
        }

        // Check that the current coordinate is not a duplicate
        // If it is, disregard and continue to the next record
        if !coords.insert(coord) {
            tracing::warn!("Duplicate coord at {contig}:{pos}. Not adding the duplicate");
            continue;
        }

        // FIX: There is a weird bug where the IndexedReader fetches records outside of the range when multithreading
        if let Some((start, stop)) = range {
            if let (Some(start), Some(stop)) = (start, stop) {
                if pos < start || pos > stop {
                    continue;
                }
            }
        }

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

        if remove_no_alt {
            let slice = &markers[(markers.len() - diff)..markers.len()];

            let mut iter = slice.iter();
            let first = iter.next().unwrap();

            // If all alleles are the same, then skip this variant
            if iter.all(|gt| gt == first) {
                let _ = markers.drain((markers.len() - diff)..);
                continue;
            }
        }

        prev_pos = pos;
    }

    Ok((markers, coords))
}

fn read_parallel(
    file_path: &Path,
    contig: &str,
    range: Option<(Option<u64>, Option<u64>)>,
    sample_indexes: &[usize],
    lookups: &[[bool; 2]],
    remove_no_alt: bool,
    include_indels: bool,
) -> Result<(Vec<u8>, BTreeSet<Coord>), Error> {
    let (first_pos, last_pos) = match range {
        Some((Some(start), Some(end))) => (start, end),
        _ => (0, contig_len_from_vcf(file_path, contig)?),
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
                include_indels,
            )
        })
        .collect::<Result<Vec<(Vec<u8>, BTreeSet<Coord>)>, Error>>()?;

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
    ploidy: Ploidy,
    variant_pos: u64,
    metadata: ReadMetadata,
) -> Result<PhasedMatrix, Error> {
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
        ploidy,
        metadata,
    );

    // vcf.variant_idx = vcf.get_first_idx_on_right_by_pos(variant_pos);
    vcf.start_coord = vcf.get_nearest_coord_by_pos(variant_pos).clone();
    vcf.variant_idx = vcf.get_coord_idx(&vcf.start_coord);

    tracing::info!(
        "Constructed a genotype matrix from {} records. Starting variant position: {}.",
        vcf.ncoords(),
        vcf.variant_idx_pos(),
    );
    Ok(vcf)
}

pub fn read_shard_of_vcf(vcf: &mut PhasedMatrix, start: u64, stop: u64) -> Result<(), Error> {
    let (indexes, lookups) = (&vcf.metadata.indexes, &vcf.metadata.lookups);

    let range = Some((Some(start), Some(stop)));

    tracing::info!(
        "Reading more from: {:?}, range: {:?}",
        vcf.metadata.file_path,
        (start, stop),
    );

    let (markers, coords) =
        match contig_len_from_vcf(&vcf.metadata.file_path, &vcf.start_coord.contig).is_err() {
            true => read_vcf_batch_to_matrix(
                &vcf.metadata.file_path,
                &vcf.start_coord.contig,
                range,
                indexes,
                lookups,
                vcf.metadata.remove_no_alt,
                vcf.metadata.include_indels,
            )?,
            false => read_parallel(
                &vcf.metadata.file_path,
                &vcf.start_coord.contig,
                range,
                indexes,
                lookups,
                vcf.metadata.remove_no_alt,
                vcf.metadata.include_indels,
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

    let start = coords.first().unwrap().clone();

    vcf.insert_matrix(start.clone(), matrix);
    tracing::debug!("Finished extending matrices.");

    let start = std::sync::Arc::new(start);

    coords.into_iter().enumerate().for_each(|(i, v)| {
        vcf.coords_mut().insert(v.clone());
        vcf.indexer.insert(v, (start.clone(), i));
    });
    tracing::debug!("Finished extending indexer and coords.");

    Ok(())
}

pub fn genotype_to_u8(g: &GenotypeAllele) -> u8 {
    match g {
        GenotypeAllele::Unphased(v) => *v as u8,
        GenotypeAllele::Phased(v) => *v as u8,
        _ => panic!(),
    }
}
