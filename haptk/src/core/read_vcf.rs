use std::path::Path;
use std::path::PathBuf;

use color_eyre::eyre::ensure;
use color_eyre::eyre::eyre;
use color_eyre::{eyre::OptionExt, Result};

use haptk_core::phased_matrix::SelectedHaplotypes;
use rust_htslib::bcf::Read;

use crate::args::Selection;
use crate::args::StandardArgs;
use crate::core::PhasedMatrix;
use crate::core::Ploidy;
use crate::error::Error;
use crate::io::get_indexes_and_sample_ids_from_vcf;
use crate::io::read_sample_ht_list_file;

pub use haptk_core::read_phased_matrix::genotype_to_u8;
pub use haptk_core::read_phased_matrix::get_vcf_reader;
pub use haptk_core::read_phased_matrix::read_shard_of_vcf;
pub use haptk_core::read_phased_matrix::read_vcf_to_matrix_by_indexes;

fn prune_by_gt(
    vcf_path: &Path,
    contig: &str,
    variant_pos: u64,
    indexes: &[usize],
    samples: Vec<String>,
    wanted_gt: u8,
) -> Result<(SelectedHaplotypes, Vec<String>)> {
    let mut reader = get_vcf_reader(
        vcf_path,
        contig,
        Some((Some(variant_pos), Some(variant_pos))),
    )?;

    let mut gt_buffer = rust_htslib::bcf::record::Buffer::new();

    let record = reader.records().next();

    // Check whether the variant exists
    let record = record.ok_or_eyre(Error::NonExistantVariant {
        contig: contig.to_string(),
        variant_pos,
    })??;

    let gts = record.genotypes_shared_buffer(&mut gt_buffer)?;

    let lookups: Vec<[bool; 2]> = indexes
        .iter()
        .map(|i| {
            let sample_gts = gts.get(*i);

            let gt1 = genotype_to_u8(&sample_gts[0]);
            let gt2 = genotype_to_u8(&sample_gts[1]);
            [gt1 == wanted_gt, gt2 == wanted_gt]
        })
        .collect();

    let lookups = SelectedHaplotypes::from_vec(lookups);

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

pub trait ReadConfiguration {
    fn file(&self) -> &PathBuf;
    fn samples(&self) -> &Option<Vec<PathBuf>>;
    fn selection(&self) -> &Selection;
    fn no_alt(&self) -> bool;
    fn include_indels(&self) -> bool;
    fn list(&self) -> &Option<PathBuf>;
}

impl ReadConfiguration for &StandardArgs {
    fn file(&self) -> &PathBuf {
        &self.file
    }

    fn samples(&self) -> &Option<Vec<PathBuf>> {
        &self.samples
    }

    fn selection(&self) -> &Selection {
        &self.selection
    }

    fn no_alt(&self) -> bool {
        self.no_alt
    }

    fn include_indels(&self) -> bool {
        self.include_indels
    }

    fn list(&self) -> &Option<PathBuf> {
        &self.list
    }
}

// TODO: Change to builder pattern to refactor out this ugly ReadConfiguration interface
pub fn read_vcf_to_matrix<T: ReadConfiguration>(
    args: T,
    contig: &str,
    variant_pos: u64,
    range: Option<(Option<u64>, Option<u64>)>,
    wanted_samples: Option<Vec<String>>,
    window: Option<u64>,
    is_genome_wide: bool,
) -> Result<PhasedMatrix> {
    let (mut indexes, samples) =
        get_indexes_and_sample_ids_from_vcf(args.file(), args.samples(), wanted_samples)?;

    let (lookups, samples) = match args.selection() {
        Selection::OnlyRefs => prune_by_gt(args.file(), contig, variant_pos, &indexes, samples, 0)?,
        Selection::OnlyAlts => prune_by_gt(args.file(), contig, variant_pos, &indexes, samples, 1)?,
        // TODO: New ugly code
        Selection::List => {
            ensure!(
                args.list().is_some(),
                eyre!("List selection enabled, but no list was given with the --list parameter")
            );
            let lookups = read_sample_ht_list_file(&args.list().clone().unwrap())?;

            indexes = samples
                .iter()
                .zip(indexes.iter())
                .filter(|(v, _i)| lookups.contains_key(*v))
                .map(|(_, i)| i)
                .cloned()
                .collect();

            (
                SelectedHaplotypes::from_vec(
                    samples
                        .iter()
                        .flat_map(|v| {
                            // if !lookups.contains_key(v) {
                            // tracing::warn!("wanted sample {} was not in the list file", v);
                            // }
                            lookups.get(v)
                        })
                        .copied()
                        .collect(),
                ),
                samples
                    .iter()
                    .filter(|&v| lookups.contains_key(v))
                    .map(|v| (v, lookups.get(v).unwrap()))
                    .flat_map(|(s, lookup)| match (lookup[0], lookup[1]) {
                        (true, true) => vec![s.clone(), s.clone()],
                        (true, false) | (false, true) => vec![s.clone()],
                        (false, false) => vec![],
                    })
                    .collect(),
            )
        }
        _ => (SelectedHaplotypes::select_all(&indexes), samples),
    };

    // The sample name to idx to ht_num system, is a catastrophy right now
    // This assert is a first step in the direction of fixing it
    // Now the sample names vector should be a 1-to-1 match with the the index vector
    // But atm it's not, because I thought being tricky with indexing using ploidy was smart. it was not.
    if !matches!(args.selection(), &Selection::All)
        && !matches!(args.selection(), &Selection::OnlyLongest)
    {
        assert_eq!(
            lookups.iter().flatten().filter(|v| **v == true).count(),
            samples.len(),
            "Found a bug while selecting haplotypes from the VCF file, please report"
        );
    }

    let ploidy: Ploidy = args.selection().as_ref().into();

    Ok(read_vcf_to_matrix_by_indexes(
        args.file(),
        variant_pos,
        contig,
        range,
        samples,
        indexes,
        lookups,
        window,
        args.no_alt(),
        ploidy,
        is_genome_wide,
        args.include_indels(),
    )?)
}
