use std::sync::Arc;
use std::path::PathBuf;

use color_eyre::Result;
use color_eyre::eyre::eyre;

use crate::args::Selection;
use crate::core::PhasedMatrix;
use crate::io::push_to_output;
use crate::io::open_csv_writer;
use crate::read_vcf::read_vcf_to_matrix;
use crate::core::Coord;
use crate::core::HapVariant;
use crate::read_vcf::ReadConfiguration;
use crate::utils::precision_f64;
use crate::traits::OnlyLongest;

#[derive(Debug, Default, Clone, PartialEq)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct Args {
    pub file: PathBuf,

    /// Haplotype for checking
    #[cfg_attr(feature = "clap", arg(long))]
    pub haplotype: PathBuf,

    /// Output directory
    #[cfg_attr(feature = "clap", arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("./"), value_hint = clap::ValueHint::DirPath))]
    pub output: PathBuf,

    /// List of samples for HST construction (one ID per row)
    #[cfg_attr(feature = "clap", arg(short = 'S', long, value_delimiter = ' ', num_args = 1.. ))]
    pub samples: Option<Vec<PathBuf>>,

    #[cfg_attr(feature = "clap", arg(short = 'a', long = "alleles", value_enum, default_value_t = Selection::All))]
    pub selection: Selection,

    /// Output filename prefix
    #[cfg_attr(feature = "clap", arg(short = 'p', long))]
    pub prefix: Option<String>,

    /// Do not include no ALTs
    #[cfg_attr(feature = "clap", arg(long))]
    pub no_alt: bool,

    /// Include only SNVs
    #[cfg_attr(feature = "clap", arg(long))]
    pub include_indels: bool,

    /// List of phase sets / haplotypes to include per sample
    #[cfg_attr(feature = "clap", arg(long))]
    pub list: Option<PathBuf>,
}

impl ReadConfiguration for &Args {
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

#[doc(hidden)]
pub fn run(args: Args) -> Result<()> {
    let mut csv_output = args.output.clone();
    push_to_output(&args.prefix, args.selection, &mut csv_output, "haplotype_check", "csv");
    let mut writer = open_csv_writer(csv_output)?;

    let ht = crate::io::read_haplotype_file(args.haplotype.clone())?;
    let start: &HapVariant = ht.first().ok_or_else(|| eyre!( 
        "Failed to get the first variant of the haplotype at {:?}. Is the haplotype file empty?", args.haplotype),
    )?;

    let end = ht.last().unwrap();

    let mut vcf = read_vcf_to_matrix(
        &args,
        &start.contig,
        start.pos,
        Some((Some(start.pos), Some(end.pos))),
        None,
        None,
        false,
    )?;

    match args.selection {
        Selection::OnlyLongest => vcf.select_only_longest_no_shard()?,
        Selection::Unphased => return Err(eyre!("Running with unphased data is not supported")),
        _ => (),
    };

    let matching_indexes = identical_haplotype_count(&vcf, &ht);
    write_matches_to_csv(&matching_indexes, &mut writer, &vcf)?;

    tracing::debug!(
        "Haplotype found in {}/{} (freq: {}) of alleles.",
        matching_indexes.len(),
        vcf.nhaplotypes(),
        precision_f64(matching_indexes.len() as f64 / vcf.nhaplotypes() as f64, 3),
    );

    Ok(())
}

fn write_matches_to_csv(
    matching_indexes: &[usize],
    writer: &mut csv::Writer<Box<dyn std::io::Write>>,
    vcf: &PhasedMatrix,
) -> Result<()> {
    writer.write_record(vec!["id", "match"])?;

    for idx in 0..vcf.nhaplotypes() {
        writer.write_record(vec![
            format!("{}", vcf.get_sample_name(idx)),
            format!("{}", matching_indexes.contains(&idx)),
        ])?;
    }

    Ok(())
}

pub fn identical_haplotype_count(vcf: &PhasedMatrix, ht: &[HapVariant]) -> Vec<usize> {
    let indexes: Vec<(&(Arc<Coord>, usize), usize)> = ht
        .iter()
        .enumerate()
        .filter_map(|(ht_idx, h)|
            match vcf.indexer.get(&h.clone().into()) {
                Some(matrix_idx) => Some((matrix_idx, ht_idx)),
                None => {
                        tracing::warn!("Haplotype variant {:?} not found in the vcf", h);
                        None
                },
            }
        )
        .collect();

    if indexes.is_empty() {
        return vec![];
    }    

    (0..vcf.nhaplotypes()).filter(|sample_idx| {
        let mut no_mismatch = true;
        for ((key_coord, var_idx), ht_idx) in &indexes {
            let matrix = vcf.matrix.get(key_coord).unwrap();

            if matrix.genotype(*sample_idx, *var_idx) != ht[*ht_idx].gt {
                no_mismatch = false;
                break;
            }
        }
        no_mismatch
    }
    )
    .collect()
}
