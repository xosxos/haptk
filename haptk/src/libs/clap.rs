use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
#[allow(unused_imports)]
use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};

use crate::args::{GraphArgs, SortOption, StandardArgs};
use crate::subcommands::{
    bhst, check_for_haplotype, compare_haplotypes, compare_to_haplotype, compare_to_hst,
    fasta_to_haplotype, haplotype_to_vcf, list_haplotypes, list_markers, list_samples, mrca, uhst,
};

// Genome-wide methods
#[cfg(feature = "experimental")]
use crate::{
    args::ConciseArgs,
    subcommands::{
        mrca_scan,
        scan::{
            hst_scan, scan_branch_mrca, scan_nodes, scan_quantitative, scan_segregate, scan_sum_hsts, scan_annotate,
        }
    }
};

#[derive(Parser, Debug)]
#[command(author, version, about, styles=get_styles())]
pub struct Arguments {
    #[command(subcommand)]
    pub cmd: SubCommand,
}

#[derive(Args, Debug, Clone)]
pub struct LogAndVerbosity {
    /// Verbosity level
    #[arg(short, long, action = clap::ArgAction::Count, default_value_t = 3)]
    pub verbosity: u8,

    /// A file path to save logs to
    #[arg(short, long)]
    pub log_file: Option<PathBuf>,

    /// Silence all warning and info messages
    #[arg(long)]
    pub silent: bool,
}

#[derive(Subcommand, Debug)]
pub enum SubCommand {
    /// Build unidirectional haplotype sharing trees at a coordinate
    Uhst {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Min amount of samples per node
        #[arg(short = 'm', long, default_value_t = 1)]
        min_size: usize,

        /// Remove sample IDs and hard cut out all nodes with samples less than or equal to `min_size`
        #[arg(long)]
        publish: bool,

        /// If you most likely dont need to read the entire VCF for a HST use a window e.g. 5000000 or 2000000 (5Mb or 2Mb) around the locus
        #[arg(long)]
        window: Option<u64>,
    },

    /// Build a bidirectional haplotype sharing tree at a coordinate
    Bhst {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Min amount of samples per node
        #[arg(short = 'm', long, default_value_t = 1)]
        min_size: usize,

        /// Remove sample IDs and hard cut out all nodes with samples less than or equal to `min_size`
        #[arg(long)]
        publish: bool,

        /// If you most likely dont need to read the entire VCF for a HST use a window e.g. 5000000 or 2000000 (5Mb or 2Mb) around the locus
        #[arg(long)]
        window: Option<u64>,
    },
    /// Analyze the MRCA based on the Gamma method at a coordinate
    Mrca {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Recombination rate file
        #[arg(short = 'r', long)]
        recombination_rates: PathBuf,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// If you most likely dont need to read the entire VCF for a HST use a window e.g. 5000000 or 2000000 (5Mb or 2Mb) around the locus
        #[arg(long)]
        window: Option<u64>,
    },

    ///  (experimental) Analyze the MRCA every x markers along a given contig
    #[cfg(feature = "experimental")]
    MrcaScan {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Recombination rate file/files
        #[arg(short = 'r', long, value_delimiter = ' ', num_args = 1.. )]
        recombination_rates: Vec<PathBuf>,

        /// Run the MRCA analysis every n markers
        #[arg(long)]
        step_size: usize,

        /// Dont output to csv
        #[arg(long)]
        no_csv: bool,

        /// List of samples for HST construction (one ID per row)
        #[arg(long, value_delimiter = ',', num_args = 1.. )]
        contigs: Option<Vec<String>>,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// The percentage of haplotypes that can overlap centromeres before tagging the position "centromeric"
        #[arg(long, default_value_t = 0.1)]
        centromere_cut_off: f32,
    },

    /// Check if samples share a given haplotype
    CheckForHaplotype {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Haplotype for checking
        #[arg(short = 'h', long)]
        haplotype: PathBuf,
    },
    /// Check differences between samples and a haplotype
    CompareToHaplotype {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Haplotype for comparing to
        #[arg(short = 'H', long)]
        haplotype: PathBuf,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Mark sample names for tagging
        #[arg(short = 'm', long)]
        mark_samples: Option<PathBuf>,

        /// Mark shorter alleles to graph
        #[arg(long)]
        mark_shorter_alleles: bool,

        #[command(flatten)]
        graph_args: GraphArgs,

        /// Output .png image
        #[arg(long)]
        png: bool,

        /// Output .npy matrix
        #[arg(long)]
        npy: bool,

        #[arg(long = "sort", value_enum, default_value_t=SortOption::Left)]
        sort_option: SortOption,
    },

    /// Check which haplotypes of the HST are present in samples
    CompareToHst {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// A HST to compare to
        #[arg(long)]
        hst: PathBuf,

        /// Filter out every match for each sample, except the longest match measured in basepairs
        #[arg(long)]
        only_longest_leafs: bool,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,
    },

    /// Compare haplotypes to each other by alignment
    CompareHaplotypes {
        haplotypes: Vec<PathBuf>,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Output directory
        #[arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("./"))]
        output: PathBuf,

        /// Output filename prefix
        #[arg(long)]
        prefix: Option<String>,

        /// Output to csv instead of stdout
        #[arg(long)]
        csv: bool,

        /// Hide missing (meaning yellow) variants from the stdout print
        #[arg(long)]
        hide_missing: bool,

        /// Tag with ok for matching, err for mismatching and mis for mismatching to null rows
        #[arg(long)]
        tag_rows: bool,

        #[arg(short = 'n', long)]
        nucleotides: bool,
    },
    /// Read the haplotypes of a given sample
    Haplotypes {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        #[arg(long)]
        selection_variant: Option<String>,

        #[arg(short = 'n', long)]
        nucleotides: bool,
    },
    /// Output the sample names from FAM / VCF / HST files
    Samples {
        file: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
    },
    /// Output the markers from a HST file
    Markers {
        file: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
    },

    /// Convert a haplotype CSV into VCF
    HaplotypeToVcf {
        file: PathBuf,

        #[arg(short = 's', long)]
        sample_name: String,

        #[arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("-"))]
        output: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
    },

    /// Convert fasta sequences to haplotype csv format
    FastaToHaplotype {
        file: PathBuf,

        /// Fasta sequences to transform into a haplotpye
        #[arg(short = 'S', long, value_delimiter = ' ', num_args = 1.. )]
        seq_name: Vec<String>,

        #[arg(short = 'o', long="outdir", default_value_os_t = PathBuf::from("-"))]
        output: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
    },

    #[cfg(feature = "experimental")]
    /// (experimental) HST scan
    HstScan {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Run the scan every n markers
        #[arg(long)]
        step_size: usize,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 1)]
        min_sample_size: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,
    },

    #[cfg(feature = "experimental")]
    /// (experimental) Scan all nodes of all trees for a node specific value
    ScanNodes {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 4)]
        min_sample_size: usize,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 10000000)]
        max_sample_size: usize,

        /// Minimum number of variants required in the haplotype
        #[arg(long, default_value_t = 1)]
        min_ht_len: usize,

        /// Maximum number of variants in the haplotype
        #[arg(long, default_value_t = 10000000)]
        max_ht_len: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Construct the HSTs and sum the leaf nodes simultaneously
        #[arg(long)]
        construct_hsts: bool,

        /// The starting coordinate, i.e. chr9:27573534
        #[arg(short = 'c', long)]
        coords: Option<String>,

        /// If constructing HSTs at the same time, scan every n markers
        #[arg(long, default_value_t = 1)]
        step_size: usize,
    },

    #[cfg(feature = "experimental")]
    /// (experimental) Scan all nodes of all trees for a difference in a quantative variable
    ScanQuantitative {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 4)]
        min_sample_size: usize,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 10000000)]
        max_sample_size: usize,

        /// Minimum number of variants required in the haplotype
        #[arg(long, default_value_t = 1)]
        min_ht_len: usize,

        /// Maximum number of variants in the haplotype
        #[arg(long, default_value_t = 10000000)]
        max_ht_len: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// List of control ids
        #[arg(long)]
        var_data: PathBuf,

        /// Path to controls vcf if controls are not included in the same file
        #[arg(long)]
        var_name: String,

        /// Coords to generate HSTs for
        #[arg(long)]
        coords: Option<PathBuf>,      
    },

    #[cfg(feature = "experimental")]
    /// (experimental) Scan all branches of all trees for the smallest MRCA value
    ScanBranchMrca {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 10)]
        min_sample_size: usize,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 10000000)]
        max_sample_size: usize,

        /// Minimum number of variants required in the haplotype
        #[arg(long, default_value_t = 1)]
        min_ht_len: usize,

        /// Maximum number of variants in the haplotype
        #[arg(long, default_value_t = 10000000)]
        max_ht_len: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Recombination rate file
        #[arg(short = 'r', long)]
        recombination_rates: PathBuf,

        /// Construct the HSTs and sum the leaf nodes simultaneously
        #[arg(long)]
        construct_hsts: bool,

        /// The starting coordinate, i.e. chr9:27573534
        #[arg(short = 'c', long)]
        coords: Option<String>,

        /// If constructing HSTs at the same time, scan every n markers
        #[arg(long, default_value_t = 1)]
        step_size: usize,

        /// Generations limit
        #[arg(long, default_value_t = 40.0)]
        limit: f64,
    },

    #[cfg(feature = "experimental")]
    /// (experimental) Find bHST based segregated haplotypes genome-wide
    ScanSegregate {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 1)]
        min_sample_size: usize,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 10000000)]
        max_sample_size: usize,

        /// Minimum number of variants required in the haplotype
        #[arg(long, default_value_t = 1)]
        min_ht_len: usize,

        /// Maximum number of variants in the haplotype
        #[arg(long, default_value_t = 10000000)]
        max_ht_len: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// Case samples
        #[arg(long)]
        case_samples: PathBuf,

        /// Ctrl samples
        #[arg(long)]
        ctrl_samples: PathBuf,

        /// Coords to generate HSTs for
        #[arg(long)]
        coords: Option<PathBuf>,      

        /// Generations limit
        #[arg(long, default_value_t = 0.005)]
        limit: f64,
    },

    #[cfg(feature = "experimental")]
    /// (experimental) Scan all leaf nodes to quantify haplotype sharing
    ScanSumHst {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        /// The percentage of haplotypes that can overlap centromeres before tagging the position "centromeric"
        #[arg(long, default_value_t = 0.1)]
        centromere_cut_off: f32,

        /// Recombination rate file
        #[arg(short = 'r', long)]
        recombination_rates: PathBuf,

        /// Construct the HSTs and sum the leaf nodes simultaneously
        #[arg(long)]
        construct_hsts: bool,

        /// The starting coordinate, i.e. chr9:27573534
        #[arg(short = 'c', long)]
        coords: Option<String>,

        /// If constructing HSTs at the same time, scan every n markers
        #[arg(long, default_value_t = 1)]
        step_size: usize,

        /// If constructing HSTs at the same time, branch sample size to end recursion
        #[arg(long, default_value_t = 1)]
        min_sample_size: usize,

        /// If constructing HSTs at the same time, do not include non-contradictory genotypes
        #[arg(long)]
        no_alt: bool,

        /// Length in basepairs instead of centimorgans
        #[arg(long)]
        length_in_bp: bool,

        /// List of samples for HST construction (one ID per row)
        #[arg(long, value_delimiter = ' ', num_args = 1.. )]
        seg_samples: Option<Vec<PathBuf>>,

    },
    
    #[cfg(feature = "experimental")]
    /// (experimental) Annotate a scan using a bed file
    AnnotateScan {
        file: PathBuf,

        annotate_file: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
    },

}


impl SubCommand {
    pub fn threads(&self) -> usize {
        match self {
            SubCommand::CompareToHaplotype { threads, .. }
            | SubCommand::CompareToHst { threads, .. }
            | SubCommand::Mrca { threads, .. }
            | SubCommand::Uhst { threads, .. }
            | SubCommand::Bhst { threads, .. } => *threads,

            #[cfg(feature = "experimental")]
            SubCommand::MrcaScan { threads, .. }
            | SubCommand::HstScan { threads, .. }
            | SubCommand::ScanSegregate { threads, .. }
            | SubCommand::ScanBranchMrca { threads, .. }
            | SubCommand::ScanQuantitative { threads, .. }
            | SubCommand::ScanSumHst { threads, .. }
            | SubCommand::ScanNodes { threads, .. } => *threads,

            _ => 1,
        }
    }

    #[rustfmt::skip]
    pub fn log_and_verbosity(&self) -> (u8, &Option<PathBuf>, bool) {
        match self {
            SubCommand::Haplotypes { log_and_verbosity, .. }
            | SubCommand::CheckForHaplotype { log_and_verbosity, .. }
            | SubCommand::CompareToHaplotype { log_and_verbosity, .. }
            | SubCommand::CompareToHst { log_and_verbosity, .. }
            | SubCommand::CompareHaplotypes { log_and_verbosity, .. }
            | SubCommand::Mrca { log_and_verbosity, .. }
            | SubCommand::Uhst { log_and_verbosity, .. }
            | SubCommand::Bhst { log_and_verbosity, .. }
            | SubCommand::Samples { log_and_verbosity, .. }
            | SubCommand::Markers { log_and_verbosity, .. }
            | SubCommand::HaplotypeToVcf { log_and_verbosity, .. }
            | SubCommand::FastaToHaplotype { log_and_verbosity, .. } 
            => (log_and_verbosity.verbosity, &log_and_verbosity.log_file, log_and_verbosity.silent),

            #[cfg(feature = "experimental")]
            SubCommand::MrcaScan { log_and_verbosity, .. }
            | SubCommand::HstScan { log_and_verbosity,  .. }
            | SubCommand::ScanSegregate { log_and_verbosity,  .. }
            | SubCommand::ScanBranchMrca { log_and_verbosity,  .. }
            | SubCommand::ScanQuantitative { log_and_verbosity,  .. }
            | SubCommand::ScanSumHst { log_and_verbosity, .. }
            | SubCommand::ScanNodes { log_and_verbosity,  .. }
            | SubCommand::AnnotateScan { log_and_verbosity, .. } 
            => (log_and_verbosity.verbosity, &log_and_verbosity.log_file, log_and_verbosity.silent),
        }
    }

    #[cfg(feature = "experimental")]
    pub fn check_sample_size(&self) -> Result<()> {
        match self {
            SubCommand::ScanNodes {
                min_sample_size, ..
            }
            | SubCommand::ScanQuantitative {
                min_sample_size, ..
            }
            | SubCommand::ScanBranchMrca {
                min_sample_size, ..
            } => {
                ensure!(
                    *min_sample_size >= 1,
                    "Min sample size needs to be atleast one"
                )
            }
            _ => return Ok(()),
        };
        Ok(())
    }

    #[rustfmt::skip]
    pub fn output(&self) -> Option<PathBuf> {
        match self {
            SubCommand::Haplotypes { args: StandardArgs { output, .. }, ..}
            | SubCommand::CheckForHaplotype { args: StandardArgs { output, .. }, ..}
            | SubCommand::CompareToHaplotype { args: StandardArgs { output, .. }, ..}
            | SubCommand::CompareToHst { args: StandardArgs { output, .. }, ..}
            | SubCommand::CompareHaplotypes { output, .. }
            | SubCommand::Mrca { args: StandardArgs { output, .. }, ..}
            | SubCommand::Uhst { args: StandardArgs { output, .. }, ..}
            | SubCommand::Bhst { args: StandardArgs { output, .. }, ..}
            => Some(output.clone()),

            SubCommand::Samples { .. }
            | SubCommand::HaplotypeToVcf { .. }
            | SubCommand::FastaToHaplotype { .. }
            | SubCommand::Markers { .. } => None,

            #[cfg(feature = "experimental")]
            SubCommand::AnnotateScan { .. } => None,

            #[cfg(feature = "experimental")]
            SubCommand::MrcaScan { args: StandardArgs { output, .. }, ..}
            | SubCommand::HstScan { args: StandardArgs { output, .. }, ..}
            | SubCommand::ScanSegregate { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanBranchMrca { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanQuantitative { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanSumHst { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanNodes { args: ConciseArgs { output, .. }, ..}
            => Some(output.clone()),

        }
    }
}

#[rustfmt::skip]
pub fn run_cmd(cmd: SubCommand) -> Result<()> {

    // Initialize threadpool
    rayon::ThreadPoolBuilder::new()
        .num_threads(cmd.threads())
        .build_global()?;

    // Initialize logger
    let (verbosity, log_file, is_silent) = cmd.log_and_verbosity();
    let _guard = init_tracing(verbosity, log_file, is_silent)?;

    // Create output directory
    if let Some(output) = cmd.output() {
        if let Err(e) = std::fs::create_dir(output.clone()) {
            match e.kind() {
                std::io::ErrorKind::AlreadyExists => (),
                _ => return Err(eyre!("Error creating directory {output:?}")),
            }
        }
    }

    match cmd {
        SubCommand::CompareToHaplotype { 
            args, haplotype, mark_samples, mark_shorter_alleles, graph_args, png, npy, sort_option, .. 
        } => compare_to_haplotype::run(
                args, haplotype, mark_samples, mark_shorter_alleles, png, npy, graph_args, sort_option,
            )?,

        SubCommand::Bhst { args,  min_size, publish, window, .. } => bhst::run(args, min_size, publish, window)?,
        SubCommand::Uhst {args,  min_size, publish, window, .. } => uhst::run(args, min_size, publish, window)?,

        SubCommand::CompareHaplotypes { haplotypes, output, prefix, csv, hide_missing, tag_rows, nucleotides, .. }
            => compare_haplotypes::run(haplotypes, output, prefix, csv, hide_missing, tag_rows, nucleotides)?,

        SubCommand::CompareToHst { args, hst, only_longest_leafs, .. }
            => compare_to_hst::run(args, hst, only_longest_leafs)?,

        SubCommand::CheckForHaplotype { args, haplotype, .. } => check_for_haplotype::run(args, haplotype)?,
        SubCommand::Mrca { args, recombination_rates, window, .. } => mrca::run(args, recombination_rates, window)?,
        SubCommand::Haplotypes { args, selection_variant, nucleotides, .. } => list_haplotypes::run(args, selection_variant, nucleotides)?,
        SubCommand::Samples { file, .. } => list_samples::run(file)?, 
        SubCommand::Markers { file, .. } => list_markers::run(file)?,
        SubCommand::HaplotypeToVcf { file, sample_name, output, .. } => haplotype_to_vcf::run(file, sample_name, output)?,
        SubCommand::FastaToHaplotype { file, seq_name, output, .. } => fasta_to_haplotype::run(file, seq_name, output)?,

        // Genome-wide methods
        #[cfg(feature = "experimental")]
        SubCommand::MrcaScan { args, recombination_rates, step_size, contigs, centromere_cut_off, .. }
        => mrca_scan::run(args, recombination_rates, step_size, contigs, centromere_cut_off)?,


        #[cfg(feature = "experimental")]
        SubCommand::HstScan { args, step_size, min_sample_size, .. } => hst_scan::run(args, step_size, min_sample_size)?,

        #[cfg(feature = "experimental")]
        SubCommand::ScanNodes { args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, construct_hsts, coords, step_size,  .. }
        => scan_nodes::run(args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), construct_hsts, coords, step_size)?,

        #[cfg(feature = "experimental")]
        SubCommand::ScanQuantitative {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, var_data, var_name, coords, ..
        } => scan_quantitative::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), var_data, var_name, coords
            )?,

        #[cfg(feature = "experimental")]
        SubCommand::ScanBranchMrca {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, recombination_rates, construct_hsts, coords, step_size, limit, ..
        } => scan_branch_mrca::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), recombination_rates, construct_hsts, coords, step_size, limit,
            )?,

        #[cfg(feature = "experimental")]
        SubCommand::ScanSumHst { args, recombination_rates, centromere_cut_off, construct_hsts, coords, step_size, min_sample_size, length_in_bp, seg_samples, .. }
        => scan_sum_hsts::run(
            args, recombination_rates, centromere_cut_off, construct_hsts, coords, step_size, min_sample_size, length_in_bp, seg_samples,
        )?,

        #[cfg(feature = "experimental")]
        SubCommand::ScanSegregate {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, case_samples, ctrl_samples, coords, limit, ..
        } => scan_segregate::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), case_samples, ctrl_samples, coords, limit,
            )?,
        #[cfg(feature = "experimental")]
        SubCommand::AnnotateScan { file, annotate_file, .. } => scan_annotate::run(file, annotate_file)?,

    };
    Ok(())
}

pub fn init_tracing(
    verbosity: u8,
    log_file: &Option<PathBuf>,
    is_silent: bool,
) -> Result<tracing_appender::non_blocking::WorkerGuard> {
    use tracing::Level;

    let level = if is_silent {
        Level::ERROR
    } else {
        match verbosity {
            0 => unreachable!(),
            1 => Level::ERROR,
            2 => Level::WARN,
            3 => Level::INFO,
            4 => Level::DEBUG,
            5..=u8::MAX => Level::TRACE,
        }
    };

    // Write logs to stderr or file
    let (wrtr, guard) = match log_file {
        Some(path) => {
            let file = std::fs::File::options()
                .create(true)
                .write(true)
                .truncate(true)
                .open(path)?;
            tracing_appender::non_blocking(file)
        }
        None => tracing_appender::non_blocking(std::io::stderr()),
    };

    let timer = time::format_description::parse("[hour]:[minute]:[second].[subsecond digits:3]")?;
    let time_offset = time::UtcOffset::current_local_offset().unwrap_or(time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, timer);

    tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(wrtr)
        .with_timer(timer)
        .init();

    Ok(guard)
}

pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

