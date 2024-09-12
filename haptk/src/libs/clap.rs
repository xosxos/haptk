use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use tracing::Level;
use tracing_appender::non_blocking::{NonBlocking, WorkerGuard};
use tracing_subscriber::fmt::time::OffsetTime;

use crate::args::{ConciseArgs, GraphArgs, SortOption, StandardArgs};
use crate::subcommands::{
    bhst, check_for_haplotype, compare_haplotypes, compare_to_haplotype, compare_to_hst, coverage,
    fasta_to_haplotype, haplotype_to_vcf, list_haplotypes, list_markers, list_samples, mrca, uhst,
};

// Genome-wide methods
use crate::subcommands::mrca_scan;
use crate::subcommands::scan::{
    hst_scan, scan_branch_mrca, scan_nodes, scan_quantitative, scan_segregate,
};

#[derive(Parser, Debug)]
#[command(author, version, about, styles=get_styles())]
pub struct Arguments {
    #[command(subcommand)]
    cmd: SubCommand,
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
    },

    ///  (experimental) Analyze the MRCA every x markers along a given contig
    MrcaScan {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
        /// Recombination rate file
        #[arg(short = 'r', long)]
        recombination_rates: PathBuf,

        /// Run the MRCA analysis every n markers
        #[arg(long)]
        step_size: usize,

        /// Dont output to csv
        #[arg(long)]
        no_csv: bool,

        /// Draw plot
        #[arg(long)]
        plot: bool,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,

        #[command(flatten)]
        graph_args: GraphArgs,

        /// Only supports hg38! Mark the centromere to the graph
        #[arg(long)]
        mark_centromere: bool,
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
        #[arg(short = 'h', long)]
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
    },
    /// Show coverage levels per contig of a VCF
    Coverage {
        file: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Minimum amount of basepairs on avg per a single snp
        #[arg(short = 'b', long, default_value_t = 10000)]
        bp_per_snp: i64,

        /// Number of pipes
        #[arg(short = 'p', long, default_value_t = 100)]
        npipes: i64,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,
    },
    /// Read the haplotypes of a given sample
    Haplotypes {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        #[arg(long)]
        selection_variant: Option<String>,
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
    /// (experimental) bHST scan
    BhstScan {
        #[command(flatten)]
        args: StandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Run the scan every n markers
        #[arg(long)]
        step_size: usize,

        /// Number of threads
        #[arg(short = 't', long, default_value_t = 8)]
        threads: usize,
    },
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
    },

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
    },

    /// (experimental) Scan all branches of all trees for the smallest MRCA value
    ScanBranchMrca {
        #[command(flatten)]
        args: ConciseArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Branch sample size to end recursion
        #[arg(long, default_value_t = 15)]
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
    },

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

        /// Wanted segregating samples
        #[arg(short = 'S', long)]
        samples: PathBuf,
    },
}

impl SubCommand {
    pub fn threads(&self) -> usize {
        match self {
            SubCommand::CompareToHaplotype { threads, .. }
            | SubCommand::CompareToHst { threads, .. }
            | SubCommand::Coverage { threads, .. }
            | SubCommand::MrcaScan { threads, .. }
            | SubCommand::Uhst { threads, .. }
            | SubCommand::Bhst { threads, .. }
            | SubCommand::BhstScan { threads, .. }
            | SubCommand::ScanSegregate { threads, .. }
            | SubCommand::ScanBranchMrca { threads, .. }
            | SubCommand::ScanQuantitative { threads, .. }
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
            | SubCommand::Coverage { log_and_verbosity, .. }
            | SubCommand::Mrca { log_and_verbosity, .. }
            | SubCommand::MrcaScan { log_and_verbosity, .. }
            | SubCommand::Uhst { log_and_verbosity, .. }
            | SubCommand::Bhst { log_and_verbosity, .. }
            | SubCommand::BhstScan { log_and_verbosity,  .. }
            | SubCommand::ScanSegregate { log_and_verbosity,  .. }
            | SubCommand::ScanBranchMrca { log_and_verbosity,  .. }
            | SubCommand::ScanQuantitative { log_and_verbosity,  .. }
            | SubCommand::ScanNodes { log_and_verbosity,  .. }
            | SubCommand::Samples { log_and_verbosity, .. }
            | SubCommand::Markers { log_and_verbosity, .. }
            | SubCommand::HaplotypeToVcf { log_and_verbosity, .. }
            | SubCommand::FastaToHaplotype { log_and_verbosity, .. }
            => (log_and_verbosity.verbosity, &log_and_verbosity.log_file, log_and_verbosity.silent),
        }
    }

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
            | SubCommand::MrcaScan { args: StandardArgs { output, .. }, ..}
            | SubCommand::Uhst { args: StandardArgs { output, .. }, ..}
            | SubCommand::Bhst { args: StandardArgs { output, .. }, ..}
            | SubCommand::BhstScan { args: StandardArgs { output, .. }, ..}
            | SubCommand::ScanSegregate { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanBranchMrca { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanQuantitative { args: ConciseArgs { output, .. }, ..}
            | SubCommand::ScanNodes { args: ConciseArgs { output, .. }, ..}
            => Some(output.clone()),
            SubCommand::Samples { .. }
            | SubCommand::HaplotypeToVcf { .. }
            | SubCommand::FastaToHaplotype { .. }
            | SubCommand::Coverage { .. }
            | SubCommand::Markers { .. } => None
        }
    }
}

pub fn run_args(args: Arguments) -> Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.cmd.threads())
        .build_global()?;

    let (verbosity, log_file, is_silent) = args.cmd.log_and_verbosity();

    let (level, wrtr, _guard) = init_tracing(verbosity, log_file, is_silent)?;

    let timer = time::format_description::parse("[hour]:[minute]:[second].[subsecond digits:3]")?;
    let time_offset = time::UtcOffset::current_local_offset().unwrap_or(time::UtcOffset::UTC);
    let timer = OffsetTime::new(time_offset, timer);

    tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(wrtr)
        .with_timer(timer)
        .init();

    if let Some(output) = args.cmd.output() {
        if let Err(e) = std::fs::create_dir(output.clone()) {
            match e.kind() {
                std::io::ErrorKind::AlreadyExists => (),
                _ => return Err(eyre!("Error creating directory {output:?}")),
            }
        }
    }

    run_cmd(args.cmd)?;

    Ok(())
}

#[rustfmt::skip]
pub fn run_cmd(cmd: SubCommand) -> Result<()> {
    match cmd {
        SubCommand::CompareToHaplotype { 
            args, haplotype, mark_samples, mark_shorter_alleles, graph_args, png, npy, sort_option, .. 
        } => compare_to_haplotype::run(
                args, haplotype, mark_samples, mark_shorter_alleles, png, npy, graph_args, sort_option,
            )?,

        SubCommand::Bhst { args,  min_size, publish, .. } => bhst::run(args, min_size, publish,)?,
        SubCommand::Uhst {args,  min_size, publish, .. } => uhst::run(args, min_size, publish,)?,

        SubCommand::CompareHaplotypes { haplotypes, output, prefix, csv, hide_missing, tag_rows, .. }
            => compare_haplotypes::run(haplotypes, output, prefix, csv, hide_missing, tag_rows)?,

        SubCommand::CompareToHst { args, hst, .. } => compare_to_hst::run(args, hst)?,
        SubCommand::CheckForHaplotype { args, haplotype, .. } => check_for_haplotype::run(args, haplotype)?,
        SubCommand::Mrca { args, recombination_rates, .. } => mrca::run(args, recombination_rates)?,
        SubCommand::Haplotypes { args, selection_variant, .. } => list_haplotypes::run(args, selection_variant)?,
        SubCommand::Samples { file, .. } => list_samples::run(file)?, 
        SubCommand::Markers { file, .. } => list_markers::run(file)?,
        SubCommand::HaplotypeToVcf { file, sample_name, output, .. } => haplotype_to_vcf::run(file, sample_name, output)?,
        SubCommand::FastaToHaplotype { file, seq_name, output, .. } => fasta_to_haplotype::run(file, seq_name, output)?,
        SubCommand::Coverage { file, bp_per_snp, npipes, .. } => coverage::run(file, bp_per_snp, npipes)?,

        // Genome-wide methods
        SubCommand::MrcaScan {
            args, recombination_rates, step_size, no_csv, plot, graph_args, mark_centromere, ..
        } => mrca_scan::run(
            args, recombination_rates, step_size, no_csv, plot, graph_args, mark_centromere
        )?,

        SubCommand::BhstScan { args, step_size, .. } => hst_scan::run(args, step_size)?,

        SubCommand::ScanNodes {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len,  ..
        } => scan_nodes::run(
            args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len),
        )?,

        SubCommand::ScanQuantitative {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, var_data, var_name, ..
        } => scan_quantitative::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), var_data, var_name,
            )?,

        SubCommand::ScanBranchMrca {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, recombination_rates, ..
        } => scan_branch_mrca::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), recombination_rates,
            )?,

        SubCommand::ScanSegregate {
            args, min_sample_size, max_sample_size, min_ht_len, max_ht_len, samples, ..
        } => scan_segregate::run(
                args, (min_sample_size, max_sample_size, min_ht_len, max_ht_len), samples,
            )?,

    };
    Ok(())
}

pub fn init_tracing(
    verbosity: u8,
    log_file: &Option<PathBuf>,
    is_silent: bool,
) -> Result<(Level, NonBlocking, WorkerGuard)> {
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
    let (wrtr, _guard) = match log_file {
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

    Ok((level, wrtr, _guard))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_init_tracing() {
        let (level, _, _) = init_tracing(1, &None, false).unwrap();
        assert_eq!(Level::ERROR, level);
        let (level, _, _) = init_tracing(2, &None, false).unwrap();
        assert_eq!(Level::WARN, level);
        let (level, _, _) = init_tracing(3, &None, false).unwrap();
        assert_eq!(Level::INFO, level);
        let (level, _, _) = init_tracing(4, &None, false).unwrap();
        assert_eq!(Level::DEBUG, level);
        let (level, _, _) = init_tracing(5, &None, false).unwrap();
        assert_eq!(Level::TRACE, level);
    }

    #[test]
    fn test_threads() {
        let subcommand = SubCommand::Samples {
            file: PathBuf::new(),
            log_and_verbosity: LogAndVerbosity {
                verbosity: 0,
                log_file: None,
                silent: false,
            },
        };

        assert_eq!(1, subcommand.threads());

        let subcommand = SubCommand::Coverage {
            file: PathBuf::new(),
            bp_per_snp: 0,
            npipes: 0,
            threads: 8,
            log_and_verbosity: LogAndVerbosity {
                verbosity: 0,
                log_file: None,
                silent: false,
            },
        };

        assert_eq!(8, subcommand.threads());
    }
}
