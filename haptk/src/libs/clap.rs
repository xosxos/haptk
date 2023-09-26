use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use tracing::Level;
use tracing_appender::non_blocking::{NonBlocking, WorkerGuard};
use tracing_subscriber::fmt::time::OffsetTime;

use crate::args::{ConciseArgs, GraphArgs, Selection, StandardArgs};
use crate::subcommands::{
    bhst, check_for_haplotype, compare_haplotypes, compare_to_haplotype, compare_to_hst, coverage,
    haplotype_to_vcf, list_haplotypes, list_markers, list_samples, mrca, uhst,
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

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum ClapSelection {
    /// Select all alleles from samples (as of now only diploid or haploid organisms are supported)
    All,
    /// Select only the allele per each sample sharing the most haplotype with the haplotypes of the other samples
    LongestHaplotype,
    /// Select only the alleles containing the REF variant at given a coordinate
    OnlyRefs,
    /// Select only the alleles containing the ALT variant at given a coordinate
    OnlyAlts,
    /// Use for unphased data (currently not supported for most commands)
    Unphased,
    /// Use for haploid genotypes
    Haploid,
}

#[derive(Args, Debug, Clone)]
pub struct ClapStandardArgs {
    pub file: PathBuf,

    /// The starting coordinate, i.e. chr9:27573534
    #[arg(short = 'c', long)]
    pub coords: String,

    /// Output directory
    #[arg(short = 'o', long, default_value_os_t = PathBuf::from("./"))]
    pub outdir: PathBuf,

    /// List of samples for HST construction (one ID per row)
    #[arg(short = 'S', long, value_delimiter = ' ', num_args = 1.. )]
    pub samples: Option<Vec<PathBuf>>,

    #[arg(short = 'a', long, value_enum)]
    pub alleles: ClapSelection,

    // #[arg(long)]
    // pub info_limit: Option<f32>,
    /// Output filename prefix
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,
}

#[derive(Args, Debug, Clone)]
pub struct ClapConciseArgs {
    pub file: PathBuf,

    /// Output directory
    #[arg(short = 'o', long, default_value_os_t = PathBuf::from("./"))]
    pub outdir: PathBuf,

    /// Output filename prefix
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,
}

#[derive(Args, Debug, Clone)]
pub struct ClapIoArgs {
    pub file: PathBuf,

    /// Output directory
    #[arg(short = 'o', long, default_value_os_t = PathBuf::from("./"))]
    pub outdir: PathBuf,

    /// Output filename prefix
    #[arg(short = 'p', long)]
    pub prefix: Option<String>,
}

impl From<ClapStandardArgs> for StandardArgs {
    fn from(value: ClapStandardArgs) -> Self {
        let prefix = match value.prefix {
            Some(v) => match v.as_str() {
                "" => None,
                "\\0" => None,
                v => Some(v.to_string()),
            },
            None => None,
        };
        Self {
            file: value.file,
            coords: value.coords,
            output: value.outdir,
            samples: value.samples,
            selection: value.alleles.into(),
            info_limit: None,
            // info_limit: value.info_limit,
            prefix,
        }
    }
}

impl From<ClapConciseArgs> for ConciseArgs {
    fn from(value: ClapConciseArgs) -> Self {
        let prefix = match value.prefix {
            Some(v) => match v.as_str() {
                "" => None,
                "\\0" => None,
                v => Some(v.to_string()),
            },
            None => None,
        };
        Self {
            file: value.file,
            output: value.outdir,
            prefix,
        }
    }
}

impl From<ClapSelection> for Selection {
    fn from(value: ClapSelection) -> Self {
        match value {
            ClapSelection::OnlyAlts => Self::OnlyAlts,
            ClapSelection::OnlyRefs => Self::OnlyRefs,
            ClapSelection::All => Self::All,
            ClapSelection::LongestHaplotype => Self::OnlyLongest,
            ClapSelection::Unphased => Self::Unphased,
            ClapSelection::Haploid => Self::Haploid,
        }
    }
}

impl From<Selection> for ClapSelection {
    fn from(value: Selection) -> Self {
        match value {
            Selection::OnlyAlts => Self::OnlyAlts,
            Selection::OnlyRefs => Self::OnlyRefs,
            Selection::All => Self::All,
            Selection::OnlyLongest => Self::LongestHaplotype,
            Selection::Unphased => Self::Unphased,
            Selection::Haploid => Self::Haploid,
        }
    }
}

#[derive(Args, Debug, Default, Clone)]
pub struct ClapGraphArgs {
    /// Graph width in px
    #[arg(long)]
    pub width: Option<f32>,

    /// Graph height in px
    #[arg(long)]
    pub height: Option<f32>,

    /// Mark the variant coordinate
    #[arg(long)]
    pub mark_locus: bool,

    // Font size
    #[arg(long)]
    pub font_size: Option<f32>,

    // Line stroke width
    #[arg(long)]
    pub stroke_width: Option<u32>,

    // Font color
    #[arg(long)]
    pub color: Option<String>,

    // Background color
    #[arg(long)]
    pub background_color: Option<String>,
}

impl From<ClapGraphArgs> for GraphArgs {
    fn from(value: ClapGraphArgs) -> Self {
        Self {
            width: value.width.unwrap_or(2560.0),
            height: value.height.unwrap_or(1440.0),
            mark_locus: value.mark_locus,
            font_size: value.font_size.unwrap_or(20.0),
            stroke_width: value.stroke_width.unwrap_or(5),
            color: value.color.unwrap_or("black".to_string()),
            background_color: value.background_color.unwrap_or("white".to_string()),
        }
    }
}

#[derive(Subcommand, Debug)]
pub enum SubCommand {
    /// Build unidirectional haplotype sharing trees at a coordinate
    Uhst {
        #[command(flatten)]
        args: ClapStandardArgs,

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
        args: ClapStandardArgs,

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
    /// bHST scan
    BhstScan {
        #[command(flatten)]
        args: ClapStandardArgs,

    /// Analyze the MRCA based on the Gamma method at a coordinate
    Mrca {
        #[command(flatten)]
        args: ClapStandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Recombination rate file
        #[arg(short = 'r', long)]
        recombination_rates: PathBuf,

        #[arg(long)]
        start: Option<u64>,

        #[arg(long)]
        stop: Option<u64>,
    },
    /// Check if samples share a given haplotype
    CheckForHaplotype {
        #[command(flatten)]
        args: ClapStandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,

        /// Haplotype for checking
        #[arg(short = 'h', long)]
        haplotype: PathBuf,
    },
    /// Check differences between samples and a haplotype
    CompareToHaplotype {
        #[command(flatten)]
        args: ClapStandardArgs,

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
        graph_args: ClapGraphArgs,

        /// Output .png image
        #[arg(long)]
        png: bool,

        /// Output .npy matrix
        #[arg(long)]
        npy: bool,
    },

    /// Check which haplotypes of the HST are present in samples
    CompareToHst {
        #[command(flatten)]
        args: ClapStandardArgs,

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
        #[arg(short = 'o', long, default_value_os_t = PathBuf::from("./"))]
        outdir: PathBuf,

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
        args: ClapStandardArgs,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
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
    ToVcf {
        file: PathBuf,

        #[arg(short = 's', long)]
        sample_name: String,

        #[arg(short = 'o', long, default_value_os_t = PathBuf::from("-"))]
        outdir: PathBuf,

        #[command(flatten)]
        log_and_verbosity: LogAndVerbosity,
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
            | SubCommand::BhstSegregation { threads, .. }
            | SubCommand::BhstScan { threads, .. }
            | SubCommand::BhstGwas { threads, .. }
            | SubCommand::BhstQtGwas { threads, .. }
            | SubCommand::BhstMrcaGwas { threads, .. } => *threads,
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
            | SubCommand::BhstSegregation { log_and_verbosity,  .. }
            | SubCommand::BhstMrcaGwas { log_and_verbosity,  .. }
            | SubCommand::BhstQtGwas { log_and_verbosity,  .. }
            | SubCommand::BhstGwas { log_and_verbosity,  .. }
            | SubCommand::Samples { log_and_verbosity, .. }
            | SubCommand::Markers { log_and_verbosity, .. }
            | SubCommand::ToVcf { log_and_verbosity, .. }
            => (log_and_verbosity.verbosity, &log_and_verbosity.log_file, log_and_verbosity.silent),
        }
    }

    pub fn check_sample_size(&self) -> Result<()> {
        match self {
            SubCommand::BhstGwas {
                min_sample_size, ..
            }
            | SubCommand::BhstQtGwas {
                min_sample_size, ..
            }
            | SubCommand::BhstMrcaGwas {
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
            SubCommand::Haplotypes { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::CheckForHaplotype { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::CompareToHaplotype { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::CompareToHst { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::CompareHaplotypes { outdir, .. }
            | SubCommand::Mrca { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::MrcaScan { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::Uhst { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::Bhst { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::BhstScan { args: ClapStandardArgs { outdir, .. }, ..}
            | SubCommand::BhstSegregation { args: ClapConciseArgs { outdir, .. }, ..}
            | SubCommand::BhstMrcaGwas { args: ClapConciseArgs { outdir, .. }, ..}
            | SubCommand::BhstQtGwas { args: ClapConciseArgs { outdir, .. }, ..}
            | SubCommand::BhstGwas { args: ClapConciseArgs { outdir, .. }, ..}
            | SubCommand::ToVcf { outdir, .. }
            => Some(outdir.clone()),
            SubCommand::Samples { .. }
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
            args, haplotype, mark_samples, mark_shorter_alleles, graph_args, png, npy, .. 
        } => compare_to_haplotype::run(
                args.into(), haplotype, mark_samples, mark_shorter_alleles, png, npy,
                GraphArgs {
                    width: graph_args.width.unwrap_or(5000.0),
                    height: graph_args.height.unwrap_or(7500.0),
                    mark_locus: graph_args.mark_locus,
                    font_size: graph_args.font_size.unwrap_or(75.0),
                    stroke_width: graph_args.stroke_width.unwrap_or(7),
                    color: graph_args.color.unwrap_or("black".into()),
                    background_color: graph_args.background_color.unwrap_or("white".into()),
                },
            )?,

        SubCommand::Bhst { args,  min_size, publish, .. } => bhst::run(args.into(), min_size, publish,)?,
        SubCommand::Uhst {args,  min_size, publish, .. } => uhst::run(args.into(), min_size, publish,)?,

        SubCommand::CompareHaplotypes { haplotypes, outdir, prefix, csv, hide_missing, tag_rows, .. }
            => compare_haplotypes::run(haplotypes, outdir, prefix, csv, hide_missing, tag_rows)?,

        SubCommand::CompareToHst { args, hst, .. } => compare_to_hst::run(args.into(), hst)?,
        SubCommand::CheckForHaplotype { args, haplotype, .. } => check_for_haplotype::run(args.into(), haplotype)?,
        SubCommand::Mrca { args, recombination_rates, start, stop, .. } => mrca::run(args.into(), recombination_rates, start, stop)?,
        SubCommand::Haplotypes { args, .. } => list_haplotypes::run(args.into())?,
        SubCommand::Samples { file, .. } => list_samples::run(file)?, 
        SubCommand::Markers { file, .. } => list_markers::run(file)?,
        SubCommand::ToVcf { file, sample_name, outdir, .. } => haplotype_to_vcf::run(file, sample_name, outdir)?,
        SubCommand::Coverage { file, bp_per_snp, npipes, .. } => coverage::run(file, bp_per_snp, npipes)?,

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
    fn test_from_traits() {
        let clap_args = ClapGraphArgs::default();
        let gr_args1: GraphArgs = clap_args.into();
        let gr_args2 = GraphArgs::default();
        assert_eq!(gr_args1, gr_args2);

        let clap_selection = ClapSelection::All;
        let selection: Selection = clap_selection.into();
        assert_eq!(selection, Selection::All)
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
