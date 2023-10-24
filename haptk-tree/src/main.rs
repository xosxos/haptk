#![feature(slice_partition_dedup)]
use std::path::PathBuf;

use clap::Parser;
use color_eyre::Result;
use eframe::Theme;

use haptk::io::read_sample_ids;
use haptk::libs::clap::{ClapStandardArgs, LogAndVerbosity};
use haptk::subcommands::bhst;
use haptk::subcommands::uhst;

mod app;
mod central_panel;
mod plot;
mod side_panel;
mod tree;
mod utils;

use crate::app::TreeApp;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    #[command(flatten)]
    pub args: ClapStandardArgs,

    #[command(flatten)]
    pub log_and_verbosity: LogAndVerbosity,

    /// Minimum number of samples per node in tree
    #[arg(short = 's', long, default_value_t = 1)]
    pub min_node_size: usize,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 8)]
    pub threads: usize,

    /// Mark samples for tagging
    #[arg(short = 'S', long)]
    pub mark_samples: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let now = std::time::Instant::now();

    let lv = args.log_and_verbosity.clone();
    let (level, wrtr, _guard) =
        haptk::libs::clap::init_tracing(lv.verbosity, &lv.log_file, lv.silent)?;

    tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(wrtr)
        .init();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    let marked_samples = read_sample_ids(&args.mark_samples)?;

    let vcf = bhst::read_vcf_with_selections(&args.args.into())?;

    let bhst = bhst::construct_bhst(&vcf, vcf.variant_idx(), 1);
    let uhst_left =
        uhst::construct_uhst(&vcf, &uhst::LocDirection::Left, vcf.variant_idx(), 1, false);
    let uhst_right = uhst::construct_uhst(
        &vcf,
        &uhst::LocDirection::Right,
        vcf.variant_idx(),
        1,
        false,
    );

    if let Err(err) = std::fs::create_dir("results_tree_haplotypes") {
        tracing::warn!("Error creating dir {err}");
    }

    let app = TreeApp::new(
        vcf,
        bhst,
        uhst_left,
        uhst_right,
        args.min_node_size,
        marked_samples,
    );

    let options = eframe::NativeOptions {
        default_theme: Theme::Light,
        ..Default::default()
    };

    tracing::debug!("App initialized at {:?}", now.elapsed());

    eframe::run_native("HAPTK Tree Viewer", options, Box::new(|_| Box::new(app))).unwrap();

    Ok(())
}
