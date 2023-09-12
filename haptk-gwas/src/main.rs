#![feature(btree_extract_if)]
use std::path::PathBuf;

use clap::Parser;
use color_eyre::Result;
use eframe::Theme;
use haptk::libs::clap::{ClapConciseArgs, LogAndVerbosity};

mod alignment;
mod app;
mod central_panel;
mod plot;
mod side_panel;
mod tree;
mod utils;

use crate::app::GwasApp;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Args {
    #[command(flatten)]
    pub args: ClapConciseArgs,

    #[command(flatten)]
    pub log_and_verbosity: LogAndVerbosity,

    /// Number of threads
    #[arg(short = 't', long, default_value_t = 8)]
    pub threads: usize,

    /// Bilateral haplotype state trees
    #[arg(long)]
    pub trees: PathBuf,

    /// List of control ids
    #[arg(long)]
    pub controls: Option<PathBuf>,

    /// Path to controls vcf if controls are not included in the same file
    #[arg(long)]
    pub control_vcf: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let now = std::time::Instant::now();

    let lv = args.log_and_verbosity.clone();
    let (level, wrtr, _guard) = haptk::libs::clap::init_tracing(lv.verbosity, &lv.log_file)?;

    tracing_subscriber::fmt()
        .with_max_level(level)
        .with_writer(wrtr)
        .init();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    let app = GwasApp::new(args)?;

    let options = eframe::NativeOptions {
        default_theme: Theme::Light,
        ..Default::default()
    };

    tracing::debug!("App initialized at {:?}", now.elapsed());

    eframe::run_native("Monitor app", options, Box::new(|_| Box::new(app))).unwrap();

    Ok(())
}
