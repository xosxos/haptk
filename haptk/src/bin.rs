#![allow(clippy::missing_errors_doc)]
use std::path::PathBuf;

use clap::Parser;
use color_eyre::Result;

use haptk::clap::{run_cmd, Arguments};

fn main() -> Result<()> {
    color_eyre::install()?;

    let args = Arguments::parse();

    // Initialize logger
    let (verbosity, log_file, is_silent) = args.cmd.log_and_verbosity();
    let _guard = init_tracing(verbosity, log_file, is_silent)?;

    run_cmd(args.cmd)?;

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
