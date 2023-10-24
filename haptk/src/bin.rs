use clap::Parser;
use color_eyre::Result;

use haptk::clap::{run_args, Arguments};

fn main() -> Result<()> {
    color_eyre::install()?;

    let args = Arguments::parse();
    run_args(args)?;

    Ok(())
}
