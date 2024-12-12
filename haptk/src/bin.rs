use clap::Parser;
use color_eyre::Result;

use haptk::clap::{run_cmd, Arguments};

fn main() -> Result<()> {
    color_eyre::install()?;

    let args = Arguments::parse();

    run_cmd(args.cmd)?;

    Ok(())
}
