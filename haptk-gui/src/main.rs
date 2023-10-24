use color_eyre::Result;
use klask::Settings;
use haptk::clap::{run_args, Arguments};

fn main() -> Result<()> {
    color_eyre::install()?;

    klask::run_derived::<Arguments, _>(Settings::default(), |args| {
        run_args(args).unwrap();
    });

    Ok(())
}
