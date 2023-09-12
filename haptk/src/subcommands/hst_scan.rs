use std::path::PathBuf;

use color_eyre::{eyre::ensure, Result};
use petgraph::Graph;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::{Selection, StandardArgs};
use crate::core::get_output;
use crate::core::parse_coords;
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::bhst::{self, Node};
use crate::utils::push_to_output;

#[doc(hidden)]
pub fn run(args: StandardArgs, step_size: usize) -> Result<()> {
    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "trees", "json.gz");

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None)?
        }
        _ => read_vcf_to_matrix(&args, contig, 0, None, None)?,
    };

    // Create a Vec because par_iter cannot be used with pure ranges and par_bridge does not
    // return in ordered fashion with .collect()
    let range: Vec<_> = (0..vcf.ncoords()).collect();

    let trees = range
        .par_iter()
        .filter(|idx| *idx % step_size == 0)
        .map(|idx| {
            if idx % (vcf.ncoords() / 10) == 0 {
                tracing::info!("A tenth has passed");
            }

            TreeRow {
                idx: *idx,
                tree: bhst::construct_bhst(&vcf, *idx, 4),
            }
        })
        .collect();

    let metadata = Metadata::new(&vcf, args);

    let trees = Trees { trees, metadata };

    write_trees(trees, output)?;

    Ok(())
}

#[derive(Serialize, Deserialize, Clone)]
pub struct TreeRow {
    pub idx: usize,
    pub tree: Graph<Node, u8>,
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Metadata {
    pub coords: String,
    pub samples: Vec<String>,
    pub selection: Selection,
    pub vcf_name: PathBuf,
    pub info_limit: Option<f32>,
}

impl Metadata {
    fn new(vcf: &PhasedMatrix, args: StandardArgs) -> Self {
        let samples = vcf.samples().to_vec();

        Self {
            coords: args.coords,
            samples,
            selection: args.selection,
            vcf_name: args.file,
            info_limit: args.info_limit,
        }
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct Trees {
    pub metadata: Metadata,
    pub trees: Vec<TreeRow>,
}

fn write_trees(trees: Trees, path: PathBuf) -> Result<()> {
    // use std::io::{BufWriter, Write};

    tracing::info!("Writing trees to file.");
    let now = std::time::Instant::now();
    let mut output = get_output(Some(path))?;

    let mut writer = bgzip::BGZFWriter::new(&mut output, bgzip::Compression::default());
    // let mut writer = BufWriter::new(output);

    serde_json::to_writer(&mut writer, &trees)?;

    writer.close()?;
    // writer.flush()?;

    tracing::debug!("Wrote trees file in {:?}", now.elapsed());
    Ok(())
}
