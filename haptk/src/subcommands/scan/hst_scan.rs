use std::{collections::BTreeMap, path::PathBuf};

use petgraph::Graph;
use rayon::prelude::*;
use serde_with::serde_as;

use crate::args::{Selection, StandardArgs};
use crate::io::push_to_output;
use crate::structs::Coord;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::utils::parse_coords;

use color_eyre::{
    eyre::{ensure, eyre, WrapErr},
    Result,
};
use serde::{Deserialize, Serialize};

use crate::io::get_output;
use crate::read_vcf::read_vcf_to_matrix;
use crate::subcommands::bhst::{HstType, Metadata, Node};

#[doc(hidden)]
pub fn run(mut args: StandardArgs, step_size: usize, min_sample_size: usize) -> Result<()> {
    ensure!(
        args.selection == Selection::All
            || args.selection == Selection::Haploid
            || args.selection == Selection::OnlyLongest,
        "Running only with phased data and all chromsomes is supported."
    );
    tracing::info!("Changed --no-alt to true");
    args.no_alt = true;

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "trees", "json.gz");

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None, None, true)?;

    tracing::info!("Starting the HST scan..");

    let hsts = Vec::from_iter(vcf.coords())
        .par_iter()
        .enumerate()
        .filter(|(n, _)| *n % step_size == 0)
        .map(|(_, &coord)| {
            let start_idxs = match args.selection {
                Selection::OnlyLongest => {
                    let res = vcf.only_longest_indexes_no_shard(coord);
                    if res.is_err() {
                        panic!("failed to find the longest-haplotype for coord {coord:?}");
                    }
                    res.ok()
                }
                _ => None,
            };

            (
                coord.clone(),
                construct_bhst_no_mut(&vcf, coord, min_sample_size, start_idxs).unwrap(),
            )
        })
        .collect();

    let metadata = Metadata::new(&vcf, &args, vcf.samples().clone(), HstType::Bhst);

    let trees = HstScan { hsts, metadata };

    write_hsts(trees, output)?;

    Ok(())
}

pub type Limits = (usize, usize, usize, usize);
pub type Hst = Graph<Node, ()>;
pub type HstMap = BTreeMap<Coord, Hst>;

#[serde_as]
#[derive(Serialize, Deserialize, Clone)]
pub struct HstScan {
    pub metadata: Metadata,
    #[serde_as(as = "Vec<(_, _)>")]
    pub hsts: BTreeMap<Coord, Hst>,
}

impl HstScan {
    pub fn nhaplotypes(&self) -> usize {
        self.metadata.samples.len() * *self.metadata.ploidy
    }

    pub fn get_sample_name(&self, index: usize) -> String {
        self.metadata
            .samples
            .get(index / *self.metadata.ploidy)
            .unwrap()
            .clone()
    }

    pub fn get_sample_names(&self, indexes: &[usize]) -> Vec<String> {
        indexes.iter().map(|i| self.get_sample_name(*i)).collect()
    }

    pub fn get_sample_idxs(&self, samples: &[String]) -> Result<Vec<usize>> {
        let idxs: Vec<_> = self
            .metadata
            .samples
            .iter()
            .enumerate()
            .filter(|(_, s)| samples.contains(s))
            .flat_map(|(i, _)| {
                ((i * *self.metadata.ploidy)..(i * *self.metadata.ploidy) + *self.metadata.ploidy)
                    .collect::<Vec<usize>>()
            })
            .collect();

        ensure!(
            !idxs.is_empty(),
            "None of the control samples are found in the vcf."
        );
        Ok(idxs)
    }

    pub fn get_sample_idx(&self, sample: &str) -> Result<usize> {
        self.metadata
            .samples
            .iter()
            .position(|s| s == sample)
            .ok_or_else(|| eyre!("sample {sample} not found in the vcf"))
    }
}

fn write_hsts(trees: HstScan, path: PathBuf) -> Result<()> {
    tracing::info!("Writing trees to file.");
    let now = std::time::Instant::now();
    let mut output = get_output(Some(path))?;

    let mut writer = bgzip::BGZFWriter::new(&mut output, bgzip::Compression::default());

    serde_json::to_writer(&mut writer, &trees)?;

    writer.close()?;

    tracing::info!("Finished writing trees to file.");
    tracing::debug!("Wrote trees file in {:?}", now.elapsed());
    Ok(())
}

pub fn read_tree_file(path: PathBuf) -> Result<HstScan> {
    let file = std::fs::File::open(path.clone()).wrap_err(eyre!("Error opening {path:?}"))?;
    let reader = bgzip::BGZFReader::new(file)?;

    let hst_scan: HstScan = serde_json::from_reader(reader).wrap_err(eyre!(
        "Failed deserializing HSTs from the json.gz. Are you sure the input file is correct?"
    ))?;

    tracing::info!(
        "Read HSTs with the following metadata: \nselection: {:?},\nnsamples:{},",
        hst_scan.metadata.selection,
        hst_scan.metadata.samples.len(),
    );

    Ok(hst_scan)
}
