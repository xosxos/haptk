use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread;

use color_eyre::eyre::ensure;
use color_eyre::Result;
use fishers_exact::fishers_exact;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::{
    filter_out_non_vcf_sample_names, open_csv_writer, push_to_output, read_coords_sample_file,
    read_sample_ids,
};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::Coord;
use crate::subcommands::bhst::Node;
use crate::subcommands::immutable_hst::construct_bhst_no_mut;
use crate::subcommands::scan::{Hst, Limits};

type Row = [String; 17];

const HEADER: [&str; 17] = [
    "contig",
    "pos",
    "ref",
    "alt",
    "marker_id",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "n",
    "case_ht_n",
    "ctrl_ht_n",
    "nhet_cases",
    "nhom_cases",
    "nhet_ctrls",
    "nhom_ctrls",
    "samples",
];

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct CoordSamples {
    pub contig: String,
    pub pos: u64,
    #[serde(rename(serialize = "ref"), alias = "ref")]
    pub reference: String,
    pub alt: String,
    #[serde(default)]
    pub marker_id: Option<String>,
}

impl From<CoordSamples> for Coord {
    fn from(hap: CoordSamples) -> Self {
        Self {
            contig: hap.contig,
            pos: hap.pos,
            reference: hap.reference,
            alt: hap.alt,
        }
    }
}

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    case_samples: PathBuf,
    ctrl_samples: PathBuf,
    coord_list_path: Option<PathBuf>,
    limit: f64,
) -> Result<()> {
    let case_samples = read_sample_ids(&Some(case_samples))?.unwrap();
    let ctrl_samples = read_sample_ids(&Some(ctrl_samples))?.unwrap();

    let args = StandardArgs {
        file: args.file.clone(),
        output: args.output,
        coords: String::from(""),
        selection: args.selection.clone(),
        prefix: args.prefix,
        samples: args.samples.clone(),
        no_alt: true,
        list: None,
        only_snv: true,
    };

    let (tx, rx): (SyncSender<Row>, _) = sync_channel(2024);

    let cargs = args.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut output = cargs.output.clone();
        push_to_output(&cargs, &mut output, "segregate_scan", "csv");
        let mut writer = open_csv_writer(output)?;
        writer.write_record(HEADER)?;

        while let Ok(row) = rx.recv() {
            writer.write_record(row)?;
        }
        Ok(())
    });

    let coord_list: Option<Vec<CoordSamples>> = coord_list_path.map(|v| {
        read_coords_sample_file(&v).unwrap_or_else(|_| panic!("Could not read coords file {v:?}"))
    });

    match coord_list {
        Some(coord_list) => {
            ensure!(!coord_list.is_empty(), "The coordinate list is empty");
            tracing::info!("Reading {} coords to HSTs.", coord_list.len());
            let mut out_map = HashMap::new();

            let contig = &coord_list[0].contig;

            let ctrl_samples = filter_out_non_vcf_sample_names(&args.file, contig, ctrl_samples)?;
            let case_samples = filter_out_non_vcf_sample_names(&args.file, contig, case_samples)?;

            for coord in coord_list {
                out_map
                    .entry(coord.contig.clone())
                    .or_insert(HashMap::new())
                    .entry(std::convert::Into::<Coord>::into(coord.clone()))
                    .or_insert(Vec::new())
                    .push(coord.marker_id);
            }

            for (contig, in_map) in out_map {
                let args = StandardArgs {
                    file: args.file.clone(),
                    output: args.output.clone(),
                    coords: contig.to_string(),
                    selection: args.selection.clone(),
                    prefix: args.prefix.clone(),
                    samples: args.samples.clone(),
                    no_alt: true,
                    list: None,
                    only_snv: true,
                };

                let vcf = read_vcf_to_matrix(&args, &contig, 0, None, None, None, true)?;

                let ctrl_list = vcf.get_idxs_for_samples(&ctrl_samples)?;
                let case_list = vcf.get_idxs_for_samples(&case_samples)?;

                tracing::info!("Starting the case-ctrl scan for {contig}.");

                in_map
                    .into_par_iter()
                    .try_for_each(|(coord, cases_list)| -> Result<()> {
                        let hst = construct_bhst_no_mut(&vcf, &coord, limits.0, None)?;

                        for marker_id in &cases_list {
                            find_optimized_value_and_create_csv_row(
                                &hst,
                                &coord,
                                limits,
                                vcf.samples(),
                                &case_list,
                                &ctrl_list,
                                *vcf.ploidy,
                                marker_id,
                                limit,
                                tx.clone(),
                            )?;
                        }
                        Ok(())
                    })?;
            }
        }
        None => {
            unreachable!()
            // let hsts = read_tree_file(args.file)?;

            // // // Filter out samples not present in the HST SCAN
            // // let case_list_names: Vec<String> = case_samples
            // //     .into_iter()
            // //     .filter(|s| {
            // //         if hsts.metadata.samples.contains(s) {
            // //             true
            // //         } else {
            // //             tracing::warn!("Sample {s:?} is on the ID list, but not in the HSTs.");
            // //             false
            // //         }
            // //     })
            // //     .collect();

            // let case_list = hsts.get_sample_idxs(&case_list_names).unwrap();

            // let ctrl_list: Vec<_> = hsts
            //     .metadata
            //     .samples
            //     .iter()
            //     .filter(|s| !case_list_names.contains(s))
            //     .flat_map(|name| hsts.get_sample_idxs(&[name.clone()]).unwrap())
            //     .collect();

            // let ploidy = hsts.metadata.ploidy;
            // let samples = hsts.metadata.samples;

            // tracing::info!("Starting the case-ctrl scan.");

            // hsts.hsts.into_par_iter().try_for_each(|(coord, hst)| {
            //     find_optimized_value_and_create_csv_row(
            //         &hst,
            //         &coord,
            //         limits,
            //         &samples,
            //         &case_list,
            //         &ctrl_list,
            //         *ploidy,
            //         limit,
            //         tx.clone(),
            //     )
            // })?;
        }
    }

    tracing::info!("Finished the case-ctrl scan.");

    drop(tx);
    let _ = writer_handle.join();

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn find_optimized_value_and_create_csv_row(
    hst: &Hst,
    coord: &Coord,
    limits: Limits,
    samples: &[String],
    case_list: &[usize],
    ctrl_list: &[usize],
    ploidy: usize,
    marker_id: &Option<String>,
    limit: f64,
    tx: SyncSender<Row>,
) -> Result<()> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;

    let mut iterator = hst.node_indices();

    // Skip first node because it contains no haplotypes
    iterator.next();

    for idx in iterator {
        let node = hst.node_weight(idx).unwrap();

        if node.indexes.len() >= nmin_samples
            && node.indexes.len() <= nmax_samples
            && node.haplotype.len() >= nmin_variants
            && node.haplotype.len() <= nmax_variants
        {
            // OPTIMIZER CODE
            let node = hst.node_weight(idx).unwrap();

            // let value = segregation_optimizer(case_list, names);
            let value = fisher_optimizer(case_list, ctrl_list, &node.indexes);

            if value < limit {
                let (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls) =
                    case_ctrl_zygosity_from_node(node, case_list, ctrl_list, samples, ploidy);
                let names = node.sample_name_list(samples, ploidy);

                if let Some(id) = marker_id {
                    if id == &node.identifier() {
                        tx.send([
                            coord.contig.clone(),
                            coord.pos.to_string(),
                            coord.reference.clone(),
                            coord.alt.clone(),
                            node.identifier(),
                            (node.stop.pos.saturating_sub(node.start.pos)).to_string(),
                            node.haplotype.len().to_string(),
                            node.start.pos.to_string(),
                            node.stop.pos.to_string(),
                            format!("{value:+.4e}"),
                            case_list.len().to_string(),
                            ctrl_list.len().to_string(),
                            nhet_cases.to_string(),
                            nhom_cases.to_string(),
                            nhet_ctrls.to_string(),
                            nhom_ctrls.to_string(),
                            names.to_string(),
                        ])?;
                        return Ok(());
                    }
                } else {
                    tx.send([
                        coord.contig.clone(),
                        coord.pos.to_string(),
                        coord.reference.clone(),
                        coord.alt.clone(),
                        node.identifier(),
                        (node.stop.pos.saturating_sub(node.start.pos)).to_string(),
                        node.haplotype.len().to_string(),
                        node.start.pos.to_string(),
                        node.stop.pos.to_string(),
                        format!("{value:+.4e}"),
                        case_list.len().to_string(),
                        ctrl_list.len().to_string(),
                        nhet_cases.to_string(),
                        nhom_cases.to_string(),
                        nhet_ctrls.to_string(),
                        nhom_ctrls.to_string(),
                        names.to_string(),
                    ])?;
                }
            }
        }
    }
    Ok(())
}

#[allow(dead_code)]
// Check that all samples are only from the wanted list and return their amount
fn segregation_optimizer(case_list: &[String], names: Vec<String>) -> f64 {
    match names.iter().all(|name| case_list.contains(name)) {
        true => names.len() as f64,
        false => 0.0,
    }
}

fn fisher_optimizer(case_list: &[usize], ctrl_list: &[usize], indexes: &[usize]) -> f64 {
    let ncases = case_list.len();
    let ncontrols = ctrl_list.len();

    let cases_true = indexes
        .iter()
        .filter(|name| case_list.contains(name))
        .count();

    let controls_true = indexes
        .iter()
        .filter(|name| ctrl_list.contains(name))
        .count();

    let cases_false = ncases - cases_true;
    let controls_false = ncontrols - controls_true;

    fishers_exact(&[
        cases_true as u32,
        controls_true as u32,
        cases_false as u32,
        controls_false as u32,
    ])
    .unwrap()
    .greater_pvalue
}

pub fn case_ctrl_zygosity_from_node(
    node: &Node,
    case_list: &[usize],
    ctrl_list: &[usize],
    samples: &[String],
    ploidy: usize,
) -> (usize, usize, usize, usize) {
    // Filter in only ctrl samples and then count the samples
    // Finally, filter out all duplicates (homozygotes) by collecting to a set
    let mut ctrl_prior = 0;
    let ctrl_post: HashSet<_> = node
        .indexes
        .iter()
        .filter(|i| ctrl_list.contains(i))
        .inspect(|_| ctrl_prior += 1)
        .map(|i| &samples[i / ploidy])
        .collect();

    let mut case_prior = 0;
    let case_post: HashSet<_> = node
        .indexes
        .iter()
        .filter(|i| case_list.contains(i))
        .inspect(|_| case_prior += 1)
        .map(|i| &samples[i / ploidy])
        .collect();

    // Calculate homozygosity statistics based on the amount of homozygotes removed above
    let nhom_cases = case_prior - case_post.len();
    let nhet_cases = (2 * case_post.len()) - case_prior;

    let nhom_ctrls = ctrl_prior - ctrl_post.len();
    let nhet_ctrls = (2 * ctrl_post.len()) - ctrl_prior;

    (nhet_cases, nhom_cases, nhet_ctrls, nhom_ctrls)
}
