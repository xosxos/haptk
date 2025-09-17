use std::{
    collections::HashMap,
    path::PathBuf,
    sync::mpsc::{sync_channel, SyncSender},
    thread,
};

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use petgraph::graph::NodeIndex;
use polars::prelude::*;
use rayon::prelude::*;

use crate::subcommands::scan::{read_tree_file, Hst};
use crate::{
    args::{ConciseArgs, Selection, StandardArgs},
    subcommands::immutable_hst::construct_bhst_no_mut,
};
use crate::{io::read_coords_file, structs::Coord};
use crate::{
    io::{open_csv_writer, push_to_output},
    read_vcf::read_vcf_to_matrix,
};

use super::hst_scan::Limits;

pub fn read_variable_data_file(path: PathBuf) -> Result<DataFrame> {
    let df = CsvReader::from_path(path)?
        .with_null_values(Some(NullValues::AllColumnsSingle("NA".to_string())))
        .infer_schema(Some(500))
        .has_header(true)
        .finish()?;

    Ok(df)
}

type Row = [String; 15];
const HEADER: [&str; 15] = [
    "contig",
    "pos",
    "ref",
    "alt",
    "marker_id",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "pvalue",
    "cases_mean",
    "ctrls_mean",
    "n_het",
    "n_hom",
    "samples",
];

#[doc(hidden)]
pub fn run(
    args: ConciseArgs,
    limits: Limits,
    var_data: PathBuf,
    var_name: String,
    coord_list_path: Option<PathBuf>,
) -> Result<()> {
    let mut var_hm = HashMap::new();

    let df = read_variable_data_file(var_data)?;

    let data: Result<Vec<f64>> = if let Ok(vec) = df[var_name.as_str()].f64() {
        Ok(vec.into_iter().flatten().collect())
    } else if let Ok(vec) = df[var_name.as_str()].i64() {
        Ok(vec.into_iter().flatten().map(|v| v as f64).collect())
    } else {
        Err(eyre!(
                "Error processing the data column {var_name:?}. The data is not floats or integers format"
            ))
    };

    let data = data?;
    let ids: Vec<&str> = df["parsed_id"].utf8()?.into_iter().flatten().collect();

    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let args = StandardArgs {
        file: args.file.clone(),
        output: args.output,
        coords: String::from(""),
        selection: args.selection.clone(),
        prefix: args.prefix,
        samples: args.samples.clone(),
        no_alt: true,
        list: None,
    };

    let (tx, rx): (SyncSender<Row>, _) = sync_channel(2024);

    let cargs = args.clone();
    let writer_handle = thread::spawn(move || -> Result<()> {
        let mut output = cargs.output.clone();
        push_to_output(&cargs, &mut output, "quantitative_scan", "csv");
        let mut writer = open_csv_writer(output)?;
        writer.write_record(HEADER)?;

        while let Ok(row) = rx.recv() {
            writer.write_record(row)?;
        }
        Ok(())
    });

    let coord_list: Option<Vec<Coord>> = coord_list_path.map(|v| {
        read_coords_file(&v).unwrap_or_else(|_| panic!("Could not read coords file {v:?}"))
    });

    match coord_list {
        Some(coord_list) => {
            tracing::info!("Reading {} coords to HSTs.", coord_list.len());
            let mut map = HashMap::new();

            for coord in coord_list {
                map.entry(coord.contig.clone())
                    .or_insert(Vec::new())
                    .push(coord);
            }

            for (contig, coords) in map {
                let args = StandardArgs {
                    file: args.file.clone(),
                    output: args.output.clone(),
                    coords: format!("{}", coords[0]),
                    selection: args.selection.clone(),
                    prefix: args.prefix.clone(),
                    samples: args.samples.clone(),
                    no_alt: true,
                    list: None,
                };

                let vcf = read_vcf_to_matrix(&args, &contig, 0, None, None, None, true)?;

                let samples = vcf.samples().clone();
                let ploidy = *vcf.ploidy;
                let nhaplotypes = vcf.nhaplotypes();

                tracing::info!("Starting the case-ctrl scan for {contig}.");

                coords
                    .into_par_iter()
                    .map(|coord| {
                        let start_idxs = match args.selection {
                            Selection::OnlyLongest => {
                                let res = vcf.only_longest_indexes_no_shard(&coord);
                                if res.is_err() {
                                    panic!(
                                        "failed to find the longest-haplotype for coord {coord:?}"
                                    );
                                }
                                res.ok()
                            }
                            _ => None,
                        };
                        (
                            construct_bhst_no_mut(&vcf, &coord, limits.0, start_idxs).unwrap(),
                            coord,
                        )
                    })
                    .filter_map(|(hst, coord)| optimizer(hst, coord, limits, &var_hm))
                    .map(|(hst, coord, top_node_idx, optimized_value)| {
                        rower(
                            hst,
                            coord,
                            top_node_idx,
                            optimized_value,
                            &var_hm,
                            nhaplotypes,
                            &samples,
                            ploidy,
                        )
                    })
                    .try_for_each(|row| tx.send(row))?;
            }
        }
        None => {
            let hsts = read_tree_file(args.file)?;
            for (id, value) in ids.iter().zip(data.iter()) {
                var_hm.insert(hsts.get_sample_idx(id)?, *value);
            }

            let samples = hsts.metadata.samples.clone();
            let ploidy = *hsts.metadata.ploidy;
            let nhaplotypes = hsts.nhaplotypes();

            hsts.hsts
                .into_par_iter()
                .filter_map(|(coord, hst)| optimizer(hst, coord, limits, &var_hm))
                .map(|(hst, coord, top_node_idx, optimized_value)| {
                    rower(
                        hst,
                        coord,
                        top_node_idx,
                        optimized_value,
                        &var_hm,
                        nhaplotypes,
                        &samples,
                        ploidy,
                    )
                })
                .try_for_each(|row| tx.send(row))?;
        }
    }

    tracing::info!("Finished the quantitative scan.");

    drop(tx);
    let _ = writer_handle.join();

    Ok(())
}

pub fn get_variable_data(indexes: &[usize], hm: &HashMap<usize, f64>) -> Vec<f64> {
    indexes.iter().filter_map(|v| hm.get(v)).copied().collect()
}

#[allow(clippy::too_many_arguments)]
fn rower(
    hst: Hst,
    coord: Coord,
    top_node_idx: NodeIndex,
    opt_value: f64,
    var_hm: &HashMap<usize, f64>,
    nhaplotypes: usize,
    samples: &[String],
    ploidy: usize,
) -> Row {
    let top_node = hst.node_weight(top_node_idx).unwrap();
    let (nhet, nhom) = top_node.zygosity(samples, ploidy);

    let other_indexes: Vec<usize> = (0..nhaplotypes)
        .filter(|i| !top_node.indexes.contains(i))
        .collect();

    // Note: v1 does not necessarily equal indexes.len() because hashmap might miss data
    let v1 = get_variable_data(&top_node.indexes, var_hm);
    let v2 = get_variable_data(&other_indexes, var_hm);

    let cases_mean = v1.iter().sum::<f64>() / v1.len() as f64;
    let ctrls_mean = v2.iter().sum::<f64>() / v2.len() as f64;

    let names = top_node.sample_name_list(samples, ploidy);

    [
        coord.contig,
        coord.pos.to_string(),
        coord.reference,
        coord.alt,
        top_node.identifier(),
        (top_node.stop.pos.saturating_sub(top_node.start.pos)).to_string(),
        top_node.haplotype.len().to_string(),
        top_node.start.pos.to_string(),
        top_node.stop.pos.to_string(),
        opt_value.to_string(),
        nhom.to_string(),
        nhet.to_string(),
        cases_mean.to_string(),
        ctrls_mean.to_string(),
        names,
    ]
}

fn optimizer(
    hst: Hst,
    coord: Coord,
    limits: Limits,
    hm: &HashMap<usize, f64>,
) -> Option<(Hst, Coord, NodeIndex, f64)> {
    let (nmin_samples, nmax_samples, nmin_variants, nmax_variants) = limits;
    let mut top_node_idx = NodeIndex::new(0);

    // START VALUE
    let start_value = f64::MAX;

    let mut optimized_value = start_value;

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
            // START OPTIMIZER CODE
            let _node = hst.node_weight(idx).unwrap();

            let node = hst.node_weight(NodeIndex::new(0)).unwrap();
            let nhaplotypes = node.indexes.len();

            let other_indexes: Vec<usize> = (0..nhaplotypes)
                .filter(|i| !node.indexes.contains(i))
                .collect();

            let mut v1 = get_variable_data(&node.indexes, hm);
            let v2 = get_variable_data(&other_indexes, hm);

            let mut genotypes = vec![1.0; node.indexes.len()];
            genotypes.extend(vec![0.0; other_indexes.len()]);

            v1.extend(v2.iter());

            // crate::stats::linear_regression(v1, genotypes)
            let value = crate::stats::two_tail_welch_t_test(&v1, &v2);
            let clause = value < optimized_value;

            // END OPTIMIZER CODE

            if clause {
                optimized_value = value;
                top_node_idx = idx;
            }
        }
    }

    if optimized_value == start_value {
        tracing::warn!(
            "No optimized value for the HST in contig {} at pos {}",
            coord.contig,
            coord.pos,
        );
        return None;
    }

    Some((hst, coord, top_node_idx, optimized_value))
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    #[test]
    fn read_variable_data() {
        let path = PathBuf::from("tests/data/clinical_data.csv");
        let df = read_variable_data_file(path).unwrap();

        let array = df["id"].utf8().unwrap();
        let vec: Vec<_> = array.into_iter().flatten().collect();
        assert_eq!(vec[0], "SAMPLE1");

        let array = df["aoo"].i64().unwrap();
        let vec: Vec<_> = array.into_iter().flatten().collect();
        assert_eq!(vec[0], 88);
        assert_eq!(vec[2], 58);

        let path = PathBuf::from("tests/data/does_not_exist.csv");
        let res = read_variable_data_file(path);
        assert!(res.is_err());
    }
}
