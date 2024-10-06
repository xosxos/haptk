use std::{
    collections::{BTreeMap, HashMap},
    path::PathBuf,
    sync::Arc,
};

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output};
use crate::subcommands::scan::{read_tree_file, return_assoc, write_assoc, AssocRow, Hst, HstScan};
use crate::utils::parse_coords;
use crate::{io::read_variable_data_file, structs::Coord};

use super::hst_scan::{top_node_from_hsts, zygosity_from_node, Limits};

const HEADER: &[&str] = &[
    "contig",
    "pos",
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
];

#[doc(hidden)]
pub fn run(args: ConciseArgs, limits: Limits, var_data: PathBuf, var_name: String) -> Result<()> {
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

    let hsts = read_tree_file(args.file)?;
    let hsts = Arc::new(hsts);

    for (id, value) in ids.iter().zip(data.iter()) {
        var_hm.insert(hsts.get_sample_idx(id)?, *value);
    }

    let args = StandardArgs {
        file: PathBuf::new(),
        output: args.output,
        info_limit: None,
        coords: hsts.metadata.contig.clone(),
        selection: hsts.metadata.selection.clone(),
        prefix: args.prefix,
        samples: None,
    };

    ensure!(
        args.selection == Selection::All || args.selection == Selection::Haploid,
        "Running only with phased data and all chromsomes is supported."
    );

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "quantitative_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    // Define the minimizer/maximiser function for the HST
    let optimizer = |hsts: Arc<HstScan>,
                     hst: &Hst,
                     coord: &Coord,
                     limits: Limits|
     -> Option<(Coord, NodeIndex, f64)> {
        optimizer_inner(hsts, hst, coord, limits, &var_hm)
    };

    // Define what information is wanted for each CSV row
    let rower = |hsts: Arc<HstScan>, coord: &Coord, top_node_idx, optimized_value| -> AssocRow {
        rower_inner(hsts, coord, top_node_idx, optimized_value, &var_hm)
    };

    let write_bam = false;
    let assoc = return_assoc(&hsts, &args, limits, write_bam, optimizer, rower);
    write_assoc(HEADER, assoc, writer)?;

    Ok(())
}

fn rower_inner(
    hsts: Arc<HstScan>,
    coord: &Coord,
    top_node_idx: NodeIndex,
    optimized_value: f64,
    var_hm: &HashMap<usize, f64>,
) -> AssocRow {
    let node = top_node_from_hsts(&hsts.hsts, coord, top_node_idx);

    let other_indexes: Vec<usize> = (0..hsts.nhaplotypes())
        .filter(|i| !node.indexes.contains(i))
        .collect();

    let v1 = get_variable_data(&node.indexes, var_hm);
    let v2 = get_variable_data(&other_indexes, var_hm);

    let cases_mean = v1.iter().sum::<f64>() / v1.len() as f64;
    let ctrls_mean = v2.iter().sum::<f64>() / v2.len() as f64;

    let (nhet, nhom) = zygosity_from_node(node);

    let mut assoc_hm = BTreeMap::new();
    assoc_hm.insert("n_hom".to_string(), nhom.to_string());
    assoc_hm.insert("n_het".to_string(), nhet.to_string());
    assoc_hm.insert("cases_mean".to_string(), cases_mean.to_string());
    assoc_hm.insert("ctrls_mean".to_string(), ctrls_mean.to_string());

    AssocRow::new(hsts.clone(), coord, top_node_idx, optimized_value, assoc_hm)
}

fn optimizer_inner(
    _hsts: Arc<HstScan>,
    hst: &Hst,
    coord: &Coord,
    limits: Limits,
    hm: &HashMap<usize, f64>,
) -> Option<(Coord, NodeIndex, f64)> {
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
            let node = hst.node_weight(idx).unwrap();

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

    Some((coord.clone(), top_node_idx, optimized_value))
}

pub fn get_variable_data(indexes: &[usize], hm: &HashMap<usize, f64>) -> Vec<f64> {
    indexes.iter().filter_map(|v| hm.get(v)).copied().collect()
}
