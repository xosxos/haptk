use std::{collections::BTreeMap, path::PathBuf, sync::Arc};

use color_eyre::Result;
use petgraph::graph::NodeIndex;

use crate::args::{ConciseArgs, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_recombination_file};
use crate::structs::Coord;
use crate::subcommands::bhst::find_majority_nodes;
use crate::subcommands::scan::{
    read_tree_file, return_assoc, top_node_from_hsts, write_assoc, zygosity_from_node, AssocRow,
    HstScan, Limits,
};
use crate::utils::parse_coords;

use super::Hst;

const HEADER: &[&str] = &[
    "contig",
    "pos",
    "marker_id",
    "bp_len",
    "marker_len",
    "start",
    "stop",
    "mrca",
    "n_het",
    "n_hom",
];

#[doc(hidden)]
pub fn run(args: ConciseArgs, limits: Limits, rec_rates: PathBuf) -> Result<()> {
    let rec_rates = read_recombination_file(rec_rates)?;
    let hsts = read_tree_file(args.file)?;
    let hsts = Arc::new(hsts);

    let args = StandardArgs {
        file: PathBuf::new(),
        output: args.output,
        info_limit: None,
        coords: hsts.metadata.contig.clone(),
        selection: hsts.metadata.selection.clone(),
        prefix: args.prefix,
        samples: None,
        no_alt: false,
    };

    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "branch_mrca_scan", "csv");
    let writer = open_csv_writer(output)?;

    let (_contig, _start, _stop) = parse_coords(&args.coords)?;

    // // Define the minimizer/maximiser function for the HST
    // let optimizer = |hsts: Arc<HstScan>,
    //                  hst_coord: &Coord,
    //                  limits: Limits|
    //  -> Option<(usize, NodeIndex, f64)> {
    //     optimizer_inner(hsts, hst_coord, limits, &rec_rates)
    // };

    // Define the minimizer/maximiser function for the HST
    let optimizer = |hsts: Arc<HstScan>,
                     hst: &Hst,
                     coord: &Coord,
                     limits: Limits|
     -> Option<(Coord, NodeIndex, f64)> {
        optimizer_inner(hsts, hst, coord, limits, &rec_rates)
    };

    // Define what information is wanted for each CSV row
    let rower = |hsts: Arc<HstScan>, coord: &Coord, top_node_idx, optimized_value| -> AssocRow {
        rower_inner(hsts, coord, top_node_idx, optimized_value)
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
) -> AssocRow {
    let node = top_node_from_hsts(&hsts.hsts, coord, top_node_idx);
    let (nhet, nhom) = zygosity_from_node(node);

    let mut hm = BTreeMap::new();
    hm.insert("n_hom".to_string(), nhom.to_string());
    hm.insert("n_het".to_string(), nhet.to_string());

    AssocRow::new(hsts.clone(), coord, top_node_idx, optimized_value, hm)
}

fn optimizer_inner(
    _hsts: Arc<HstScan>,
    hst: &Hst,
    coord: &Coord,
    limits: Limits,
    rec_rates: &BTreeMap<u64, f32>,
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
            // OPTIMIZER CODE

            // Lowest MRCA
            let start_node = idx;
            let nodes = find_majority_nodes(hst, start_node);

            // Iterate nodes from the leaf towards the root
            // NOTE: Sometimes genotyping data runs out and leaf nodes have more than 1 sample
            // NOTE: This is not really accounted for in the Gamma method without chromosome length corrections
            let mut nodes_iter = nodes.iter().map(|(_, idx)| idx).rev();
            let mut prev_idx = nodes_iter.next().unwrap();

            let mut used_indexes = vec![];
            let mut cm_sum = 0.;
            let mut n = 0;

            for idx in nodes_iter {
                let node = hst.node_weight(*idx).unwrap();
                let prev_node = hst.node_weight(*prev_idx).unwrap();

                for idx in &node.indexes {
                    if !used_indexes.contains(&idx) {
                        let (start, stop) = (prev_node.start.pos, prev_node.stop.pos);

                        let (_, start_cm) = rec_rates
                            .range(..=start)
                            .next_back()
                            .unwrap_or_else(|| rec_rates.range(start..).nth(1).unwrap());

                        let (_, stop_cm) = rec_rates
                            .range(stop..)
                            .next()
                            .unwrap_or_else(|| rec_rates.range(..stop).next_back().unwrap());

                        let distance = ((stop_cm - start_cm) / 100.).max(0.);

                        cm_sum += distance;
                        n += 1;

                        used_indexes.push(idx);
                    }
                }
                prev_idx = idx;
            }

            // Calculate MRCA estimate
            let n = n as f32;
            let b_c = (2.0 * n - 1.0) / (2.0 * n);
            let i_tau_hat = (b_c * 2.0 * n) / cm_sum;
            let value = i_tau_hat as f64;

            // Define the clause
            let clause = value < optimized_value;

            // OPTIMIZER CODE END

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
            coord.pos
        );
        return None;
    }

    Some((coord.clone(), top_node_idx, optimized_value))
}
