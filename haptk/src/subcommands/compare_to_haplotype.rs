#![allow(clippy::comparison_chain)]
use std::{collections::BTreeSet, path::PathBuf};

use color_eyre::{
    eyre::{eyre, WrapErr},
    Result,
};
use ndarray::parallel::prelude::*;
use ndarray::{s, Axis};
use ndarray::{Array2, ShapeBuilder};

use crate::{
    args::{GraphArgs, Selection, SortOption, StandardArgs},
    graphs::MatrixGraph,
    io::{open_csv_writer, push_to_output, read_haplotype_file, read_sample_ids},
    read_vcf::read_vcf_to_matrix,
    structs::{HapVariant, PhasedMatrix},
    utils::parse_snp_coord,
};
// use crate::graphs::matrix_graph::matrix_graph_png;

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn run(
    args: StandardArgs,
    haplotype_path: PathBuf,
    decoy_samples: Option<PathBuf>,
    mark_shorter_alleles: bool,
    want_png: bool,
    want_npy: bool,
    graph_args: GraphArgs,
    sort_option: SortOption,
) -> Result<()> {
    if want_png {
        return Err(eyre!("PNG is not supported for matrix graph."));
    }
    if args.selection == Selection::Unphased {
        return Err(eyre!("Running with unphased data is not supported."));
    }

    // File reads
    let decoy_samples = read_sample_ids(&decoy_samples)?;
    let ht = read_haplotype_file(haplotype_path)?;

    // IMG output
    let mut img_output = args.output.clone();
    push_to_output(
        &args,
        &mut img_output,
        match mark_shorter_alleles {
            true => "ht_comparison_shorter_alleles_marked",
            false => "ht_comparison",
        },
        match want_png {
            true => "png",
            false => "svg",
        },
    );
    svg::save(&img_output, &svg::Document::new())
        .wrap_err(eyre!("Failed writing to {:?}", img_output))?;

    // CSV output
    let mut output = args.output.clone();
    push_to_output(&args, &mut output, "ht_shared_segments", "csv");
    let mut writer = open_csv_writer(output)?;

    // VCF read
    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;
    let (start, end) = (ht.first().unwrap(), ht.last().unwrap());

    let mut only_longest = None;
    let vcf = match args.selection {
        Selection::All => {
            let vcf = read_vcf_to_matrix(
                &args,
                contig,
                variant_pos,
                Some((Some(start.pos), Some(end.pos))),
                None,
                false,
            )?;
            let mut vcf = transform_gt_matrix_to_match_matrix(vcf, &ht, variant_pos)?;
            only_longest = Some(vcf.only_longest_indexes()?);
            vcf
        }
        Selection::OnlyAlts | Selection::OnlyRefs => {
            let vcf = read_vcf_to_matrix(
                &args,
                contig,
                variant_pos,
                Some((Some(start.pos), Some(end.pos))),
                None,
                false,
            )?;

            transform_gt_matrix_to_match_matrix(vcf, &ht, variant_pos)?
        }
        Selection::OnlyLongest => {
            let vcf = read_vcf_to_matrix(
                &args,
                contig,
                variant_pos,
                Some((Some(start.pos), Some(end.pos))),
                None,
                false,
            )?;

            let mut vcf = transform_gt_matrix_to_match_matrix(vcf, &ht, variant_pos)?;

            // Do only longest selection after swithing to the given haplotype as the refernce
            vcf.select_only_longest()?;

            // Use start and end from the haplotype to select columns from the matrix by range

            let start = vcf.get_first_idx_on_left_by_pos(start.pos);
            let stop = vcf.get_first_idx_on_right_by_pos(end.pos);

            vcf.select_columns_by_range_idx(start..=stop);
            vcf
        }
        Selection::Haploid => {
            let vcf = read_vcf_to_matrix(
                &args,
                contig,
                variant_pos,
                Some((Some(start.pos), Some(end.pos))),
                None,
                false,
            )?;

            transform_gt_matrix_to_match_matrix(vcf, &ht, variant_pos)?
        }
        Selection::Unphased => unreachable!(),
    };

    if ht.len() > vcf.ncoords() {
        tracing::warn!(
            "Haplotype has more variants than the VCF {} vs {}",
            ht.len(),
            vcf.ncoords()
        );
    } else if ht.len() < vcf.ncoords() {
        tracing::warn!(
            "Haplotype has less variants than the VCF {} vs {}",
            ht.len(),
            vcf.ncoords()
        );
    }

    if want_npy {
        let mut npy_output = args.output.clone();
        push_to_output(&args, &mut npy_output, "ht_comparison", "npy");
        ndarray_npy::write_npy(npy_output, vcf.matrix())?;
    }

    let index_order = sort_indexes_for_diff_graph(
        &vcf,
        decoy_samples.as_ref(),
        mark_shorter_alleles,
        only_longest.as_ref(),
        sort_option,
    );

    let shared_ranges = find_shared_haplotype_ranges(&vcf);

    print_ranges_to_csv(
        &vcf,
        &shared_ranges,
        decoy_samples.as_ref(),
        only_longest.as_ref(),
        &mut writer,
    )?;

    tracing::info!(
        "Sample haplotypes: {}, average length: {}, median length: {}",
        shared_ranges.len(),
        range_length_avg(&shared_ranges),
        range_length_median(&shared_ranges)
    );

    let mut dg = MatrixGraph::new(
        &vcf,
        graph_args,
        decoy_samples,
        mark_shorter_alleles,
        only_longest.as_ref(),
    );
    dg.draw_graph(&index_order);
    svg::save(img_output, &dg.document)?;
    tracing::debug!("Finished drawing graph");

    Ok(())
}

pub fn transform_gt_matrix_to_match_matrix(
    mut vcf: PhasedMatrix,
    ht: &[HapVariant],
    variant_pos: u64,
) -> Result<PhasedMatrix> {
    let mut match_matrix: Vec<u8> = vec![];
    let mut coords_na = 0;
    let mut coords = BTreeSet::new();
    vcf.matrix_axis_iter(1)
        .enumerate()
        .for_each(|(var_idx, col)| {
            if let Some(hv) = ht.iter().find(|hv| *hv == vcf.get_coord(var_idx)) {
                coords.insert(vcf.get_coord(var_idx).clone());
                match_matrix.extend(col.iter().map(|gt| match hv.gt == *gt {
                    true => 1,
                    false => 0,
                }));
            } else {
                coords_na += 1;
                tracing::trace!(
                    "Coord not found in the haplotype: {}",
                    &vcf.get_coord(var_idx)
                );
            }
        });

    if coords_na != 0 {
        tracing::warn!(
            "In the haplotype range {coords_na} markers were present in the VCF, but not in the haplotype. These were disregarded."
        );
    }

    tracing::info!(
        "The haplotype and the VCF had {} markers in common.",
        coords.len()
    );

    let vcf_coords = vcf.coords_mut();
    *vcf_coords = coords;
    let offset_coord_start = vcf_coords.first().unwrap().clone();

    vcf.set_matrix(
        offset_coord_start,
        Array2::from_shape_vec((vcf.nhaplotypes(), vcf.ncoords()).f(), match_matrix)?,
    );

    vcf.start_coord = vcf.get_nearest_coord_by_pos(variant_pos).clone();
    vcf.variant_idx = vcf.get_coord_idx(vcf.start_coord());
    // vcf.variant_idx = vcf.get_first_idx_on_right_by_pos(variant_pos);

    tracing::debug!("Finished transforming the genotype matrix to a match matrix.");
    tracing::info!(
        "Constructed a match matrix from {} records. Starting variant position: {}.",
        vcf.ncoords(),
        vcf.variant_idx_pos(),
    );
    Ok(vcf)
}

fn sort_indexes_for_diff_graph(
    vcf: &PhasedMatrix,
    decoy_samples: Option<&Vec<String>>,
    mark_shorter_alleles: bool,
    only_longest: Option<&Vec<usize>>,
    sort_option: SortOption,
) -> Vec<usize> {
    match sort_option {
        SortOption::Left => tracing::info!(
            "Sorting the samples based on haplotype sharing to the left of the starting locus",
        ),
        SortOption::Right => tracing::info!(
            "Sorting the samples based on haplotype sharing to the right of the starting locus",
        ),
        SortOption::Total => tracing::info!(
            "Sorting the samples based on haplotype sharing on both sides of the starting locus",
        ),
    }

    let mut values: Vec<(usize, i32)> = match sort_option {
        // Take all from start to variant index, reverse and calculate 1 count in parallel
        // to get the amount of markers shared to the left
        SortOption::Left => vcf
            .slice_cols(0..vcf.variant_idx() + 1)
            .axis_iter(Axis(0))
            .into_par_iter()
            .enumerate()
            .map(|(y, row)| {
                let mut count = 0;
                for i in row.iter().rev() {
                    match i {
                        0 => break,
                        1 => count += 1,
                        _ => panic!(),
                    }
                }
                (y, count)
            })
            .collect(),
        SortOption::Right => vcf
            .slice_cols(vcf.variant_idx()..vcf.ncoords())
            // .slice(s![.., vcf.variant_idx()..vcf.ncoords()])
            .axis_iter(Axis(0))
            .into_par_iter()
            .enumerate()
            .map(|(y, row)| {
                let mut count = 0;
                for i in row.iter() {
                    match i {
                        0 => break,
                        1 => count += 1,
                        _ => panic!(),
                    }
                }
                (y, count)
            })
            .collect(),
        SortOption::Total => vcf
            .slice_cols(0..vcf.ncoords())
            .axis_iter(Axis(0))
            .into_par_iter()
            .enumerate()
            .map(|(y, row)| {
                let mut count = 0;
                let right = row.into_iter().skip(vcf.variant_idx());
                let left = row.slice(s![0..vcf.variant_idx() + 1]);
                for i in right {
                    match i {
                        0 => break,
                        1 => count += 1,
                        _ => panic!(),
                    }
                }
                for i in left.iter().rev() {
                    match i {
                        0 => break,
                        1 => count += 1,
                        _ => panic!(),
                    }
                }
                (y, count)
            })
            .collect(),
    };

    // Sort shortest alleles to the top
    if mark_shorter_alleles {
        let only_longest =
            only_longest.expect("Failed to select the only-longest haplotypes for tagging. For this to succeed, the option `--select all` should be used.");
        values.sort_by(|a, b| {
            only_longest
                .contains(&a.0)
                .cmp(&only_longest.contains(&b.0))
        });
    }

    // Sort decoy alleles to the top
    if let Some(samples) = decoy_samples {
        values.sort_by(|b, a| {
            samples
                .contains(&vcf.get_sample_name(a.0))
                .cmp(&samples.contains(&vcf.get_sample_name(b.0)))
        });
    }

    // Sort by length
    values.sort_by(|a, b| a.1.cmp(&b.1));

    tracing::debug!("Finished sorting.");

    values.iter().map(|v| v.0).collect::<Vec<usize>>()
}

pub fn find_shared_haplotype_ranges(vcf: &PhasedMatrix) -> Vec<(usize, usize, usize)> {
    vcf.matrix_axis_iter(0)
        .into_par_iter()
        .enumerate()
        .map(|(idx, row)| {
            let (mut last_start, mut last_end) = (None, None);
            for (i, gt) in row.iter().enumerate() {
                match gt {
                    1 => {
                        if last_start.is_some() {
                            last_end = Some(i);
                        } else {
                            last_start = Some(i);
                            last_end = Some(i);
                        }
                    }
                    0 => {
                        if i == vcf.variant_idx() {
                            return (idx, i, i);
                        }
                        if let (Some(end), Some(start)) = (last_end, last_start) {
                            if vcf.variant_idx() <= end && vcf.variant_idx() >= start {
                                return (idx, start, end);
                            }
                        }

                        (last_start, last_end) = (None, None);
                    }
                    _ => panic!(),
                }
            }
            //
            (idx, last_start.unwrap(), vcf.ncoords() - 1)
            // panic!("No haplotye found in check-diff algo")
        })
        .collect()
}

pub fn range_length_avg(ranges: &[(usize, usize, usize)]) -> f32 {
    // ranges.iter().for_each(|n| println!("{}", n.2 - n.1));
    let sum: usize = ranges.iter().map(|n| n.2 - n.1).sum();
    sum as f32 / ranges.len() as f32
}

pub fn range_length_median(ranges: &[(usize, usize, usize)]) -> usize {
    let mut vec: Vec<_> = ranges.iter().map(|n| n.2 - n.1).collect();
    vec.sort();
    vec[vec.len() / 2]
}

pub fn print_ranges_to_csv(
    vcf: &PhasedMatrix,
    ranges: &[(usize, usize, usize)],
    decoy_samples: Option<&Vec<String>>,
    only_longest: Option<&Vec<usize>>,
    writer: &mut csv::Writer<Box<dyn std::io::Write>>,
) -> Result<()> {
    writer.write_record(vec![
        "id",
        "start",
        "stop",
        "length",
        "markers",
        "is_marked",
        "is_longest",
    ])?;
    for (idx, start, stop) in ranges {
        let start_pos = vcf.get_pos(*start);
        let stop_pos = vcf.get_pos(*stop);

        let mut row = vec![
            format!("{}", vcf.get_sample_name(*idx)),
            start_pos.to_string(),
            stop_pos.to_string(),
            (stop_pos - start_pos + 1).to_string(),
            (stop - start + 1).to_string(),
        ];

        match decoy_samples {
            Some(decoys) => match decoys.contains(&vcf.get_sample_name(*idx)) {
                true => row.push("true".to_string()),
                false => row.push("false".to_string()),
            },
            None => row.push("na".to_string()),
        }
        match only_longest {
            Some(longest) => match longest.contains(idx) {
                true => row.push("true".to_string()),
                false => row.push("false".to_string()),
            },
            None => row.push("na".to_string()),
        }
        writer.write_record(row)?;
    }
    Ok(())
}
