use std::collections::BTreeMap;
use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use itertools::Itertools;
use plotters::prelude::*;
use rayon::prelude::*;

use crate::args::{GraphArgs, Selection, StandardArgs};
use crate::io::{open_csv_writer, push_to_output, read_recombination_file};
use crate::read_vcf::read_vcf_to_matrix;
use crate::structs::PhasedMatrix;
use crate::subcommands::mrca::mrca_gamma_method;
use crate::utils::{centromeres_hg38, parse_coords};

#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn run(
    args: StandardArgs,
    rec_rates: PathBuf,
    step_size: usize,
    no_csv: bool,
    yes_plot: bool,
    graph_args: GraphArgs,
    mark_centromere: bool,
) -> Result<()> {
    match args.selection {
        Selection::OnlyAlts | Selection::OnlyRefs | Selection::Unphased => {
            return Err(eyre!(
                "Running with selection {:?} is not supported.",
                args.selection
            ))
        }
        _ => (),
    };

    let (contig, start, stop) = parse_coords(&args.coords)?;

    let range = match (start, stop) {
        (Some(start), Some(stop)) => Some((start, stop)),
        (None, None) => None,
        _ => panic!("Open ended ranges are not supported for now"),
    };

    let vcf = read_vcf_to_matrix(&args, contig, 0, range, None)?;

    let rates = read_recombination_file(rec_rates)?;

    let ages = find_ages(&vcf, &args, rates, step_size)?;

    if !no_csv {
        write_ages_to_csv(&args, &ages, args.output.clone(), &vcf)?;
    }

    if yes_plot {
        draw_plot(
            &vcf,
            &args,
            graph_args,
            args.output.clone(),
            ages,
            contig,
            mark_centromere,
        )?;
    }

    Ok(())
}

fn write_ages_to_csv(
    args: &StandardArgs,
    ages: &Vec<(usize, f64)>,
    mut output: PathBuf,
    vcf: &PhasedMatrix,
) -> Result<()> {
    push_to_output(args, &mut output, "mrca_scan", "csv");

    let mut writer = open_csv_writer(output)?;
    writer.write_record(vec!["POS", "mrca"])?;
    for (idx, age) in ages {
        writer.write_record(vec![format!("{}", vcf.get_pos(*idx)), format!("{age:.5}")])?;
    }
    Ok(())
}

fn find_ages(
    vcf: &PhasedMatrix,
    args: &StandardArgs,
    rates: BTreeMap<u64, f32>,
    step_size: usize,
) -> Result<Vec<(usize, f64)>> {
    // Create a Vec because par_iter cannot be used with pure ranges and par_bridge does not
    // return in ordered fashion with .collect()
    let range: Vec<usize> = (0..vcf.ncoords()).collect();

    match args.selection == Selection::OnlyLongest {
        true => range
            .par_iter()
            .filter(|n| *n % step_size == 0)
            .map(|i| -> Result<(usize, f64)> {
                let only_longest_lengths = vcf.only_longest_lengths(*i);
                let ((i_tau_hat, _, _), _) =
                    mrca_gamma_method(vcf, only_longest_lengths, vcf.get_pos(*i), &rates)?;
                Ok((*i, i_tau_hat))
            })
            .collect(),
        false => range
            .par_iter()
            .filter(|n| *n % step_size == 0)
            .map(|i| -> Result<(usize, f64)> {
                let shared_lengths = vcf.get_lengths_from_uhst(*i);
                let ((i_tau_hat, _, _), _) =
                    mrca_gamma_method(vcf, shared_lengths, vcf.get_pos(*i), &rates)?;
                Ok((*i, i_tau_hat))
            })
            .collect(),
    }
}

fn draw_plot(
    vcf: &PhasedMatrix,
    args: &StandardArgs,
    graph_args: GraphArgs,
    mut output: PathBuf,
    data: Vec<(usize, f64)>,
    contig: &str,
    mark_centromere: bool,
) -> Result<()> {
    push_to_output(args, &mut output, "mrca_scan", "png");

    let root = BitMapBackend::new(&output, (graph_args.width as u32, graph_args.height as u32))
        .into_drawing_area();

    let y_start = data.iter().min_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;
    let y_stop = data.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().1;

    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .set_label_area_size(LabelAreaPosition::Left, 60)
        .set_label_area_size(LabelAreaPosition::Right, 60)
        .set_label_area_size(LabelAreaPosition::Bottom, 40)
        .build_cartesian_2d(0..vcf.ncoords(), y_start - 0.2..y_stop + 0.2)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_labels(30)
        .max_light_lines(4)
        .y_desc("MRCA")
        .x_desc("POS")
        .draw()?;

    // Raw log10 data line
    // chart.draw_series(LineSeries::new(
    // data.iter().map(|(pos, mrca)|(*pos,*mrca)),
    // &BLUE,
    // ))?;

    // Last 10 elements moving average
    let mut moving_average = vec![];
    for (a, b, c, d, e, f, g, h, i, j) in data.iter().tuples() {
        let sum = a.1 + b.1 + c.1 + d.1 + e.1 + f.1 + g.1 + h.1 + i.1 + j.1;
        let mean = sum / 10.0;
        moving_average.push((j.0, mean));
    }

    chart.draw_series(LineSeries::new(
        moving_average.iter().map(|(pos, mrca)| (*pos, *mrca)),
        ShapeStyle {
            color: BLACK.mix(0.9),
            filled: true,
            stroke_width: graph_args.stroke_width,
        },
    ))?;

    // Centromere location
    if mark_centromere {
        let (centro_start, centro_stop) = centromeres_hg38(contig);
        let centro_start = vcf.get_nearest_idx_by_pos(centro_start);
        let centro_stop = vcf.get_nearest_idx_by_pos(centro_stop);

        chart.draw_series(LineSeries::new(
            [(centro_start, y_start), (centro_start, y_stop)].into_iter(),
            ShapeStyle {
                color: GREEN.mix(0.6),
                filled: true,
                stroke_width: 2,
            },
        ))?;

        chart.draw_series(LineSeries::new(
            [(centro_stop, y_start), (centro_stop, y_stop)].into_iter(),
            ShapeStyle {
                color: GREEN.mix(0.6),
                filled: true,
                stroke_width: 2,
            },
        ))?;
    }

    // Variant location
    if graph_args.mark_locus {
        chart.draw_series(LineSeries::new(
            [(vcf.variant_idx(), y_start), (vcf.variant_idx(), y_stop)].into_iter(),
            ShapeStyle {
                color: BLUE.mix(0.6),
                filled: true,
                stroke_width: 2,
            },
        ))?;
    }

    // Mean line
    let sum: f64 = data.iter().map(|(_pos, mrca)| mrca).sum();
    let mean = sum / data.len() as f64;
    chart.draw_series(LineSeries::new(
        [(0, mean), (vcf.ncoords(), mean)].into_iter(),
        ShapeStyle {
            color: RED.mix(0.6),
            filled: true,
            stroke_width: graph_args.stroke_width / 2,
        },
    ))?;

    root.present().expect("Unable to write result to file");

    Ok(())
}
