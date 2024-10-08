use std::path::PathBuf;

use color_eyre::{
    eyre::eyre,
    Result,
};

use crate::{
    args::{Selection, StandardArgs},
    io::{open_csv_writer, push_to_output},
    read_vcf::read_vcf_to_matrix,
    structs::{HapVariant, PhasedMatrix},
    utils::{parse_snp_coord, precision_f64},
};

#[doc(hidden)]
pub fn run(args: StandardArgs, haplotype_path: PathBuf) -> Result<()> {
    let (contig, variant_pos) = parse_snp_coord(&args.coords)?;

    let mut csv_output = args.output.clone();
    push_to_output(&args, &mut csv_output, "haplotype_check", "csv");
    let mut writer = open_csv_writer(csv_output)?;

    let ht = crate::io::read_haplotype_file(haplotype_path.clone())?;
    let start = ht.first().ok_or_else(|| eyre!( 
        "Failed to get the first variant of the haplotype at {haplotype_path:?}. Is the haplotype file empty?"),
    )?;

    let end = ht.last().unwrap();

    let mut vcf = read_vcf_to_matrix(
        &args,
        contig,
        variant_pos,
        Some((Some(start.pos), Some(end.pos))),
        None,
        false,
    )?;
    match args.selection {
        Selection::OnlyLongest => vcf.select_only_longest()?,
        Selection::Unphased => return Err(eyre!("Running with unphased data is not supported")),
        _ => (),
    };

    let matching_indexes = identical_haplotype_count(&vcf, &ht);
    write_matches_to_csv(&matching_indexes, &mut writer, &vcf)?;

    tracing::debug!(
        "Haplotype found in {}/{} (freq: {}) of alleles.",
        matching_indexes.len(),
        vcf.nhaplotypes(),
        precision_f64(matching_indexes.len() as f64 / vcf.nhaplotypes() as f64, 3),
    );

    Ok(())
}

fn write_matches_to_csv(
    matching_indexes: &[usize],
    writer: &mut csv::Writer<Box<dyn std::io::Write>>,
    vcf: &PhasedMatrix,
) -> Result<()> {
    writer.write_record(vec!["id", "match"])?;
    for idx in 0..vcf.nhaplotypes() {
        writer.write_record(vec![
            format!("{}", vcf.get_sample_name(idx)),
            format!("{}", matching_indexes.contains(&idx)),
        ])?;
    }
    Ok(())
}

pub fn identical_haplotype_count(vcf: &PhasedMatrix, ht: &[HapVariant]) -> Vec<usize> {
    let indexes: Vec<(usize, usize)> = ht
        .iter()
        .enumerate()
        .filter_map(|(ht_idx, h)|
            match vcf.try_coord_by_hapvariant(h) {
                    Some(matrix_idx) => Some((matrix_idx, ht_idx)),
                    None => {
                        tracing::warn!("Haplotype variant {:?} not found in the vcf", h);
                        None
                    }   
                }
        )
        .collect();

    if indexes.is_empty() {
        return vec![];
    }

    let mut matching = vec![true; vcf.nhaplotypes()];

    for (n, value) in matching.iter_mut().enumerate() {
        'inner: for (matrix_idx, ht_idx) in &indexes {
            if vcf.matrix_point(n, *matrix_idx) != ht[*ht_idx].gt {
                *value = false;
                break 'inner;
            }
        }
    }

    matching
        .iter()
        .enumerate()
        .filter_map(|(i, m)| if *m { Some(i) } else { None })
        .collect()
}
