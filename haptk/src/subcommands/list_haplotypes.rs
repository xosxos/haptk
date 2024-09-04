#![allow(clippy::unwrap_or_default)]

use color_eyre::{eyre::ensure, Result};
use indexmap::{IndexMap, IndexSet};

use crate::{
    args::{Selection, StandardArgs},
    io::{open_csv_writer, push_to_output},
    read_vcf::read_vcf_to_matrix,
    structs::{HapVariant, PhasedMatrix},
    utils::{parse_coords, parse_snp_coord},
};

#[doc(hidden)]
pub fn run(args: StandardArgs, selection_variant: Option<String>) -> Result<()> {
    if args.selection != Selection::All {
        ensure!(
            selection_variant.is_some(),
            "--selection-variant is required if the selection is not All"
        );
    }

    let pos = match selection_variant {
        Some(coords) => parse_snp_coord(&coords)?.1,
        None => 0,
    };

    let (contig, start, stop) = parse_coords(&args.coords)?;

    // let vcf = match (start, stop) {
    //     (Some(start), Some(stop)) => {
    //         read_vcf_to_matrix(&args, contig, 0, Some((Some(start), Some(stop))), None)?
    //     }
    //     _ => read_vcf_to_matrix(&args, contig, 0, None, None)?,
    // };

    let mut vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(&args, contig, pos, Some((start, stop)), None)?
        }
        _ => read_vcf_to_matrix(&args, contig, pos, None, None)?,
    };

    match &args.selection {
        Selection::OnlyAlts | Selection::OnlyRefs => vcf.select_carriers(pos, &args.selection)?,
        Selection::OnlyLongest => vcf.select_only_longest(),
        _ => (),
    }

    ensure!(!vcf.samples().is_empty(), "No samples found in VCF");

    let unique_haplotypes_map = get_unique_haplotypes_map(&vcf)?;

    write_haplotype_file(&args, &unique_haplotypes_map)?;
    write_genotype_file(&args, &unique_haplotypes_map)?;

    Ok(())
}

pub fn get_unique_haplotypes_map(
    vcf: &PhasedMatrix,
) -> Result<IndexMap<Vec<HapVariant>, Vec<String>>> {
    let mut map: IndexMap<Vec<HapVariant>, Vec<String>> = IndexMap::new();

    let samples = vcf.samples().iter().collect::<IndexSet<_>>();

    for sample_name in samples {
        let idxs = vcf.get_sample_idxs(&[sample_name.clone()])?;

        for idx in idxs.iter() {
            let haplotype = vcf.find_haplotype_for_sample(vcf.get_contig(), 0..vcf.ncoords(), *idx);

            map.entry(haplotype)
                .or_insert_with(Vec::new)
                .push(sample_name.clone());
        }
    }
    Ok(map)
}

pub fn write_haplotype_file(
    args: &StandardArgs,
    map: &IndexMap<Vec<HapVariant>, Vec<String>>,
) -> Result<()> {
    let mut sh_output = args.output.clone();
    push_to_output(args, &mut sh_output, "haplotype_file", "csv");
    let mut writer = open_csv_writer(sh_output)?;

    // Header
    let mut header: Vec<String> = vec!["contig".into(), "pos".into(), "ref".into(), "alt".into()];
    let ht_headers = (0..map.len()).map(|i| format!("ht_{i}"));
    header.extend(ht_headers);

    // Init csv
    let mut haplotype_file: Vec<Vec<String>> = vec![];
    haplotype_file.push(header);

    // Push rows to csv
    // Loop a map of vectors so that the first element of all vectors is added to the same row
    for (i, _) in map.get_index(0).unwrap().0.iter().enumerate() {
        let coords = map.get_index(0).unwrap().0[i].clone();

        let mut row = vec![
            coords.contig.to_string(),
            coords.pos.to_string(),
            coords.reference.to_string(),
            coords.alt.to_string(),
        ];

        for n in 0..map.len() {
            let gt = map.get_index(n).unwrap().0[i].gt;
            row.push(gt.to_string());
        }

        haplotype_file.push(row);
    }

    // Write csv
    for line in haplotype_file {
        writer.write_record(&line)?;
    }

    Ok(())
}

pub fn write_genotype_file(
    args: &StandardArgs,
    map: &IndexMap<Vec<HapVariant>, Vec<String>>,
) -> Result<()> {
    let mut sh_output = args.output.clone();
    push_to_output(args, &mut sh_output, "genotype_file", "csv");
    let mut writer = open_csv_writer(sh_output)?;

    let mut gts: IndexMap<String, Vec<String>> = IndexMap::new();

    // Iterate all haplotypes and store into a hashmap. Key: sample name, Value: haplotype ids
    for (i, (_, samples)) in map.iter().enumerate() {
        for sample in samples.iter().cloned() {
            gts.entry(sample)
                .or_insert_with(Vec::new)
                .push(i.to_string());
        }
    }

    let mut gt_file: Vec<Vec<String>> = vec![vec!["id".into(), "ht".into()]];

    for (sample, genotypes) in gts {
        let genotypes = genotypes.join("/");
        gt_file.push(vec![format!("{sample}"), format!("{genotypes}")]);
    }

    for line in gt_file {
        writer.write_record(&line)?;
    }
    Ok(())
}
