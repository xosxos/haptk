#![allow(clippy::unwrap_or_default)]

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use indexmap::{IndexMap, IndexSet};

use crate::{
    args::{Selection, StandardArgs},
    io::{get_htslib_contig_len, open_csv_writer, push_to_output},
    read_vcf::{read_vcf_to_matrix, read_vcf_to_matrix_by_indexes},
    structs::{HapVariant, PhasedMatrix},
    subcommands::check_for_haplotype::identical_haplotype_count,
    utils::{parse_coords, parse_snp_coord, precision_f64},
};

type HaplotypeMap = IndexMap<(Vec<HapVariant>, usize), Vec<String>>;

#[doc(hidden)]
pub fn run(args: StandardArgs, selection_variant: Option<String>, nucleotides: bool) -> Result<()> {
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

    let range = match (start, stop) {
        (Some(_), Some(_)) => Some((start, stop)),
        _ => None,
    };

    let vcf = if args.selection == Selection::OnlyLongest {
        let (only_longest_lookups, vcf) = if get_htslib_contig_len(&args.file, contig).is_ok() {
            let mut vcf =
                read_vcf_to_matrix(&args, contig, pos, None, None, Some(5_000_000), false)?;
            let lookups = vcf.get_only_longest_lookups()?;
            (lookups, vcf)
        } else {
            let vcf = read_vcf_to_matrix(&args, contig, pos, None, None, None, false)?;
            let lookups = vcf.get_only_longest_lookups_no_shard()?;
            (lookups, vcf)
        };

        read_vcf_to_matrix_by_indexes(
            &args.file,
            pos,
            contig,
            range,
            vcf.samples().clone(),
            vcf.metadata.indexes,
            only_longest_lookups,
            None,
            args.no_alt,
            &Selection::Haploid,
            false,
        )?
    } else {
        read_vcf_to_matrix(&args, contig, pos, range, None, None, false)?
    };

    ensure!(!vcf.samples().is_empty(), "No samples found in VCF");

    let unique_haplotypes_map = get_unique_haplotypes_map(&vcf)?;

    write_haplotype_file(&vcf, &args, &unique_haplotypes_map, nucleotides)?;
    write_genotype_file(&args, &unique_haplotypes_map)?;

    Ok(())
}

pub fn get_unique_haplotypes_map(vcf: &PhasedMatrix) -> Result<HaplotypeMap> {
    let mut map: HaplotypeMap = IndexMap::new();

    let samples = vcf.samples().iter().collect::<IndexSet<_>>();

    for sample_name in samples {
        let idxs = vcf.get_idxs_for_samples(&[sample_name.clone()])?;

        for idx in idxs.iter() {
            let haplotype = vcf.find_haplotype_for_sample(
                vcf.coords().first().unwrap()..=vcf.coords().last().unwrap(),
                *idx,
            );

            let matching_indexes = identical_haplotype_count(vcf, &haplotype);

            let n = matching_indexes.len();

            map.entry((haplotype, n))
                .or_insert_with(Vec::new)
                .push(sample_name.clone());
        }
    }

    map.sort_by(|ak, _, bk, _| bk.1.cmp(&ak.1));

    Ok(map)
}

pub fn write_haplotype_file(
    vcf: &PhasedMatrix,
    args: &StandardArgs,
    map: &HaplotypeMap,
    nucleotides: bool,
) -> Result<()> {
    let mut sh_output = args.output.clone();

    let core_name = match nucleotides {
        true => "haplotype_file_nucleotides",
        false => "haplotype_file",
    };

    push_to_output(args, &mut sh_output, core_name, "csv");
    let mut writer = open_csv_writer(sh_output)?;

    let mut header = match nucleotides {
        true => ["contig", "pos", "ref"].map(String::from).to_vec(),
        false => ["contig", "pos", "ref", "alt"].map(String::from).to_vec(),
    };

    // Header
    let ht_headers = (0..map.len()).map(|i| format!("ht_{i}"));
    header.extend(ht_headers);

    // Init csv
    let mut haplotype_file: Vec<Vec<String>> = vec![];
    haplotype_file.push(header);

    // Push rows to csv
    // Loop a map of vectors so that the first element of all vectors is added to the same row
    for (i, _) in map.get_index(0).unwrap().0 .0.iter().enumerate() {
        let coord = map.get_index(0).unwrap().0 .0[i].clone();

        let mut row = match nucleotides {
            true => vec![
                coord.contig.to_string(),
                coord.pos.to_string(),
                coord.reference.to_string(),
            ],

            false => vec![
                coord.contig.to_string(),
                coord.pos.to_string(),
                coord.reference.to_string(),
                coord.alt.to_string(),
            ],
        };

        for n in 0..map.len() {
            let coord = &map.get_index(n).unwrap().0 .0[i];

            let gt = match nucleotides {
                true => match coord.gt {
                    0 => coord.reference.clone(),
                    1 => coord.alt.clone(),
                    _ => return Err(eyre!("multiallelic haplotypes not yet allowed")),
                },
                false => coord.gt.to_string(),
            };

            row.push(gt.to_string());
        }

        haplotype_file.push(row);
    }

    // Write csv
    for line in haplotype_file {
        writer.write_record(&line)?;
    }

    let mut freq_line = match nucleotides {
        true => ["freq", "-", "-"].map(String::from).to_vec(),
        false => ["freq", "-", "-", "-"].map(String::from).to_vec(),
    };
    let mut n_line = match nucleotides {
        true => ["n", "-", "-"].map(String::from).to_vec(),
        false => ["n", "-", "-", "-"].map(String::from).to_vec(),
    };

    for ((_ht, n), _sample) in map {
        n_line.push(n.to_string());
        freq_line.push(format!(
            "{}",
            precision_f64(*n as f64 / vcf.nhaplotypes() as f64, 3)
        ));
    }
    writer.write_record(freq_line)?;
    writer.write_record(n_line)?;

    Ok(())
}

pub fn write_genotype_file(args: &StandardArgs, map: &HaplotypeMap) -> Result<()> {
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
