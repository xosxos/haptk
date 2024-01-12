use color_eyre::Result;

use crate::{
    args::StandardArgs,
    io::{open_csv_writer, push_to_output},
    read_vcf::read_vcf_to_matrix,
    structs::HapVariant,
    utils::parse_coords,
};

#[doc(hidden)]
#[tracing::instrument]
pub fn run(args: StandardArgs) -> Result<()> {
    let (contig, start, stop) = parse_coords(&args.coords)?;

    let vcf = match (start, stop) {
        (Some(start), Some(stop)) => {
            read_vcf_to_matrix(&args, contig, 0, Some((start, stop)), None)?
        }
        _ => read_vcf_to_matrix(&args, contig, 0, None, None)?,
    };

    for sample_name in vcf.samples() {
        let idxs = vcf.get_sample_idxs(&[sample_name.clone()])?;

        for (i, idx) in idxs.iter().enumerate() {
            let mut sh_output = args.output.clone();

            push_to_output(
                &args,
                &mut sh_output,
                &format!("{sample_name}_haplotype_{}", i + 1),
                "csv",
            );

            let mut writer = open_csv_writer(sh_output)?;

            let ht = vcf.find_haplotype_for_sample(contig, 0..vcf.ncoords(), *idx);
            write_haplotype(&mut writer, ht)?;
        }
    }

    Ok(())
}

#[doc(hidden)]
fn write_haplotype(
    writer: &mut csv::Writer<Box<dyn std::io::Write + 'static>>,
    haplotype: Vec<HapVariant>,
) -> Result<()> {
    for variant in haplotype {
        writer.write_record(vec![
            &variant.contig.to_string(),
            &variant.pos.to_string(),
            &variant.reference,
            &variant.alt,
            &variant.gt.to_string(),
        ])?;
    }
    Ok(())
}
