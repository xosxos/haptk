use std::collections::BTreeMap;
use std::path::PathBuf;

use color_eyre::Result;
use termion::color;

use crate::io::{open_csv_writer, read_haplotype_file};
use crate::structs::HapVariant;
use crate::utils::strip_prefix;

pub fn run(
    haplotypes: Vec<PathBuf>,
    mut output: PathBuf,
    prefix: Option<String>,
    csv: bool,
    hide_missing: bool,
    tag_rows: bool,
) -> Result<()> {
    match strip_prefix(prefix) {
        None => output.push("ht_comparison.csv"),
        Some(prefix) => output.push(format!("{prefix}_ht_comparison.csv")),
    }

    let mut names = vec![];
    let mut hts = vec![];
    for path in haplotypes {
        names.push(path.display().to_string());
        hts.push(read_haplotype_file(path)?);
    }

    let aligner = HaplotypeAligner::new(hts, names);

    match csv {
        true => aligner.to_csv(output)?,
        false => aligner.print(hide_missing, tag_rows),
    }
    Ok(())
}

pub struct HaplotypeAligner {
    names: Vec<String>,
    alignment: BTreeMap<u64, Vec<Option<String>>>,
}

impl HaplotypeAligner {
    pub fn new(haplotypes: Vec<Vec<HapVariant>>, names: Vec<String>) -> Self {
        let mut last_contig = &haplotypes[0][0].contig;

        // Register all positions first
        let mut positions = vec![];
        for ht in &haplotypes {
            let curr_contig = &ht[0].contig;
            if last_contig != curr_contig {
                panic!("Haplotypes are from different contigs {last_contig} vs {curr_contig}")
            }
            last_contig = curr_contig;
            for variant in ht {
                positions.push((variant.pos, vec![None; names.len()]))
            }
        }
        positions.sort();
        positions.dedup();

        // Register present or missing genotypes for each position
        let mut alignment = BTreeMap::from_iter(positions);
        for (k, v) in alignment.iter_mut() {
            for (idx, ht) in haplotypes.iter().enumerate() {
                for var in ht {
                    if k == &var.pos {
                        v[idx] = Some(var.genotype().to_owned());
                    }
                }
            }
        }

        Self { names, alignment }
    }

    pub fn print(&self, hide_missing: bool, tag_rows: bool) {
        // Prepare header
        let names = self
            .names
            .iter()
            .fold(" pos".to_string(), |s, r| format!("{s} {r}"));

        println!(
            "{}{names} status",
            color::Fg(color::AnsiValue::rgb(3, 3, 5))
        );

        for (k, v) in &self.alignment {
            // If missing genotypes present, yellow is true
            // If contradictory genotypes present, red is also true
            let (yellow, red) = self.find_mismatch(v);
            let color = if red {
                color::Fg(color::AnsiValue::rgb(5, 0, 1))
            } else if yellow {
                color::Fg(color::AnsiValue::rgb(5, 5, 1))
            } else {
                color::Fg(color::AnsiValue::rgb(0, 5, 1))
            };

            // Switch Nones to - and fold the vec into a string
            let line = v
                .iter()
                .map(|v| match v {
                    None => "-".to_string(),
                    Some(v) => v.to_owned(),
                })
                .fold(format!(" {k}"), |s, r| format!("{s} {r}"));

            if !hide_missing || !yellow {
                if tag_rows {
                    match (yellow, red) {
                        (_, true) => println!("{color}{line} err"),
                        (true, false) => println!("{color}{line} mis"),
                        (false, false) => println!("{color}{line} ok"),
                    }
                } else {
                    println!("{color}{line}");
                }
            }
        }
    }

    fn find_mismatch(&self, v: &[Option<String>]) -> (bool, bool) {
        let mut values = v.iter().flatten().collect::<Vec<&String>>();
        let yellow = values.len() != v.len();
        values.sort();
        values.dedup();
        let red = values.len() > 1;
        (yellow, red)
    }

    pub fn to_csv(&self, path: PathBuf) -> Result<()> {
        let mut wrtr = open_csv_writer(path)?;
        let mut header = vec!["pos".to_string()];
        header.extend(self.names.clone());

        wrtr.write_record(&header)?;
        for (k, v) in &self.alignment {
            let mut record = vec![format!("{k}")];
            v.iter().for_each(|v| match v {
                None => record.push("-".to_string()),
                Some(v) => record.push(v.to_owned()),
            });
            wrtr.write_record(&record)?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[rustfmt::skip]
mod tests {
    use super::*;

    fn get_ref_alt() -> (String, String) {
        ("A".to_string(), "T".to_string())
    }

    fn create_haplotype(missing: u64, gt: u64) -> Vec<HapVariant> {
        let mut haps = vec![];
        for i in 0..50 {
            if i % missing != 0 {
                let (reference, alt) = get_ref_alt();
                let gt = match i % gt == 0 {
                    true => 1,
                    false => 0,
                };
                let hap = HapVariant {
                    contig: "".to_string(),
                    pos: i,
                    reference,
                    alt,
                    gt,
                };
                haps.push(hap);
            }
        }
        haps
    }

    #[test]
    fn test_aligner_alignment() {
        let hap1 = create_haplotype(11, 9);
        let hap2 = create_haplotype(2, 10);
        let hap3 = create_haplotype(5, 11);
        let aligner = HaplotypeAligner::new(
            vec![hap1, hap2, hap3],
            vec!["hap1".to_string(), "hap2".to_string(), "hap3".to_string()],
        );
        let pos: Vec<u64> = (1..50).collect();
        let colors = [(false, false),
            (true, false),
            (false, false),
            (true, false),
            (true, false),
            (true, false),
            (false, false),
            (true, false),
            (false, true),
            (true, false)];

        for (i, (k, v)) in aligner.alignment.iter().enumerate() {
            assert_eq!(&pos[i], k);
            assert_eq!(colors[i], aligner.find_mismatch(v));

            if i == 9 {
                break;
            }
        }
    }

    #[test]
    fn test_aligner_mismatch() {
        let aligner = HaplotypeAligner { names : vec![], alignment: BTreeMap::new() };

        let foo: Vec<Option<String>> = vec![Some("A".into()), None, Some("T".into())];
        let (yellow, red) = aligner.find_mismatch(&foo);
        assert_eq!((yellow, red), (true, true));

        let foo: Vec<Option<String>> = vec![Some("A".into()), Some("A".into()), Some("T".into())];
        let (yellow, red) = aligner.find_mismatch(&foo);
        assert_eq!((yellow, red), (false, true));

        let foo: Vec<Option<String>> = vec![Some("A".into()), None, None];
        let (yellow, red) = aligner.find_mismatch(&foo);
        assert_eq!((yellow, red), (true, false));

        let foo: Vec<Option<String>> = vec![None, None, None];
        let (yellow, red) = aligner.find_mismatch(&foo);
        assert_eq!((yellow, red), (true, false));

        let foo: Vec<Option<String>> = vec![Some("A".into()), Some("A".into()), Some("A".into())];
        let (yellow, red) = aligner.find_mismatch(&foo);
        assert_eq!((yellow, red), (false, false));
    }
}
