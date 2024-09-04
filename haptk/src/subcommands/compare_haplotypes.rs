use std::collections::btree_map::Entry;

use std::collections::{BTreeMap, HashSet};
use std::io::{self};
use std::path::PathBuf;

use color_eyre::{eyre::eyre, Result};
use csv::{Reader, ReaderBuilder};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use termion::color;

use crate::io::{get_input, open_csv_writer};
use crate::structs::Coord;
use crate::utils::strip_prefix;

pub fn get_csv_reader<R: io::Read>(input: R) -> Reader<R> {
    ReaderBuilder::new()
        .delimiter(b',')
        .has_headers(true)
        .flexible(false)
        .from_reader(input)
}

#[derive(Debug, Default, Clone, PartialOrd, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CompareHaplotype {
    pub contig: String,
    pub pos: u64,
    pub ann: Option<String>,
    pub reference: String,
    pub alt: String,
    pub gts: Vec<String>,
}

impl From<&CompareHaplotype> for Coord {
    fn from(variant: &CompareHaplotype) -> Self {
        Coord {
            contig: variant.contig.clone(),
            pos: variant.pos,
            reference: variant.reference.clone(),
            alt: variant.alt.clone(),
        }
    }
}

impl PartialEq<CompareHaplotype> for Coord {
    fn eq(&self, other: &CompareHaplotype) -> bool {
        let clause1 = self.contig == other.contig
            && self.pos == other.pos
            && self.reference == other.reference;

        let clause2 = self.alt == other.alt || self.alt == "-" || other.alt == "-";

        clause1 && clause2
    }
}

impl PartialEq<Coord> for CompareHaplotype {
    fn eq(&self, other: &Coord) -> bool {
        let clause1 = self.contig == other.contig
            && self.pos == other.pos
            && self.reference == other.reference;

        let clause2 = self.alt == other.alt || self.alt == "-" || other.alt == "-";

        clause1 && clause2
    }
}

impl From<IndexMap<String, String>> for CompareHaplotype {
    fn from(hm: IndexMap<String, String>) -> Self {
        CompareHaplotype {
            contig: hm["contig"].clone(),
            pos: hm["pos"]
                .parse::<u64>()
                .unwrap_or_else(|_| panic!("Pos {:?} is not an integer", hm["pos"])),
            ann: hm.get("ann").cloned(),
            reference: hm["ref"].clone(),
            alt: hm["alt"].clone(),
            gts: hm
                .keys()
                .filter(|k| k.starts_with("ht") || k.starts_with("gt"))
                .map(|k| hm[k].clone())
                .collect(),
        }
    }
}

pub fn read_haplotype_file(ht_path: PathBuf) -> Result<(Vec<CompareHaplotype>, Vec<String>)> {
    let mut rdr = get_csv_reader(get_input(Some(ht_path.clone()))?);

    let mut ht_names: Vec<String> = vec![];
    let line = rdr.deserialize().next().unwrap();
    let record: IndexMap<String, String> = line?;
    ht_names.extend(
        record
            .keys()
            .filter(|k| k.starts_with("ht") || k.starts_with("gt"))
            // .inspect(|&k| println!("{k:?}"))
            .cloned(),
    );

    let mut rdr = get_csv_reader(get_input(Some(ht_path.clone()))?);

    let mut variants: Vec<CompareHaplotype> = vec![];

    for line in rdr.deserialize() {
        let record: IndexMap<String, String> = line?;
        let variant: CompareHaplotype = record.into();

        if variant.gts.is_empty() {
            tracing::info!(
                "No genotypes at pos: {}, is the column correctly named starting with ht or gt",
                variant.pos
            );
        }

        if let Some(latest) = variants.last() {
            if latest.pos > variant.pos {
                return Err(eyre!(
                    "The haplotype file is not sorted by position. {} is larger than {}",
                    latest.pos,
                    variant.pos
                ));
            }
        }
        variants.push(variant);
    }

    Ok((variants, ht_names))
}

pub fn run(
    haplotypes: Vec<PathBuf>,
    mut output: PathBuf,
    prefix: Option<String>,
    csv: bool,
    hide_missing: bool,
    tag_rows: bool,
) -> Result<()> {
    tracing::debug!("Reading in haplotypes: {:?}", haplotypes);
    match strip_prefix(prefix) {
        None => output.push("ht_comparison.csv"),
        Some(prefix) => output.push(format!("{prefix}_ht_comparison.csv")),
    }

    let mut names = vec![];
    let mut hts = vec![];
    let mut ht_names = vec![];
    for path in haplotypes {
        names.push(path.display().to_string());
        let (ht, ht_name_vec) = read_haplotype_file(path)?;
        ht_names.push(ht_name_vec);

        hts.push(ht);
    }

    tracing::debug!("Finished reading haplotype files");

    let aligner = HaplotypeAligner::new(hts, names, ht_names);

    tracing::debug!("Finished aligning haplotypes");

    match csv {
        true => aligner.to_csv(output)?,
        false => aligner.print(hide_missing, tag_rows),
    }
    Ok(())
}

pub struct HaplotypeAligner {
    filenames: Vec<String>,
    ht_names: Vec<Vec<String>>,
    alignment: BTreeMap<Coord, (Vec<Option<String>>, String)>,
}

impl HaplotypeAligner {
    pub fn new(
        haplotypes: Vec<Vec<CompareHaplotype>>,
        filenames: Vec<String>,
        ht_names: Vec<Vec<String>>,
    ) -> Self {
        let n_gts_count: usize = ht_names.iter().map(|v| v.len()).sum();
        let mut last_contig = &haplotypes[0][0].contig;

        // Register all positions first and initialize a vec of None's corresponding to the number of haplotypes to study
        let mut tree_map: BTreeMap<Coord, (Vec<Option<String>>, String)> = BTreeMap::new();

        for ht in &haplotypes {
            let curr_contig = &ht[0].contig;
            if last_contig != curr_contig {
                panic!("Haplotypes are from different contigs {last_contig} vs {curr_contig}")
            }
            last_contig = curr_contig;

            for variant in ht {
                match tree_map.entry(variant.into()) {
                    Entry::Occupied(entry) => {
                        if entry.key().alt == "-" {
                            entry.remove_entry();
                            tree_map
                                .insert(variant.into(), (vec![None; n_gts_count], String::new()));
                        }
                    }
                    Entry::Vacant(entry) => {
                        entry.insert((vec![None; n_gts_count], String::new()));
                    }
                };
            }
        }

        // Register present or missing genotypes for each position into the vec for each position
        for (coord, (genotypes, annotations)) in tree_map.iter_mut() {
            let mut cum_sum = 0;
            for ht in haplotypes.iter() {
                let mut found = false;
                for hap_variant in ht {
                    if coord == hap_variant {
                        found = true;
                        // Change a None into Some if the haplotype has a genotype at this position
                        for gt in &hap_variant.gts {
                            genotypes[cum_sum] = Some(gt.clone());
                            cum_sum += 1;
                        }

                        if let Some(annotation) = &hap_variant.ann {
                            annotations.push(' ');
                            annotations.push_str(annotation);
                        }
                        break;
                    }
                }
                if !found {
                    cum_sum += ht[0].gts.len();
                }
            }
        }

        Self {
            ht_names,
            filenames,
            alignment: tree_map,
        }
    }

    pub fn print(&self, hide_missing: bool, tag_rows: bool) {
        // Prepare header
        let mut header_names = vec![];
        for (i, name) in self.filenames.iter().enumerate() {
            for ht_name in &self.ht_names[i] {
                header_names.push(format!("{name}/{ht_name}"))
            }
        }
        let header_names = header_names.join(" ");

        println!(" pos {header_names} annotation");

        for (coord, (genotypes, ann)) in &self.alignment {
            // If missing genotypes present, yellow is true
            // If contradictory genotypes present, red is also true

            let (yellow, red) = self.find_mismatch(genotypes);
            let color = if red {
                color::Fg(color::AnsiValue::rgb(5, 0, 1))
            } else if yellow {
                // color::Fg(color::AnsiValue::rgb(5, 5, 1))
                color::Fg(color::AnsiValue::rgb(0, 0, 0))
            } else {
                color::Fg(color::AnsiValue::rgb(0, 5, 1))
            };

            // Switch Nones to - and fold the vec into a string
            let gts = genotypes
                .iter()
                .map(|gt| match gt {
                    None => "-".to_string(),
                    Some(gt) => gt.to_owned(),
                })
                .fold(String::new(), |s, r| format!("{s} {r}"));

            let line = format!("{color} {coord} {gts} {ann}");

            if !hide_missing || !yellow {
                if tag_rows {
                    match (yellow, red) {
                        (_, true) => println!("{line} err"),
                        (true, false) => println!("{line} mis"),
                        (false, false) => println!("{line} ok"),
                    }
                } else {
                    println!("{line}");
                }
            }
        }
    }

    fn find_mismatch(&self, genotypes: &[Option<String>]) -> (bool, bool) {
        let mut is_some_count = 0;
        let hash_set = genotypes
            .iter()
            .flatten()
            .inspect(|_| is_some_count += 1)
            .collect::<HashSet<_>>();

        let has_missing = is_some_count != genotypes.len();
        let has_mismatch = hash_set.len() > 1;

        (has_missing, has_mismatch)
    }

    pub fn to_csv(&self, path: PathBuf) -> Result<()> {
        let mut wrtr = open_csv_writer(path)?;
        let mut header = vec!["pos".to_string()];
        header.extend(self.filenames.clone());

        wrtr.write_record(&header)?;
        for (pos, (genotypes, _ann)) in &self.alignment {
            let mut record = vec![format!("{pos}")];
            genotypes.iter().for_each(|gt| match gt {
                None => record.push("-".to_string()),
                Some(gt) => record.push(gt.to_owned()),
            });
            wrtr.write_record(&record)?;
        }
        Ok(())
    }
}
