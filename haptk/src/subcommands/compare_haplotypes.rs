use std::collections::btree_map::Entry;

use std::collections::{BTreeMap, HashSet};
use std::io::{self};
use std::path::PathBuf;

use color_eyre::{
    eyre::{ensure, eyre},
    Result,
};
use csv::{Reader, ReaderBuilder};
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use termion::color;

use crate::io::{get_csv_writer, get_input};
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
    pub gts: Vec<u8>,
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
                .map(|k| hm[k].parse::<u8>().expect("cannot parse genotype to u8"))
                .collect(),
        }
    }
}

#[allow(clippy::type_complexity)]
pub fn read_haplotype_file(
    ht_path: PathBuf,
) -> Result<(Vec<CompareHaplotype>, Vec<String>, Vec<String>, Vec<String>)> {
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
    let mut freq: Vec<String> = vec![];
    let mut n: Vec<String> = vec![];

    for line in rdr.deserialize() {
        let record: IndexMap<String, String> = line?;

        if record["contig"] == "freq" {
            freq.extend(
                record
                    .keys()
                    .filter(|k| k.starts_with("ht") || k.starts_with("gt"))
                    .map(|k| record[k].to_string()),
            );
        } else if record["contig"] == "n" {
            n.extend(
                record
                    .keys()
                    .filter(|k| k.starts_with("ht") || k.starts_with("gt"))
                    .map(|k| record[k].to_string()),
            );
        } else {
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
    }
    if n.is_empty() {
        n.extend(ht_names.iter().map(|_| "NA".to_string()));
    }

    if freq.is_empty() {
        freq.extend(ht_names.iter().map(|_| "NA".to_string()));
    }

    Ok((variants, ht_names, n, freq))
}

pub fn run(
    haplotypes: Vec<PathBuf>,
    mut output: PathBuf,
    prefix: Option<String>,
    csv: bool,
    hide_missing: bool,
    tag_rows: bool,
    nucleotides: bool,
) -> Result<()> {
    tracing::debug!("Reading in haplotypes: {:?}", haplotypes);

    ensure!(!haplotypes.is_empty(), "No haplotypes were given");

    match strip_prefix(prefix) {
        None => output.push("ht_comparison.csv"),
        Some(prefix) => output.push(format!("{prefix}_ht_comparison.csv")),
    }

    let mut names = vec![];
    let mut hts = vec![];
    let mut ht_names = vec![];
    let mut freq_data = vec![];
    for path in haplotypes {
        names.push(path.display().to_string());
        let (ht, ht_name_vec, n, freq) = read_haplotype_file(path)?;
        ht_names.push(ht_name_vec);

        hts.push(ht);
        freq_data.push((n, freq));
    }

    tracing::debug!("Finished reading haplotype files");

    let aligner = HaplotypeAligner::new(hts, names, ht_names, freq_data);

    tracing::debug!("Finished aligning haplotypes");

    match csv {
        true => aligner.to_csv(output, nucleotides)?,
        false => aligner.print(hide_missing, tag_rows, nucleotides),
    }
    Ok(())
}

pub struct HaplotypeAligner {
    filenames: Vec<String>,
    ht_names: Vec<Vec<String>>,
    freq_data: Vec<(Vec<String>, Vec<String>)>,
    alignment: BTreeMap<Coord, (Vec<Option<u8>>, String)>,
}

fn create_btreemap(
    haplotypes: &[Vec<CompareHaplotype>],
    n_gts_count: usize,
) -> BTreeMap<Coord, (Vec<Option<u8>>, String)> {
    let mut tree_map: BTreeMap<Coord, (Vec<Option<u8>>, String)> = BTreeMap::new();
    let mut last_contig = &haplotypes[0][0].contig;
    let mut cum_sum = 0;

    // Loop all haplotypes and variants in haplotypes
    for ht in haplotypes {
        let ngenotypes = ht[0].gts.len();
        let curr_contig = &ht[0].contig;
        if last_contig != curr_contig {
            panic!("Haplotypes are from different contigs {last_contig} vs {curr_contig}")
        }
        last_contig = curr_contig;

        for variant in ht {
            // Insert variant to map
            match tree_map.entry(variant.into()) {
                Entry::Occupied(entry) => {
                    if entry.key().alt == "-" {
                        let (_gts, ann) = entry.get();
                        let new_ann = match &variant.ann {
                            Some(new_ann) => format!("{ann} {new_ann}"),
                            None => ann.to_string(),
                        };
                        entry.remove_entry();

                        let init_gts = initiate_genotypes(cum_sum, n_gts_count, variant);
                        tree_map.insert(variant.into(), (init_gts, new_ann));
                    } else {
                        let (gts, ann) = entry.into_mut();
                        gts.iter_mut()
                            .skip(cum_sum)
                            .enumerate()
                            .take(ngenotypes)
                            .for_each(|(i, v)| *v = Some(variant.gts[i]));

                        if let Some(new_ann) = &variant.ann {
                            *ann = format!("{ann} {new_ann}");
                        }
                    }
                }
                Entry::Vacant(entry) => {
                    let init_gts = initiate_genotypes(cum_sum, n_gts_count, variant);
                    entry.insert((init_gts, variant.ann.clone().unwrap_or(String::from(""))));
                }
            };
        }
        cum_sum += ngenotypes;
    }
    tree_map
}

fn initiate_genotypes(
    cum_sum: usize,
    n_gts_count: usize,
    variant: &CompareHaplotype,
) -> Vec<Option<u8>> {
    let genotypes: Vec<Option<u8>> = (0..cum_sum)
        // Add cum_sum amount of Nones
        .map(|_| None)
        // Add genotypes
        .chain(variant.gts.iter().map(|gt| Some(*gt)))
        // Fill the rest with nones
        .chain((0..(n_gts_count - cum_sum - variant.gts.len())).map(|_| None))
        .collect();

    assert_eq!(genotypes.len(), n_gts_count);
    genotypes
}

impl HaplotypeAligner {
    pub fn new(
        haplotypes: Vec<Vec<CompareHaplotype>>,
        filenames: Vec<String>,
        ht_names: Vec<Vec<String>>,
        freq_data: Vec<(Vec<String>, Vec<String>)>,
    ) -> Self {
        let n_gts_count: usize = ht_names.iter().map(|v| v.len()).sum();
        Self {
            ht_names,
            freq_data,
            filenames,
            alignment: create_btreemap(&haplotypes, n_gts_count),
        }
    }

    fn header_names(&self) -> Vec<String> {
        let mut header_names = vec![];
        for (i, name) in self.filenames.iter().enumerate() {
            for ht_name in &self.ht_names[i] {
                header_names.push(format!("{name}/{ht_name}"))
            }
        }
        header_names
    }

    pub fn print(&self, hide_missing: bool, tag_rows: bool, yes_nucleotides: bool) {
        // Prepare header
        let mut header = self.header(yes_nucleotides);
        header.extend(self.header_names());
        header.push("annotation".to_string());
        let header = header.join(" ");

        // Print out header
        println!("{header}");

        // If missing genotypes present, yellow is true
        // If contradictory genotypes present, red is also true
        for (coord, (genotypes, ann)) in &self.alignment {
            let (is_mismatch, is_missing) = self.find_mismatch(genotypes);

            let color = match (is_mismatch, is_missing) {
                // Red
                (true, _) => color::Fg(color::AnsiValue::rgb(5, 0, 1)),
                // Yellow
                (false, true) => color::Fg(color::AnsiValue::rgb(0, 0, 0)),
                // Green
                (false, false) => color::Fg(color::AnsiValue::rgb(0, 5, 1)),
            };

            // Switch Nones to - and fold the vec into a string
            let gts = genotypes
                .iter()
                .map(|gt| match gt {
                    None => "-".to_string(),
                    Some(gt) => gt.to_string(),
                })
                .fold(String::new(), |s, r| format!("{s} {r}"));

            let line = format!(
                "{color} {} {} {} {} {gts} {ann}",
                coord.contig, coord.pos, coord.reference, coord.alt
            );

            if !hide_missing || !is_missing {
                if tag_rows {
                    match (is_mismatch, is_missing) {
                        (true, _) => println!("{line} err"),
                        (false, true) => println!("{line} mis"),
                        (false, false) => println!("{line} ok"),
                    }
                } else {
                    println!("{line}");
                }
            }
        }

        let (mut freq_line, mut n_line) = self.freq_and_n_line(yes_nucleotides);

        for (ns, freqs) in &self.freq_data {
            for freq in freqs {
                freq_line.push(freq.clone());
            }
            for n in ns {
                n_line.push(n.clone());
            }
        }
        freq_line.push("-".to_string());
        n_line.push("-".to_string());
        let freq_line = freq_line.join(" ");
        let n_line = n_line.join(" ");
        println!("{freq_line}");
        println!("{n_line}");
    }

    fn find_mismatch(&self, genotypes: &[Option<u8>]) -> (bool, bool) {
        let mut is_some_count = 0;
        let hash_set = genotypes
            .iter()
            .flatten()
            .inspect(|_| is_some_count += 1)
            .collect::<HashSet<_>>();

        let has_missing = is_some_count != genotypes.len();
        let has_mismatch = hash_set.len() > 1;

        (has_mismatch, has_missing)
    }

    fn header(&self, yes_nucleotides: bool) -> Vec<String> {
        match yes_nucleotides {
            true => ["contig", "pos", "ref"].map(String::from).to_vec(),
            false => ["contig", "pos", "ref", "alt"].map(String::from).to_vec(),
        }
    }

    fn freq_and_n_line(&self, yes_nucleotides: bool) -> (Vec<String>, Vec<String>) {
        let freq_line = match yes_nucleotides {
            true => ["freq", "", ""].map(String::from).to_vec(),
            false => ["freq", "", "", ""].map(String::from).to_vec(),
        };
        let n_line = match yes_nucleotides {
            true => ["n", "", ""].map(String::from).to_vec(),
            false => ["n", "", "", ""].map(String::from).to_vec(),
        };

        (freq_line, n_line)
    }

    pub fn to_csv(&self, _path: PathBuf, yes_nucleotides: bool) -> Result<()> {
        // let mut wrtr = open_csv_writer(path)?;
        let mut wrtr = get_csv_writer(Box::new(std::io::stdout()));

        let mut header = self.header(yes_nucleotides);
        header.extend(self.header_names());
        header.push("annotation".to_string());

        wrtr.write_record(&header)?;

        let (mut freq_line, mut n_line) = self.freq_and_n_line(yes_nucleotides);

        for (ns, freqs) in &self.freq_data {
            for freq in freqs {
                freq_line.push(freq.clone());
            }
            for n in ns {
                n_line.push(n.clone());
            }
        }
        freq_line.push("".to_string());
        n_line.push("".to_string());

        wrtr.write_record(&freq_line)?;
        wrtr.write_record(&n_line)?;

        for (coord, (genotypes, ann)) in &self.alignment {
            let mut record = match yes_nucleotides {
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

            genotypes.iter().for_each(|gt| match gt {
                None => record.push("-".to_string()),
                Some(gt) => match yes_nucleotides {
                    true => match gt {
                        0 => record.push(coord.reference.clone()),
                        1 => record.push(coord.alt.clone()),
                        _ => panic!("multiallelic haplotypes not yet allowed"),
                    },
                    false => record.push(gt.to_string()),
                },
            });

            record.push(ann.to_string());
            wrtr.write_record(&record)?;
        }
        Ok(())
    }
}
