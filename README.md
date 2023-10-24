## Features

```
Commands:
  uhst                  Build unidirectional haplotype sharing trees at a coordinate
  bhst                  Build a bidirectional haplotype sharing tree at a coordinate
  mrca                  Analyze the MRCA based on the Gamma method at a coordinate
  check-for-haplotype   Check if samples share a given haplotype
  compare-to-haplotype  Check differences between samples and a haplotype
  compare-to-hst        Check what haplotypes of the HST are present in samples
  compare-haplotypes    Compare haplotypes to each other by alignment
  coverage              Show coverage levels per contig in a VCF
  haplotypes            Read the haplotypes of a given sample
  samples               Output the sample names from FAM / VCF files
  to-vcf                Convert a haplotype CSV into VCF
  help                  Print this message or the help of the given subcommand(s)
```

## How to install
```bash
# HAPTK has been tested with Rust 1.69 and higher

# Install Rust and the Cargo package manager
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# The Rust binaries directory is added automatically to the shells $PATH variable when it is restarted
# If you face problems, add the directory to $PATH.
# For example for `bash`, add this line to your configuration file:
export PATH="$HOME/.cargo/bin:$PATH"

# Install HAPTK
cargo install haptk --locked
```

## Installing from source
```bash
# Clone the repository
git clone https://github.com/xosxos/haptk

# Change directory 
cd haptk

# Build. The binary will be available in `./target/release/` 
cargo build --release -p haptk

```

## Running commands
To get started, try the `run_haptk.sh` script, which uses the 1000 genomes reference panel and Finnish ALS haplotypes. The script can be found in the `examples` directory.

Often a common requirement is to select only certain alleles of the samples based on a condition, this is done with the `--select` argument.
Options:
- `all` (default) [Select all alleles]
- `only-refs` [Select only the samples carrying the ref variant at the coordinate]
- `only-alts` [Select only the samples carrying the alt variant at the coordinate]
- `only-longest` [Select only the longest haplotype sharing alleles per sample]

```bash

# Download the 1000 genomes reference panel for chromosome 9
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

# Define variables
file=1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
biallelic=1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.biallelic.vcf.gz
ids=${file%%.*}.ids

# Extract only biallelic SNPs from the panel
bcftools view -m2 -M2 -v snps $file -Ou | 
  bcftools annotate -x "INFO" -Oz -o $biallelic

# Index the file
tabix $biallelic

# Save all the ids to a file
haptk samples $biallelic > $ids

# Save a known over 20 GGGGCC repeat carrying sample ID into a file
echo HG01109 > over_20_repeats.ids

# The major bottleneck for HST creation is reading in the genotype data
# Reading the 1k genomes reference panel matrix should take around 2 minutes using 8-cores
# For now to speed up HST creation, you can try to filter out variants
# from the reference panel before running the analysis

# Create unidirectional HSTs of the samples and tag over 20 repeat carriers to find interesting branches
haptk uhst $biallelic \
  --select all \
  --coords chr9:27573534 \
  --samples $ids \
  --mark-samples over_20_repeats.ids \
  --min-node-size 50 \
  -t 8 \
  -vvv \
  -o results \
  --prefix 20_repeat_tagged

# The example files are in Github at `examples/data/`
# Check which samples share the Finnish C9orf72 HRE core haplotype
haptk check-for-haplotype $biallelic \
  --haplotype examples/data/finnish_als_core_ht.csv \
  --coords chr9:27573534 \
  --select all \
  -vvv \
  -o results \
  --prefix finnish_als

# Extract these samples into a list
cat results/finnish_als_haplotype_check.csv \
  | grep true \
  | sed 's/,true//g' \
  > results/finnish_als_core_ht_carriers.ids

# Compare these samples to a Finnish C9orf72 HRE ancestral haplotype
haptk compare-to-haplotype $biallelic \
  --select only-longest \
  --coords chr9:27573534 \
  --samples results/finnish_als_core_ht_carriers.ids \
  --mark-samples over_20_repeats.ids \
  --haplotype examples/data/finnish_als_uhst_mbah.csv \
  -vvv \
  -o results \
  --prefix finnish_als

# Extract the samples sharing over 20 kb segments of the ancestral haplotype
cat results/finnish_als_ht_shared_segments_only_longest.csv \
  | awk 'BEGIN { FS = "," } ;{if ($4 > 20000) {print $1}}' \
  | grep -v id \
  > results/over_20kb_sharing_samples.ids

# Tag the over 20 kb Finnish ALS ancestral haplotype sharing samples
# in the bidirectional HST of the C9orf72 HRE core haplotype carriers
haptk bhst $biallelic \
  --select only-longest \
  --coords chr9:27573534 \
  --samples results/finnish_als_core_ht_carriers.ids \
  --mark-samples results/over_20kb_sharing_samples.ids \
  -vvv \
  -o results \
  --prefix 20kb_tagged


# To visualize graphs for publication, using ETE3 is recommended
# Install ETE from http://etetoolkit.org/download/ and the required python dependencies

# Try adding `conda-forge` and `bioconda` channels to `mamba` and running:
# mamba install ete3 rustworkx

# Example scripts for ETE3 visualization
python python-scripts/article_tagged_tree.py results/20kb_tagged_bhst_only_longest.hst.gz \
  --min-size 1 \
  --ids results/over_20kb_sharing_samples.ids \
  --output 20kb_tagged_bhst_only_longest.png

python python-scripts/article_normal_tree.py results/20_repeat_tagged_uhst_left.hst.gz \
  --min-size 10 \
  --output 20_repeat_tagged_uhst_left.png

# You can also view the trees interactively using the HAPTK Tree Viewer
# Clicking nodes saves the haplotype and a list of carrier IDs into `results_tree_haplotypes` directory
# NOTE: Requires building `haptk-tree` binary from source using Rust nightly 
# To install Rust nightly: rustup default nightly
# After installing nightly: cargo build --release -p haptk-tree
haptk-tree 1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.biallelic.vcf.gz \
  --select all \
  --coords chr9:27573534 \
  --samples finnish_als_core_ht_carriers.ids \
  --mark-samples over_20_repeats.ids

```
## Example graphs from the original [article](https://www.biorxiv.org/content/10.1101/2023.07.28.550820v3)

### Unidirectional haplotype sharing trees (uHST)
The left side uHST starting from the C9orf72 expansion.

![ALS_203_left_only_longest](https://github.com/xosxos/haptk/assets/44613540/6be3bcfc-e7f7-432b-926b-006d07aa2498)

The right side uHST starting from the C9orf72 expansion.

![ALS_203_right_only_longest](https://github.com/xosxos/haptk/assets/44613540/b963bac7-1407-40ae-bb60-3afa7ec6f7f1)

### Bidirectional haplotype sharing tree (bHST) with tagged samples
A Finnish C9orf72 expansion carrier bHST with FTD cases tagged in green and ALS cases tagged in red.

![ALS_203_ftd_tagged](https://github.com/xosxos/haptk/assets/44613540/f1193d27-cd78-43b1-b217-f8150dd6c6dd)

### Haplotype comparison graph
All the samples are compared to a haplotype, in this case the C9orf72 expansion ancestral haplotype.
The matching variants are colored in white (or blue) and deviating variants in pink.
Because the actual sample haplotypes are seen in relation to the ancestral haplotype, potential switch and genotyping errors can be evaluated from this graph.
The selection of C9orf72 haplotypes is done by discarding the haplotypes tagged in blue (the ones sharing less ancestral haplotype sequence per sample).

![diff_only_longest_marked](https://github.com/xosxos/haptk/assets/44613540/f9a9b3b5-24b5-467f-9785-76ceba00754c)

## Citations
_The shared ancestry between the C9orf72 hexanucleotide repeat expansion and intermediate-length alleles using haplotype sharing trees and HAPTK_

[https://www.biorxiv.org/content/10.1101/2023.07.28.550820v3](https://www.biorxiv.org/content/10.1101/2023.07.28.550820v3)
