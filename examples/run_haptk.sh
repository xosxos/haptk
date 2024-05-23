#!/bin/bash

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

file=1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
biallelic=1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.biallelic.vcf.gz
ids=${file%%.*}.ids

# Extract only biallelic SNPs from the panel
bcftools view -m2 -M2 -v snps $file -Ou | 
  bcftools annotate -x "INFO" -Oz -o $biallelic

# Index the file
tabix $biallelic

# Save a known over 20 GGGGCC repeat carrying sample ID into a file
haptk samples $biallelic > $ids

# Save a known over 20 GGGGCC repeat carrying sample ID into a file
echo HG01109 > over_20_repeats.ids

# At the moment reading the 1k genomes reference panel into a genotype matrix takes around 2 minutes on a 8-core computer

# Create unidirectional HSTs of all the panel samples and tag the over 20 repeat carrier to find interesting branches
haptk uhst $biallelic \
  --select all \
  --coords chr9:27573534 \
  --samples $ids \
  --mark-samples over_20_repeats.ids \
  --min-size 50 \
  -t 8 \
  -vvv \
  -o results \
  --prefix 20_repeat_tagged

# Check which samples share the Finnish C9orf72 HRE core haplotype
haptk check-for-haplotype $biallelic \
  --haplotype data/finnish_als_core_ht.csv \
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
RUST_BACKTRACE=full haptk compare-to-haplotype $biallelic \
  --select only-longest \
  --coords chr9:27573534 \
  --samples results/finnish_als_core_ht_carriers.ids \
  --mark-samples over_20_repeats.ids \
  --haplotype data/finnish_als_core_ht.csv \
  -vvv \
  -o results \
  --prefix finnish_als

# Extract the samples sharing over 20 kb segments of the ancestral haplotype
cat results/finnish_als_ht_shared_segments_only_longest.csv \
  | awk 'BEGIN { FS = "," } ;{if ($4 > 20000) {print $1}}' \
  | grep -v id \
  > results/over_20kb_sharing_samples.ids

# Tag the over 20 kb Finnish ALS ancestral haplotype sharing samples in the bHST of the 1k genomes C9orf72 HRE core haplotype carriers
haptk bhst $biallelic \
  --select only-longest \
  --coords chr9:27573534 \
  --samples results/finnish_als_core_ht_carriers.ids \
  --mark-samples results/over_20kb_sharing_samples.ids \
  -vvv \
  -o results \
  --prefix 20kb_tagged

# Tag the over 20 repeat carrier in the bHST of the 1k genomes C9orf72 HRE core haplotype carriers
# View it interactively using HAPTK Tree Viewer
# Clicking nodes saves the haplotype and a list of the samples in the node into a `results_tree_haplotypes` directory
haptk-tree 1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.biallelic.vcf.gz --select all --coords chr9:27573534 --samples finnish_als_core_ht_carriers.ids --mark-samples over_20_repeats.ids
