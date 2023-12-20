#!/bin/bash

### Remember to activate your mamba environment before running

file=$1
samples=$2
coords=$3
recombination_rates=$4
outdir=$5
gene=$6

if [[ -z $file || -z $samples || -z $coords || -z $recombination_rates || -z $outdir || -z $gene ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

ancestral_haplotype=./$outdir/uhst_mbah.csv
ancestral_haplotype_bhst=./$outdir/bhst_mbah.csv
core_haplotype=./$outdir/uhst_shared_core_haplotype.csv

# Construct the unidirectional HSTs
haptk uhst $file -c $coords -S $samples -s all -o $outdir

# Draw the HST of the left side
python ./python-scripts/article_normal_tree.py $outdir/uhst_left.hst.gz \
  --min-size 10 \
  --output $outdir/uhst_left.png

# Draw the HST of the right side
python ./python-scripts/article_normal_tree.py $outdir/uhst_right.hst.gz \
  --min-size 10 \
  --output $outdir/uhst_right.png

# Construct the bidirectional HST
haptk bhst $file -c $coords -S $samples -s all -o $outdir

# Draw the bidirectional HST
python ./python-scripts/article_normal_tree.py $outdir/bhst.hst.gz \
  --min-size 10 \
  --output $outdir/bhst.png

# Compare all the haplotypes to the ancestral uHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s all \
  --haplotype $ancestral_haplotype \
  --mark-shorter-alleles \
  --prefix "uhst" \
  -o $outdir

# Compare all the haplotypes to the ancestral uHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s all \
  --haplotype $ancestral_haplotype \
  --prefix "uhst" \
  -o $outdir

# Compare all the sample haplotypes to the ancestral bHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s all \
  --haplotype $ancestral_haplotype_bhst \
  --mark-shorter-alleles \
  --prefix "bhst" \
  -o $outdir

# Compare all the sample haplotypes to the ancestral bHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s all \
  --haplotype $ancestral_haplotype_bhst \
  --prefix "bhst" \
  -o $outdir


# Plot uHST ancestral segments with recombination rates
python ./python-scripts/article_plot_ancestral_segments.py $outdir/uhst_ht_shared_segments.csv \
  -r $recombination_rates \
  -c $coords \
  --gene $gene \
  -o $outdir/uhst_mbah_segments.png

# Plot bHST ancestral segments with recombination rates
python ./python-scripts/article_plot_ancestral_segments.py $outdir/bhst_ht_shared_segments.csv \
  -r $recombination_rates \
  -c $coords \
  --gene $gene \
  -o $outdir/bhst_mbah_segments.png


# Calculate MRCA
haptk mrca $file -c $coords -S $samples -s all -r $recombination_rates -o $outdir

# Construct the publishable HSTs
# This cuts out all nodes with less than 10 samples and removes all IDs
haptk uhst $file -c $coords -S $samples --select all -o $outdir --min-size 10 --publish -p "pub"

python ./python-scripts/article_normal_tree.py $outdir/pub_uhst_right.hst.gz \
  --min-size 1 \
  --output $outdir/pub_uhst_right.png

python ./python-scripts/article_normal_tree.py $outdir/pub_uhst_left.hst.gz \
  --min-size 1 \
  --output $outdir/pub_uhst_left.png

haptk bhst $file -c $coords -S $samples -s all -o $outdir --min-size 10 --publish -p "pub"

python ./python-scripts/article_normal_tree.py $outdir/pub_bhst.hst.gz \
  --min-size 1 \
  --output $outdir/pub_bhst.png

# Print majority based ancestral haplotype sharing summaries
cd $outdir
printf "\nbHST all shared segments\n"
cat bhst_ht_shared_segments.csv | zsv stats -s markers,length --median | zsv table

printf "\nuHST all shared segments\n"
cat uhst_ht_shared_segments.csv | zsv stats -s markers,length --median | zsv table

echo ""

# Archive and compress all figures
tar -cvf figures.tar.gz pub_* *.png *ht_comparison*

