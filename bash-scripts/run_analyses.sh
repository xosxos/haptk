#!/bin/bash

### Remember to activate your mamba environment before running

file=$1
cohort=$2
samples=$3
coords=$4
recombination_rates=$5
outdir=$6
gene=$7

if [[ -z $file || -z $cohort || -z $samples || -z $coords || -z $recombination_rates || -z $outdir || -z $gene ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

ancestral_haplotype=./$outdir/uhst_mbah_only_longest.csv
ancestral_haplotype_bhst=./$outdir/bhst_mbah_only_longest.csv
core_haplotype=./$outdir/uhst_shared_core_haplotype_only_longest.csv

# Construct the unidirectional HSTs
haptk uhst $file -c $coords -S $samples --select only-longest -o $outdir

# Draw the HST of the left side
python ./python-scripts/article_normal_tree.py $outdir/uhst_left_only_longest.hst.gz \
  --min-size 10 \
  --output $outdir/uhst_left_only_longest.png

# Draw the HST of the right side
python ./python-scripts/article_normal_tree.py $outdir/uhst_right_only_longest.hst.gz \
  --min-size 10 \
  --output $outdir/uhst_right_only_longest.png

# Construct the bidirectional HST
haptk bhst $file -c $coords -S $samples -s only-longest -o $outdir

# Draw the bidirectional HST
python ./python-scripts/article_normal_tree.py $outdir/bhst_only_longest.hst.gz \
  --min-size 10 \
  --output $outdir/bhst_only_longest.png

# Compare all the sample haplotypes to the ancestral uHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s all \
  --haplotype $ancestral_haplotype \
  --mark-shorter-alleles \
  --prefix "uhst" \
  -o $outdir

# Compare only-longest sample haplotypes to the ancestral uHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s only-longest \
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

# Compare only-longest sample haplotypes to the ancestral bHST haplotype
haptk compare-to-haplotype $file \
  -c $coords \
  -S $samples \
  -s only-longest \
  --haplotype $ancestral_haplotype_bhst \
  --prefix "bhst" \
  -o $outdir

# Plot uHST ancestral segments with recombination rates
python ./python-scripts/article_plot_ancestral_segments.py $outdir/uhst_ht_shared_segments_only_longest.csv \
  -r $recombination_rates \
  -c $coords \
  --gene $gene \
  -o $outdir/uhst_mbah_segments_only_longest.png

# Plot bHST ancestral segments with recombination rates
python ./python-scripts/article_plot_ancestral_segments.py $outdir/bhst_ht_shared_segments_only_longest.csv \
  -r $recombination_rates \
  -c $coords \
  --gene $gene \
  -o $outdir/bhst_mbah_segments_only_longest.png


# Construct the unidirectional HSTs
haptk mrca $file -c $coords -S $samples -s only-longest -r $recombination_rates -o $outdir

# Construct the publishable HSTs
# This cuts out all nodes with less than 10 samples and removes all IDs
haptk uhst $file -c $coords -S $samples --select only-longest -o $outdir --min-size 10 --publish -p "pub"

python ./python-scripts/article_normal_tree.py $outdir/pub_uhst_right_only_longest.hst.gz \
  --min-size 1 \
  --output $outdir/pub_uhst_right_only_longest.png

python ./python-scripts/article_normal_tree.py $outdir/pub_uhst_left_only_longest.hst.gz \
  --min-size 1 \
  --output $outdir/pub_uhst_left_only_longest.png

haptk bhst $file -c $coords -S $samples -s only-longest -o $outdir --min-size 10 --publish -p "pub"

python ./python-scripts/article_normal_tree.py $outdir/pub_bhst_only_longest.hst.gz \
  --min-size 1 \
  --output $outdir/pub_bhst_only_longest.png

if [ $cohort = "tampere" ]; then
  # Draw the HST of the right side with FTD cases tagged
  python ./python-scripts/article_tagged_tree.py $outdir/uhst_right_only_longest.hst.gz \
    --min-size 1 \
    --ids test_data/ftd.ids \
    --output $outdir/uhst_right_only_longest_ftd.png

  # Draw the HST of the left side with FTD cases tagged
  python ./python-scripts/article_tagged_tree.py $outdir/uhst_left_only_longest.hst.gz \
    --min-size 1 \
    --ids test_data/ftd.ids \
    --output $outdir/uhst_left_only_longest_ftd.png

  # Draw the circular bHST with FTD cases tagged
  python ./python-scripts/article_tagged_tree.py $outdir/bhst_only_longest.hst.gz \
    --min-size 1 \
    --ids test_data/ftd.ids \
    --output $outdir/bhst_only_longest_ftd.png

  # Draw a normal bHST with FTD cases tagged
  python ./python-scripts/article_normal_tree.py $outdir/bhst_only_longest.hst.gz \
    --min-size 1 \
    --ids test_data/ftd.ids \
    --output $outdir/bhst_only_longest_ftd.png
fi

# Print majority based ancestral haplotype sharing summaries
cd $outdir
printf "\nbHST all shared segments\n"
cat bhst_ht_shared_segments.csv | qsv stats -s markers,length --median | qsv table

printf "\nbHST only-longest shared segments\n"
cat bhst_ht_shared_segments_only_longest.csv | qsv stats -s markers,length --median | qsv table

printf "\nuHST all shared segments\n"
cat uhst_ht_shared_segments.csv | qsv stats -s markers,length --median | qsv table

printf "\nuHST only-longest shared segments\n"
cat uhst_ht_shared_segments_only_longest.csv | qsv stats -s markers,length --median | qsv table

echo ""

# Archive and compress all figures
tar -cvf figures.tar.gz pub_* *.png *ht_comparison*

