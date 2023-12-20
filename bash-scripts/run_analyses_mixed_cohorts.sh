#!/bin/bash

### Remember to activate your mamba environment before running

file=$1
samples_no_tag=$2
samples_yes_tag=$3
coords=$4
outdir=$5

if [[ -z $file || -z $samples_no_tag || -z $samples_yes_tag || -z $coords || -z $outdir ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

samples=./temp.ids

cat $samples_no_tag $samples_yes_tag > $samples

# Construct the unidirectional HSTs
haptk uhst $file -c $coords -S $samples -s all -o $outdir

# Draw the HST of the left side
python ./python-scripts/article_normal_tree.py $outdir/uhst_left.hst.gz \
  --min-size 1 \
  --ids $samples_yes_tag \
  --output $outdir/uhst_left.png

# Draw the HST of the right side
python ./python-scripts/article_normal_tree.py $outdir/uhst_right.hst.gz \
  --min-size 1 \
  --ids $samples_yes_tag \
  --output $outdir/uhst_right.png

# Construct the bidirectional HST
haptk bhst $file -c $coords -S $samples -s all -o $outdir

# Draw the bidirectional HST
python ./python-scripts/article_normal_tree.py $outdir/bhst.hst.gz \
  --min-size 1 \
  --ids $samples_yes_tag \
  --output $outdir/bhst.png

echo ""

rm $samples

# Archive and compress all figures
# tar -cvf figures.tar.gz pub_* *.png *ht_comparison*

