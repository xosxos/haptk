#!/bin/bash

### Remember to activate your mamba environment before running

coords=chr9:27573534
outdir=results_misc
helsinki=test_data/ALS_2020.all.phased.vcf.gz
tampere=test_data/DEB_CDP_HBS_ALS.all.phased.vcf.gz

# Pub HSTs comparison
 haptk compare-to-hst $helsinki \
   -S test_data/c9pos.no_dupl_to_tampere.ids \
   -t 63 \
   -s only-longest \
   -c $coords \
   -o $outdir \
   -p uhst_right \
   --hst results_tampere_c9orf72/pub_uhst_right_only_longest.hst.gz

 python ./python-scripts/article_compare_hsts.py \
   --hst results_tampere_c9orf72/pub_uhst_right_only_longest.hst.gz \
   --match-hst $outdir/uhst_right_match_hst_only_longest.hst.gz \
   --min-size 1 \
   --output $outdir/uhst_right_comparison_only_longest.png


 haptk compare-to-hst $helsinki \
   -S test_data/c9pos.no_dupl_to_tampere.ids \
   -t 63 \
   -s only-longest \
   -c $coords \
   -o $outdir \
   -p uhst_left \
   --hst results_tampere_c9orf72/pub_uhst_left_only_longest.hst.gz

 python ./python-scripts/article_compare_hsts.py \
   --hst results_tampere_c9orf72/pub_uhst_left_only_longest.hst.gz \
   --match-hst $outdir/uhst_left_match_hst_only_longest.hst.gz \
   --min-size 1 \
   --output $outdir/uhst_left_comparison_only_longest.png


 haptk compare-to-hst $helsinki \
   -S test_data/c9pos.no_dupl_to_tampere.ids \
   -t 63 \
   -s only-longest \
   -c $coords \
   -o $outdir \
   -p bhst \
   --hst results_tampere_c9orf72/pub_bhst_only_longest.hst.gz

 python ./python-scripts/article_compare_hsts.py \
   --hst results_tampere_c9orf72/pub_bhst_only_longest.hst.gz \
   --match-hst $outdir/bhst_match_hst_only_longest.hst.gz \
   --min-size 1 \
   --output $outdir/bhst_comparison_only_longest.png


haptk compare-to-haplotype $helsinki \
  -S test_data/c9pos.no_dupl_to_tampere.ids \
  -t 63 \
  -s only-longest \
  -c $coords \
  -o $outdir \
  --haplotype results_tampere_c9orf72/pub_uhst_mbah_only_longest.csv


# Mix helsinki HSTs with 20 IAs
haptk bhst $helsinki \
  -S ./test_data/affy_20_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "20_HEL" \
  --mark-samples ./test_data/affy_20_50.ids \
  -c $coords

python ./python-scripts/article_tagged_tree.py $outdir/20_HEL_bhst_only_longest.hst.gz \
  --min-size 1 \
  --ids ./test_data/affy_20_50.ids \
  --output $outdir/20_HEL_bhst_only_longest.png


# Mix helsinki HSTs with 15-19 IAs
haptk bhst $helsinki \
  -S ./test_data/affy_15_19_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "15_HEL" \
  --mark-samples ./test_data/affy_15_19.ids \
  -c $coords

python ./python-scripts/article_tagged_tree.py $outdir/15_HEL_bhst_only_longest.hst.gz \
  --min-size 1 \
  --ids ./test_data/affy_15_19.ids \
  --output $outdir/15_HEL_bhst_only_longest.png

# Mix Tampere HSTs with 20 IAs
haptk bhst $tampere \
  -S ./test_data/V1_20_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "20_TAM" \
  --mark-samples ./test_data/V1_20_lengths.ids \
  -c $coords

python ./python-scripts/article_tagged_tree.py $outdir/20_TAM_bhst_only_longest.hst.gz \
  --min-size 1 \
  --ids ./test_data/V1_20_lengths.ids \
  --output $outdir/20_TAM_bhst_only_longest.png

# Mix Tampere HSTs with 15-19 IAs
haptk bhst $tampere \
  -S ./test_data/V1_15_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "15_TAM" \
  --mark-samples ./test_data/V1_15_19_lengths.ids \
  -c $coords

python ./python-scripts/article_tagged_tree.py $outdir/15_TAM_bhst_only_longest.hst.gz \
  --min-size 1 \
  --ids ./test_data/V1_15_19_lengths.ids \
  --output $outdir/15_TAM_bhst_only_longest.png



haptk compare-to-haplotype $tampere \
  -S ./test_data/V1_20_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "20_TAM" \
  --mark-samples ./test_data/V1_20_lengths.ids \
  --haplotype ./results_tampere_c9orf72/uhst_mbah_only_longest.csv \
  -c $coords


haptk compare-to-haplotype $tampere \
  -S ./test_data/V1_15_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "15_TAM" \
  --mark-samples ./test_data/V1_15_19_lengths.ids \
  --haplotype ./results_tampere_c9orf72/uhst_mbah_only_longest.csv \
  -c $coords


# Mix helsinki HSTs with 20 IAs
haptk compare-to-haplotype $helsinki \
  -S ./test_data/affy_20_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "20_HEL" \
  --mark-samples ./test_data/affy_20_50.ids \
  --haplotype ./results_helsinki_c9orf72/uhst_mbah_only_longest.csv \
  -c $coords

haptk compare-to-haplotype $helsinki \
  -S ./test_data/affy_15_19_and_c9pos.ids \
  -t 63 \
  -s only-longest \
  -o $outdir \
  -p "15_HEL" \
  --mark-samples ./test_data/affy_15_19.ids \
  --haplotype ./results_helsinki_c9orf72/uhst_mbah_only_longest.csv \
  -c $coords

echo ""

cd $outdir
# Archive and compress all figures
tar -cvf figures.tar.gz *.png *ht_comparison*


