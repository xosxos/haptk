#!/bin/bash

file=$1
min_branch_size=$2
outdir=$3
prefix=$4

if [[ -z $file || -z $min_branch_size || -z $outdir || -z $prefix ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

for i in {1..22} X;
do
  CHR=chr${i}
  rec_rates=test_data/recombination_rates_${CHR}_GRCh38.tsv

  cargo run --release --bin haptk -- \
    bhst-mrca-gwas \
    $file \
    --trees $outdir/${prefix}_${CHR}_trees.json.gz \
    --recombination-rates $rec_rates \
    --prefix ${prefix}_${CHR} \
    --min-sample-size $min_branch_size \
    -o $outdir \
    -t 64 \
    -vvv &

done

wait

for i in {1..22} X;
do
  CHR=chr${i}

  python python-scripts/plot_mrca_scan.py \
    $outdir/${prefix}_${CHR}_hst_mrca_scan.csv \
    --output $outdir/${prefix}_${CHR}_hst_mrca_scan.png \
    --ymax 3 &

done

wait
mkdir -p $outdir/pdfs
convert $outdir/${prefix}_chr*_hst_mrca_scan.png $outdir/pdfs/${prefix}_hst_mrca_scans.pdf
