#!/bin/bash

file=$1
outdir=$2
samples=$3

if [[ -z $file || -z $outdir || -z $samples ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

for i in {1..22} X;
do
  CHR=chr${i}
  rec_rates=test_data/recombination_rates_${CHR}_GRCh38.tsv

  cargo run --release --bin haptk -- \
    mrca-scan \
    $file \
    --coords $CHR \
    --select only-longest \
    --step-size 1 \
    --samples $samples \
    --recombination-rates $rec_rates \
    --prefix $CHR \
    -o $outdir \
    -t 64 \
    -vvv &

  python python-scripts/plot_mrca_scan.py \
    $outdir/${CHR}_mrca_scan_only_longest.csv \
    --output $outdir/${CHR}_mrca_scan_only_longest.png &

done

wait
mkdir -p $outdir/pdfs
convert $outdir/chr*_mrca_scan_only_longest.png $outdir/pdfs/mrca_scans_only_longest.pdf
