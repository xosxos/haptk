#!/bin/bash

file=$1
samples=$2
outdir=$3
prefix=$4

if [[ -z $file || -z $outdir || -z $samples || -z $prefix ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

for i in {1..22} X;
do
  CHR=chr${i}

  cargo run --release --bin haptk -- \
    bhst-segregation \
    $file \
    --trees $outdir/${prefix}_${CHR}_trees.json.gz \
    --samples $samples \
    --prefix ${prefix}_$CHR \
    -o $outdir \
    -t 64 \
    -vvv &

done

wait

for i in {1..22} X;
do
  CHR=chr${i}

  python python-scripts/plot_hst_segregation.py \
    $outdir/${prefix}_${CHR}_hst_segregation_scan.csv \
    --output $outdir/${prefix}_${CHR}_hst_segregation.png &
done

wait
mkdir -p $outdir/pdfs
convert $outdir/${prefix}_chr*_hst_segregation.png $outdir/pdfs/${prefix}_hst_segregation.pdf
