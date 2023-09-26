#!/bin/bash

file=$1
controls=$2
outdir=$3
prefix=$4

if [[ -z $file || -z $outdir || -z $controls || -z $prefix ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

for i in {1..22} X;
# for i in 22;
do
  CHR=chr${i}

  cargo run --release --bin haptk -- \
    bhst-gwas \
    $file \
    --trees $outdir/${prefix}_${CHR}_trees.json.gz \
    --controls $controls \
    --prefix ${prefix}_$CHR \
    -o $outdir \
    -t 64 \
    -vvv &

done

wait

for i in {1..22} X;
# for i in 22;
do
  CHR=chr${i}

  python python-scripts/plot_gwas.py \
    $outdir/${prefix}_${CHR}_hst_gwas.csv \
    --output $outdir/${prefix}_${CHR}_hst_gwas.png &
done

wait
mkdir -p $outdir/pdfs
convert $outdir/${prefix}_chr*_hst_gwas.png $outdir/pdfs/${prefix}_hst_gwas.pdf
