#!/bin/bash

file=$1
outdir=$2
stepsize=$3
samples=$4
prefix=$5


if [[ -z $file || -z $outdir || -z $stepsize || -z $samples || -z $prefix ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

for i in {1..22} X;
do
  CHR=chr${i}
  cargo run --release --bin haptk -- \
    bhst-scan \
    $file \
    --coords $CHR \
    --select all \
    --step-size $stepsize \
    --samples $samples \
    --prefix ${prefix}_${CHR} \
    -o $outdir \
    -t 64 \
    -vvv &
done
