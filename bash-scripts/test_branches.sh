#!/bin/bash

branch1=$1
branch2=$2

if [[ -z $branch1 || -z $branch2 ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

git checkout $branch1
cargo run --release --bin haptk bhst haptk/tests/data/1k_genomes_test.vcf.gz -c chr9:27573534 -s all -p tmp_$branch1

hyperfine -w 2 "cargo run --release --bin haptk bhst haptk/tests/data/1k_genomes_test.vcf.gz -c chr9:27573534 -s all -p tmp_${branch1}"

git checkout $branch2
cargo run --release --bin haptk bhst haptk/tests/data/1k_genomes_test.vcf.gz -c chr9:27573534 -s all -p tmp_$branch2

hyperfine -w 2 "cargo run --release --bin haptk bhst haptk/tests/data/1k_genomes_test.vcf.gz -c chr9:27573534 -s all -p tmp_${branch2}"


python ./python-scripts/assert_equal_hst.py tmp_${branch1}_bhst.hst.gz tmp_${branch2}_bhst.hst.gz
