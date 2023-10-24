#!/bin/bash

file=$1
lengths=$2
cohort=$3

if [[ -z $lengths || -z $file || -z $cohort ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi

coords=chr9:27573534
# anc_haplotype=./results_${cohort}_c9orf72_no_imp/ancestral_haplotype.csv
# long_haplotype=./results_${cohort}_c9orf72_no_imp/longest_haplotype.csv
recombination_rates=./test_data/recombination_rates_chr9_GRCh38.tsv

core_haplotype=./results_${cohort}_als_uhst_shared_core_haplotype_only_longest.csv
anc_haplotype=./results_${cohort}_als/uhst_mbah_only_longest.csv

# for i in 0,6 7,9 10,14 15,19 20,50 100,100;
# for i in 0,2 3,4 5,6 7,8, 9,10 11,12 13,14 15,16 17,18 19,20 21,22 23,24 25,26 27,28 29,30 31,32 33,34 35,36 37,38 39,40 41,42 43,44 45,46 100,100;
# for i in 0,2 3,4 5,6 7,8, 9,10 11,12 13,14 15,16 17,18 19,20 21,22 23,24 25,26 27,28 29,30 31,32 33,34 35,36 37,38 39,40 41,42 43,44 45,46 100,100;
for i in 2,2 3,3 4,4 5,5 6,6 7,7 8,8 9,9 10,10 11,11 12,12 13,13 14,14 15,15 16,16 17,17 18,18 19,19 20,20 21,21 22,22 23,23 24,24 25,25 26,26 27,27 28,28 29,29 30,30 31,31 31,31 32,32 33,33 34,34 35,35 36,36 37,37 38,38 39,39 40,40 41,41 42,42 43,43 44,44 45,45 46,46 100,100;
 
do IFS=",";
  set -- $i
  output="step_${1}_${2}.ids"
  cat $lengths \
    | qsv py filter "${2} >= int(n) >= ${1}" \
    | qsv select id \
    | rg -v id > ${1}_temp1.ids;
  hatk samples $file > ${1}_temp2.ids
  cat ${1}_temp1.ids ${1}_temp2.ids | sort | uniq -d > $output

  # haptk compare-to-haplotype $file --coords $coords --samples $output --haplotype $anc_haplotype --select only-longest -vvv

  echo "${1}_${2}:"
  haptk compare-to-haplotype $file --coords $coords --samples $output --haplotype $anc_haplotype --select all --prefix "${1}_${2}" -o results_avg_lengths | sed "s/$/,${1}/g" > results_avg_lengths/${1}_sharing.csv
  # haptk compare-to-haplotype $file --coords $coords --samples $output --haplotype $anc_haplotype --select all --prefix "${1}_${2}" -o results_avg_lengths

  # haptk mrca $file --coords $coords --samples $output --only-longest -r $recombination_rates --prefix "${1}_${2}" -o results_mrca
  # haptk check-for-haplotype $file --coords $coords --samples $output --haplotype $core_haplotype -vvv --select all

  # rm $output
  rm ${1}_temp1.ids
  rm ${1}_temp2.ids
done;
