#!/bin/bash

file=$1
lengths=$2
cohort=$3
coords=$4
outdir=$5

if [[ -z $lengths || -z $file || -z $cohort || -z $coords || -z $outdir ]]
then
  echo "error: missing input variables, check script for info"
  exit
fi
mkdir -p $outdir

recombination_rates=./test_data/recombination_rates_chr9_GRCh38.tsv

core_haplotype=./results_${cohort}_c9orf72/uhst_shared_core_haplotype_only_longest.csv
ancestral_haplotype=./results_${cohort}_c9orf72/uhst_mbah_only_longest.csv

printf "group,n,hap_n,has_core_ht,avg_markers,median,max,mrca,mrca_95_low,mrca_95_high\n"

for i in 0,6 7,9 10,14 15,19 20,50 100,100;
# for i in 2,2 3,3 4,4 5,5 6,6 7,7 8,8 9,9 10,10 11,11 12,12 13,13 14,14 15,15 16,16 17,17 18,18 19,19 20,20 21,21 22,22 23,23 24,24 25,25 26,26 27,27 28,28 29,29 30,30 31,31 31,31 32,32 33,33 34,34 35,35 36,36 37,37 38,38 39,39 40,40 41,41 42,42 43,43 44,44 45,45 46,46 100,100;
 
do IFS=",";
  set -- $i

  group="step_${1}_${2}.ids"

  cat $lengths \
    | qsv py filter "${2} >= int(n) >= ${1}" \
    | qsv select id \
    | rg -v id > $outdir/${1}_temp1.ids;

  haptk samples $file > $outdir/${1}_temp2.ids

  cat $outdir/${1}_temp1.ids $outdir/${1}_temp2.ids | sort | uniq -d > $outdir/$group

  if [ -s $outdir/$group ];then
      nsamples=$(cat $outdir/$group | wc -l)

      haptk compare-to-haplotype $file -c $coords --samples $outdir/$group --haplotype $ancestral_haplotype -s all --prefix "${1}_${2}" -o $outdir --silent

      haptk mrca $file --coords $coords -S $outdir/$group -s all -r $recombination_rates -p "${1}_${2}" -o $outdir --silent

      haptk check-for-haplotype $file -c $coords -S $outdir/$group --haplotype $core_haplotype -s all -p "${1}_${2}" -o $outdir --silent

      avg_markers=$(cat $outdir/${1}_${2}_ht_shared_segments.csv | qsv stats -s markers --median | qsv select mean | qsv table | rg -v mean)
      median=$(cat $outdir/${1}_${2}_ht_shared_segments.csv | qsv stats -s markers --median | qsv select median | qsv table | rg -v median)
      max=$(cat $outdir/${1}_${2}_ht_shared_segments.csv | qsv stats -s markers --median | qsv select max | qsv table | rg -v max)


      if [ "1" == $nsamples ]; then
        mrca="na"
        mrca_ci="na,na"
      else
        mrca_ci=$(cat $outdir/${1}_${2}_mrca_gamma_method_only_longest.txt | rg "Corr" -A 1 | rg -v "Corr" | sed 's/^.*CI //g' | sed -E 's/\(|\)| //g')
        # mrca_ci=$(cat $outdir/${1}_${2}_mrca_gamma_method.txt | rg "Corr" -A 1 | rg -v "Corr" | sed 's/^.*CI //g' | sed -E 's/\(|\)| //g')
        mrca=$(cat $outdir/${1}_${2}_mrca_gamma_method_only_longest.txt | rg "Corr" -A 1 | rg -v "Corr" | sed 's/age: //g' | sed 's/ CI.*$//g')
        # mrca=$(cat $outdir/${1}_${2}_mrca_gamma_method.txt | rg "Corr" -A 1 | rg -v "Corr" | sed 's/age: //g' | sed 's/ CI.*$//g')
      fi

      core_ht=$(cat $outdir/${1}_${2}_haplotype_check.csv | rg true | wc -l)
      nhaplos=$(( $nsamples * 2 ))
      printf "${1}_${2},${nsamples},${nhaplos},${core_ht},${avg_markers},${median},${max},${mrca},${mrca_ci}\n"
  fi

  rm $outdir/${1}_temp1.ids
  rm $outdir/${1}_temp2.ids
done;
