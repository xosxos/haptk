# Fish script
# 
set -l triplets $argv[1]
set -l data_csv $argv[2]

set -l trio_csv "$(cat $triplets | rg -wf  $(cat $data_csv | sed 's/-/,/g' | sed 's/id/id,ht/g' | rg -wf $(cat $triplets | zsv select 1 | psub) | zsv select 1 | psub))";

set -l member 1;
set -l header 'c_id,c_a1,c_a2,c_d1,c_d2';
set -l cohort_gt0 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',0,' '\-0,')";
set -l cohort_gt1 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',1,' '\-1,')";
set -l alleles "$(zsv join --full 1 $(echo $cohort_gt0 | psub) 1 $(echo $cohort_gt1 | psub) | zsv py map -n 'col[0] if col[0] != "" else col[4]' | zsv select 9,3,7,4,8 | sed 's/\.0//g')"
set -l child_joined "$(zsv join --full $member $(echo $trio_csv | psub) 1 $(begin; echo $header; echo $alleles; end | psub))"

set -l member 2;
set -l header 'f_id,f_a1,f_a2,f_d1,f_d2';
set -l cohort_gt0 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',0,' '\-0,')";
set -l cohort_gt1 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',1,' '\-1,')";
set -l alleles "$(zsv join --full 1 $(echo $cohort_gt0 | psub) 1 $(echo $cohort_gt1 | psub) | zsv py map -n 'col[0] if col[0] != "" else col[4]' | zsv select 9,3,7,4,8 | sed 's/\.0//g')"
set -l father_joined "$(zsv join --full $member $(echo $child_joined | psub) 1 $(begin; echo $header; echo $alleles; end | psub))"

set -l member 3;
set -l header 'm_id,m_a1,m_a2,m_d1,m_d2';
set -l cohort_gt0 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',0,' '\-0,')";
set -l cohort_gt1 "$(cat $data_csv | rg -wf $(echo $trio_csv | zsv select $member | psub) | rg -r ',1,' '\-1,')";
set -l alleles "$(zsv join --full 1 $(echo $cohort_gt0 | psub) 1 $(echo $cohort_gt1 | psub) | zsv py map -n 'col[0] if col[0] != "" else col[4]' | zsv select 9,3,7,4,8 | sed 's/\.0//g')"
set -l mother_joined "$(zsv join --full $member $(echo $father_joined | psub) 1 $(begin; echo $header; echo $alleles; end | psub))"

echo $mother_joined | zsv py filter 'f_id != "" and m_id != ""' | zsv select !1,2,3 | table
