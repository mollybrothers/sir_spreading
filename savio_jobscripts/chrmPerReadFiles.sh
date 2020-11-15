#!/bin/bash
#
#SBATCH --job-name=readPerChrom
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_htc
#SBATCH --output=perreadChr.out
#SBATCH --error=perreadChr.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=koenvdberge@berkeley.edu
#SBATCH --time=19:00:00

cd /global/scratch/molly_brothers/201012_Doudna/megalodon_outputv2

awk '$2 ~ /^I$/' per_read_modified_base_calls.txt > chrI.txt
awk '$2 ~ /^II$/' per_read_modified_base_calls.txt > chrII.txt
awk '$2 ~ /^III$/' per_read_modified_base_calls.txt > chrIII.txt
awk '$2 ~ /^IV$/' per_read_modified_base_calls.txt > chrIV.txt
awk '$2 ~ /^V$/' per_read_modified_base_calls.txt > chrV.txt
awk '$2 ~ /^VI$/' per_read_modified_base_calls.txt > chrVI.txt
awk '$2 ~ /^VII$/' per_read_modified_base_calls.txt > chrVII.txt
awk '$2 ~ /^VIII$/' per_read_modified_base_calls.txt > chrVII.txt
awk '$2 ~ /^IX$/' per_read_modified_base_calls.txt > chrIX.txt
awk '$2 ~ /^X$/' per_read_modified_base_calls.txt > chrX.txt
awk '$2 ~ /^XI$/' per_read_modified_base_calls.txt > chrXI.txt
awk '$2 ~ /^XII$/' per_read_modified_base_calls.txt > chrXII.txt
awk '$2 ~ /^XIII$/' per_read_modified_base_calls.txt > chrXII.txt
awk '$2 ~ /^XIV$/' per_read_modified_base_calls.txt > chrXIV.txt
awk '$2 ~ /^XV$/' per_read_modified_base_calls.txt > chrXV.txt
awk '$2 ~ /^XVI$/' per_read_modified_base_calls.txt > chrXVI.txt




#awk -F\t '{print>"$2"}' per_read_modified_base_calls.txt


# get second column (chromosomes), sort, uniq, rm header, and save.
#cut -f 2 per_read_modified_base_calls.txt | sort | uniq | sed '/chrm/d' > allChr
#
#for chr in allChr
#do
#    cat per_read_modified_base_calls.txt | awk '$2 ~ /^"$chr"$/' > "chr${chr}_plus_2.txt"
#done
#
## cat per_read_modified_base_calls.txt | mawk '$2 ~ /^III$/ && $3==1' {print $0}' > chrIII_plus.txt
## cat per_read_modified_base_calls.txt | awk '$2 ~ /^III$/ && $3==1' {print $0}' > hlpKoen.txt
#
#
#
#cut -f 2 chrHlp.txt | sort | uniq > allChr
#while read chr
#do
#    #cat chrHlp.txt | awk '$2 ~ /^"${chr}"$/' > "chr${chr}_plus_2.txt"
#    # cat chrHlp.txt | awk '$2 ~ /^III$/ && $3==1'
#    #awk '$2 ~ /^III$/' chrHlp.txt
#    awk 'BEGIN {print "'"$chr"'"}'
#    awk 'BEGIN {print }'
#    #awk 'BEGIN {$2 ~ /^"'"$chr"'"$/}' chrHlp.txt
#    #awk -v var="$chr"  'BEGIN {$2 ~ /^{$chr}$/}' chrHlp.txt
#done < allChr

