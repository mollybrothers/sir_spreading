#!/bin/bash

#input a megalodon_output directory that contains a per_read_modified_basecalls.txt file
#input a chromosome number in Roman numeral format
#OUTPUT: a chr.txt file containing the per_read modified base calls for that chromosome

if [ ! $# == 2 ] ; then
    echo "Incorrect number of arguments"
    echo "Usage: bash extract_chr.sh path/to/megalodon_output_xx/per_read_modified_base_calls.txt chr-number"
    exit
fi

dir=$1
input="$dir/per_read_modified_base_calls.txt"
chr=$2
barcode="${dir##*_}" #get the barcode from megalodon file name
pattern="^$chr$"

if [ -f "$input" ] ; then
    echo "starting chr$chr extraction from $dir's per_read file"
else
    echo "$input does not exist"
    exit
fi

cat $input | awk -v pat="$pattern" '$2 ~ pat {print $0}' > "$dir/chr${chr}_${barcode}.txt"
