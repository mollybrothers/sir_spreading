#!/bin/bash

exp="/home/mbrothers/210607_Dotty"
chr="III"

ls -d "$exp/megalodon_output_"* | parallel -u -I{} /home/mbrothers/nanopore/scripts/extract_chr.sh {} $chr
