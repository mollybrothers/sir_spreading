#!/bin/bash

source activate megalodon-env

EXP="/home/mbrothers/nanopore/210607_Dotty"

ls -d "$exp/megalodon_output_"* | parallel -u -I{} megalodon_extras per_read_text modified_bases {}
