#!/bin/bash

source activate megalodon-env

MEGALODON_DIRECTORY="/home/mbrothers/nanopore/210607_Dotty/megalodon_output_01"
OUTPUT_TYPE="mods"
OUTPUT_SUFFIX="aggregate01"
PROCESSES=4

megalodon_extras aggregate run \
--megalodon-directory $MEGALODON_DIRECTORY \
--outputs $OUTPUT_TYPE \
--output-suffix $OUTPUT_SUFFIX \
--processes $PROCESSES

##remove the empty 5mC bed file
rm $MEGALODON_DIRECTORY/*5mC*
