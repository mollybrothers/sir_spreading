#!/bin/bash

source activate megalodon-env

MEGALODON_DIRECTORY="/home/mbrothers/nanopore/210607_Dotty/megalodon_output_01"

megalodon_extras per_read_text modified_bases $MEGALODON_DIRECTORY
