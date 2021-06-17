#!/bin/bash

GUPPY="/home/mbrothers/tools/ont-guppy/bin/guppy_basecall_server"
GUPPY_PARAMS="-d /home/mbrothers/tools/rerio/basecall_models/"
GUPPY_CONFIG="res_dna_r941_min_modbases-all-context_v001.cfg"
INPUT="/home/mbrothers/210607_Dotty/raw_data_multifast5"
OUTPUT="/home/mbrothers/210607_Dotty/megalodon_output_01"
MOD_MOTIF="Y A 0"
FILES_OUT="basecalls mod_mappings per_read_mods"
GENOME="/home/mbrothers/genomes/genome_mat_to_N.fa"
PROCESSES=4
NUM_READS="1000"
READ_IDS="/home/mbrothers/210607_Dotty/barcode01_readIDs.txt"
GPU=0

megalodon $INPUT \
--output-directory $OUTPUT \
--guppy-server-path $GUPPY \
--guppy-params "$GUPPY_PARAMS" \
--guppy-config $GUPPY_CONFIG \
--mod-motif $MOD_MOTIF \
--reference $GENOME \
--outputs $FILES_OUT \
--processes $PROCESSES \
--devices $GPU \
--read-ids-filename $READ_IDS \
--num-reads $NUM_READS
