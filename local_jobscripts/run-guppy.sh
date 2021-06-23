#!/bin/bash

export PATH=$PATH:/home/mbrothers/tools/ont-guppy/bin

source activate megalodon-env

INPUT1="/mnt/cerberus_ftp/DATA/MollyBrothers/nanopore/210607_Dotty/JRY13114/test"
OUTPUT1="/mnt/cerberus_ftp/DATA/MollyBrothers/nanopore/210607_Dotty/JRY13114/test/guppy_basecall"
CELL="FLO-MIN106"
KIT="SQK-LSK109"

guppy_basecaller --input_path $INPUT1 --save_path $OUTPUT1 --flowcell $CELL --kit $KIT --compress_fastq --device cuda:0,1

echo 'making guppy_logs directory'
mkdir $OUTPUT1/guppy_logs
mv $OUTPUT1/*.log $OUTPUT1/guppy_logs

INPUT2="/mnt/cerberus_ftp/DATA/MollyBrothers/nanopore/210607_Dotty/JRY13114/test/guppy_basecall/pass"
OUTPUT2="/mnt/cerberus_ftp/DATA/MollyBrothers/nanopore/210607_Dotty/JRY13114/test/guppy_basecall/barcode"
KIT="EXP-NBD104"

guppy_barcoder --input_path $INPUT2 --save_path $OUTPUT2 --barcode_kits $KIT --recursive --device cuda:0,1
