#!/bin/bash

#SBATCH --job-name=guppy
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_1080ti
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --time=40:00:00
#SBATCH --output=guppy.out
#SBATCH --error=guppy.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=molly_brothers@berkeley.edu

##Commands to run:

export PATH=$PATH:/global/home/users/molly_brothers/sources/ont-guppy/bin

source activate guppy
module load cuda/10.0

INPUT1="/global/scratch/molly_brothers/201218_Mariah/raw_data_multifast5"
OUTPUT1="/global/scratch/molly_brothers/201218_Mariah/guppy_basecall"
CELL="FLO-MIN106"
KIT="SQK-RBK004"

guppy_basecaller --input_path $INPUT1 --save_path $OUTPUT1 --flowcell $CELL --kit $KIT --compress_fastq --device cuda:0

echo 'making fastq directory'
mkdir $OUTPUT1/fastq
mv $OUTPUT1/*.fastq $OUTPUT1/fastq

echo 'making guppy_logs directory'
mkdir $OUTPUT1/guppy_logs
mv $OUTPUT1/*.log $OUTPUT1/guppy_logs

INPUT2="/global/scratch/molly_brothers/201218_Mariah/guppy_basecall/fastq"
OUTPUT2="/global/scratch/molly_brothers/201218_Mariah/guppy_basecall/barcode"
KIT="SQK-RBK004"

guppy_barcoder --input_path $INPUT2 --save_path $OUTPUT2 --barcode_kits $KIT --recursive --device cuda:0
