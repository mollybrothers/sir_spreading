#!/bin/bash

#SBATCH --job-name=megalodon
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_1080ti
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2
#SBATCH --output=megalodontests.out
#SBATCH --error=megalodontests.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=molly_brothers@berkeley.edu
#SBATCH --time=48:00:00

##commands to run:
##still to do: separate the GPU and CPU functions of megalodon for faster computing time.
##Right now this script runs everything on the GPU, but there is a way to separate the aggregation functions.
##This will require two different job scripts. One for the basecalling and one for the aggregation.

source activate guppy
module load cuda/10.0
module load samtools

GUPPY="/global/home/users/molly_brothers/sources/ont-guppy/bin/guppy_basecall_server"
GUPPY_PARAMS="-d /global/home/users/molly_brothers/tools/rerio/basecall_models/"
GUPPY_CONFIG="res_dna_r941_min_modbases-all-context_v001.cfg"
INPUT="/global/scratch/molly_brothers/201012_Doudna/raw_data_multifast5"
OUTPUT="/global/scratch/molly_brothers/201012_Doudna/megalodon_tests"
MOD_MOTIF="Y A 0"
FILES_OUT="basecalls mod_mappings"
GENOME="/global/scratch/molly_brothers/genomes/genome_mat_to_N.fa"
PROCESSES="$SLURM_CPUS_ON_NODE"
##NUM_READS=100
##READ_IDS="/global/scratch/molly_brothers/200814_McClintock/try2/basecall/barcode/barcode03_readIDs.txt"

megalodon $INPUT \
--output-directory $OUTPUT \
--guppy-server-path $GUPPY \
--guppy-params "$GUPPY_PARAMS" \
--guppy-config $GUPPY_CONFIG \
--mod-motif $MOD_MOTIF \
--reference $GENOME \
--outputs $FILES_OUT \
--write-mods-text \
--processes $PROCESSES \
--devices 0 1 \
##--num-reads $NUM_READS
##--read-ids-filename $READ_IDS

echo 'done'
