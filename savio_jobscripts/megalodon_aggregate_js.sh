#!/bin/bash

#SBATCH --job-name=megalodon_aggregate
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio
#SBATCH --nodes=1
#SBATCH --output=aggregate_test.out
#SBATCH --error=aggregate_test.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=molly_brothers@berkeley.edu
#SBATCH --time=48:00:00

##Commands to run:

source activate guppy #conda environment configuration is text file guppy_conda.yml
module load samtools

MEGALODON_DIRECTORY="/global/scratch/molly_brothers/201012_Doudna/megalodon_tests"
OUTPUT_TYPE="mods"
OUTPUT_SUFFIX="aggregate"
PROCESSES=$SLURM_CPUS_ON_NODE

megalodon_extras aggregate run \
--megalodon-directory $MEGALODON_DIRECTORY \
--outputs $OUTPUT_TYPE \
--output-suffix $OUTPUT_SUFFIX \
--processes $PROCESSES

##remove the empty 5mC bed file
rm $MEGALODON_DIRECTORY/*5mC*
