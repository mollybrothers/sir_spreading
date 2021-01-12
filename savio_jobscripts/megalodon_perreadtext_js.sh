#!/bin/bash
#
#SBATCH --job-name=megalodon_perreadtext
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio1
#SBATCH --output=perreadtext.out
#SBATCH --error=perreadtext.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=molly_brothers@berkeley.edu
#SBATCH --time=48:00:00

##Commands to run:

source activate guppy #conda environment, config guppy_conda.yml

MEGALODON_DIRECTORY="/global/scratch/molly_brothers/201012_Doudna/megalodon_tests/"

megalodon_extras per_read_text modified_bases $MEGALODON_DIRECTORY
