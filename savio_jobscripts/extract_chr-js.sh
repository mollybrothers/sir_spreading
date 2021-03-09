#!/bin/bash

#SBATCH --job-name=extract_chr
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --output=extract_chr.out
#SBATCH --error=extract_chr.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=molly_brothers@berkeley.edu
#SBATCH --time=02:00:00


##Commands to run:
module load gnu-parallel

exp="/global/scratch/molly_brothers/210304_Amanita"
chr="III"

ls -d "$exp/megalodon_output_"* | parallel -u -I{} /global/home/users/molly_brothers/jobscripts/extract_chr.sh {} $chr
