#!/bin/bash
#
#SBATCH --job-name=condaEnv
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_htc
#SBATCH --output=condaEnvGuppy.out
#SBATCH --error=condaEnvGuppy.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=koenvdberge@berkeley.edu
#SBATCH --time=48:00:00

module load python #to get conda
cd /global/scratch/molly_brothers/
conda env create -f guppy_conda.yml
# activate with conda activate guppy

#  createCondaEnvironment.sh
#  
#
#  Created by Koen Van den Berge on 1/1/21.
#  
