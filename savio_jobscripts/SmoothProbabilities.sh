#!/bin/bash
#
#SBATCH --job-name=smoothProbR
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_htc
#SBATCH --output=smoothProb.out
#SBATCH --error=smoothProb.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=koenvdberge@berkeley.edu
#SBATCH --time=19:00:00


cd /global/scratch/molly_brothers/201012_Doudna/megalodon_outputv2

module load r
module load r-packages

# added new pandoc to path
#PATH=$PATH:/global/scratch/molly_brothers/201012_Doudna/megalodon_outputv2/pandoc-2.11.1.1/bin/

R -e "rmarkdown::render('20201114_smoothedProbChromsome.Rmd')"



#  SmoothProbabilities.sh
#  
#
#  Created by Koen Van den Berge on 11/14/20.
#  
