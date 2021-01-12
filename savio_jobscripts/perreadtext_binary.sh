#!/bin/bash
#
#SBATCH --job-name=perReadBinary
#SBATCH --account=fc_nanosir
#SBATCH --partition=savio2_htc
#SBATCH --output=perReadBinary.out
#SBATCH --error=perREadBinary.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=koenvdberge@berkeley.edu
#SBATCH --time=44:00:00

cd /global/scratch/molly_brothers/201125_Turkey/megalodon_output_06/
awk '(NR>1) && ($5 > -0.2231436) || ($5 < -1.609438) ' per_read_modified_base_calls.txt > /global/scratch/koenvdberge/201125_Turkey/megalodon_output_06/per_read_modified_base_calls_binary.txt

cd /global/scratch/molly_brothers/201125_Turkey/megalodon_output_07/
awk '(NR>1) && ($5 > -0.2231436) || ($5 < -1.609438) ' per_read_modified_base_calls.txt > /global/scratch/koenvdberge/201125_Turkey/megalodon_output_07/per_read_modified_base_calls_binary.txt


#  perreadtext_binary.sh
#
#
#  Created by Koen Van den Berge on 1/1/21.
#
