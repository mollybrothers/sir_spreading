#!/bin/awk

##############################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 2020-12-14
#############################

# ex run:
# awk -f single_bp_Chereji-etal.awk input.bed

#this script takes in the bed/bedGraph file of nucleosome occupancy from Chereji et al 2018 (Henikoff lab)
#and turns it into a file where each line is a single base pair

BEGIN {
	printf "chr\tstart\tend\toccupancy\n"
} 
{
	for(i=$2; i<$3; ++i) 
		printf $1 "\t" i "\t" i + 1 "\t" $4 "\n"
}
