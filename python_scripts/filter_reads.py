#!/bin/python3

##########################
# Author: Molly Brothers
# Github: mollybrothers
# Date: 11/18/2020
##########################

# This script uses an input sequencing_summary.txt file from a guppy basecalling run to find the read_ids
# that have a mean_qscore_template < 7 and generate a list of those read_ids.
# Then, it reads through the per_read_modified_base_calls.txt file from a megalodon output and filters out
# any lines corresponding to those read_ids

# import sequencing_summary.txt as a numpy array

# read each line. If mean_qscore_template is < 7, add the read_id to a txt file (one read_id per line)

# import per_read_modified_base_calls.txt as a numpy array

# read each line. If read_id matches any read_id in the read_ids to filter out (text file made above),
# get rid of that line

# save the new array as per_read_modified_base_calls_filteredQ7.txt
