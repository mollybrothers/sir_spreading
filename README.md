# SIR Spreading

# Intro

# Running Megalodon
https://github.com/nanoporetech/megalodon
our current running version: 2.2.3

## necessary configurations
Megalodon is run from a conda environment specified in the `guppy_conda.yml` file

Modified basecalling within megalodon uses the rerio basecalling config file `res_dna_r941_min_modbases-all-context_v001.cfg` (https://github.com/nanoporetech/rerio). Due to an oddity, you must copy the guppy barcoding directory into rerio or megalodon will not run correctly (https://github.com/nanoporetech/rerio/issues/10). This might be fixed in the future.

As of typing this up, Megalodon does not support demultiplexing. It might in the future (https://github.com/nanoporetech/megalodon/issues/43).
For now, a workaround is to basecall using guppy_basecaller including barcoding, which will give you a `barcoding_summary.txt` file. Extract the readIDs corresponding to each barcode using `get_barcodes.py`. Then, run megalodon using the `--read-ids-filename` option. You can then run megalodon separately for each barcode.

## output from each megalodon step (in `savio_jobscripts`):
1. `megalodon_js.sh`
+ `mod_mappings.bam`: BAM file including the tags Mm and Ml, which contain information on modification probabilities for each A in each read.
+ `per_read_modified_base_calls.db`: SQLite database
  
2. `megalodon_aggregate_js.sh`
+ `modified_bases.aggregate.6mA.bed`: bedMethyl file (https://www.encodeproject.org/data-standards/wgbs/)

(does not work yet. For now include the --write-mods-text flag in the first megalodon script. Eventually want this to be separate for computing speed.):
3. `meglaodon_perreadtext_js.sh`
+ `per_read_modified_base_calls.txt`
