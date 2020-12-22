# SIR Spreading

# Intro

# Running Megalodon
https://github.com/nanoporetech/megalodon
our current running version: 2.2.8

## necessary configurations
Megalodon is run from a conda environment specified in the `guppy_conda.yml` file

Modified basecalling within megalodon uses the rerio basecalling config file `res_dna_r941_min_modbases-all-context_v001.cfg` (https://github.com/nanoporetech/rerio). Due to an oddity, you must copy the guppy barcoding directory into rerio or megalodon will not run correctly (https://github.com/nanoporetech/rerio/issues/10). This might be fixed in the future.

As of typing this up, Megalodon does not support demultiplexing. It might in the future (https://github.com/nanoporetech/megalodon/issues/43).
For now, a workaround is to basecall using guppy_basecaller and guppy_barcoding, which will give you a `barcoding_summary.txt` file (see `guppy-js.sh` script). Extract the readIDs corresponding to each barcode using `get_barcodes.py`. Then, run megalodon using the `--read-ids-filename` option. You can then run megalodon separately for each barcode.

## output from each megalodon step (in `savio_jobscripts`):
1. `megalodon_js.sh`
    + `mod_mappings.bam`: BAM file including the tags `Mm` and `Ml`, which contain information on modification probabilities for each A in each read.
    + `per_read_modified_base_calls.db`: SQLite database
    + `sequencing_summary.txt`: summary of the guppy run. Used for filtering reads on Qscores in post-megalodon analyses (does not contain modification information)
    + `basecalls.fastq`: basecalls from guppy (canonical bases only) in fastq form
  
2. `megalodon_aggregate_js.sh`
    + `modified_bases.aggregate.6mA.bed`: bedMethyl file (https://www.encodeproject.org/data-standards/wgbs/). Used for `percentage_methylation.R`

3. `meglaodon_perreadtext_js.sh`
    + `per_read_modified_base_calls.txt`: each row is a single adenine position on a single read. Used for `per_read_probabilities.R`

# Post-Megalodon Analysis

1. `percentage_methylation.R`: make % methylation plots based on data in the bedMethyl output file from `megalodon_aggregate_js.sh`

2. `per_read_probabilities.R`: make single-read plots of methylation probabilities based on data in the `per_read_modified_base_calls.txt` file from `megalodon_per_read_text_js.sh`. This script also filters out the reads with a basecalling qscore < 9 using data from the `sequencing_summary.txt` file output by `megalodon_js.sh`

# Info on experiments
## 200814_McClintock (megalodon v 2.2.5)
### Strains:
1. JRY11699 (no EcoGII): barcode01
2. JRY12838 (_sir3∆::EcoGII_): barcode02
3. JRY12840 (_SIR3-EcoGII_): barcode03

## 201012_Doudna (megalodon v 2.2.4)
### Strains:
1. JRY13027 (_SIR3-EcoGII_)

## 201125_Turkey (megalodon v 2.2.8)
### Strains:
1. JRY11699 (no EcoGII): barcode04
2. JRY12838 (_sir3∆::EcoGII_): barcode05
3. JRY12840 (_SIR3-EcoGII_): barcode06
4. JRY13027 (_SIR3-EcoGII_): barcode07

## 201218_Mariah (megalodon v 2.2.8)
### Strains:
1. JRY13110 (_sir3-8_) – 25C: barcode04
2. JRY13110 (_sir3-8_) – 37C: barcode05
3. JRY13114 (_sir3-8-EcoGII_) – 25C: barcode06
4. JRY13114 (_sir3-8-EcoGII_) – 37C: barcode07
5. JRY13115 (_sir3-8-EcoGII_) – 25C: barcode08
6. JRY13115 (_sir3-8-EcoGII_) – 37C: barcode09
7. JRY13027 (_SIR3-EcoGII_) – 25C: barcode10
8. JRY13027 (_SIR3-EcoGII_) – 37C: barcode11
