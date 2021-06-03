# SIR Spreading

# Intro

# Running Megalodon
https://github.com/nanoporetech/megalodon
our current running version: 2.2.10

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
    + `modified_bases.aggregate.6mA.bed`: bedMethyl file (https://www.encodeproject.org/data-standards/wgbs/). Used for `percentage_methylation.R` and `loess-timecourse.R`

3. `meglaodon_perreadtext_js.sh`
    + `per_read_modified_base_calls.txt`: each row is a single adenine position on a single read. Used for `per_read_probabilities.R`
    
    A [conversation](https://github.com/nanoporetech/rerio/issues/12) with the megalodon developer made clear that the `per_read_modified_base_calls.txt` file contains **all** the data, i.e., data on all adenine nucleotides on all reads considered. In comparison, the `.bed` files are summarized/aggregated methylation calls for each adenine nucleotide. Here, the coverage per adenine is lower than in the `per_read_modified_base_calls.txt` file since by default a non-methylated position is considered a position with methylation probability <0.2, while a methylated position has a methylation probability >0.8.

# Post-Megalodon Analysis

1. `percentage_methylation.R`: make % methylation plots based on data in the bedMethyl output file from `megalodon_aggregate_js.sh`

2. `per_read_probabilities.R`: make single-read plots of methylation probabilities based on data in the `per_read_modified_base_calls.txt` file from `megalodon_per_read_text_js.sh`. This script also filters out the reads with a basecalling qscore < 9 using data from the `sequencing_summary.txt` file output by `megalodon_js.sh`

# Info on experiments
## 200814_McClintock (megalodon v2.2.5)
### Strains:
+ JRY11699 (no EcoGII): barcode01
+ JRY12838 (_sir3∆::EcoGII_): barcode02
+ JRY12840 (_SIR3-EcoGII_): barcode03

Flowcell: FLO-MIN106 / Kit: SQK-RBK004 / N50: 21kb / Yield: 6Gb

## 201012_Doudna (megalodon v2.2.4)
### Strains:
+ JRY13027 (_SIR3-EcoGII_)

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 / N50: 3kb / Yield: 13Gb

## 201125_Turkey (megalodon v2.2.8)
### Strains:
+ JRY11699 (no EcoGII): barcode04
+ JRY12838 (_sir3∆::EcoGII_): barcode05
+ JRY12840 (_SIR3-EcoGII_): barcode06
+ JRY13027 (_SIR3-EcoGII_): barcode07

Flowcell: FLO-MIN106 / Kit: SQK-RBK004 / N50: 22kb / Yield: 12Gb

## 201218_Mariah (megalodon v2.2.8)
Unfortunately bad run. Reads too long (N50 33kb) and pore occupancy started low and quickly depleted. Unfortunately unusable data.
### Strains:
+ JRY13110 (_sir3-8_) – 25C: barcode04
+ JRY13110 (_sir3-8_) – 37C: barcode05
+ JRY13114 (_sir3-8-EcoGII_) – 25C: barcode06
+ JRY13114 (_sir3-8-EcoGII_) – 37C: barcode07
+ JRY13115 (_sir3-8-EcoGII_) – 25C: barcode08
+ JRY13115 (_sir3-8-EcoGII_) – 37C: barcode09
+ JRY13027 (_SIR3-EcoGII_) – 25C: barcode10
+ JRY13027 (_SIR3-EcoGII_) – 37C: barcode11

Flowcell: FLO-MIN106 / Kit: SQK-RBK004 / N50: 33kb / Yield: 6Gb

## 201223_Elf (megalodon v2.2.8)
Data here is still not great in terms of coverage and pore occupancy, but qualitatively you can tell that the sir3-8-EcoGII strains have no methylation at 37C and do have methylation at 25C
### Strains:
+ JRY13114 (_sir3-8-EcoGII_) – 25C: barcode01
+ JRY13114 (_sir3-8-EcoGII_) – 37C: barcode02
+ JRY13115 (_sir3-8-EcoGII_) – 25C: barcode03
+ JRY13115 (_sir3-8-EcoGII_) – 37C: barcode04
+ JRY12840 (_SIR3-EcoGII_) – 25C: barcode05
+ JRY12840 (_SIR3-EcoGII_) – 37C: barcode06

Flowcell: FLO-MIN106 / Kit: SQK-RBK004 / N50: 3kb / Yield: 3Gb

## 210205_Ocular (megalodon v2.2.8)
Fragment size is a little smaller than I would have wanted (15-20kb), due to spinning the Covaris g-TUBE faster than I was supposed to. BUT the nanopore ran great and I got plenty of data.
### Strains:
+ JRY9316 (no EcoGII) - 25C: barcode07
+ JRY13114 (_sir3-8-EcoGII_) - 25C: barcode08
+ JRY13027 (_SIR3-EcoGII_) - 25C: barcode09
+ JRY9316 (no EcoGII) - 37C: barcode10
+ JRY13114 (_sir3-8-EcoGII_) - 37C: barcode11
+ JRY13027 (_SIR3-EcoGII_) - 37C: barcode12

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 7.5kb / Yield: 15.7Gb

## 210208_Rescue (megalodon v2.2.10)
### Strains:
+ JRY13114 (_sir3-8-EcoGII_) - 0min: barcode01
+ JRY13114 (_sir3-8-EcoGII_) - 5min: barcode02
+ JRY13114 (_sir3-8-EcoGII_) - 15min: barcode03
+ JRY13114 (_sir3-8-EcoGII_) - 45min: barcode04
+ JRY13114 (_sir3-8-EcoGII_) - 90min: barcode05
+ JRY13114 (_sir3-8-EcoGII_) - 150min: barcode06

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 10.6kb / Yield: 16Gb

## 210304_Amanita (megalodon v2.2.10)
Definitely see methylation here, but appears much more slowly (barely anything by time 90min) than with 210208_Rescue. Perhaps because I'm using a different strain...perhaps because of variability, perhaps because of culture density (these were at a higher OD when I started the experiment than 210208_Rescue). 210310_Russula, which uses JRY13114 like 210208_Rescue, was also taken at a higher culture density than 210208_Rescue and also seems to have slower kinetics.
### Strains:
+ JRY13134 (_sir3-8-EcoGII_) - 0min: barcode07
+ JRY13134 (_sir3-8-EcoGII_) - 15min: barcode08
+ JRY13134 (_sir3-8-EcoGII_) - 30min: barcode09
+ JRY13134 (_sir3-8-EcoGII_) - 45min: barcode10
+ JRY13134 (_sir3-8-EcoGII_) - 60min: barcode11
+ JRY13134 (_sir3-8-EcoGII_) - 90min: barcode12

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 11.3kb / Yield: 15.9Gb

## 210306_Ecology (have not called with megalodon yet)
Accidentally did the experiment with JRY13140 (_sir3-8, hmr-i∆_) instead of JRY13141 (_sir3-8-EcoGII, hmr-i∆_). So I am not analyzing it yet but could potentially use it as a (way overkill) negative control in the future.
### Strains:
+ JRY13140 (_sir3-8, HMR-i∆_) - 0min: barcode01
+ JRY13140 (_sir3-8, HMR-i∆_) - 15min: barcode02
+ JRY13140 (_sir3-8, HMR-i∆_) - 30min: barcode03
+ JRY13140 (_sir3-8, HMR-i∆_) - 45min: barcode04
+ JRY13140 (_sir3-8, HMR-i∆_) - 60min: barcode05
+ JRY13140 (_sir3-8, HMR-i∆_) - 90min: barcode06

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 11kb / Yield: 24Gb

## 210310_Russula (megalodon v2.2.10)
Same strain as 210208_Rescue, but taken at a higher culture density and seems to have slower kinetics (similar to 210304_Amanita).
### Strains:
+ JRY13114 (_sir3-8-EcoGII_) - 0min: barcode01
+ JRY13114 (_sir3-8-EcoGII_) - 15min: barcode02
+ JRY13114 (_sir3-8-EcoGII_) - 30min: barcode03
+ JRY13114 (_sir3-8-EcoGII_) - 45min: barcode04
+ JRY13114 (_sir3-8-EcoGII_) - 60min: barcode05
+ JRY13114 (_sir3-8-EcoGII_) - 90min: barcode06

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 9.3kb / Yield: 16.5Gb

## 210311_Hygrocybe (megalodon v2.2.10)
Kinetics are slower than Rescue, just like Amanita and Russula runs.
### Strains:
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 0min: barcode07
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 15min: barcode08
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 30min: barcode09
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 45min: barcode10
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 60min: barcode11
+ JRY13141 (_sir3-8-EcoGII, HMR-i∆_) - 90min: barcode12

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 11.5kb / Yield: 16Gb

## 210403_Hello (megalodon v2.2.10)
### Strains:
+ JRY9316 (no EcoGII) - 25C: barcode01
+ JRY13114 (_sir3-8-EcoGII_) - 25C: barcode02
+ JRY13220 (_sir3-8-Hia5_) - 25C: barcode03

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 11kb / Yield: 12.5Gb

## 210531_Texas (megalodon ??)
Grew these strains up to log phase overnight. Unfortunately, 12 samples is too much for the nanopore, so I don't have enough data to be happy with this run.
### Strains:
+ JRY13134 (_sir3-8-EcoGII_) - 0min: barcode01
+ JRY13134 (_sir3-8-EcoGII_) - 30min: barcode02
+ JRY13134 (_sir3-8-EcoGII_) - 60min: barcode03
+ JRY13134 (_sir3-8-EcoGII_) - 90min: barcode04
+ JRY13134 (_sir3-8-EcoGII_) - 120min: barcode05
+ JRY13134 (_sir3-8-EcoGII_) - 150min: barcode06
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 0min: barcode07
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 30min: barcode08
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 60min: barcode09
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 90min: barcode10
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 120min: barcode11
+ JRY13213 (_sir3-8-EcoGII, dot1∆_) - 150min: barcode12

Flowcell: FLO-MIN106 / Kit: SQK-LSK109 + EXP-NBD104 / N50: 9.6kb / Yield: 16.7Gb
