# SIR SPREADING

# INTRO

# RUNNING MEGALODON
https://github.com/nanoporetech/megalodon
our current running version: 2.2.3

Megalodon is run from a conda environment specified in the `guppy_conda.yml` file

Modified basecalling within megalodon uses the rerio basecalling config file `res_dna_r941_min_modbases-all-context_v001.cfg` (https://github.com/nanoporetech/rerio). Due to an oddity, you must copy the guppy barcoding directory into rerio or megalodon will not run correctly (https://github.com/nanoporetech/rerio/issues/10). This might be fixed in the future.

As of typing this up, Megalodon does not support demultiplexing. It might in the future (https://github.com/nanoporetech/megalodon/issues/43).
For now, a workaround is to basecall using guppy_basecaller including barcoding, which will give you a `barcoding_summary.txt` file. Extract the readIDs corresponding to each barcode using `get_barcodes.py`. Then, run megalodon using the `--read-ids-filename` option. You can then run megalodon separately for each barcode.
