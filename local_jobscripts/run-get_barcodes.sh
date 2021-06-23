#!/bin/bash

for x in "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" "13"
    do
        python /home/mbrothers/scripts/get_barcodes.py -i /home/mbrothers/210607_Dotty/guppy_basecall/barcode/barcoding_summary.txt -o /home/mbrothers/210607_Dotty/ -b "barcode$x"
    done
