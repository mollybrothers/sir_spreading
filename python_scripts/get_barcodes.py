#!/bin/env python3

import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_file', required=True, type=str, help='input "barcoding_summary.txt" file')
parser.add_argument('-b', '--barcode', required=True, type=str, help="name of barcode in format 'barcodeXX'")
parser.add_argument('-o', '--out_dir', required=True, type=str, help='output directory path')
args = parser.parse_args()

input = args.in_file
in_file = open(input, 'r')
out_dir = args.out_dir
barcode = args.barcode

#create a text file of readIDs that correspond to the barcode
with open(out_dir + barcode +'_readIDs' + '.txt', 'w') as out_file:
    for line in in_file:
        values = line.split()
        if values[1]==barcode:
            out_file.write(values[0] + '\n')
