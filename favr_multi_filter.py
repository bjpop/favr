#!/bin/env python

# Tool to filter out 'multis' from e.g. exome variant lists
#
# Revision History:
#
# 6th August 2012, Daniel J Park, initial version
# 7th August 2012, Bernie Pope, removed loop on varlist

import csv
from argparse import ArgumentParser

def parseArgs():
    parser = ArgumentParser(
        description = 'Tool to filter out multis from vcf files without files')
    parser.add_argument(
        '--VCF', metavar = 'VCFFILE', required = True,
        help = 'name of input vcf file')
    parser.add_argument(
        '--number', metavar = 'NUMBERVARS', type = int, required = True,
        help = 'number of other variants acceptable within window, not counting the queried variant')
    parser.add_argument(
        '--halfwindow', metavar = 'HALFWINDOW', type = int, required = True,
        help = 'size of half the window within which to assess acceptable number')
    return parser.parse_args()

def main():
    args = parseArgs()
    with open(args.VCF, 'U') as vcf_file:
        reader = csv.reader(open(args.VCF, 'rU'), delimiter='\t', quotechar='"')
        varlist = list(reader)
        for var in varlist:
            chr = var[0]
            coord = int(var[1])
            if checker(chr, coord, args.number, args.halfwindow, varlist):
                print("\t".join(var))

def checker(chr, coord, number, halfwindow, varlist):
    counter = 0 
    for item in varlist:
        if item[0] == chr and abs(int(item[1]) - coord) <= halfwindow:
            counter += 1
    return counter <= number

if __name__ == '__main__':
    main()
