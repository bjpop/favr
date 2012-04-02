#!/bin/env python

'''
Variant filter script. Bins variants which only appear on 35 base reads in sample.

Authors: Bernie Pope, Danny Park, Tu Nguyen-Dumont.

Revision history:

29 Aug 2011.   Initial version.
'''

import os
import pysam
import sys
import csv
import getopt
from favr_common import (parsePolymorphism, lookupPileup, makeSafeFilename)

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --variants=<variant list as TSV file>
    --bam=<bam file of reads for the same sample as variants>
    --bin=<bin filename>
    --keep=<keep filename>
    --log=<log filename>""" % sys.argv[0])

longOptionsFlags = ["help", "variants=", "bam=", "bin=", "keep=", "log="]
shortOptionsFlags = "h"

class Options(object):
    def __init__(self):
        self.variants = None
        self.bin = None
        self.keep = None
        self.log = None
        self.bam = None
    def check(self):
        return all([self.variants, self.bam, self.bin, self.keep, self.log])

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptionsFlags, longOptionsFlags)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    options = Options()
    for o, a in opts:
        if o == "--variants":
            options.variants = a
        elif o == "--bam":
            options.bam = a
        elif o == "--bin":
            options.bin = a
        elif o == "--keep":
            options.keep = a
        elif o == "--log":
            options.log = a
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
    if not options.check():
        print('Incorrect arguments')
        usage()
        exit(2)
    filterVariants(options)

def filterVariants(options):
    variantFile = open(options.variants)
    with pysam.Samfile(options.bam, "rb") as bam:
        with open(options.bin, 'w') as binFile:
            with open(options.keep,'wb') as keepFile:
                with open(options.log, 'wb') as logFile:
                    for variant in csv.reader(variantFile, delimiter=',', quotechar='|'):
                        variantStr = ','.join(variant)
                        thirty_fives,fifties = count_read_sizes(variant, bam)
                        logFile.write('%s: 35s=%d, 50s=%d' % (variantStr, thirty_fives, fifties))
                        if thirty_fives > 0 and fifties == 0:
                            binFile.write('%s\n' % variantStr)
                            logFile.write(', bin\n')
                        else:
                            keepFile.write('%s\n' % variantStr)
                            logFile.write(', keep\n')
    variantFile.close()

def count_read_sizes(variant, bamFile):
    thirty_fives = 0 # number of variants on 35 read in pair
    fifties = 0 # number of variants on 50 read in pair
    info = parseVariantRow(variant)
    # only process valid rows in the variant TSV file
    if info:
        # get the pileup information for this particular variant coordinates
        # the pileup tells us what base was called in each of the reads at
        # this particular coordinate
        pileupCol = lookupPileup(bamFile, info.chromosome, info.position)
        if pileupCol:
            #reads = pileupCol
            pos,coverage,reads = pileupCol
            # Count the number of reads that are the same as the variant
            for pileupread in reads:
                # find the base at the same position as the variant
                readBase = pileupread.alignment.seq[pileupread.qpos]
                if not pileupread.is_del and readBase == info.variantBase:
                    # get the size of the read for the base using the cigar of the read
                    cigar = pileupread.alignment.cigar
                    if cigar:
                        read_size = 0
                        for code,length in cigar:
                            read_size += length
                        # count the number of 35, 50 and unknown length reads
                        if read_size >= 20 and read_size <= 35:
                            thirty_fives += 1
                        elif read_size >= 40 and read_size <= 50:
                            fifties += 1
                    # we couldn't figure out the size from the cigar
                    # skip and print a warning
                    else:
                        print('Warning: could not find cigar info for variant: %s, skipping' % str(variant))
        else:
            print('Warning: could not find pileup for variant: %s' % str(variant))
    return thirty_fives, fifties

class VariantInfo(object):
    '''Interesting information about a particular variant.'''
    def __init__(self, id, chromosome, position, refBase, variantBase, inputRow):
        self.id = id
        self.chromosome = chromosome
        self.position = position
        self.refBase = refBase
        self.variantBase = variantBase
        self.inputRow = inputRow

def parseVariantRow(row):
    '''Extract the interesting information about a variant from a SIFT TSV row.'''
    if len(row) >= 4:
        fromTo = parsePolymorphism(row[3])
        if fromTo:
           chrName = "chr" + row[0]
           return VariantInfo(id = "%s:%s" % (chrName,row[1]),
                              chromosome = chrName,
                              position = int(row[1]), # XXX should really check it is all digits
                              refBase = fromTo[0],
                              variantBase = fromTo[1],
                              inputRow = row)
    return None

if __name__ == '__main__':
    main()
