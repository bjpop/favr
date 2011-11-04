#!/bin/env python

'''
Variant filter script. Bins variants which only appear on 35 base reads in sample.

Authors: Bernie Pope, Danny Park, Tu Ng.

Revision history:

29 Aug 2011.   Initial version.
'''

import os
import pysam
import sys
import csv

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: %s variantfile.csv sample.bam ..." % sys.argv[0])
    variantFilename = sys.argv[1]
    bamFilename = sys.argv[2]
    filterVariants(variantFilename, bamFilename)

def filterVariants(variantFilename, bamFilename):
    variantFile = open(variantFilename)
    binFilename = makeSafeFilename('35s_binfile')
    keepFilename = makeSafeFilename('35s_keepfile')
    logFilename = makeSafeFilename('35s_logfile')
    with pysam.Samfile(bamFilename, "rb") as bam:
        with open(binFilename,'w') as binFile:
            with open(keepFilename,'wb') as keepFile:
                with open(logFilename,'wb') as logFile:
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
    # only process valid rows in the variant CSV file
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
    '''Extract the interesting information about a variant from a SIFT CSV row.'''
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

def parsePolymorphism(s):
    '''Parse the X/Y polymorphism notation, extracting X and Y.'''
    fromTo = s.split('/')
    if len(fromTo) == 2:
        fromBase = fromTo[0]
        toBase = fromTo[1]
        if validBase(fromBase) and validBase(toBase):
            return fromBase,toBase
    return None

def validBase(s):
    '''Check that a letter is a vaid code for a base.'''
    return s in ['G', 'A', 'T', 'C']

def lookupPileup(bam, chr, col):
    '''Retrieve the pileup for a particular chromosome:position coordinate.'''
    # assume argument is in 1-based numbering.
    # samtools gives back a range of pilups that cover the requested coordinates,
    # (no idea why) so we have to search for the particular one we want
    for pileupcolumn in bam.pileup(chr, col-1, col):
        if pileupcolumn.pos == col-1:
            # XXX can't return pileupcolumn here because: segfault!
            # this appears to be a bug in either pysam or samtools (or both)
            return (pileupcolumn.pos , pileupcolumn.n, pileupcolumn.pileups)
    return None

def makeSafeFilename(name):
    if not os.path.exists(name):
        return name
    else:
        count = 2
        while (count < sys.maxint):
            newName = name + str(count)
            if not os.path.exists(newName):
                return newName
            count += 1
    # a safety condition in case we can't return a safe name
    return None

if __name__ == '__main__':
    main()
