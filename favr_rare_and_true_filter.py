#!/bin/env python

'''
Variant filter script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Nguyen-Dumont.

Reads a list of variants from a TSV file for a sample and compares
them to the sequence reads in BAM files for other samples.

Revision history:

2 May 2011.    Initial incomplete version. Can read TSV file and BAMS and
               generate a table of which variants from the original sample
               are found in the other samples.

30 May 2011.   Filter out reads which are considered 'is_del'. These appear
               to come from deletions, and don't belong in the pileup for
               a given position.

31 May 2011.   Change the supported input format from a custom TSV file to
               a SIFT (csv) file.

25 Aug 2011.   Changed the binning rule to bin 35-length only reads, in
               addition to the previous binning rule.

29 Aug 2011.   Changed back to not bin 35-length only reads. This is now
               done in another tool.

19 Sep 2011.   Added classification for family samples.

20 Sep 2011.   Allowed the input to be tab separated, with coordinates comma separated.

14 Nov 2011.   Added classify arguments to the command line.

'''

import os
import pysam
import sys
import csv
import yaml
import getopt
from favr_common import (safeReadInt, getEvidence, makeSafeFilename, sortByCoord)
from favr_rare_and_true_classify import classify

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --variants=<variant list as TSV file>
    --bin=<bin filename>
    --keep=<keep filename>
    --log=<log filename>
    --varLikeThresh=<variant read threshold>
    --samplesPercent=<percent of total samples which pass the threshold>
    reads1.bam reads2.bam ...""") % sys.argv[0]

longOptionsFlags = ["help", "variants=", "bin=", "keep=", "log=", "varLikeThresh=", "samplesPercent="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.variants = None
        self.bin = None
        self.keep = None
        self.log = None
        self.varLikeThresh = None
        self.samplesPercent = None
    def check(self):
        return all([self.variants, self.bin, self.keep, self.log, self.varLikeThresh, self.samplesPercent])

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
        elif o == "--bin":
            options.bin = a
        elif o == "--keep":
            options.keep = a
        elif o == "--log":
            options.log = a
        elif o == "--varLikeThresh":
            options.varLikeThresh = safeReadInt(a)
        elif o == "--samplesPercent":
            options.samplesPercent = safeReadInt(a)
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
    if not options.check():
        print('Incorrect arguments')
        usage()
        exit(2)
    bamFilenames = args
    # Read the rows of the variants TSV file into a list.
    with open(options.variants) as variants:
        variantList = list(csv.reader(variants, delimiter='\t', quotechar='|'))
    # compute the presence/absence of each variant in the bam files
    evidence = getEvidence(variantList, bamFilenames)
    # filter the variants
    filter(options, evidence)

def filter(options, evidence):
    '''Decide which variants to keep and which to bin.'''
    binFilename =  options.bin
    keepFilename = options.keep
    logFilename = options.log
    with open (logFilename,'w') as logFile:
        with open(binFilename,'w') as binFile:
            with open(keepFilename,'wb') as keepFile:
                csvWriter = csv.writer(keepFile, delimiter='\t', quotechar='|')
                # sort the variants by coordinate
                for key,info in sortByCoord(evidence):
                    classification = classify(options, info.counts)
                    # record the classification of this variant in the logfile
                    logFile.write("%s: %s: %s\n" % (key, classification.action, classification.reason))
                    if classification.action == 'bin':
                        # bin the variant
                        binFile.write('%s\n' % key)
                        for readCount,depth in info.counts:
                            binFile.write('    <vars/coverage: %d/%d>\n' % (readCount,depth))
                    elif classification.action == 'keep':
                        # keep the variant
                        csvWriter.writerow(info.inputRow)

if __name__ == "__main__":
   main()
