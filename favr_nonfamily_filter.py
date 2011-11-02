#!/bin/env python

'''
Variant filter script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Ng.

Reads a list of variants from a CSV file for a sample and compares
them to the sequence reads in BAM files for other samples.

Revision history:

2 May 2011.    Initial incomplete version. Can read CSV file and BAMS and
               generate a table of which variants from the original sample
               are found in the other samples.

30 May 2011.   Filter out reads which are considered 'is_del'. These appear
               to come from deletions, and don't belong in the pileup for
               a given position.

31 May 2011.   Change the supported input format from a custom CSV file to
               a SIFT (csv) file.

25 Aug 2011.   Changed the binning rule to bin 35-length only reads, in
               addition to the previous binning rule.

29 Aug 2011.   Changed back to not bin 35-length only reads. This is now
               done in another tool.

19 Sep 2011.   Added classification for family samples.

20 Sep 2011.   Allowed the input to be tab separated, with coordinates comma separated.

'''

import os
import pysam
import sys
import csv
import yaml
import getopt
from favr_common import (getEvidence, makeSafeFilename, sortByCoord)

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --variants=<variant list as CSV file>
    reads1.bam reads2.bam ...""") % sys.argv[0]

longOptionsFlags = ["help", "variants="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.variants = None
    def check(self):
        return self.variants != None

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
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
    if not options.check():
        print('Incorrect arguments')
        usage()
        exit(2)
    bamFilenames = args
    # Read the rows of the variants CSV file into a list.
    with open(options.variants) as variants:
        variantList = list(csv.reader(variants, delimiter='\t', quotechar='|'))
    # compute the presence/absence of each variant in the bam files
    evidence = getEvidence(variantList, bamFilenames)
    # filter the variants
    filter(evidence)

def filter(evidence):
    '''Decide which variants to keep and which to bin.'''
    binFilename = makeSafeFilename('rare_binfile')
    keepFilename = makeSafeFilename('rare_keepfile')
    logFilename = makeSafeFilename('rare_logfile')
    with open (logFilename,'w') as logFile:
        with open(binFilename,'w') as binFile:
            with open(keepFilename,'wb') as keepFile:
                csvWriter = csv.writer(keepFile, delimiter=',', quotechar='|')
                # sort the variants by coordinate
                for key,info in sortByCoord(evidence):
                    classification = classify(info.counts)
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
    print('Kept reads saved into file: %s' % keepFilename)
    print('Binned reads saved into file: %s' % binFilename)
    print('Log record saved into file: %s' % logFilename)

# The classify function decides if a particular variant from the
# variant list should be kept or binned (discarded). It returns
# a classification of the variant which can be used to
# determine if it should be kept or discarded and a reason
# why.
#
# The decision to bin or keep is based on how many times
# we see the same variant in other samples.
#
# If you want to use a different binning criteria, then
# rewrite the classify function. The one provided below is just
# an example.

readCountThreshold = 1 # how many reads of the variant do we need to see in a single sample?
samplesPercent = 30    # what percentage of the total samples need to pass the above threashold?

class Classify(object):
    def __init__(self, action, reason):
        self.action = action # 'bin' or 'keep'
        self.reason = reason # some text explaining why

def classify(variantInfo):
    totalSamples = 0     # number of sample files
    binableSamples = 0   # number of samples which are considered binable
    totalSameAsVariants = 0    # total number of reads which had a base the same as the variant in the same position
    # we ignore the depth in this particular example
    for readCount,_depth in variantInfo:
        totalSamples += 1
        if readCount >= readCountThreshold:
            binableSamples += 1
        totalSameAsVariants += readCount
    if totalSamples > 0:
        if (binableSamples * 100 / totalSamples) >= samplesPercent:
            return Classify('bin', binThresholdMessage % (binableSamples, totalSamples, samplesPercent))
        else:
            return Classify('keep', keepMessage % (binableSamples, totalSamples, samplesPercent))
    else:
        return Classify('bin', binZeroSamplesMessage)

binThresholdMessage = '(binableSamples(=%d) * 100 / totalSamples(=%d)) >= samplesPercent(=%d)'
binZeroSamplesMessage = 'there were zero samples to compare with'
keepMessage = '(binableSamples(=%d) * 100 / totalSamples(=%d)) < samplesPercent(=%d)'

if __name__ == '__main__':
    main()
