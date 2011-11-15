#!/bin/env python

'''
Variant annotation script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Ng.
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
    --variants=<variant list as TSV file>
    --annotations=<output TSV file with annotations added>
    reads1.bam reads2.bam ...""") % sys.argv[0]

longOptionsFlags = ["help", "variants=", "annotations="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.variants = None
        self.annotations = None
    def check(self):
        return all([self.variants, self.annotations])

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
        elif o == "--annotations":
            options.annotations = a
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
    # annotate the variants
    annotate(options, evidence)

def annotate(options, evidence):
    '''Annotate variants which appear in a sample of a family member'''
    annotateFilename = options.annotations
    with open(annotateFilename,'w') as annotateFile:
        inFamily = []
        notInFamily = []
        # sort the variants by coordinate
        # XXX does this include the variant multiple times?
        for key,info in sortByCoord(evidence):
            for readCount,depth in info.counts:
                if readCount > 0:
                    inFamily.append(info.inputRow)
                else:
                    notInFamily.append(info.inputRow)
        csvWriter = csv.writer(annotateFile, delimiter='\t', quotechar='|')
        for row in inFamily:
             csvWriter.writerow(row + ['IN RELATIVE'])
        for row in notInFamily:
             csvWriter.writerow(row + ['NOT IN RELATIVE'])
    print('Annotated reads saved into file: %s' % annotateFilename)

if __name__ == '__main__':
    main()
