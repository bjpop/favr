#!/bin/env python

'''
RefGene annotation script.

Authors: Bernie Pope, Danny Park, Tu Ng.

Reads a list of variants from a CSV file for a sample and compares
them to features in the refGene.txt database from UCSC.

Revision history:

21 Sep 2011.    Initial version.
24 Oct 2011.    Fixed coordinate issue.

'''

import os
import sys
import csv
import getopt

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --variants=<variant list as CSV file>
    --startslack=<distance from start of coding region>
    --spliceslack=<distance from exon start/end sites>
    --refGene=<refGene.txt file>
    --output=<output file name>""") % sys.argv[0]

longOptionsFlags = ["help", "variants=", "refGene=", "startslack=", "spliceslack=", "output="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.variants = None
        self.refGene = None
        self.spliceslack = None
        self.startslack = None
        self.output = None
    def check(self):
        return (self.refGene != None and
                self.variants != None and
                self.spliceslack != None and
                self.startslack != None and
                self.output != None)

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptionsFlags, longOptionsFlags)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    options = Options()
    for o, a in opts:
        if o == "--refGene":
            options.refGene = a
        elif o == "--variants":
            options.variants = a
        elif o == "--spliceslack":
            options.spliceslack = safeReadInt(a)
        elif o == "--startslack":
            options.startslack = safeReadInt(a)
        elif o == "--output":
            options.output = a
        elif o in ('-h', '--help'):
            usage()
            sys.exit(0)
    if not options.check():
        print('Incorrect arguments')
        usage()
        exit(2)
    refGene = readRefGene(options)
    #showRefGene(refGene)
    annotate(options, refGene)

def annotate(options, refGene):
    with open(options.output, 'w') as output:
        csvWriter = csv.writer(output, delimiter='\t', quotechar='|')
        # Read the rows of the variants CSV file into a list.
        with open(options.variants) as variants:
            for row in csv.reader(variants, delimiter='\t', quotechar='|'):
                if len(row) >= 1:
                    coords = row[0].split(',')
                    if len(coords) >= 4 and coords[1].isdigit():
                        chrName = "chr" + coords[0]
                        pos = safeReadInt(coords[1])
                        searchResult = search(chrName, pos, refGene)
                        if searchResult != None:
                            csvWriter.writerow(row + [searchResult.annotation])
                        else:
                            csvWriter.writerow(row)

def search(chr, pos, refGene):
    features = refGene.get(chr, [])
    for f in features:
        if pos >= f.lowerBound and pos <= f.upperBound:
            return f
    return None

def showRefGene(refGene):
    for chr, vals in refGene.items():
        print("%s" % chr)
        for v in vals:
            print("\t%s" % str(v))

class CodingRegionStart(object):
    def __init__(self, direction, slack, start, end):
        self.annotation = 'Within %d before coding region start' % slack
        if direction == '+':
            self.lowerBound = start - slack
            self.upperBound = start - 1 # don't include the start in the region
        else:
            self.upperBound = end + slack
            self.lowerBound = end + 1 # don't include the end in the region
    def __str__(self):
        return "CodingRegionStart %d %d" % (self.upperBound, self.lowerBound)

# all kinds of exon boundaries are sub-types of ExonBoundary
# they differ in annotation and in the way they are printed
class ExonBoundary(object):
    def __init__(self, slack, coord, type):
        if type == 'start':
            self.lowerBound = coord - slack
            self.upperBound = coord + (slack - 1)
        else: # type == 'end'
            self.lowerBound = coord - (slack - 1)
            self.upperBound = coord + slack

class CodingExonBoundary(ExonBoundary):
    def __init__(self, slack, coord, type):
        super(CodingExonBoundary, self).__init__(slack, coord, type)
        self.annotation = 'Within +/- %d of coding exon %s boundary' % (slack, type)
    def __str__(self):
        return "CodingExonBoundary %d %d" % (self.upperBound, self.lowerBound)

class NonCodingExonBoundary(ExonBoundary):
    def __init__(self, slack, coord, type):
        super(NonCodingExonBoundary, self).__init__(slack, coord, type)
        self.annotation = 'Within +/- %d of NON-coding exon %s boundary' % (slack, type)
    def __str__(self):
        return "NonCodingExonBoundary %d %d" % (self.upperBound, self.lowerBound)

class PartialCodingExonBoundary(ExonBoundary):
    def __init__(self, slack, coord, type):
        super(PartialCodingExonBoundary, self).__init__(slack, coord, type)
        self.annotation = 'Within +/- %d of PARTIAL-coding exon %s boundary' % (slack, type)
    def __str__(self):
        return "NonCodingExonBoundary %d %d" % (self.upperBound, self.lowerBound)

def readRefGene(options):
    with open(options.refGene) as refs:
       refGene = {}
       for row in csv.reader(refs, delimiter='\t'):
           if len(row) >= 11:
               chr = row[2]
               if chr not in refGene:
                   refGene[chr] = []
               direction = row[3]
               transcriptStart = safeReadInt(row[4]) + 1
               transcriptEnd = safeReadInt(row[5])
               # All coordinates for all fields in the refGene file
               # are 0-based start and 1-based end, including the transcript, CDS and exons.
               # see: http://genome.ucsc.edu/FAQ/FAQtracks.html#tracks1
               codingRegionStart = safeReadInt(row[6]) + 1 # add one to fix up zero-based start coordinate
               codingRegionEnd = safeReadInt(row[7])
               # For non-coding genes (who, by definition, have a coding region size of
               # 0), cdsStart will always equal cdsEnd in the genePred format.
               # As a convention to help with standardization, we
               # have made the cdsStart equal the txtStart for non-coding genes.

               # check if this is a coding gene
               if codingRegionStart < (codingRegionEnd + 1):
                   codingStart = CodingRegionStart(direction, options.startslack, codingRegionStart, codingRegionEnd)
                   refGene[chr].append(codingStart)
               exonStarts = map(lambda x: safeReadInt(x) + 1, row[9].rstrip(',').split(',')) # add one to fix up zero-based start coordinate
               exonEnds = map(lambda x: safeReadInt(x), row[10].rstrip(',').split(','))
               # we classify the start and ends of exons depending on whether they fall
               # inside or outside the coding region, including the partial case where
               # one end is inside and the other end is outside.
               for start,end in zip(exonStarts, exonEnds):
                   isStartCoding = isCoding(codingRegionStart, codingRegionEnd, start)
                   isEndCoding = isCoding(codingRegionStart, codingRegionEnd, end)
                   if isStartCoding and isEndCoding:
                       refGene[chr].append(CodingExonBoundary(options.spliceslack, start, 'start'))
                       refGene[chr].append(CodingExonBoundary(options.spliceslack, end, 'end'))
                   elif not isStartCoding and not isEndCoding:
                       refGene[chr].append(NonCodingExonBoundary(options.spliceslack, start, 'start'))
                       refGene[chr].append(NonCodingExonBoundary(options.spliceslack, end, 'end'))
                   elif isStartCoding and not isEndCoding:
                       refGene[chr].append(CodingExonBoundary(options.spliceslack, start, 'start'))
                       refGene[chr].append(PartialCodingExonBoundary(options.spliceslack, end, 'end'))
                   elif not isStartCoding and isEndCoding:
                       refGene[chr].append(PartialCodingExonBoundary(options.spliceslack, start, 'start'))
                       refGene[chr].append(CodingExonBoundary(options.spliceslack, end, 'end'))
    return refGene

def isCoding(codingStart, codingEnd, coord):
    return coord >= codingStart and coord <= codingEnd

def safeReadInt(str):
    if str.isdigit():
        return int(str)
    else:
        raise Exception, 'not an integer: ' + str

# sort the variants by coordinate.
def sortByCoord(evidence):
    return sorted(evidence.items(), cmp=compareCoord)

# compare two chromosome coordinates for ordering.
def compareCoord(coord1, coord2):
    chr1,pos1 = coord1[0].split(':')
    chr2,pos2 = coord2[0].split(':')
    if chr1 == chr2:
        return cmp(int(pos1), int(pos2))
    else:
        return compareChrCode(chr1[3:], chr2[3:])

def compareChrCode(code1, code2):
    isNumerical1 = code1.isdigit()
    isNumerical2 = code2.isdigit()
    if isNumerical1:
        if isNumerical2:
            return cmp(int(code1), int(code2))
        else:
            return -1
    else:
        return cmp(code1, code2)

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
