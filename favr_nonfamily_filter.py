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

# print a usage message
def usage():
    print("""Usage: %s
    [-h | --help]
    --family=[True|False]
    --variants=<variant list as CSV file>
    reads1.bam reads2.bam ...""") % sys.argv[0]

longOptionsFlags = ["help", "family=", "variants="]
shortOptionsFlags = "h"

# A place to store command line arguments.
class Options(object):
    def __init__(self):
        self.family = None
        self.variants = None
    def check(self):
        return self.family != None and self.variants != None

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptionsFlags, longOptionsFlags)
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)
    options = Options()
    for o, a in opts:
        if o == "--family" and a in ['True','False']:
            options.family = eval(a)
        elif o == "--variants":
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
    # Initialise all the evidence to zero.
    # The evidence records the number of samples which contain a particular variant.
    evidence = initEvidence(variantList)
    # Iterate over sample BAM files.
    for bamFile in bamFilenames:
       with pysam.Samfile(bamFile, "rb") as bam:
           # Count how many samples have each particular variant.
           countVariants(evidence, variantList, bam)
    # either annotate the variants or filter them, depending on whether
    # we are dealing with familial samples or not.
    if options.family:
        annotate(evidence)
    else:
        filter(evidence)

def annotate(evidence):
    '''Annotate variants which appear in a sample of a family member'''
    annotateFilename = makeSafeFilename('annotated_variants')
    with open(annotateFilename,'w') as annotateFile:
        inFamily = []
        notInFamily = []
        # sort the variants by coordinate
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

class EvidenceInfo(object):
    def __init__(self, inputRow, counts):
        self.inputRow = inputRow
        self.counts = counts

def initEvidence(variantList):
    '''Initialise the frequency counter for each variant to be zero.'''
    evidence = {}
    for variant in variantList:
        info = parseVariantRow(variant)
        if info:
            evidence[info.id] = EvidenceInfo(inputRow = info.inputRow, counts = [])
    return evidence

def showEvidence(evidence):
    '''Print out the frequency counter for each variant.'''
    for chrPos,info in evidence.items():
        print("%s %s" % (info.inputRow,str(info.counts)))

def countVariants(evidence, variantList, bam):
    '''For each variant in the list, check if it is evident in this particular sample BAM.'''
    for variant in variantList:
        info = parseVariantRow(variant)
        # only process valid rows in the variant CSV file
        if info:
            # get the pileup information for this particular variant coordinates
            # the pileup tells us what base was called in each of the reads at
            # this particular coordinate
            pileupCol = lookupPileup(bam, info.chromosome, info.position)
            sameAsVariant = 0
            coverage = 0
            unknown_size_reads = 0 # number of variants on a read of unknown size
            if pileupCol:
                pos,coverage,reads = pileupCol
                # Count the number of reads that are the same as the variant
                for pileupread in reads:
                    readBase = pileupread.alignment.seq[pileupread.qpos]
                    # check if the sample base (at the same position) is the same as the variant base
                    # but skip reads which are marked as "is_del", deletions
                    if not pileupread.is_del and readBase == info.variantBase:
                        sameAsVariant += 1
            evidence[info.id].counts.append((sameAsVariant,coverage))

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
    if len(row) >= 1:
        coordinates = row[0].split(',')
        if len(coordinates) >= 4:
            fromTo = parsePolymorphism(coordinates[3])
            if fromTo:
               chrName = "chr" + coordinates[0]
               return VariantInfo(id = "%s:%s" % (chrName,coordinates[1]),
                                  chromosome = chrName,
                                  position = int(coordinates[1]), # XXX should really check it is all digits
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
