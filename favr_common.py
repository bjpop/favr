'''
Common utilities for FAVR programs.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Ng.
'''

import os
import pysam
import sys

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

def getEvidence(variantList, bamFilenames):
    evidence = initEvidence(variantList)
    # Iterate over sample BAM files.
    for bamFile in bamFilenames:
       with pysam.Samfile(bamFile, "rb") as bam:
           # Count how many samples have each particular variant.
           countVariants(evidence, variantList, bam)
    return evidence

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
