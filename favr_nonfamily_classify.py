'''
Classifier for the Variant filter script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Ng.
'''

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

# the argument variantInfo is a list of pairs:
#    (readCount, depth)
# one per sample, where
#    readCount is the number of reads in a sample which match the variant (at the same position)
#    depth is the number of reads in a sample covering the same position as the variant (aka coverage)
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
