#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import random

# External imports
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def createRandomSeqRecord(seqPair, windowSize, limits=None, isFasta=True):
    seq  = ""

    if (isFasta):
        ## Get a random sequence that doesn't contain unexpected characters
        while (len(seq) == 0 or seq.count("N") != 0):
            f   = _getRandomPosition
            (randStartPos, endPos) = f(seqPair, windowSize, limits)
            seq = seqPair[1][randStartPos:endPos]

        seqId  = "{}:{}:{}".format(seqPair[0], str(randStartPos), str((endPos - 1)))
        seqRec = SeqRecord(Seq(seq), id=seqId, description='')

    else:
        f       = createRandomSeqRecord
        seqRec  = f(seqPair, windowSize, limits=limits, isFasta=True)
        qValues = [91] * len(seqRec.seq)    ## Random number for a quality
        seqRec.letter_annotations['phred_quality'] = qValues

    return seqRec

def createSpecificSeqRecord(seqPair, pos):
    (startPos, endPos) = pos
    seqId = "{}:{}:{}".format(seqPair[0], str(startPos), str((endPos - 1)))
    seq   = seqPair[1][startPos:endPos]
    return SeqRecord(Seq(seq), id=seqId, description='')

#------------------- Private Classes & Functions ------------#

def _getRandomPosition(seqPair, windowSize, limits=None):
    if (limits is not None):
        randStartPos = random.randint(limits[0], limits[1])

    else:
        idx          = len(seqPair[1]) - windowSize - 1
        randStartPos = random.randint(0, idx)

    endPos = randStartPos + windowSize
    return (randStartPos, endPos)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
