#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import random

# External imports

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getSeqPair(fastxId, seqPairs):
    f       = lambda x: x[0] == fastxId
    seqPair = list(filter(f, seqPairs))
    if (len(seqPair) == 0):
        return None

    return seqPair[0]

def getRandomSeqPair(seqPairs):
    randFaIdx = random.randint(0, len(seqPairs) - 1)
    seqPair   = [seqPairs[randFaIdx]]
    return seqPair[0]

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
