#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np

# Internal imports
from .. import bio

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getSeqPairs(seqPairRdd, nWindows, nSteps):
    ## Create a (ID, Seq) pair for each (S, E) pair
    ## (ID, Seq) => ((ID, Seq), SeqLen)
    ##           => ((ID, Seq), (SPos, EPos))
    ##           => (SeqRec)
    ##           => (ID, Seq)
    f = lambda x: (x, len(x[1]))
    g = lambda x: (getPositions(x, nWindows, nSteps))
    h = lambda x: bio.fixed.createSpecificSeqRecord(x[0], x[1])
    i = lambda x: (x.id, str(x.seq))
    seqPairRdd = seqPairRdd.map(f).flatMapValues(g) \
                           .map(h).map(i)

    return seqPairRdd

#------------------- Private Classes & Functions ------------#

def getPositions(seqLen, nWindows, nSteps):
    windowSize = seqLen / nWindows
    stepSize   = windowSize / nSteps
    nSeqPairs  = seqLen - windowSize + 1
    positions  = [(int(round(s)), int(round(s + windowSize)))
                  for s in np.arange(0, nSeqPairs, stepSize)]
    return positions

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
