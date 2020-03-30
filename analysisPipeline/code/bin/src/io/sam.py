#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import gzip
import os

# External imports
import pysam

# Internal imports
from .constants import *
from .common import isZFile

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*filepaths):
    alignRecs  = [_readXam(f) for f in filepaths]
    alignRecs  = [(r.query_name, r.query_sequence)
                  for sublist in alignRecs for r in sublist]
    return alignRecs

#------------------- Private Classes & Functions ------------#

def _readXam(filepath):
    if (isZFile(filepath)):
        with gzip.open(filepath, 'rt') as filehandle:
            tmp     = os.path.splitext(filepath)[0]
            xamType = _getXamType(filepath)

            f       = pysam.AlignmentFile(filepath, 'r' + xamType)
            alignRecs = [alignRec
                         for alignRec in f.fetch(until_eof=True)]

    else:
        xamType = _getXamType(filepath)
        f = pysam.AlignmentFile(filepath, 'r' + xamType)
        alignRecs = [alignRec
                     for alignRec in f.fetch(until_eof=True)]

    return alignRecs

def _getXamType(filepath):
    if (filepath.endswith('.bam')):
        return BAM

    elif (filepath.endswith('.sam')):
        return SAM

    else:
        raise NotImplementedError("Unknown Xam file")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
