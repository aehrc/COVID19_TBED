#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import gzip
import os
import re

# External imports
from Bio import SeqIO

# Internal imports
from .constants import *
from .common import isZFile
from .common import createDirIfNone
from .common import removeFileIfExists

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*filepaths):
    seqRecs  = (_readFastX(f) for f in filepaths)
    seqPairs = [(r.id, str(r.seq)) for sublist in seqRecs for r in sublist]
    return seqPairs

def write(filepath, seqRecords, fastXtype):
    ## Create DIR if it doesnt exist
    outputDir = os.path.dirname(filepath)
    createDirIfNone(outputDir)
    removeFileIfExists(filepath)
    SeqIO.write(seqRecords, filepath, fastXtype)

#------------------- Private Classes & Functions ------------#

def _readFastX(filepath):

    """
    Description:
        Reads a FASTX file/s and generates a list of sequence records.

    Args:
        filepath (str):
            Filepath string.

    Returns:
        seqRecs (list<Bio.SeqRecord.SeqRecord>):
            List of sequence records.
    """

    if (isZFile(filepath)):
        with gzip.open(filepath, 'rt') as filehandle:
            tmp       = os.path.splitext(filepath)[0]
            fastxType = _getType(tmp)
            seqRecs   = [seqRec.upper()
                         for seqRec in SeqIO.parse(filehandle, fastxType)]

            ## Ensembl .fa.gz files contain both
            ## chromosome and scaffold / haplotype sequences
            ## Probably better to ignore scaffold / haplotyle sequences
            ## If we want to consider them, just remove the condition
            if (fastxType == FASTA):
                seqRecs = _filterFasta(seqRecs)

    else:
        fastxType = _getType(filepath)
        seqRecs   = [seqRec.upper()
                     for seqRec in SeqIO.parse(filepath, fastxType)]

    return seqRecs

def _getType(filepath):
    if (filepath.endswith('.fa') \
        or filepath.endswith('.fasta')):
        return FASTA

    elif (filepath.endswith('.fq') \
          or filepath.endswith('.fastq')):
        return FASTQ

    else:
        raise NotImplementedError("Unknown fastX file")

def _filterFasta(seqRecs):
    f = lambda x: re.search(" REF$", x.description) \
              and re.search(":chromosome chromosome:", x.description)
    newSeqRecs = list(filter(f, seqRecs))
    return newSeqRecs

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
