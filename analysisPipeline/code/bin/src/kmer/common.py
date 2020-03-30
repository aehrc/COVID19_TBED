#!/bin/python

#------------------- Description & Notes --------------------#

'''
Kmer frequencies only need to be calculated in
the forward or reverse direction. We don't need to
do both, since they basically contain the same information.
In fact, the forward proportion is roughly equivalent to the
forward & reverse proportion
'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
from itertools import product
from pyspark.sql.functions import *
from pyspark.sql.types import *

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def repartition(kmerDf, nParts):
    if (kmerDf.rdd.getNumPartitions() > nParts):
        kmerDf = kmerDf.coalesce(nParts)

    elif (kmerDf.rdd.getNumPartitions() < nParts):
        kmerDf = kmerDf.repartition(nParts)

    return kmerDf

def groupRows(kmerDf):
    ## (ID, Kmer, Percent, *) => (ID, Kmer -> Percent, *)
    ##                        => (ID, [Kmer -> Percent], *)
    gCols = kmerDf.schema.names
    gCols.remove(KMER_COL_NAME)
    gCols.remove(PERCENTAGE_COL_NAME)

    schema = struct([KMER_COL_NAME, PERCENTAGE_COL_NAME])
    f      = map_from_entries(collect_list(schema))
    kmerDf = kmerDf.groupBy(*gCols).agg(f.alias(KMER_COL_NAME))
    return kmerDf

def ungroupRows(kmerDf):
    gCols = kmerDf.schema.names
    gCols.remove(KMER_COL_NAME)

    ## (ID, [Kmer -> Percent], *) => (ID, Kmer, Percent, *)
    f        = explode(KMER_COL_NAME)
    colNames = f.alias(KMER_COL_NAME, PERCENTAGE_COL_NAME)
    kmerDf   = kmerDf.select(*gCols, colNames)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def getExpectedSequences(kmerLength):

    """
    Description:
        Generates a list of all possible Kmer sequences of length K.

    Args:
        kmerLength (int):
            Length of Kmers. Must be a positive integer.

    Returns:
        kmerSeqs (generator):
            List of all possible Kmer sequences.

    Raises:
        ValueError:
            If kmerLength is not a positive integer.
    """

    f = product(NUCLEOTIDES.keys(), repeat=kmerLength)
    return (''.join(c) for c in f)

def getExpectedTotal(kmerLength):
    return (4 ** kmerLength)

def getObservedSequences(seq, kmerLength):

    """
    Description:
        Generates a list of Kmer sequences of length K in a sequence.

    Args:
        seq (str):
            Sequence to be examined.

        kmerLength (int):
            Length of Kmers. Must be a positive integer.

    Returns:
        kmerSeqs (list):
            List of str containing the Kmer sequences of length K
            in the sequence record
    """

    obsTotal = getObservedTotal(seq, kmerLength)
    return [seq[i:i + kmerLength] for i in range(0, obsTotal)]

def getObservedTotal(seq, kmerLength):
    return len(seq) - kmerLength + 1

@udf(returnType=FloatType())
def getPercentage(counts, total):
    return ((counts / total) * 100)

def getValidRows(kmerDf):

    """
    Description:
        Generates a Spark DataFrame containing only valid rows. Rows are
        considered invalid if the Kmer sequence contains unexpected characters.

    Args:
        kmerDf (pyspark.DataFrame):
            Spark DataFrame containing the counts / percentages
            of Kmer sequences.

    Returns:
        kmerDf (pyspark.DataFrame):
            Filtered Spark DataFrame containing the counts / percentages
            of Kmer sequences.
    """

    ## Ensure that we don't have Kmers containing ambiguous bases
    ## i.e., N's, R's, S's, etc...
    nStr = ''.join(NUCLEOTIDES.keys())
    p    = '^[' + nStr + ']+$'
    return kmerDf.filter(kmerDf.kmer.rlike(p))

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
