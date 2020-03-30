#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
from operator import add

# External imports
from pyspark.sql.functions import broadcast

# Internal imports
from .constants import *
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getCounts(rdd, kmerLength, includeN):
    ## **********
    ## *** Not sure if this is the "best" approach for split counts.
    ## *** But seems to run the fastest!
    ## ***
    ## *** Using RDD ReduceByKey once to minimise shuffling can encounter
    ## *** a few errors (not sure why).
    ## ***
    ## *** I've tried:
    ## ***  * Using fewer lambda functions to reduce double serialisation.
    ## ***  * Repartitioning before reducing. However, this doesn't seem to
    ## ***    result in any major improvements since we're shuffling data twice.
    ## **********
    kmerRdd  = _getKmerRdd(rdd, kmerLength)
    kmerDf   = kmerRdd.toDF([ID_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME])
    kmerDf   = getValidRows(kmerDf)

    totalRdd = _getTotalRddWithN(rdd, kmerLength) if includeN else _getTotalRddWithoutN(kmerDf.rdd)
    totalDf  = totalRdd.toDF([ID_COL_NAME, TOTAL_COL_NAME])
    kmerDf   = _countsToPercentage(kmerDf, totalDf)
    kmerDf   = groupRows(kmerDf)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def _getKmerRdd(rdd, kmerLength):
    ## Get the Kmer counts for each record
    ## Counts are (K, V) pairs, where K = sequence and V = frequency
    ## (ID, Seq) => (ID, kmerSeq)
    ##           => ((ID, kmerSeq), 1)
    ##           => ((ID, kmerSeq), count)
    f = lambda x: getObservedSequences(x, kmerLength)
    kmerRdd = rdd.flatMapValues(f) \
                 .map(lambda x: (x, 1)) \
                 .reduceByKey(add) \
                 .map(lambda x: (*x[0], x[1]))
    ## We only count Kmers that occur at least once. Kmers that
    ## do not occur (zero counts) are ignored, but need to be
    ## added when we're analysing results. This will save us quite
    ## a bit of space.
    kmerRdd.cache()
    return kmerRdd

def _getTotalRddWithN(rdd, kmerLength):
    ## Get the total number of Kmer for each records
    ## (ID, Seq) => (ID, total)
    f = lambda x: getObservedTotal(x, kmerLength)
    totalRdd = rdd.mapValues(f)
    return totalRdd

def _getTotalRddWithoutN(rdd):
    ## Get the total number of Kmer for each records
    ## (ID, kmerSeq, count) => (ID, count)
    ##                      => (ID, total)
    f = lambda x: (x[0], x[2])
    totalRdd = rdd.map(f).reduceByKey(add)
    return totalRdd

def _countsToPercentage(kmerDf, totalDf):
    cond     = kmerDf.id == totalDf.id
    f        = getPercentage(COUNT_COL_NAME, TOTAL_COL_NAME)
    colNames = [kmerDf.id, KMER_COL_NAME, f.alias(PERCENTAGE_COL_NAME)]
    kmerDf   = kmerDf.join(broadcast(totalDf), cond, 'left') \
                     .select(colNames)
    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
