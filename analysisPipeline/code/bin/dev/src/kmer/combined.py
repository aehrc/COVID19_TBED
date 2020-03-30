#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
from operator import add
import random

# External imports
from pyspark.sql.functions import *

# Internal imports
from .constants import *
from .common import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getCounts(rdd, kmerLength, includeN):
    ## **********
    ## *** This is (probably) the "best" approach for combined counts.
    ## *** RDD ReduceByKey runs better than DataFrame GroupBy.
    ## *** So DO NOT change!!!
    ## ***
    ## *** Using RDD ReduceByKey once to minimise shuffling can encounter
    ## *** a few errors (not sure why).
    ## ***
    ## *** My understanding is that DataFrames GroupBy is faster, but
    ## *** I don't understand why this isn't the case. I've tried
    ## ***  * Explicitly defining the schema
    ## ***  * Repartitioning the DataFrame before grouping
    ## **********
    kmerRdd = _getKmerRdd(rdd, kmerLength)
    kmerDf  = kmerRdd.toDF([ID_COL_NAME, KMER_COL_NAME, COUNT_COL_NAME])
    kmerDf  = getValidRows(kmerDf)

    total   = _getTotalWithN(rdd, kmerLength) if includeN else _getTotalWithoutN(kmerDf.rdd)
    kmerDf  = _countsToPercentage(kmerDf, total)
    kmerDf.show()

    ## Writing large Kmers (i.e., > 10) to file encounters
    ## a few warnings. So, we'll change the format of the
    ## table so that the columns are divided across multiple
    ## rows. This should make it easier for us to read and
    ## write, but need to make sure we combine them when
    ## we're analysing results.
    kmerDf = _pivot(kmerDf, kmerLength)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def _getKmerRdd(rdd, kmerLength):
    ## Get the Kmer counts across all records
    ## Counts are (K, V) pairs, where K = sequence and V = frequency
    ## (ID, Seq) => (ID, kmerSeq)
    ##           => (kmerSeq, 1)
    ##           => (kmerSeq, count)
    f = lambda x: getObservedSequences(x, kmerLength)
    kmerRdd = rdd.flatMapValues(f) \
                 .map(lambda x: (x[1], 1)) \
                 .reduceByKey(add) \
                 .map(lambda x: (0, *x))
    kmerRdd.cache()
    ## We only count Kmers that occur at least once. Kmers that
    ## do not occur (zero counts) are ignored, but need to be
    ## added when we're analysing results. This will save us quite
    ## a bit of space.
    return kmerRdd

def _getTotalWithN(rdd, kmerLength):
    ## Get the total number of Kmers across all records
    ## (ID, Seq) => (ID, count)
    ##           => total
    f = lambda x: getObservedTotal(x[1], kmerLength)
    total = rdd.map(f).reduce(add)
    return total

def _getTotalWithoutN(rdd):
    ## Get the total number of Kmers across all records
    ## (ID, kmerSeq, count) => count
    ##                      => total
    f = lambda x: x[2]
    total = rdd.map(f).reduce(add)
    return total

def _countsToPercentage(kmerDf, total):
    f        = getPercentage(COUNT_COL_NAME, lit(total))
    colNames = [ID_COL_NAME, KMER_COL_NAME, f.alias(PERCENTAGE_COL_NAME)]
    kmerDf   = kmerDf.select(colNames)
    return kmerDf

def _pivot(kmerDf, kmerLength):
    ## **********
    ## *** Select is suppose to be faster than withColumn
    ## **********

    ## Replace the ID column to something
    ## that can be used for grouping
    ## (ID, Kmer, Percent) => (ID, [Kmer -> Percent])
    f        = round(rand()*(kmerLength-1))
    colNames = [f.alias(ID_COL_NAME), KMER_COL_NAME, PERCENTAGE_COL_NAME]
    kmerDf   = kmerDf.select(colNames)
    kmerDf   = groupRows(kmerDf)

    ## Replace the ID column to something
    ## more suitable for Identification purposes
    f        = lit(random.randint(0, 100000))
    colNames = [f.alias(ID_COL_NAME).cast(StringType()), KMER_COL_NAME]
    kmerDf   = kmerDf.select(colNames)
    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
