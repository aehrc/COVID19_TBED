#!/bin/python

#------------------- Description & Notes --------------------#

'''
Description:
    Given a list of FASTA/FASTQ or SAM/BAM files, output (to file) the
    frequencies for each Kmer sequence of length K in a sequence/s.
    Frequencies are calculated either:
        * Across all sequence records (combined).
        * For each sequence record (split).
        * Along the length of a given sequence record (sequential).

Args:
    fastXFiles (filepath):
        List containing the filepath of each file. Files can be
        compressed (.gz) and should contain at least one
        sequence record (FASTA/FASTQ or SAM/BAM). Our processing limit
        seems to be around 3.5 million (3,500,000) sequence records.

    kmerLength (int):
        Length of Kmer sequences. Must be a positive integer. Ideally,
        this should be <= 13 since the total number of Kmer sequences
        exponentially increases (4^K).
            * 4^13 Kmers =    67,108,864  ## Possible
            * 4^14 Kmers =   268,435,456  ## Sometimes possible
            * 4^15 Kmers = 1,073,741,824  ## Probably not possible

Returns:
    oFile (dir):
        Directory containing a list of files. Each file is compressed
        in Parquet format and contains the frequencies for each
        Kmer sequence of length K in a sequence/s.
'''

'''
Design considerations
    Two general approaches:
        * Create table of all possible Kmers, iterate through each Kmer and
          count its frequency in the sequence.
        * Sliding window across the sequence and count the occurrence of
          each Kmer.

    Problems encountered:
        * Table of Kmers does not scale well (4^K Kmers).
        * Sliding window encounters unexpected characters (N, M, etc...)
        * Frequency tables can be big (4^K rows).
'''

'''
Tool evaluation:
    Jellyfish
        * Seems to be the literature standard.
        * Fast, inexact (approximate) Kmer counting.

    KmerCounter (implemented by Brendan's brother)
        * Fast, exact Kmer counting.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import math
import sys

sys.path.append('/home/wil9cq/scratch1/genomicsignatures/scripts/dev/')
sys.path.append('/home/wil9cq/scratch1/genomicsignatures/scripts/')

# External imports

# Internal imports
from src.util import params
from src.util import spark
from src import bio
from src import io
from src import kmer

#------------------- Constants ------------------------------#

## Maximum number of output files
MAX_N_FILES = 16

## Maximum number of FASTA/FASTQ entries per Partition
CHUNK_SIZE = 1000000

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

def combineKmerFrequencies(fastXFiles, kmerLength, includeN, oFile):
    seqPairs = io.fastx.read(*fastXFiles)

    print("Combining counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Number of partitions should be proportional to the
            ## number of sequence records we are querying
            eParts    = math.floor(len(seqPairs) / CHUNK_SIZE)
            nParts    = sc.defaultParallelism + (eParts * spark.N_PARTS_MULTI)
            seqPairRdd = sc.parallelize(seqPairs, nParts)

            ## Get a table containing the Kmer percentages across all records
            kmerDf = kmer.combined.getCounts(seqPairRdd, kmerLength, includeN)

            ## Prepare the table for writing
            nFiles = 1
            if (kmerLength > 10):
                ## Start writing to more files
                ## once we start reaching > 10-mers
                r = kmerLength % 10
                nFiles = 2 ** math.ceil((r / 2))
            nParts = min(MAX_N_FILES, nFiles)
            kmerDf = kmer.repartition(kmerDf, nParts)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def splitKmerFrequencies(fastXFiles, kmerLength, includeN, oFile):
    seqPairs = io.fastx.read(*fastXFiles)

    print("Splitting counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Number of partitions should be proportional to the
            ## number of sequence records we are querying
            eParts    = math.floor(len(seqPairs) / CHUNK_SIZE)
            nParts    = sc.defaultParallelism + (eParts * spark.N_PARTS_MULTI)
            nParts    = nParts * kmerLength
            seqPairRdd = sc.parallelize(seqPairs, nParts)

            ## Get a table containing the kmer percentages across all records
            kmerDf = kmer.split.getCounts(seqPairRdd, kmerLength, includeN)

            ## Prepare the table for writing
            ## Number of output files should be proportional
            ## to the number of sequence records we are querying
            ## and the Kmer length
            nFiles = math.ceil(len(seqPairs) / CHUNK_SIZE)
            nFiles = nFiles + math.floor((kmerLength / 5))
            nParts = min(MAX_N_FILES, nFiles)
            kmerDf = kmer.repartition(kmerDf, nParts)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def sequentialKmerFrequencies(fastXFiles,
    fastxId, nWindows, nSteps,
    kmerLength, includeN, oFile):

    seqPairs = io.fastx.read(*fastXFiles)

    ## Look for a sequence record if there is a
    ## specific one we want to query
    if (fastxId is not None):
        seqPair = bio.getSeqPair(fastxId, seqPairs)
        if (seqPair is None):
            raise IndexError("Invalid FASTX ID")
        seqPairs = [seqPair]

    print("Sequential counts")
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            ## Number of partitions should be proportional to the
            ## number of sequence records we are querying
            eParts     = math.floor(len(seqPairs) / CHUNK_SIZE)
            nParts     = sc.defaultParallelism + (eParts * spark.N_PARTS_MULTI)
            nParts     = nParts * kmerLength
            seqPairRdd = sc.parallelize(seqPairs, nParts)

            ## Get a table containing the kmer percentages across all records
            seqPairRdd = kmer.sequential.getSeqPairs(seqPairRdd, nWindows, nSteps)
            kmerDf     = kmer.split.getCounts(seqPairRdd, kmerLength, includeN)

            ## Prepare the table for writing
            ## Number of output files should be proportional
            ## to the number of sequence records we are querying
            ## and the Kmer length
            nFiles = math.ceil(len(seqPairs) / CHUNK_SIZE)
            nFiles = nFiles + math.floor((kmerLength / 5))
            nParts = min(MAX_N_FILES, nFiles)
            kmerDf = kmer.repartition(kmerDf, nParts)

            ## Write the table to disk
            print("Writing output")
            io.kmer.write(oFile, kmerDf)
            ## Instead of outputing, we could extend the pipeline here
            ## so that it goes directly to analysing Kmer frequencies

def main():
    argParser = params.CalculateKmerArgParser()
    argParser.parse()

    if (argParser.cmd == 'combined'):
        argParser.printArgs()

        combineKmerFrequencies(argParser.iFiles,
            argParser.kmerLength, argParser.includeN, argParser.oFile)

    elif (argParser.cmd == 'split'):
        argParser.printArgs()

        splitKmerFrequencies(argParser.iFiles,
            argParser.kmerLength, argParser.includeN, argParser.oFile)

    elif (argParser.cmd == 'sequential'):
        argParser.parseSequentialArgs()
        argParser.printArgs()

        sequentialKmerFrequencies(argParser.iFiles,
            argParser.fastxId, argParser.nWindows, argParser.nSteps,
            argParser.kmerLength, argParser.includeN, argParser.oFile)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
