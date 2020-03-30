#!/bin/python

#------------------- Description & Notes --------------------#

'''
Description:
    Given a list of FASTA, output (to file) FASTA/FASTQ sequences
    of a given length, or a given feature (i.e., exon or intron).
    Sequences can be either:
        * Fixed-sized sequences (fixed).
        * Feature-sized sequences (feature).
        * Hybrid sequence made up of two FASTA sequences.

Args:
    fastXFiles (filepath):
        List containing the filepath of each file. Files can be
        compressed (.gz) and should contain at least one
        sequence record (FASTA).

    nSeqRecs (int):
        Number of sequence records. Must be a positive integer.

Returns:
    oFile (file):
        Single FASTA/FASTQ file containing multiple sequence records.
        The number of sequence records == nSeqRecs.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import sys
import random

# External imports

# Internal imports
from src.util import params
from src import bio
from src import io

#------------------- Constants ------------------------------#

MAX_ITER_COUNT  = 100000

FASTQ_TARGET_SPLIT = 0.2

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

def generateFixedFasta(fastxFiles, nSeqRecs, windowSize, fastxId, toFastq, oFile):
    seqPairs  = io.fastx.read(*fastxFiles)
    seqRecs   = []
    iterCount = 0   ## Ensure that we don't end up in an infinite loop

    print("Generating fixed-sized sequence records")
    while (len(seqRecs) < nSeqRecs):
        if (iterCount == MAX_ITER_COUNT):
            msg = "Only found {} out of {} features".format(len(seqRecs), nSeqRecs)
            print("Ran out of features. " + msg)
            break

        if (fastxId is None):
            seqPair = bio.getRandomSeqPair(seqPairs)

        else:
            seqPair = bio.getSeqPair(fastxId, seqPairs)
            if (seqPair is None):
                raise IndexError("Invalid FASTX ID")

        if (not toFastq):
            seqRec = bio.fixed.createRandomSeqRecord(seqPair, windowSize)

            ## If we encounter a duplicate record, try again
            f = lambda x: x.id == seqRec.id
            if (len(list(filter(f, seqRecs))) != 0):
                iterCount += 1

            print(seqRec.id)
            seqRecs.append(seqRec)

        else:
            ## Check whether we are generating reads for a hybrid sequence
            if (bio.hybrid.isValid(seqPair)):
                print(seqPair[0])
                ## Ensure that at least X% of the FASTQ reads align
                ## (or partially align) to the source
                nSRecs = int(nSeqRecs * FASTQ_TARGET_SPLIT)
                nTRecs = int(nSeqRecs - nSRecs)

                ## Get the position of the source
                ## within in the hybrid sequence
                limits = bio.hybrid.getSourcePositions(seqPair, windowSize)

                ## Generate reads. However, theres an issue whereby reads
                ## will only be generated for the first hybrid sequence
                ## we encounter
                f       = bio.fixed.createRandomSeqRecord
                sRecs   = [f(seqPair, windowSize, limits, isFasta=False)
                           for n in range(0, nSRecs)]
                tRecs   = [f(seqPair, windowSize, isFasta=False)
                           for n in range(0, nTRecs)]
                seqRecs = sRecs + tRecs

                ## Shuffle the reads so that we don't
                ## know which ones are source or target
                random.shuffle(seqRecs)

            else:
                f      = bio.fixed.createRandomSeqRecord
                seqRec = f(seqPair, windowSize, isFasta=False)
                print(seqRec.id)
                seqRecs.append(seqRec)

    if (not toFastq):
        io.fastx.write(oFile, seqRecs, io.FASTA)
    else:
        io.fastx.write(oFile, seqRecs, io.FASTQ)

def generateFeatureFasta(fastxFiles, nSeqRecs, annoFiles, featureId, oFile):
    seqPairs   = io.fastx.read(*fastxFiles)
    featureDBs = io.gff.read(*annoFiles, asPdf=False)
    seqRecs    = []
    iterCount  = 0   ## Ensure that we don't end up in an infinite loop

    ## Combine features across multiple files into a single featureDB
    featureDB = None
    for db in featureDBs:
        if (featureDB is None):
            featureDB = db
        else:
            featureDB.update(db)

    print("Generating FASTA sequences from features")
    while (len(seqRecs) < nSeqRecs):
        if (iterCount == MAX_ITER_COUNT):
            msg = "Only found {} out of {} features".format(len(seqRecs), nSeqRecs)
            print("Ran out of features. " + msg)
            break

        ## Get a seqRec that isn't already in seqRecords
        ## There is a possibility that we will run out of seqRecs
        featureRec = bio.feature.getRandomRecord(featureDB, featureId)
        fastxId    = featureRec.chrom

        seqPair    = bio.getSeqPair(fastxId, seqPairs)
        if (seqPair is None):
            continue

        seqRec = bio.feature.createSeqRecord(seqPair, featureRec)

        ## If we encounter a duplicate record, try again
        f = lambda x: x.id == seqRec.id
        if (len(list(filter(f, seqRecs))) != 0):
            iterCount += 1
            continue

        print(seqRec.id)
        seqRecs.append(seqRec)

    io.fastx.write(oFile, seqRecs, io.FASTA)

def generateHybridFasta(fastxFiles, nSeqRecs, sId, tId, oFile):
    seqPairs  = io.fastx.read(*fastxFiles)
    seqRecs   = []
    iterCount = 0   ## Ensure that we don't end up in an infinite loop

    print("Generating hybrid FASTA sequence")
    while (len(seqRecs) < nSeqRecs):
        if (iterCount == MAX_ITER_COUNT):
            msg = "Only found {} out of {} features".format(len(seqRecs), nSeqRecs)
            print("Ran out of features. " + msg)
            break

        sSeqPair = bio.getSeqPair(sId, seqPairs)
        tSeqPair = bio.getSeqPair(tId, seqPairs)
        if (sSeqPair is None and tSeqPair is None):
            raise IndexError("Invalid FASTX ID")

        seqRec = bio.hybrid.createRandomSeqRecord(sSeqPair, tSeqPair)

        ## Get a hybrid sequence that doesn't hybridise with itself
        if (sSeqPair[0] == tSeqPair[0]):
            iterCount += 1
            continue

        ## If we encounter a duplicate record, try another interval
        f = lambda x: x.id == seqRec.id
        if (len(list(filter(f, seqRecs))) != 0):
            iterCount += 1
            continue

        print(seqRec.id)
        seqRecs.append(seqRec)

    io.fastx.write(oFile, seqRecs, io.FASTA)

def main():
    argParser = params.GenerateFastxArgParser()
    argParser.parse()

    if (argParser.cmd == 'fixed'):
        argParser.parseFixedArgs()
        argParser.printArgs()

        generateFixedFasta(argParser.iFiles, argParser.nSeqRecs,
            argParser.windowSize, argParser.fastxId,
            argParser.toFastq, argParser.oFile)

    elif (argParser.cmd == 'feature'):
        argParser.parseFeatureArgs()
        argParser.printArgs()

        generateFeatureFasta(argParser.iFiles, argParser.nSeqRecs,
            argParser.annoFiles, argParser.featureId, argParser.oFile)

    elif (argParser.cmd == 'hybrid'):
        argParser.parseHybridArgs()
        argParser.printArgs()

        generateHybridFasta(argParser.iFiles, argParser.nSeqRecs,
            argParser.sId, argParser.tId, argParser.oFile)

    print("DONE")

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
