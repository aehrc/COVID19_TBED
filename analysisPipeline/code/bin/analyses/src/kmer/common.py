#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd

# Internal imports
from dev.src import kmer

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def getDirs(dirs, **kwargs):
    dataStr  = getDataStr(**kwargs)
    kmerDirs = [str(d) + "_" + dataStr for d in dirs]
    return kmerDirs

def rotateAndSplit(kmerDf):
    kmerDf              = rotate(kmerDf)
    (kmerId, kmerCount) = split(kmerDf)
    return (kmerId, kmerCount)

def filterByFilename(kmerId, kmerCount, *dirPaths):
    if (len(dirPaths) != 0):
        filenames = [d.name for d in dirPaths]
        pattern   = '|'.join(filenames)

        cond      = ~kmerId[kmer.FILE_COL_NAME].str.contains(pattern)
        kmerId    = kmerId.copy()[cond]
        kmerCount = kmerCount.iloc[kmerId.index,]

        kmerId.reset_index(drop=True, inplace=True)
        kmerCount.reset_index(drop=True, inplace=True)

    return (kmerId, kmerCount)

#------------------- Private Classes & Functions ------------#

def getDataStr(**kwargs):
    dataStr   = []

    fStr = kwargs.get('f')
    if (fStr is not None):
        dataStr.append(fStr)

    wStr = kwargs.get('w')
    if (wStr is not None):
        wStr = str(wStr) + 'w'
        dataStr.append(wStr)
    
    sStr = kwargs.get('s')
    if (sStr is not None):
        sStr = str(sStr) + 's'
        dataStr.append(sStr)

    bStr = kwargs.get('b')
    if (bStr is not None):
        bStr = str(bStr) + 'bp'
        dataStr.append(bStr)

    kmerStr = str(kwargs.get('k')) + 'mer'
    dataStr.append(kmerStr)
    return '_'.join(dataStr)

def rotate(kmerDf):
    kmerDf = pd.pivot_table(kmerDf, index=[kmer.ID_COL_NAME, kmer.FILE_COL_NAME],
        columns=kmer.KMER_COL_NAME, values=kmer.PERCENTAGE_COL_NAME, fill_value=0)
    kmerDf.reset_index(inplace=True)
    kmerDf.columns.name = None
    return kmerDf

def split(kmerDf):
    idColNames    = [kmer.ID_COL_NAME, kmer.FILE_COL_NAME]
    countColNames = list(kmerDf.columns)
    countColNames.remove(kmer.ID_COL_NAME)
    countColNames.remove(kmer.FILE_COL_NAME)

    idCols    = kmerDf.loc[:, idColNames]
    countCols = kmerDf.loc[:, countColNames]
    return (idCols, countCols)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
