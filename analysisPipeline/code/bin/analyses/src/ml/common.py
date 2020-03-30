#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
from itertools import product

# External imports
import numpy as np
import pandas as pd

# Internal imports
from .. import kmer
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def filterByFilename(pca, *dirPaths):
    filenames = [d.name for d in dirPaths]
    pattern   = '|'.join(filenames)

    cond = pca[kmer.FILE_COL_NAME].str.contains(pattern)
    pca.drop(pca[cond].index, inplace=True)
    pca.reset_index(drop=True, inplace=True)
    return pca

def filterByOutlier(kmerPca):
    cond    = (kmerPca[OLABEL_COL_NAME] == 1)
    kmerPca = kmerPca.loc[cond].copy()
    kmerPca.reset_index(drop=True, inplace=True)
    return kmerPca

def prototypeAnalysis(tKmerPca, tKmerCount):
    isOutlier = tKmerPca[OLABEL_COL_NAME] == -1
    outliers  = tKmerPca.loc[isOutlier]
    oKmerId   = outliers.iloc[:, 0:2]
    oFreq     = tKmerCount.iloc[outliers.index, :]

    oKmerId.reset_index(inplace=True, drop=True)
    oFreq.reset_index(inplace=True, drop=True)

    (oPca, oPcaColNames) = runPca(oFreq)
    oCols                = [oPca]
    oColNames            = [oPcaColNames]

    (oCluster, oClusterColNames) = runClustering(oPca)
    oCols.append(oCluster)
    oColNames.append(oClusterColNames)

    oKmerPca = formatAnalysis(oKmerId, oCols, oColNames)
    return oKmerPca

def getParameterGrid(*params):
    pName = [f[0] for f in params]
    pVals = [f[1] for f in params]
    pVals = list(product(*pVals))

    params = [dict(zip(pName, v)) for v in pVals]
    return params

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
