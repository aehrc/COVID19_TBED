#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import os
from pathlib import Path

# External imports
import pandas as pd

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def readTable(filepath):
    compression = 'gzip' if isZFile(filepath) else None
    t           = pd.read_csv(filepath, sep='\t', compression=compression)
    return t

def writePca(filepath, pca):
    ## Create DIR if it doesnt exist
    outputDir = os.path.dirname(filepath)
    createDirIfNone(outputDir)
    removeFileIfExists(filepath)
    pca.to_csv(filepath,
        index=False, sep='\t')

#------------------- Private Classes & Functions ------------#

def isZFile(filepath):
    if (filepath.endswith('.gz') \
        or filepath.endswith('.gzip')):
        return True

    return False

def createDirIfNone(filepath):
    p = Path(filepath)
    p.mkdir(parents=True, exist_ok=True)

def removeFileIfExists(filepath):
    p = Path(filepath)
    if (os.path.exists(filepath)):
        p.unlink()

def getTempDir(filepath):
    tmpName = ".tmp_" + os.path.basename(filepath)
    f       = os.path.abspath(__file__)
    d       = os.path.dirname(f)
    tmpDir  = os.path.join(d, tmpName)
    return tmpDir

def isValidDir(filepath):
    ## Check that we can safely write the output
    p = Path(filepath)
    if (p.exists()):
        errMsg = """Directory ({}) exists. Might be important,
                    so we better not continue. Please remove
                 """.format(p)
        raise NameError(errMsg)

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
