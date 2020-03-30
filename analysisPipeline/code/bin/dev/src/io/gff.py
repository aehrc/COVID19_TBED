#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import os

# External imports
import pandas as pd

# Internal imports
from .common import isZFile

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*filepaths, asPdf=True):
    featureDBs = (_readGFF(f, asPdf) for f in filepaths)
    return featureDBs

#------------------- Private Classes & Functions ------------#

def _readGFF(filepath, asPdf):
    if (isZFile(filepath)):
        tmp = os.path.splitext(filepath)[0]
        if (not _isGFFFile(tmp)):
            raise NotImplementedError("Unknown GFF file")

        featureDB = _getFeatureDB(filepath, asPdf)

    else:
        if (not _isGFFFile(filepath)):
            raise NotImplementedError("Unknown GFF file")

        featureDB = _getFeatureDB(filepath, asPdf)

    return featureDB

def _getFeatureDB(filepath, asPdf):
    featureDB = None
    if (asPdf):
        featureDB = pd.read_csv(filepath, sep='\t', 
            comment='#', header=None);

    else:
        import gffutils             ## Requires python 3.5; not 3.7
        featureDB = gffutils.create_db(filepath,
            ':memory:', force=True, keep_order=True,
            merge_strategy='merge', sort_attribute_values=True)

        # featureDB.update(featureDB.create_introns())
    return featureDB

def _isGFFFile(filepath):
    if (filepath.endswith('.gff') \
        or filepath.endswith('.gff3')):
        return True

    return False

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
