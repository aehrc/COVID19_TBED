#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports
import os
import shutil
from pathlib import Path

# External imports
from pyspark.sql.functions import *
from pyspark.sql.types import *

# Internal imports
from ..util import spark
from .. import kmer
from .common import createDirIfNone
from .common import getTempDir

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def read(*dirs):
    with spark.getSparkSession() as ss:
        with ss.sparkContext as sc:
            filepaths = [str(f) for d in dirs for f in Path(d).glob("*.parquet")]
            kmerSdf   = ss.read.parquet(*filepaths)
            if (not _isKmerFrequencyFile(kmerSdf)):
                raise NotImplementedError("Unknown kmer frequency format")

            kmerSdf = kmerSdf.withColumn(kmer.FILE_COL_NAME, _getDirname(input_file_name()))
            kmerSdf = kmer.ungroupRows(kmerSdf)
            kmerPdf = kmerSdf.toPandas()

    return kmerPdf

def write(filepath, kmerDf):
    ## Create DIR if it does not exist
    createDirIfNone(filepath)

    try:
        tmpDir = getTempDir(filepath)

        ## Write the table to file in PARQUET format
        ## We can change the format in the future if it doesn't work out...
        kmerDf.write.parquet(tmpDir, mode='overwrite', compression="snappy")
        # kmerDf.write.csv(tmpDir, mode='overwrite')

        ## Remove the hidden files, we don't need them
        [p.unlink() for p in Path(tmpDir).glob("*.crc")]

        ## Move the result files to the filepath
        for source in Path(tmpDir).glob("*.parquet"):
            dest = Path(filepath, source.name)
            source.replace(dest)

    except:
        raise NotImplementedError("Not Implemented!")
        sys.exit(2)

    finally:
        ## Cleanup the temp directory
        shutil.rmtree(tmpDir, ignore_errors=True)

#------------------- Private Classes & Functions ------------#

@udf(returnType=StringType())
def _getDirname(filepath):
    return Path(filepath).parent.name

def _isKmerFrequencyFile(kmerDf):
    nCol = len(kmerDf.columns)
    nCol = nCol - 1 if (kmer.ID_COL_NAME in kmerDf.columns) else nCol
    if (nCol == 1):
        return True
    return False

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
