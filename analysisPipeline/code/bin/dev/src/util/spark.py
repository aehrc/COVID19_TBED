#!/bin/python

#------------------- Description & Notes --------------------#

'''
If we encounter ModuleNotFound, check that the file
paths for the zipped files are relative to the main scripts
'''

#------------------- Dependencies ---------------------------#

# Standard library imports
import os

# External imports
from pyspark import SparkConf
from pyspark import SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.dataframe import DataFrame

# Internal imports

#------------------- Constants ------------------------------#

THIS_FILE = os.path.abspath(__file__)
THIS_DIR  = os.path.dirname(THIS_FILE)
SRC_DIR   = os.path.dirname(THIS_DIR)
SRC_ZIP   = os.path.join(SRC_DIR, "modules.zip")

## Number of Partitions multiplier for Spark
N_PARTS_MULTI = 6

#------------------- Public Classes & Functions -------------#

def getSparkContext():

    """
    Description:
        Generate a SparkContext

    Returns:
        sc (SparkContext)
            SparkContext environment
    """

    conf = getSparkConfiguration()
    sc   = SparkContext(pyFiles=[SRC_ZIP], conf=conf)
    return sc

def getSparkSession():

    """
    Description:
        Generate a SparkSession

    Returns:
        sc (SparkContext)
            SparkSession environment
    """

    conf = getSparkConfiguration()
    ss   = SparkSession.builder.config(conf=conf).getOrCreate()
    ss.sparkContext.addPyFile(SRC_ZIP)
    return ss

#------------------- Private Classes & Functions ------------#

def getSparkConfiguration():
    conf = SparkConf()
    conf.set('spark.driver.memory', '64G')
    conf.set('spark.driver.maxResultSize', '0')
    conf.set('spark.executor.memory', '32G')
    conf.set('spark.local.dir', THIS_DIR)
    conf.set('spark.network.timeout', '600s')
    conf.set('spark.executor.heartbeatInterval', '60s')
    conf.set('spark.sql.execution.arrow.enabled', 'true')
    return conf

## For chaining dataframe transformations
## https://adatis.co.uk/pyspark-dataframe-transformation-chaining/
def transform(self, f):
    return f(self)

DataFrame.transform = transform

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
