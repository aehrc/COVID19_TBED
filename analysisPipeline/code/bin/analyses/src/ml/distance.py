#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

# Internal imports
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def euclidean(df):
    df   = df.values.astype(np.float32)
    dist = euclidean_distances(df, df)
    dist = pd.DataFrame(dist)
    dist = dist.stack().reset_index()
    dist.rename(columns={0:'distance'}, inplace=True)
    return dist

def addIdColumns(ids, dist):
    ## Append ID information
    dist = dist.merge(ids.reset_index(), left_on='level_0', right_on='index')
    dist = dist.merge(ids.reset_index(), left_on='level_1', right_on='index')

    ## Remove redundant columns
    colsToRemove = ['level_0', 'level_1', 'index_x', 'index_y']
    dist.drop(columns=colsToRemove, inplace=True)
    return dist

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
