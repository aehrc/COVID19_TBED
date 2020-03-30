#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports

# Internal imports

#------------------- Constants ------------------------------#

""" PCA related constants """
PCA1_COL_NAME   = 'PCA1'
PCA2_COL_NAME   = 'PCA2'
PCA3_COL_NAME   = 'PCA3'

## Data column names
DATA_COL_NAMES_2D = [PCA1_COL_NAME, PCA2_COL_NAME]
DATA_COL_NAMES_3D = [PCA1_COL_NAME, PCA2_COL_NAME, PCA3_COL_NAME]

## [OPTIONAL] Update to DATA_COL_NAMES_3D for 3D PCA
PCA_DATA_COL_NAMES = DATA_COL_NAMES_2D
N_PCA_COMPONENTS   = len(PCA_DATA_COL_NAMES)

""" Clustering related constants """
CLABEL_COL_NAME = 'cLabel'
EPS_COL_NAME    = 'eps'

""" Density related constants """
DPROB_COL_NAME  = 'dProb'
DLABEL_COL_NAME = 'dLabel'

""" Outlier related constants """
OLABEL_COL_NAME = 'oLabel'

#------------------- Public Classes & Functions -------------#

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
