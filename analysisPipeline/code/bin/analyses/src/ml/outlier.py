#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import numpy as np
import pandas as pd
import holoviews as hv
from holoviews import opts
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM

# Internal imports
from .. import plot
from .constants import *

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def detect(kmerPca, contamination=0.01):
    ## We expect 1% of the points are outliers
    pca     = kmerPca.loc[:, PCA_DATA_COL_NAMES]
    algo    = EllipticEnvelope(contamination=contamination)
    labels  = _getLabels(algo, pca)
    labels  = pd.DataFrame(labels, columns=[OLABEL_COL_NAME])
    kmerPca = pd.concat([kmerPca, labels], axis=1)
    return kmerPca

def investigateOptimalAlgorithms(kmerId, kmerPca):
    plot.setLibrary('bokeh')

    pca   = kmerPca.loc[:, PCA_DATA_COL_NAMES]
    plots = {}
    algos = (
        ('Elliptic', EllipticEnvelope()),
        ('SVM', OneClassSVM()),
        ('Forest', IsolationForest()),
        ('Local', LocalOutlierFactor()))

    ## Visualise data and manually determine which algorithm will be good
    for i, (name, algo) in enumerate(algos, 1):
        labels   = _getLabels(algo, pca)
        labels   = pd.DataFrame(labels, columns=[OLABEL_COL_NAME])
        kmerDf   = pd.concat([kmerId, pca, labels], axis=1)

        dataset  = hv.Dataset(kmerDf, PCA_DATA_COL_NAMES)
        scatter  = dataset.to(hv.Scatter, PCA_DATA_COL_NAMES, groupby=OLABEL_COL_NAME).overlay()
        scatter.opts(opts.Scatter(size=10, show_legend=True))
        plots[name] = scatter

    plots = hv.HoloMap(plots, kdims='algo')
    plots = plots.collate()
    return plots

#------------------- Private Classes & Functions ------------#

def _getLabels(algo, pca):
    labels = algo.fit_predict(pca) 
    labels = np.reshape(labels, (-1, 1))
    return labels

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
