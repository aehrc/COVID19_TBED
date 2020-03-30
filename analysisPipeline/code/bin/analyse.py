#!/usr/bin/env python
# coding: utf-8

# # Coronavirus Kmer Analyses

# This notebook contains analyses of Kmers arising from Coronavirus sequences.

# ## Notebook Initialisation

# In[ ]:


## Date of analysis
import datetime

currDate = datetime.datetime.now().date()
currTime = datetime.datetime.now().time()

print('This analysis was run on:')
print(currDate.strftime('%Y-%m-%d' + ' ')
      + currTime.strftime('%H:%M:%S'))


# In[ ]:


## Library imports
## Standard library imports
import os
import sys
import warnings
warnings.simplefilter('ignore')

from itertools import product
from pathlib import Path
path = Path('../..').resolve()
if (str(path) not in sys.path):
    ## Add the parent path to the sys path
    ## This let's us import modules from analyse and dev
    sys.path.append(str(path))

## External imports
import numpy as np
import pandas as pd
import holoviews as hv

## Internal imports
from analyses.src import cov
from analyses.src import kmer
from analyses.src import ml
from analyses.src import plot
from dev.src import io


# In[ ]:


## IO related

## **********
## *** Use the following if we are running locally
## **********
# OUTPUT_DIR    = Path('kmerSigs')
# COVID_WG_DIR  = Path(OUTPUT_DIR, 'freqs')
# DUP_INFO_FILE = Path('files/200317.trimmed.noDup.fasta.identicalIDs.txt')

## **********
## *** Use the following if we are running on HPC
## **********
# FREQUENCY_OUTPUT_DIR = Path(os.environ['FREQUENCY_OUTPUT_DIR'])
# OUTPUT_DIR           = Path('.')
# WG_DIR               = Path(FREQUENCY_OUTPUT_DIR, 'default/fasta/cov_genomes')
# COVID_WG_DIR         = Path(WG_DIR, 'COVID-19_10mer')

## **********
## *** Use the following if we are running on through the pipelin
## **********
OUTPUT_DIR    = Path(sys.argv[3])
COVID_WG_DIR  = Path(OUTPUT_DIR, sys.argv[1])
DUP_INFO_FILE = Path(sys.argv[2])

DIRS          = [COVID_WG_DIR]


# In[ ]:


## Plotting related
PLOT_LIBRARY = 'matplotlib'
ZOOMED_1X    = 0.01
ZOOMED_2X    = 0.2
ZOOMED_3X    = 0.4

def getPlots(kmerPcas, oFilename):
    scatters = [cov.plot.scatter(PLOT_LIBRARY, kmerPca)
                for kmerPca in kmerPcas]

    ## Create the map of plots
    plots    = {a:b for a, b in zip(ARGS, scatters)}
    plots    = hv.HoloMap(plots, kdims='kmer')
    plots    = plots.collate()

    if (PLOT_LIBRARY != 'bokeh'):
        oFile    = Path(OUTPUT_DIR, oFileName)
        hv.save(plots, oFile)

    return plots


# # Genome Kmer frequency analyses

# ## Whole genome

# In[ ]:


## Table formatting-related
def getKmerPca(kmerId, kmerCount, oFilter=0):
    kmerPca = ml.feature.reduce(kmerCount)
    if (oFilter != 0):
        kmerPca = ml.outlier.detect(kmerPca, contamination=oFilter)

    kmerPca = pd.concat([kmerId, kmerPca], axis=1)
    return kmerPca

def getPcaDistance(kmerPca):
    ## Scale the coordinates so that its a bit easier
    kmerPca['PCA1'] = kmerPca['PCA1'] * 10000
    kmerPca['PCA2'] = kmerPca['PCA2'] * 10000

    idDf     = kmerPca[[kmer.ID_COL_NAME, kmer.FILE_COL_NAME]]
    pcaDf    = kmerPca[ml.PCA_DATA_COL_NAMES]
    kmerDist = ml.distance.euclidean(pcaDf)
    kmerDist = ml.distance.addIdColumns(idDf, kmerDist)
    return kmerDist


# In[ ]:


## Read Kmer frequencies
KMERS      = [i for i in range(10, 11)]    ## Doesn't actually do anything 
ARGS       = KMERS

kmerDirs   = [DIRS for k in ARGS]
kmerDfs    = [io.kmer.read(*kmerDir) for kmerDir in kmerDirs]
kmerDfs    = [kmer.rotateAndSplit(kmerDf) for kmerDf in kmerDfs]


# ### SARS-CoV-2

# #### All

# In[ ]:


## Pre-filtering
filtKmerDfs    = [kmer.filterByFilename(kmerId, kmerCount)
                  for kmerId, kmerCount in kmerDfs]
filtKmerIds    = [kmerDf[0] for kmerDf in filtKmerDfs]
filtKmerCounts = [kmerDf[1] for kmerDf in filtKmerDfs]

## Table formatting
dupInfo        = cov.io.parseDupInfo(DUP_INFO_FILE)
filtKmerIds    = [kmerId.merge(dupInfo, how='left', on='id')
                  for kmerId in filtKmerIds]
filtKmerIds    = [cov.style.default(kmerId) for kmerId in filtKmerIds]

## Run PCA
kmerPcas0x     = [getKmerPca(kmerId, kmerCount)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas1x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_1X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas2x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_2X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas3x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_3X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]

## Post-filtering
kmerPcas1x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas1x]
kmerPcas2x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas2x]
kmerPcas3x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas3x]


# In[ ]:


## Create PCA distance matrix
kmerDists = [getPcaDistance(kmerPca.copy()) for kmerPca in kmerPcas0x]

oFile = Path(OUTPUT_DIR, 'SARS-CoV-2_distance.txt')  
[kmerDist.to_csv(oFile, sep='\t', index=False) for kmerDist in kmerDists]

datasets  = [hv.Dataset(kmerDist, ['id_x', 'id_y'])
             for kmerDist in kmerDists]
heatmaps  = [d.to(hv.HeatMap, ['id_x', 'id_y'])
             for d in datasets]

plot.setLibrary(PLOT_LIBRARY)
if (PLOT_LIBRARY == 'bokeh'):
    ## Style plots

    ## Create the map of plots
    plots    = {a:b for a, b in zip(ARGS, heatmaps)}
    plots    = hv.HoloMap(plots, kdims='kmer')
    plots    = plots.collate()
    plots

else:
    ## Matplotlib doesnt work the same way as Bokeh
    ## Labels cant be removed (for some reason)
    pass


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.0x_all.svg'
plots     = getPlots(kmerPcas0x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.1x_all.svg'
plots     = getPlots(kmerPcas1x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.2x_all.svg'
plots     = getPlots(kmerPcas2x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.3x_all.svg'
plots     = getPlots(kmerPcas3x, oFileName)


# #### Wuhan-Hu-1 and AUS (Isolates)

# In[ ]:


## Pre-filtering
filtKmerDfs    = [kmer.filterByFilename(kmerId, kmerCount)
                  for kmerId, kmerCount in kmerDfs]
filtKmerIds    = [kmerDf[0] for kmerDf in filtKmerDfs]
filtKmerCounts = [kmerDf[1] for kmerDf in filtKmerDfs]

## Table formatting
dupInfo     = cov.io.parseDupInfo(DUP_INFO_FILE)
filtKmerIds = [kmerId.merge(dupInfo, how='left', on='id')
               for kmerId in filtKmerIds]
filtKmerIds = [cov.style.refAusIsolates(kmerId) for kmerId in filtKmerIds]

## Run PCA
kmerPcas0x     = [getKmerPca(kmerId, kmerCount)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas1x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_1X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas2x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_2X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas3x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_3X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]

## Post-filtering
kmerPcas1x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas1x]
kmerPcas2x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas2x]
kmerPcas3x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas3x]


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.0x_Wuhan.Aus.Animal_isolates.svg'
plots     = getPlots(kmerPcas0x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.1x_Wuhan.Aus.Animal_isolates.svg'
plots     = getPlots(kmerPcas1x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.2x_Wuhan.Aus.Animal_isolates.svg'
plots     = getPlots(kmerPcas2x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.3x_Wuhan.Aus.Animal_isolates.svg'
plots     = getPlots(kmerPcas3x, oFileName)


# #### Wuhan-Hu-1 and AUS (Isolate groups)

# In[ ]:


## Pre-filtering
filtKmerDfs    = [kmer.filterByFilename(kmerId, kmerCount)
                  for kmerId, kmerCount in kmerDfs]
filtKmerIds    = [kmerDf[0] for kmerDf in filtKmerDfs]
filtKmerCounts = [kmerDf[1] for kmerDf in filtKmerDfs]

## Table formatting
dupInfo     = cov.io.parseDupInfo(DUP_INFO_FILE)
filtKmerIds = [kmerId.merge(dupInfo, how='left', on='id')
               for kmerId in filtKmerIds]
filtKmerIds = [cov.style.refAusIsolateGroups(kmerId)
               for kmerId in filtKmerIds]

## Run PCA
kmerPcas0x     = [getKmerPca(kmerId, kmerCount)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas1x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_1X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas2x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_2X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]
kmerPcas3x     = [getKmerPca(kmerId, kmerCount, oFilter=ZOOMED_3X)
                  for kmerId, kmerCount in zip(filtKmerIds, filtKmerCounts)]

## Post-filtering
kmerPcas1x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas1x]
kmerPcas2x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas2x]
kmerPcas3x     = [ml.filterByOutlier(kmerPca) for kmerPca in kmerPcas3x]


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.0x_Wuhan.Aus.Animal_isolate.groups.svg'
plots     = getPlots(kmerPcas0x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.1x_Wuhan.Aus.Animal_isolate.groups.svg'
plots     = getPlots(kmerPcas1x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.2x_Wuhan.Aus.Animal_isolate.groups.svg'
plots     = getPlots(kmerPcas2x, oFileName)


# In[ ]:


## Create plot
plot.setLibrary(PLOT_LIBRARY)
oFileName = 'SARS-CoV-2.3x_Wuhan.Aus.Animal_isolate.groups.svg'
plots     = getPlots(kmerPcas3x, oFileName)

