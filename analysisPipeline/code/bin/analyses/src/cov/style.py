#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
from holoviews.plotting.util import process_cmap

# Internal imports

#------------------- Constants ------------------------------#

AUS_ISOLATES = 'AUS isolate'
AM_ISOLATES  = 'Animal model isolate'
RC_ISOLATES  = 'Isolate with possible recombination'
CC_ISOLATES  = 'Consensus sequence'

#------------------- Public Classes & Functions -------------#

def default(kmerDf):
    kmerDf = _addColumns(kmerDf)
    return kmerDf

def refAusIsolates(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows)
    return kmerDf

def refAusIsolateGroups(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows, newClass=AUS_ISOLATES)
    return kmerDf

def refAusRecombinationIsolates(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows) \
                                .pipe(_updateRecombinationRows)
    return kmerDf

def refAusRecombinationIsolateGroups(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows, newClass=AUS_ISOLATES) \
                                .pipe(_updateRecombinationRows, newClass=RC_ISOLATES)
    return kmerDf

def refAusAnimalModelIsolates(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows) \
                                .pipe(_updateAnimalModelRows)
    return kmerDf

def refAusAnimalModelIsolateGroups(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows, newClass=AUS_ISOLATES) \
                                .pipe(_updateAnimalModelRows, newClass=AM_ISOLATES)
    return kmerDf

def refAusAnimalModelConsensusIsolates(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows) \
                                .pipe(_updateAnimalModelRows) \
                                .pipe(_updateConsensusRows)
    return kmerDf

def refAusAnimalModelConsensusIsolateGroups(kmerDf):
    kmerDf = _addColumns(kmerDf).pipe(_updateRefRow) \
                                .pipe(_updateAusRows, newClass=AUS_ISOLATES) \
                                .pipe(_updateAnimalModelRows, newClass=AM_ISOLATES) \
                                .pipe(_updateConsensusRows, newClass=CC_ISOLATES)
    return kmerDf

#------------------- Private Classes & Functions ------------#

def _addColumns(kmerDf):
    kmerDf['class'] = 'SARS-CoV-2'
    kmerDf['color'] = '#bcbd22'
    kmerDf['bSize'] = 5
    kmerDf['mSize'] = 50
    return kmerDf

def _updateRefRow(kmerDf):
    cond = ((kmerDf['id'].str.contains('Wuhan-Hu-1'))
            | (kmerDf['dupId'].str.contains('Wuhan-Hu-1')))

    kmerDf = _updateClass(kmerDf, cond, 'Wuhan-Hu-1')
    kmerDf.loc[cond, 'bSize'] = 10
    kmerDf.loc[cond, 'mSize'] = 100
    kmerDf.loc[cond, 'color'] = '#000000'
    return kmerDf

def _updateAusRows(kmerDf, newClass=None):
    cond = ((kmerDf['id'].str.contains('Australia'))
            | (kmerDf['id'].str.contains('Sydney'))
            | (kmerDf['dupId'].str.contains('Australia'))
            | (kmerDf['dupId'].str.contains('Sydney')))
    
    numRows = len(kmerDf.loc[cond])
    colors  = process_cmap('glasbey_dark')[:numRows]
    
    kmerDf = _updateClass(kmerDf, cond, newClass)
    kmerDf.loc[cond, 'bSize'] = 10
    kmerDf.loc[cond, 'mSize'] = 100
    if (newClass is None):
        kmerDf.loc[cond, 'color'] = colors

    else:
        kmerDf.loc[cond, 'color'] = '#1f77b4'

    return kmerDf

def _updateAnimalModelRows(kmerDf, newClass=None):
    cond = ((kmerDf['id'].str.contains('Australia/VIC01'))
            | (kmerDf['id'].str.contains('BetaCoV_France_IDF0372_2020_C2'))
            | (kmerDf['id'].str.contains('Canada/ON-VIDO-01'))
            | (kmerDf['id'].str.contains('Germany/BavPat1'))
            | (kmerDf['id'].str.contains('USA/WA1'))
            | (kmerDf['dupId'].str.contains('Australia/VIC01'))
            | (kmerDf['dupId'].str.contains('BetaCoV_France_IDF0372_2020_C2'))
            | (kmerDf['dupId'].str.contains('Canada/ON-VIDO-01'))
            | (kmerDf['dupId'].str.contains('Germany/BavPat1'))
            | (kmerDf['dupId'].str.contains('USA/WA1')))
    
    numRows = len(kmerDf.loc[cond])
    colors  = set(process_cmap('glasbey_cool'))[:numRows]
    
    kmerDf = _updateClass(kmerDf, cond, newClass)
    kmerDf.loc[cond, 'bSize'] = 10
    kmerDf.loc[cond, 'mSize'] = 100
    if (newClass is None):
        kmerDf.loc[cond, 'color'] = colors
    else:
        kmerDf.loc[cond, 'color'] = '#2ca02c'

    return kmerDf

def _updateRecombinationRows(kmerDf, newClass=None):
    cond = ((kmerDf['id'].str.contains('Korea/KCDC05'))
            | (kmerDf['id'].str.contains('Chongqing/IVDC-CQ-001'))
            | (kmerDf['id'].str.contains('Italy/INMI1-cs'))
            | (kmerDf['dupId'].str.contains('Korea/KCDC05'))
            | (kmerDf['dupId'].str.contains('Chongqing/IVDC-CQ-001'))
            | (kmerDf['dupId'].str.contains('Italy/INMI1-cs')))

    numRows = len(kmerDf.loc[cond])
    colors  = set(process_cmap('glasbey_warm'))[:numRows]
    
    kmerDf = _updateClass(kmerDf, cond, newClass)
    kmerDf.loc[cond, 'bSize'] = 10
    kmerDf.loc[cond, 'mSize'] = 100
    if (newClass is None):
        kmerDf.loc[cond, 'color'] = colors

    else:
        kmerDf.loc[cond, 'color'] = '#ce6dbd'

    return kmerDf

def _updateConsensusRows(kmerDf, newClass=None):
    cond = (kmerDf['id'].str.match('C[1-6]'))
    
    numRows = len(kmerDf.loc[cond])
    colors  = process_cmap('glasbey_dark')[:numRows]
    
    kmerDf = _updateClass(kmerDf, cond, newClass)
    kmerDf.loc[cond, 'bSize'] = 10
    kmerDf.loc[cond, 'mSize'] = 100
    if (newClass is None):
        kmerDf.loc[cond, 'color'] = colors
    else:
        kmerDf.loc[cond, 'color'] = '#d62728'

    return kmerDf

def _updateClass(kmerDf, cond, newClass=None):
    if (newClass is None):
        kmerDf.loc[cond, 'class'] = kmerDf.loc[cond, 'id']

    else:
        kmerDf.loc[cond, 'class'] = newClass

    return kmerDf

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
