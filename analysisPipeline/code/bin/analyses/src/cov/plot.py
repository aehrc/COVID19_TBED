#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import holoviews as hv
from holoviews import dim         ## Requires python 3.7; not 3.5
from holoviews import opts

# Internal imports
from .. import ml

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def scatter(library, kmerPca):
    ## Create plots
    d = hv.Dataset(kmerPca, ml.PCA_DATA_COL_NAMES)
    s = d.to(hv.Scatter, ml.PCA_DATA_COL_NAMES, groupby='class').overlay()

    ## Style plots
    options = _getOptions(library)
    s.opts(options)
    return s

#------------------- Private Classes & Functions ------------#

def _getOptions(library):
    if (library == 'bokeh'):
        options  = [opts.Scatter(size='bSize', color='color',
                                 show_legend=False)]

    else:
        options  = [opts.Scatter(s='mSize', color='color',
                                 show_legend=True)]
    return options

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
