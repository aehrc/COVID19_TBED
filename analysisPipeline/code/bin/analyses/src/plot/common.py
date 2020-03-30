#!/bin/python

#------------------- Description & Notes --------------------#

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import holoviews as hv
from holoviews import dim         ## Requires python 3.7; not 3.5
from holoviews import opts

# Internal imports

#------------------- Constants ------------------------------#

#------------------- Public Classes & Functions -------------#

def setLibrary(library='bokeh'):
    ## Use Bokeh by default
    if (library == 'bokeh'):
        hv.extension('bokeh')
        hv.archive.auto(filename_formatter="{obj:.7}")    ## For notebooks

        opts.defaults(
            opts.Scatter(tools=['hover'], width=700, height=700, padding=0.05),
            opts.HeatMap(tools=['hover'], width=700, height=700, labelled=[],
                         xrotation=45, colorbar=True, cmap=('Blues')))

        #     opts.HeatMap(tools=['hover'], width=700, height=700, labelled=[],
        #                  xrotation=45, colorbar=True, cmap=('Blues')))
        
        ## The library that Bokeh uses to export to SVG is not longer supported
        ## and so cannot be exported to SVG

    elif (library == 'matplotlib'):
        hv.extension('matplotlib')
        hv.output(fig='svg')

        opts.defaults(
            opts.Scatter(fig_size=300, padding=0.05),
            opts.HeatMap(fig_size=300, labelled=[], xrotation=45, 
                         colorbar=True, cmap=('Blues')))

    else:
        raise NotImplementedError("Unknown plotting library.")

#------------------- Private Classes & Functions ------------#

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()

#------------------------------------------------------------------------------
