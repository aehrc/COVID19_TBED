#!/bin/python

#------------------- Description & Notes --------------------#

'''
tSNE seems to represent both global and local differences well
Both MDS and PCA appear to represent global differences well. However,
they don't seem to represent the local differences very well.
'''

#------------------- Dependencies ---------------------------#

# Standard library imports

# External imports
import pandas as pd
import holoviews as hv
from holoviews import opts
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.pipeline import make_pipeline
from sklearn import decomposition
from sklearn import preprocessing
from sklearn import manifold

# Internal imports
from .. import kmer
from .. import plot
from ..util import spark
from .constants import *
from .common import getParameterGrid

#------------------- Constants ------------------------------#


#------------------- Public Classes & Functions -------------#

def select(kmerId, kmerCount):
    ## Is there a way to see which columns are the most important?
    ## TODO
    method = ExtraTreesClassifier(n_estimators=50)
    fit    = method.fit(kmerCount, kmerId)

    model  = SelectFromModel(fit, prefit=True)
    kmerCount = model.transform(kmerCount)
    return kmerCount

def reduce(kmerCount, algo='PCA', random_state=42):
    algo = _getAlgorithm(algo, random_state)
    com  = _getComponents(algo, kmerCount)
    com  = pd.DataFrame(com , columns=PCA_DATA_COL_NAMES)
    return com 

def investigateOptimalAlgorithms(kmerId, kmerCount):
    plot.setLibrary('bokeh')

    plots  = {}
    params = {'n_components':N_PCA_COMPONENTS, 'random_state':42}
    algos  = (
        ('PCA', decomposition.PCA(**params)),
        ('LLE', manifold.LocallyLinearEmbedding(method='standard', **params)),
        ('LTSA', manifold.LocallyLinearEmbedding(method='ltsa', **params)),
        ('Hessian LLE', manifold.LocallyLinearEmbedding(method='hessian',
            n_neighbors=10, **params)),
        ('Modified LLE', manifold.LocallyLinearEmbedding(method='modified',
            **params)),
        ('tSNE', manifold.TSNE(**params)),
        ('Isomap', manifold.Isomap(n_components=N_PCA_COMPONENTS)),
        ('MDS', manifold.MDS(**params)),
        ('SE', manifold.SpectralEmbedding(**params)))

    ## Visualise data and manually determine which algorithm will be good
    for i, (name, algo) in enumerate(algos, 1):
        com     = _getComponents(algo, kmerCount)
        com     = pd.DataFrame(com, columns=PCA_DATA_COL_NAMES)
        kmerDf  = pd.concat([kmerId, com], axis=1)

        dataset = hv.Dataset(kmerDf, PCA_DATA_COL_NAMES)
        scatter = dataset.to(hv.Scatter, PCA_DATA_COL_NAMES)
        scatter.opts(opts.Scatter(size=10, show_legend=True))
        plots[name] = scatter

    plots = hv.HoloMap(plots, kdims='algo')
    plots = plots.collate()
    return plots

def investigateOptimalParameters(kmerId, kmerCount):
    ## TO FIX
    import holoviews as hv
    from holoviews import dim         ## Requires python 3.7; not 3.5
    from holoviews import opts
    hv.extension('bokeh')

    kmerId = updateDuplicates(kmerId, kmer.ID_COL_NAME)
    cond = ((kmerId[kmer.ID_COL_NAME].str.contains('Wuhan-Hu-1'))
             | (kmerId[kmer.ID_COL_NAME].str.match('Australia'))
             | (kmerId[kmer.ID_COL_NAME].str.match('Sydney'))
             | (kmerId[kmer.ID_COL_NAME].str.match('C6')))
    kmerId.loc[cond, kmer.FILE_COL_NAME] = kmerId.loc[cond, kmer.ID_COL_NAME]

    labels = []
    # algos = investigatePCAParameters(kmerId, kmerCount)
    # algos = investigateMDSParameters(kmerId, kmerCount)
    algos = investigateTSNEParameters(kmerId, kmerCount)

    plots = {}
    for p, a in algos:
        print(p)
        pca     = getPca(a, kmerCount)
        kmerPca = joinColumns(kmerId, PCA_DATA_COL_NAMES, pca)
        dataset = hv.Dataset(kmerPca, PCA_DATA_COL_NAMES)
        scatter = dataset.to(hv.Scatter, PCA_DATA_COL_NAMES, 
                             groupby=kmer.FILE_COL_NAME) \
                         .overlay()
        scatter.opts(opts.Scatter(tools=['hover'], height=700, width=700,
                                  size=10, show_legend=True))

        plots[tuple(p.values())] = scatter
        labels = p.keys()

    ## Create the map of plots
    plots = hv.HoloMap(plots, kdims=[*labels])
    plots = plots.collate()
    hv.save(plots, 'plot.html')

#------------------- Private Classes & Functions ------------#

def _getAlgorithm(algo='PCA', random_state=42):
    if (algo == 'PCA'):
        algo = decomposition.PCA(n_components=N_PCA_COMPONENTS, 
                                 random_state=random_state)

    elif (algo == 'TSNE'):
        algo = manifold.TSNE(n_components=N_PCA_COMPONENTS, 
                             random_state=random_state, 
                             perplexity=30, n_iter=1000)

    elif (algo == 'MDS'):
        algo = manifold.MDS(n_components=N_PCA_COMPONENTS, 
                            random_state=random_state, 
                            max_iter=500, n_init=10)

    return algo

def _getComponents(algo, kmerCount):
    pipe   = make_pipeline(preprocessing.StandardScaler(), algo)
    comps  = pipe.fit_transform(kmerCount)
    ## With 2D PCA, we are only explaining
    ## about 70% of the variance (not great, but OK)
    # print(pipe.steps[1][1].explained_variance_ratio_)

    ## MinMaxScaler doesn't like being in a pipeline...
    scaler = preprocessing.MinMaxScaler()
    comps    = scaler.fit_transform(comps)
    return comps

def investigatePCAParameters(kmerId, kmerCount):
    params = getParameterGrid(
        ('n_components', [N_PCA_COMPONENTS]),
        ('random_state', [42]))

    algos  = ((p, decomposition.PCA(**p)) for p in params)
    return algos

def investigateMDSParameters(kmerId, kmerCount):
    params = getParameterGrid(
        ('n_components', [N_PCA_COMPONENTS]),
        ('n_init', list(range(4, 20, 4))),
        ('max_iter', list(range(100, 1000, 100))))

    algos  = ((p, manifold.MDS(**p)) for p in params)
    return algos

def investigateTSNEParameters(kmerId, kmerCount):
    params = getParameterGrid(
        ('n_components', [N_PCA_COMPONENTS]),
        ('perplexity', list(range(10, 100, 10))),
        ('n_iter', list(range(1000, 50000, 5000))))

    ## I think we're looking at an optimal perplexity somewhere between 50-70
    algos  = ((p, manifold.TSNE(**p)) for p in params)
    return algos

#------------------- Main -----------------------------------#

if (__name__ == "__main__"):
    main()
