import numpy as np
from prody import LOGGER, SETTINGS
from prody.utilities import showFigure

__all__ = ['KMeans', 'Hierarchy', 'showLinkage']

def KMeans(V, **kwargs):
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        raise ImportError('Use of this function (KMeans) requires the '
                          'installation of sklearn.')
    
    n_clusters = kwargs.get('n_clusters', None)
    if n_clusters is None:
        raise ValueError('KMeans requires to desiginate the number of clusters.')
    
    n_init = kwargs.pop('n_init', 100)
    
    kmeans = KMeans(n_init=n_init, **kwargs).fit(V)
    return kmeans.labels_

def Hierarchy(V, **kwargs):
    try:
        from scipy.cluster.hierarchy import linkage, fcluster, inconsistent
    except ImportError:
        raise ImportError('Use of this function (Hierarchy) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'complete')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)
    I = inconsistent(Z)
    i = np.percentile(I[:,3], 99.9)

    t = kwargs.pop('t', i)
    criterion = kwargs.pop('criterion', 'inconsistent')
    depth = kwargs.pop('depth', 2)
    R = kwargs.pop('R', None)
    monocrit = kwargs.pop('monocrit', None)

    n_clusters = kwargs.pop('n_clusters', None)
    if n_clusters is not None:
        criterion = 'maxclust'
        t = n_clusters
    labels = fcluster(Z, t, criterion=criterion, depth=depth, R=R, monocrit=monocrit)
    return labels.flatten()

def showLinkage(V, **kwargs):
    from .functions import _getEigvecs

    V = _getEigvecs(V, row_norm=True)
    try:
        from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    except ImportError:
        raise ImportError('Use of this function (showLinkage) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'complete')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)

    no_labels = kwargs.pop('no_labels', True)
    dendrogram(Z, no_labels=no_labels, **kwargs)
    if SETTINGS['auto_show']:
        showFigure()
    return Z