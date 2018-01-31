import numpy as np
from prody import LOGGER, SETTINGS
from prody.utilities import showFigure
from functions import _getEigvecs

__all__ = ['getGNMDomains', 'KMeans', 'Hierarchy', 'Discretize', 'showLinkage']

def KMeans(V, **kwargs):
    """Performs k-means clustering on *V*. The function uses :func:`sklearn.cluster.KMeans`. See sklearn documents 
    for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`numpy.ndarray`

    :arg n_clusters: specifies the number of clusters. 
    :type n_clusters: int
    """

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
    """Performs hierarchical clustering on *V*. The function essentially uses two scipy functions: ``linkage`` and 
    ``fcluster``. See :func:`scipy.cluster.hierarchy.linkage` and :func:`scipy.cluster.hierarchy.fcluster` for the 
    explaination of the arguments. Here lists arguments that are different from those of scipy.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`numpy.ndarray`

    :arg inconsistent_percentile: if the clustering *criterion* for :func:`scipy.cluster.hierarchy.fcluster`
    is ``inconsistent`` and threshold *t* is not given (default), then the function will use the percentile specified 
    by this argument as the threshold.
    :type inconsistent_percentile: double

    :arg n_clusters: specifies the maximal number of clusters. If this argument is given, then the function will 
    automatically set *criterion* to ``maxclust`` and *t* equal to *n_clusters*.
    :type n_clusters: int
    """

    try:
        from scipy.cluster.hierarchy import linkage, fcluster, inconsistent
    except ImportError:
        raise ImportError('Use of this function (Hierarchy) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'single')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)
    
    criterion = kwargs.pop('criterion', 'inconsistent')
    t = kwargs.get('t', None)
    ip = kwargs.pop('inconsistent_percentile', 99.9)
    if t is None and criterion == 'inconsistent':
        I = inconsistent(Z)
        i = np.percentile(I[:,3], ip)

    t = kwargs.pop('t', i)
    depth = kwargs.pop('depth', 2)
    R = kwargs.pop('R', None)
    monocrit = kwargs.pop('monocrit', None)

    n_clusters = kwargs.pop('n_clusters', None)
    if n_clusters is not None:
        criterion = 'maxclust'
        t = n_clusters
    labels = fcluster(Z, t, criterion=criterion, depth=depth, R=R, monocrit=monocrit)
    return labels.flatten()

def Discretize(V, **kwargs):
    try:
        from sklearn.cluster.spectral import discretize
    except ImportError:
        raise ImportError('Use of this function (Discretize) requires the '
                          'installation of sklearn.')

    n_clusters = kwargs.pop('n_clusters', None)
    if n_clusters is not None:
        print('Discretization does not need to desiginate the number of clusters.')

    labels = discretize(V, **kwargs)
    return labels

def showLinkage(V, **kwargs):
    """Shows the dendrogram of hierarchical clustering on *V*. See :func:`scipy.cluster.hierarchy.dendrogram` for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`numpy.ndarray`

    """

    V, _ = _getEigvecs(V, row_norm=True, remove_zero_rows=True)
    try:
        from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    except ImportError:
        raise ImportError('Use of this function (showLinkage) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'single')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)

    no_labels = kwargs.pop('no_labels', True)
    dendrogram(Z, no_labels=no_labels, **kwargs)
    if SETTINGS['auto_show']:
        showFigure()
    return Z
    
def getGNMDomains(modes, method=Hierarchy, **kwargs):
    """Uses spectral clustering to separate structural domains in the chromosome.
    
    :arg modes: GNM modes used for segmentation
    :type modes: :class:`ModeSet`

    :arg method: Label assignment algorithm used after Laplacian embedding of loci.
    :type method: func
    """

    V, mask = _getEigvecs(modes, row_norm=True, remove_zero_rows=True)

    labels_ = method(V, **kwargs)

    labels = np.empty(len(mask))
    labels.fill(np.nan)
    labels[mask] = labels_

    currlbl = labels_[np.argmax(~np.isnan(labels_))]

    for i in range(len(labels)):
        l = labels[i]
        if np.isnan(l):
            labels[i] = currlbl
        elif currlbl != l:
            currlbl = l

    return labels