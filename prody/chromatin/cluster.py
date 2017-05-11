import numpy as np

__all__ = ['KMeans', 'Hierarchy']

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
        raise ImportError('Use of this function (Hierarchical) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'single')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)
    I = inconsistent(Z)
    i = np.percentile(I[:,3], 99.9)

    t = kwargs.pop('t', i)
    criterion = kwargs.pop('criterion', 'inconsistent')
    depth = kwargs.pop('depth', 2)
    R = kwargs.pop('R', None)
    monocrit = kwargs.pop('monocrit', None)
    labels = fcluster(Z, t, criterion=criterion, depth=depth, R=R, monocrit=monocrit)
    return labels.flatten()