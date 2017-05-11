__all__ = ['KMeans']

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
        from scipy.cluster.hierarchy import linkage
        from scipy.cluster.hierarchy import fcluster
    except ImportError:
        raise ImportError('Use of this function (Hierarchical) requires the '
                          'installation of scipy.')
    
    method = kwargs.pop('method', 'ward')
    metric = kwargs.pop('metric', 'euclidean')
    Z = linkage(V, method=method, metric=metric)

    t = kwargs.pop('t', 30)
    criterion = kwargs.pop('criterion', 'inconsistent')
    depth = kwargs.pop('depth', 2)
    R = kwargs.pop('R', None)
    monocrit = kwargs.pop('monocrit', None)
    labels = fcluster(Z, t, criterion=criterion, depth=depth, R=R, monocrit=monocrit)
    return labels.flatten()