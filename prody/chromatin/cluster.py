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
    
    n_init = kwargs.get('n_init')
    if n_init is None: n_init = 100
    
    kmeans = KMeans(n_init=n_init, **kwargs).fit(V)
    return kmeans.labels_