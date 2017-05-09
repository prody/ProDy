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

    kmeans = KMeans(**kwargs).fit(V)
    return kmeans.labels_