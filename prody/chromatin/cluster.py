__all__ = ['KMeans']

def KMeans(V, **kwargs):
    try:
        from sklearn.cluster import KMeans
    except ImportError:
        raise ImportError('Use of this function (KMeans) requires the '
                          'installation of sklearn.')
    
    n_clusters = kwargs.get('n_clusters', 5)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(V)