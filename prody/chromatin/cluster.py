import numpy as np
from prody import LOGGER, SETTINGS
from prody.dynamics import MaskedGNM
from prody.utilities import showFigure, bin2dec
from numbers import Integral

from .functions import _getEigvecs

__all__ = ['calcGNMDomains', 'Hingeplane', 'KMeans', 'Hierarchy', 'Discretize', 'showLinkage', 'GaussianMixture', 'BayesianGaussianMixture']

def Hingeplane(V, **kwargs):
    S = np.sign(np.sign(V) + 1)
    n, m = S.shape

    labels = np.zeros(n)
    for i, s in enumerate(S):
        labels[i] = bin2dec(s)

    uniq_labels = np.unique(labels)

    for i, l in enumerate(uniq_labels):
        labels[labels==l] = i
    return labels

def KMeans(V, **kwargs):
    """Performs k-means clustering on *V*. The function uses :func:`sklearn.cluster.KMeans`. See sklearn documents 
    for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`~numpy.ndarray`

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
    :type V: :class:`~numpy.ndarray`

    :arg inconsistent_percentile: if the clustering *criterion* for :func:`scipy.cluster.hierarchy.fcluster`
    is ``inconsistent`` and threshold *t* is not given (default), then the function will use the percentile specified 
    by this argument as the threshold.
    :type inconsistent_percentile: double

    :arg n_clusters: specifies the maximal number of clusters. If this argument is given, then the function will 
    automatically set *criterion* to ``maxclust`` and *t* equal to *n_clusters*.
    :type n_clusters: int
    """

    from scipy.cluster.hierarchy import linkage, fcluster, inconsistent
    
    method = kwargs.pop('linkage', 'single')
    metric = kwargs.pop('metric', 'cosine')
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
    """Adapted from :func:`~sklearn.cluster.spectral.discretize`. Copyright please 
    see LICENSE.rst.
    """

    from scipy.sparse import csc_matrix
    from scipy.linalg import LinAlgError

    copy = kwargs.pop('copy', False)
    max_svd_restarts = kwargs.pop('max_svd_restarts', 30)
    n_iter_max = kwargs.pop('n_iter_max', 20)
    random_state = kwargs.pop('random_state', None)
    info = kwargs.pop('info', None)

    def check_random_state(seed):
        """Adapted from :func:`~sklearn.utils.validation.check_random_state`."""
        if seed is None or seed is np.random:
            return np.random.mtrand._rand
        if isinstance(seed, Integral):
            return np.random.RandomState(seed)
        if isinstance(seed, np.random.RandomState):
            return seed
        raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                        ' instance' % seed)

    random_state = check_random_state(random_state)

    eps = np.finfo(float).eps
    vectors = np.array(V, dtype=float, copy=copy)
    n_samples, n_components = vectors.shape

    # Normalize the eigenvectors to an equal length of a vector of ones.
    # Reorient the eigenvectors to point in the negative direction with respect
    # to the first element.  This may have to do with constraining the
    # eigenvectors to lie in a specific quadrant to make the discretization
    # search easier.
    norm_ones = np.sqrt(n_samples)
    for i in range(vectors.shape[1]):
        vectors[:, i] *= norm_ones
        if vectors[0, i] != 0:
            vectors[:, i] = -1 * vectors[:, i] * np.sign(vectors[0, i])

    # Normalize the rows of the eigenvectors.  Samples should lie on the unit
    # hypersphere centered at the origin.  This transforms the samples in the
    # embedding space to the space of partition matrices.
    vectors = vectors / np.sqrt((vectors ** 2).sum(axis=1))[:, np.newaxis]

    svd_restarts = 0
    has_converged = False

    # If there is an exception we try to randomize and rerun SVD again
    # do this max_svd_restarts times.
    while (svd_restarts < max_svd_restarts) and not has_converged:

        # Initialize first column of rotation matrix with a row of the
        # eigenvectors
        rotation = np.zeros((n_components, n_components))
        rotation[:, 0] = vectors[random_state.randint(n_samples), :].T

        # To initialize the rest of the rotation matrix, find the rows
        # of the eigenvectors that are as orthogonal to each other as
        # possible
        c = np.zeros(n_samples)
        for j in range(1, n_components):
            # Accumulate c to ensure row is as orthogonal as possible to
            # previous picks as well as current one
            c += np.abs(np.dot(vectors, rotation[:, j - 1]))
            rotation[:, j] = vectors[c.argmin(), :].T

        last_objective_value = 0.0
        n_iter = 0

        while not has_converged:
            n_iter += 1

            t_discrete = np.dot(vectors, rotation)

            labels = t_discrete.argmax(axis=1)
            vectors_discrete = csc_matrix(
                (np.ones(len(labels)), (np.arange(0, n_samples), labels)),
                shape=(n_samples, n_components))

            t_svd = vectors_discrete.T * vectors

            try:
                U, S, Vh = np.linalg.svd(t_svd)
                svd_restarts += 1
            except LinAlgError:
                print("SVD did not converge, randomizing and trying again")
                break

            ncut_value = 2.0 * (n_samples - S.sum())
            if ((abs(ncut_value - last_objective_value) < eps) or
                    (n_iter > n_iter_max)):
                has_converged = True
            else:
                # otherwise calculate rotation and continue
                last_objective_value = ncut_value
                rotation = np.dot(Vh.T, U.T)

    if info is not None:
        if not isinstance(info, dict):
            raise TypeError('info must be a dict')
        info['indicators'] = t_discrete
        info['vectors'] = vectors

    if not has_converged:
        raise LinAlgError('SVD did not converge')
    return labels


def GaussianMixture(V, **kwargs):
    """Performs clustering on *V* by using Gaussian mixture models. The function uses :func:`sklearn.micture.GaussianMixture`. See sklearn documents 
    for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`~numpy.ndarray`

    :arg n_clusters: specifies the number of clusters. 
    :type n_clusters: int
    """

    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        raise ImportError('Use of this function (GaussianMixture) requires the '
                          'installation of sklearn.')
    
    n_components = kwargs.pop('n_components', None)
    if n_components == None:
        n_components = kwargs.pop('n_clusters',None)
        if n_components == None:
            n_components = 1
    
    n_init = kwargs.pop('n_init', 1)
    
    mixture = GaussianMixture(n_init=n_init, n_components=n_components, **kwargs).fit(V)

    return mixture.fit_predict(V)

def BayesianGaussianMixture(V, **kwargs):
    """Performs clustering on *V* by using Gaussian mixture models with variational inference. The function uses :func:`sklearn.micture.GaussianMixture`. See sklearn documents 
    for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`~numpy.ndarray`

    :arg n_clusters: specifies the number of clusters. 
    :type n_clusters: int
    """

    try:
        from sklearn.mixture import BayesianGaussianMixture
    except ImportError:
        raise ImportError('Use of this function (BayesianGaussianMixture) requires the '
                          'installation of sklearn.')
    
    n_components = kwargs.pop('n_components', None)
    if n_components == None:
        n_components = kwargs.pop('n_clusters',None)
        if n_components == None:
            n_components = 1
    
    n_init = kwargs.pop('n_init', 1)
    
    mixture = BayesianGaussianMixture(n_init=n_init, **kwargs).fit(V)

    return mixture.fit_predict(V)

def showLinkage(V, **kwargs):
    """Shows the dendrogram of hierarchical clustering on *V*. See :func:`scipy.cluster.hierarchy.dendrogram` for details.

    :arg V: row-normalized eigenvectors for the purpose of clustering.
    :type V: :class:`~numpy.ndarray`

    """

    V = _getEigvecs(V, row_norm=True)
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
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
    
def calcGNMDomains(modes, method=Discretize, **kwargs):
    """Uses spectral clustering to separate structural domains in chromosomes and proteins.
    
    :arg modes: GNM modes used for segmentation
    :type modes: :class:`ModeSet`

    :arg method: Label assignment algorithm used after Laplacian embedding of loci.
    :type method: func
    """

    dummy_mode = kwargs.pop('dummy_mode', True)
    row_norm = kwargs.pop('row_norm', True)
    linear = kwargs.pop('linear', False)

    V = _getEigvecs(modes, row_norm=row_norm, dummy_mode=dummy_mode)

    labels_ = method(V, **kwargs)

    if linear:
        split_labels = lambda l: np.split(l, np.where(np.diff(l) != 0)[0]+1)

        labels = split_labels(labels_)
        for i in range(len(labels)):
            l = np.empty_like(labels[i])
            l.fill(i)
            labels[i] = l

        labels = np.hstack(labels)
    else:
        labels = labels_

    if hasattr(modes, '_model'):
        model = modes._model
    else:
        model = modes
        
    if isinstance(model, MaskedGNM):
        currlbl = labels[0]
        labels = model._extend(labels, -1)
        
        for i, l in enumerate(labels):
            if l < 0:
                labels[i] = currlbl
            elif currlbl != l:
                currlbl = l

    return labels
    