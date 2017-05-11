import numpy as np

__all__ = ['showMap', 'showDomains', 'showEmbedding']

## normalization methods ##
def div0(a, b):
    """ Performs ``true_divide`` but ignores the error when division by zero 
    (result is set to zero instead). """

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        if np.isscalar(c):
            if not np.isfinite(c):
                c = 0
        else:
            c[~np.isfinite(c)] = 0.  # -inf inf NaN
    return c

def showMap(map, spec='', **kwargs):
    """A convenient function that can be used to visualize Hi-C contact map. 
    *kwargs* will be passed to :func:`matplotlib.pyplot.imshow`.

    :arg map: a Hi-C contact map.
    :type map: :class:`numpy.ndarray`

    :arg spec: a string specifies how to preprocess the matrix. Blank for no preprocessing,
    'p' for showing only data from *p*-th to *100-p*-th percentile.
    :type spec: str

    :arg p: specifies the percentile threshold.
    :type p: double
    """

    assert isinstance(map, np.ndarray), 'map must be a numpy.ndarray.'
    
    from matplotlib.pyplot import figure, imshow

    if not '_' in spec:
        figure()
    
    if 'p' in spec:
        p = kwargs.pop('p', 5)
        vmin = np.percentile(map, p)
        vmax = np.percentile(map, 100-p)
    else:
        vmin = vmax = None
  
    return imshow(map, vmin=vmin, vmax=vmax, **kwargs)

def showDomains(domains, linespec='r-', **kwargs):
    """A convenient function that can be used to visualize Hi-C structural domains. 
    *kwargs* will be passed to :func:`matplotlib.pyplot.plot`.

    :arg domains: a 2D array of Hi-C domains, such as [[start1, end1], [start2, end2], ...].
    :type domains: :class:`numpy.ndarray`
    """

    domains = np.array(domains)
    shape = domains.shape

    if len(shape) < 2:
        # convert to domain list if labels are provided
        indicators = np.diff(domains)
        indicators = np.append(1., indicators)
        indicators[-1] = 1
        sites = np.where(indicators != 0)[0]
        starts = sites[:-1]
        ends = sites[1:]
        domains = np.array([starts, ends]).T

    from matplotlib.pyplot import figure, plot

    x = []; y = []
    lwd = kwargs.pop('linewidth', 1)
    linewidth = np.abs(lwd)
    for i in range(len(domains)):
        domain = domains[i]
        start = domain[0]; end = domain[1]
        if lwd > 0:
            x.extend([start, end, end])
            y.extend([start, start, end])
        else:
            x.extend([start, start, end])
            y.extend([start, end, end])

    return plot(x, y, linespec, linewidth=linewidth, **kwargs)

def showEmbedding(modes, labels=None):
    if isinstance(modes, ModeSet):
        V = modes.getEigvecs()
    elif isinstance(modes, Mode):
        V = modes.getEigvec()
    elif isinstance(modes, np.ndarray):
        V = modes
    else:
        try:
            mode0 = modes[0]
            if isinstance(mode0, Mode):
                V = np.empty((len(mode0),0))
                for mode in modes:
                    assert isinstance(mode, Mode), 'Modes should be a list of modes.'
                    v = mode.getEigvec()
                    v = np.expand_dims(v, axis=1)
                    V = np.hstack((V, v))
            else:
                V = np.array(modes)
        except TypeError:
            TypeError('Modes should be a list of modes.')
    if V.ndim == 1:
        V = np.expand_dims(V, axis=1)

    if labels is not None:
        if len(labels) != V.shape[0]:
            raise ValueError('Modes (%d) and the Hi-C map (%d) should have the same number'
                                ' of atoms. Turn off "useTrimed" if you intended to apply the'
                                ' modes to the full map.'
                                %(V.shape[0], len(labels)))
    
    # normalize the rows so that feature vectors are unit vectors
    la = importLA()
    norms = la.norm(V, axis=1)
    N = np.diag(div0(1., norms))
    V = np.dot(N, V)

    X, Y, Z = V[:,:3].T

    from matplotlib.pyplot import figure
    from mpl_toolkits.mplot3d import Axes3D
    f = figure()
    ax = Axes3D(f)
    ax.plot(X, Y, Z, ':', color=[0.3, 0.3, 0.3])
    if labels is None:
        C = 'b'
    else:
        C = labels
    ax.scatter(X, Y, Z, s=30, c=C, depthshade=True)