import numpy as np

from prody.dynamics import NMA, MaskedGNM
from prody.dynamics.mode import Mode
from prody.dynamics.modeset import ModeSet

from prody.utilities import importLA, copy, showFigure, div0
from prody import LOGGER, SETTINGS

__all__ = ['showDomains', 'showEmbedding', 'getDomainList']

## normalization methods ##

def showDomains(domains, linespec='-', **kwargs):
    """A convenient function that can be used to visualize Hi-C structural domains. 
    *kwargs* will be passed to :func:`matplotlib.pyplot.plot`.

    :arg domains: a 2D array of Hi-C domains, such as [[start1, end1], [start2, end2], ...].
    :type domains: :class:`numpy.ndarray`
    """

    fill_ends = kwargs.pop('fill_ends', 'close')
    domains = np.array(domains)
    shape = domains.shape

    if len(shape) == 1:
        # convert to domain list if labels are provided
        indicators = np.diff(domains)
        length = len(domains)
        if fill_ends in ['open', 'close']:
            indicators = np.append(1., indicators)
            indicators[-1] = 1
        elif fill_ends == 'skip':
            indicators = np.append(0., indicators)
        else:
            raise ValueError('invalid fill_ends mode: %s'%str(fill_ends))
        sites = np.where(indicators != 0)[0]
        starts = sites[:-1]
        ends = sites[1:]
        domains = np.array([starts, ends]).T
        consecutive = True
    elif len(shape) == 2:
        if domains.dtype == bool:
            length = domains.shape[1]
            domains_ = []
            for h in domains:
                start = None
                for i, b in enumerate(h):
                    if b:
                        if start is None:  # start
                            start = i
                    else:
                        if start is not None: # end
                            domains_.append([start, i-1])
                            start = None
                if start is not None:
                    domains_.append([start, i])
            domains = np.array(domains_)
        else:
            length = domains.max()
        consecutive = False
    else:
        raise ValueError('domains must be either one or two dimensions')

    from matplotlib.pyplot import figure, plot

    a = []; b = []
    lwd = kwargs.pop('linewidth', 1)
    lwd = kwargs.pop('lw', lwd)
    linewidth = np.abs(lwd)
    if fill_ends == 'open' and len(domains) == 1:
        domains = []

    for i in range(len(domains)):
        domain = domains[i]
        start = domain[0]; end = domain[1]
        if fill_ends == 'open' and start == 0:
            a.extend([end, end])
            b.extend([start, end])
        elif fill_ends == 'open' and end == length-1:
            a.extend([start, end])
            b.extend([start, start])
        else:
            a.extend([start, end, end])
            b.extend([start, start, end])

        if not consecutive:
            a.append(np.nan)
            b.append(np.nan)

    if lwd > 0:
        x = a; y = b
    else:
        x = b; y = a

    plt = plot(x, y, linespec, linewidth=linewidth, **kwargs)
    if SETTINGS['auto_show']:
        showFigure()
    return plt

def _getEigvecs(modes, row_norm=False, dummy_mode=False):
    la = importLA()

    if isinstance(modes, (Mode, ModeSet, NMA)):
        if hasattr(modes, '_model'):
            model = modes._model
        else:
            model = modes
        if isinstance(model, MaskedGNM):
            masked = model.masked
            model.masked = True
            V = modes.getArray()
            model.masked = masked
        else:
            V = modes.getArray()
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
            raise TypeError('Modes should be a list of modes.')
    if V.ndim == 1:
        V = np.expand_dims(V, axis=1)

    # add a dummy zero mode to the modeset
    if dummy_mode:
        v0 = V[:, 0]
        if np.allclose(v0, np.mean(v0)):
            dummy_mode = False
            LOGGER.warn('at least one zero mode is detected therefore dummy mode will NOT be added')

    if dummy_mode:
        n, _ = V.shape
        v0 = np.ones((n, 1), dtype=V.dtype)
        v0 /= la.norm(v0)
        V = np.hstack((v0, V))
        LOGGER.debug('a dummy zero mode is added')

    # normalize the rows so that feature vectors are unit vectors
    if row_norm:
        norms = la.norm(V, axis=1)
        N = np.diag(div0(1., norms))
        V = np.dot(N, V)

    return V

def showEmbedding(modes, labels=None, trace=True, headtail=True, cmap='prism'):
    """Visualizes Laplacian embedding of Hi-C data. 

    :arg modes: modes in which loci are embedded. It can only have 2 or 3 modes for the purpose 
    of visualization.
    :type modes: :class:`.ModeSet`

    :arg labels: a list of integers indicating the segmentation of the sequence.
    :type labels: list

    :arg trace: if **True** then the trace of the sequence will be indicated by a grey dashed line.
    :type trace: bool

    :arg headtail: if **True** then a star and a closed circle will indicate the head and the tail 
    of the sequence respectively.
    :type headtail: bool

    :arg cmap: the color map used to render the *labels*.
    :type cmap: str
    """
    V, _ = _getEigvecs(modes, True)
    m, n = V.shape

    if labels is not None:
        if len(labels) != m:
            raise ValueError('Modes (%d) and the Hi-C map (%d) should have the same number'
                                ' of atoms. Turn off "masked" if you intended to apply the'
                                ' modes to the full map.'
                                %(m, len(labels)))
    if n > 3:
        raise ValueError('This function can only visualize the embedding of 2 or 3 modes.')
    
    from matplotlib.pyplot import figure, plot, scatter
    from mpl_toolkits.mplot3d import Axes3D

    if n == 2:
        la = importLA()

        X, Y = V[:,:2].T
        R = np.array(range(len(X)))
        R = R / la.norm(R)
        X *= R; Y *= R
        
        f = figure()
        if trace:
            plot(X, Y, ':', color=[0.3, 0.3, 0.3])
        if labels is None:
            C = 'b'
        else:
            C = labels
        scatter(X, Y, s=30, c=C, cmap=cmap)
        if headtail:
            plot(X[:1], Y[:1], 'k*', markersize=12)
            plot(X[-1:], Y[-1:], 'ko', markersize=12)
    elif n == 3:
        X, Y, Z = V[:,:3].T
        
        f = figure()
        ax = Axes3D(f)
        if trace:
            ax.plot(X, Y, Z, ':', color=[0.3, 0.3, 0.3])
        if labels is None:
            C = 'b'
        else:
            C = labels
        ax.scatter(X, Y, Z, s=30, c=C, depthshade=True, cmap=cmap)
        if headtail:
            ax.plot(X[:1], Y[:1], Z[:1], 'k*', markersize=12)
            ax.plot(X[-1:], Y[-1:], Z[-1:], 'ko', markersize=12)

    if SETTINGS['auto_show']:
        showFigure()
    return f

def getDomainList(labels):
    """Returns a list of domain separations. The list has two columns: the first is for 
    the domain starts and the second is for the domain ends."""

    indicators = np.diff(labels)
    indicators = np.append(1., indicators)
    indicators[-1] = 1
    sites = np.where(indicators != 0)[0]
    starts = sites[:-1]
    ends = sites[1:]
    domains = np.array([starts, ends]).T

    return domains