"""This module defines miscellaneous utility functions that is public to users."""

import numpy as np

from prody import PY3K
from .misctools import addEnds, interpY, index, isListLike
from .checkers import checkCoords
from .logger import LOGGER


__all__ = ['calcTree', 'clusterMatrix', 'showLines', 'showMatrix', 
           'reorderMatrix', 'findSubgroups', 'getCoords',  
           'getLinkage', 'getTreeFromLinkage', 'clusterSubfamilies']

class LinkageError(Exception):
    pass

def clusterSubfamilies(similarities, n_clusters=0, linkage='all', method='tsne', cutoff=0.0, **kwargs):
    """Perform clustering based on members of the *ensemble* projected into lower a reduced
    dimension.
    
    :arg similarities: a matrix of similarities for each structure in the ensemble, such as
                        RMSD-matrix, dynamics-based spectral overlap, sequence similarity
    :type similarities: :class:`~numpy.ndarray`

    :arg n_clusters: the number of clusters to generate. If **0**, will scan a range of 
                        number of clusters and return the best one based on highest
                        silhouette score. Default is **0**.
    :type n_clusters: int

    :arg linkage: if **all**, will test all linkage types (ward, average, complete,
                    single). Otherwise will use only the one(s) given as input. Default is
                    **all**.
    :type linkage: str, list, tuple, :class:`~numpy.ndarray`

    :arg method: if set to **spectral**, will generate a Kirchoff matrix based on the 
                    cutoff value given and use that as input as clustering instead of
                    the values themselves. Default is **tsne**.
    :type method: str

    :arg cutoff: only used if *method* is set to **spectral**. This value is used for 
                    generating the Kirchoff matrix to use for generating clusters when
                    doing spectral clustering. Default is **0.0**.
    :type cutoff: float
    """

    # Import necessary packages
    try:
        from sklearn.manifold import SpectralEmbedding
        from sklearn.cluster import AgglomerativeClustering
        from sklearn.metrics import silhouette_score
        from sklearn.manifold import TSNE
    except ImportError:
        raise ImportError('need sklearn module')
        '''
        try: 
            import Bio 
        except ImportError:
            raise ImportError('Phylo module could not be imported. '
                'Reinstall ProDy or install Biopython '
                'to solve the problem.')
        '''
        

    # Check inputs to make sure are of valid types/values
    if not isinstance(similarities, np.ndarray):
        raise TypeError('similarities should be a numpy ndarray')

    dim = similarities.shape
    if dim[0] != dim[1]:
        raise ValueError('similarities must be a square matrix')

    if n_clusters != 0:
        if not isinstance(n_clusters, int):
            raise TypeError('clusters must be an instance of int')
        if n_clusters < 1:
            raise ValueError('clusters must be a positive integer')
        elif n_clusters > similarities.shape[0]:
            raise ValueError('clusters can\'t be longer than similarities matrix')
        nclusts = range(n_clusters,n_clusters+1)
    else:
        nclusts = range(2,10,1)

    if linkage != 'all':
        # Check if given input for linkage is list-like
        if isListLike(linkage):
            for val in linkage:
                if val.lower() not in ['ward', 'average', 'complete', 'single']:
                    raise ValueError('linkage must be one or more of: \'ward\', \'average\', \'complete\', or \'single\'')
            if len(linkage) > 4:
                raise ValueError('linkage must be one or more of: \'ward\', \'average\', \'complete\', or \'single\'')
            linkages = [ x.lower() for x in linkage ]

        # If not, check if it is a valid string and method name
        else:
            if not isinstance(linkage, str):
                raise TypeError('linkage must be an instance of str or list-like of strs')

            if linkage not in ['ward', 'average', 'complete', 'single']:
                raise ValueError('linkage must one or more of: \'ward\', \'average\', \'complete\', or \'single\'')

            linkages = [linkage]
    else:
        linkages = ['ward', 'average', 'complete', 'single']

    if method != 'tsne':
        if not isinstance(method, str):
            raise TypeError('method must be an instance of str')
        if method != 'spectral':
            raise ValueError('method must be either \'tsne\' or \'spectral\'')

        if not isinstance(cutoff, float):
            raise TypeError('cutoff must be an instance of float')

    best_score = -1
    best_nclust = 0
    best_link = ''
    best_labels = []

    # Scan over range of clusters
    for x in nclusts:
        if method == 'tsne':
            embedding = TSNE(n_components=2)
            transform = embedding.fit_transform(similarities)

        else:
            kirchhoff = np.where(similarities > cutoff, 0, -1)
            embedding = SpectralEmbedding(n_components=2)
            transform = embedding.fit_transform(kirchhoff)

        for link in linkages:
            clustering = AgglomerativeClustering(linkage=link, n_clusters=x)
            clustering.fit(transform)

            silhouette_avg = silhouette_score(transform, clustering.labels_)
            
            if silhouette_avg > best_score:
                best_score = silhouette_avg
                best_nclust = x
                best_link = link
                best_labels = clustering.labels_


    return best_labels

def getCoords(data):

    try:
        data = (data._getCoords() if hasattr(data, '_getCoords') else
                data.getCoords())
    except AttributeError:
        try:
            checkCoords(data)
        except TypeError:
            raise TypeError('data must be a Numpy array or an object '
                            'with `getCoords` method')

    return data

def getLinkage(names, tree):
    """ Obtain the :func:`~scipy.cluster.hierarchy.linkage` matrix encoding 
    ``tree``. 
    
    :arg names: a list of names, the order determines the values in the 
                linkage matrix
    :type names: list, :class:`~numpy.ndarray`

    :arg tree: tree to be converted
    :type tree: :class:`~Bio.Phylo.BaseTree.Tree`
    """

    tree_terminals = tree.get_terminals()

    if len(tree_terminals) != len(names):
        raise ValueError('inconsistent number of terminals in tree and names')
    
    terminals = [None] * len(names)
    for clade in tree_terminals:
        i = index(names, clade.name)
        terminals[i] = clade

    n = len(terminals)
    nonterminals = [c for c in reversed(tree.get_nonterminals())]
    if len(nonterminals) != n-1:
        raise LinkageError('wrong number of terminal clades')

    Z = np.zeros((n-1, 4))

    root = tree.root

    def _indexOfClade(clade):
        if clade.is_terminal():
            i = index(terminals, clade)
        else:
            i = index(nonterminals, clade) + n
        return i

    def _height_of(clade):
        if clade.is_terminal():
            height = 0 
        else:
            height = max(_height_of(c) + c.branch_length for c in clade.clades)

        return height

    def _dfs(clade):
        if clade.is_terminal():
            return

        i = _indexOfClade(clade)
        clade_a = clade.clades[0]
        clade_b = clade.clades[1]

        a = _indexOfClade(clade_a)
        b = _indexOfClade(clade_b) 

        l = min(a, b)
        r = max(a, b)

        Z[i-n, 0] = l
        Z[i-n, 1] = r
        Z[i-n, 2] = _height_of(clade) * 2.
        Z[i-n, 3] = clade.count_terminals()

        _dfs(clade_a)
        _dfs(clade_b)
    
    _dfs(root)

    return Z

def getTreeFromLinkage(names, linkage):
    """ Obtain the tree encoded by ``linkage``. 
    
    :arg names: a list of names, the order should correspond to the values in  
                linkage
    :type names: list, :class:`~numpy.ndarray`

    :arg linkage: linkage matrix
    :type linkage: :class:`~numpy.ndarray`
    """
    try: 
        import Bio 
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')

    from Bio.Phylo.BaseTree import Tree, Clade
    
    if not isinstance(linkage, np.ndarray):
        raise TypeError('linkage must be a numpy.ndarray instance')

    if linkage.ndim != 2:
        raise LinkageError('linkage must be a 2-dimensional matrix')

    if linkage.shape[1] != 4:
        raise LinkageError('linkage must have exactly 4 columns')

    n_terms = len(names)
    if linkage.shape[0] != n_terms-1:
        raise LinkageError('linkage must have exactly len(names)-1 rows')
    
    clades = []
    heights = []
    for name in names:
        clade = Clade(None, name)
        clades.append(clade)
        heights.append(0.)

    for link in linkage:
        l = int(link[0])
        r = int(link[1])
        height = link[2]

        left = clades[l]
        right = clades[r]

        lh = heights[l]
        rh = heights[r]

        left.branch_length = height - lh
        right.branch_length = height - rh

        clade = Clade(None, None)
        clade.clades.append(left)
        clade.clades.append(right)

        clades.append(clade)
        heights.append(height)

    return Tree(clade)

def calcTree(names, distance_matrix, method='upgma', linkage=False):
    """ Given a distance matrix, it creates an returns a tree structure.

    :arg names: a list of names
    :type names: list, :class:`~numpy.ndarray`

    :arg distance_matrix: a square matrix with length of ensemble. If numbers does not match *names*
                          it will raise an error
    :type distance_matrix: :class:`~numpy.ndarray`

    :arg method: method used for constructing the tree. Acceptable options are ``"upgma"``, ``"nj"``, 
                 or methods supported by :func:`~scipy.cluster.hierarchy.linkage` such as ``"single"``, 
                 ``"average"``, ``"ward"``, etc. Default is ``"upgma"``
    :type method: str

    :arg linkage: whether the linkage matrix is returned. Note that NJ trees do not support linkage
    :type linkage: bool
    """
    try: 
        import Bio 
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')
            
    from .TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    
    if len(names) != distance_matrix.shape[0] or len(names) != distance_matrix.shape[1]:
        raise ValueError("Mismatch between the sizes of matrix and names.")
    
    method = method.lower().strip()

    if method in ['ward', 'single', 'average', 'weighted', 'centroid', 'median']:
        from scipy.cluster.hierarchy import linkage as hlinkage
        from scipy.spatial.distance import squareform
        
        Z = hlinkage(squareform(distance_matrix), method=method)
        tree = getTreeFromLinkage(names, Z)
    else:
        matrix = []
        k = 1
        Z = None
        for row in distance_matrix:
            matrix.append(list(row[:k]))
            k = k + 1
        
        if isinstance(names, np.ndarray):
            names = names.tolist()
        dm = DistanceMatrix(names, matrix)
        constructor = DistanceTreeConstructor()

        method = method.strip().lower()
        if method == 'nj':
            tree = constructor.nj(dm)
        elif method == 'upgma':
            tree = constructor.upgma(dm)
            if linkage:
                Z = getLinkage(names, tree)
        else:
            raise ValueError('Method can be only either "nj", "upgma" or '
                             'hierarchical clustering such as "single", "average", etc.')

        for node in tree.get_nonterminals():
            node.name = None

    if linkage:
        return tree, Z
    else:
        return tree

def writeTree(filename, tree, format_str='newick'):
    """ Write a tree to file using Biopython.

    :arg filename: name for output file
    :type filename: str

    :arg tree: a square matrix with length of ensemble. If numbers does not match *names*
                          it will raise an error
    :type tree: :class:`~Bio.Phylo.BaseTree.Tree`

    :arg format_str: a string specifying the format for the tree
    :type format_str: str
    """
    try: 
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')

    if not isinstance(filename, str):
        raise TypeError('filename should be a string')

    if not isinstance(tree, Phylo.BaseTree.Tree):
        raise TypeError('tree should be a Biopython.Phylo Tree object')

    if not isinstance(format_str, str):
        raise TypeError('format_str should be a string')

    Phylo.write(tree, filename, format_str)


def clusterMatrix(distance_matrix=None, similarity_matrix=None, labels=None, return_linkage=None, **kwargs):
    """
    Cluster a distance matrix using scipy.cluster.hierarchy and 
    return the sorted matrix, indices used for sorting, sorted labels (if **labels** are passed),  
    and linkage matrix (if **return_linkage** is **True**). 
    
    :arg distance_matrix: an N-by-N matrix containing some measure of distance 
         such as 1. - seqid_matrix (Hamming distance), rmsds, or distances in PCA space
    :type distance_matrix: :class:`~numpy.ndarray`

    :arg similarity_matrix: an N-by-N matrix containing some measure of similarity 
         such as sequence identity, mode-mode overlap, or spectral overlap.
         Each element will be subtracted from 1. to get distance, so make sure this 
         is reasonable.
    :type similarity_matrix: :class:`~numpy.ndarray`
    
    :arg labels: labels for each matrix row that can be returned sorted
    :type labels: list

    :arg no_plot: if **True**, don't plot the dendrogram.
         default is **True**
    :type no_plot: bool
    
    :arg reversed: if set to **True**, then the sorting indices will be reversed.
    :type reversed: bool

    Other arguments for :func:`~scipy.hierarchy.linkage` and :func:`~scipy.hierarchy.dendrogram`
    can also be provided and will be taken as **kwargs**.
    """

    import scipy.cluster.hierarchy as sch
    from scipy import spatial
    if similarity_matrix is None and distance_matrix is None:
        raise ValueError('Please provide a distance matrix or a similarity matrix')
    
    orientation = kwargs.pop('orientiation', 'right')
    reversed = kwargs.pop('reversed', False)
    no_plot = kwargs.pop('no_plot', True)

    if distance_matrix is None:
        matrix = similarity_matrix
        distance_matrix = 1. - similarity_matrix
    else:
        matrix = distance_matrix
        
    formatted_distance_matrix = spatial.distance.squareform(distance_matrix)
    linkage_matrix = sch.linkage(formatted_distance_matrix, **kwargs)
    sorting_dendrogram = sch.dendrogram(linkage_matrix, orientation=orientation, labels=labels, no_plot=no_plot)

    indices = sorting_dendrogram['leaves']
    sorted_labels = sorting_dendrogram['ivl']

    if reversed:
        indices = indices[::-1]
        sorted_labels = sorted_labels[::-1]
    
    sorted_matrix = matrix[indices, :]
    sorted_matrix = sorted_matrix[:, indices]
    
    return_vals = [sorted_matrix, indices]

    if labels is not None:
        return_vals.append(sorted_labels)
    if return_linkage:
        return_vals.append(linkage_matrix)
    return tuple(return_vals) # convert to tuple to avoid [pylint] E0632:Possible unbalanced tuple unpacking

def showLines(*args, **kwargs):
    """
    Show 1-D data using :func:`~matplotlib.axes.Axes.plot`. 
    
    :arg x: (optional) x coordinates. *x* can be an 1-D array or a 2-D matrix of 
            column vectors.
    :type x: :class:`~numpy.ndarray`

    :arg y: data array. *y* can be an 1-D array or a 2-D matrix of 
            column vectors.
    :type y: :class:`~numpy.ndarray`

    :arg dy: an array of variances of *y* which will be plotted as a 
             band along *y*. It should have the same shape with *y*.
    :type dy: :class:`~numpy.ndarray`

    :arg lower: an array of lower bounds which will be plotted as a 
                band along *y*. It should have the same shape with *y* and should be 
                paired with *upper*.
    :type lower: :class:`~numpy.ndarray`

    :arg upper: an array of upper bounds which will be plotted as a 
                band along *y*. It should have the same shape with *y* and should be 
                paired with *lower*.
    :type upper: :class:`~numpy.ndarray`

    :arg alpha: the transparency of the band(s) for plotting *dy*.
    :type alpha: float

    :arg beta: the transparency of the band(s) for plotting *miny* and *maxy*.
    :type beta: float

    :arg ticklabels: user-defined tick labels for x-axis.
    :type ticklabels: list
    """
    
    # note for developers: this function serves as a low-level 
    # plotting function which provides basic utilities for other 
    # plotting functions. Therefore showFigure is not handled 
    # in this function as it should be already handled in the caller.

    ticklabels = kwargs.pop('ticklabels', None)
    dy = kwargs.pop('dy', None)
    miny = kwargs.pop('lower', None)
    maxy = kwargs.pop('upper', None)
    alpha = kwargs.pop('alpha', 0.5)
    beta = kwargs.pop('beta', 0.25)
    gap = kwargs.pop('gap', False)
    labels = kwargs.pop('label', None)

    from matplotlib import cm, ticker
    from matplotlib.pyplot import figure, gca, xlim
    from .drawtools import IndexFormatter

    ax = gca()
    lines = ax.plot(*args, **kwargs)

    polys = []
        
    for i, line in enumerate(lines):
        color = line.get_color()
        x, y = line.get_data()
        
        if gap:
            x_new, y_new = addEnds(x, y)
            line.set_data(x_new, y_new)
        else:
            x_new, y_new = x, y
        
        if labels is not None:
            if np.isscalar(labels):
                line.set_label(labels)
            else:
                try:
                    line.set_label(labels[i])
                except IndexError:
                    raise ValueError('The number of labels ({0}) and that of y ({1}) do not match.'
                                     .format(len(labels), len(line)))
        
        # the following function needs to be here so that line exists
        def sub_array(a, i, tag='a'):
            ndim = 0
            if a is not None:
                if np.isscalar(a[0]):
                    ndim = 1   # a plain list (array)
                else:
                    ndim = 2   # a nested list (array)
            else:
                return None

            if ndim == 1:
                _a = a
            else:
                try:
                    _a = a[i]
                except IndexError:
                    raise ValueError('The number of {2} ({0}) and that of y ({1}) do not match.'
                                     .format(len(miny), len(line), tag))

            if len(_a) != len(y):
                raise ValueError('The shapes of {2} ({0}) and y ({1}) do not match.'
                                 .format(len(_miny), len(y), tag))
            return _a

        if miny is not None and maxy is not None:
            _miny = sub_array(miny, i)
            _maxy = sub_array(maxy, i)

            if gap:
                _, _miny = addEnds(x, _miny)
                _, _maxy = addEnds(x, _maxy)
                
            poly = ax.fill_between(x_new, _miny, _maxy,
                                    alpha=beta, facecolor=color, edgecolor=None,
                                    linewidth=1, antialiased=True)
            polys.append(poly)

        if dy is not None:
            _dy = sub_array(dy, i)

            if gap:
                _, _dy = addEnds(x, _dy)
                
            poly = ax.fill_between(x_new, y_new-_dy, y_new+_dy,
                                    alpha=alpha, facecolor=color, edgecolor=None,
                                    linewidth=1, antialiased=True)
            polys.append(poly)

    ax.margins(x=0)
    if ticklabels is not None:
        if callable(ticklabels):
            ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(ticklabels))
        else:
            ax.get_xaxis().set_major_formatter(IndexFormatter(ticklabels))
    
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    return lines, polys

def showMatrix(matrix, x_array=None, y_array=None, **kwargs):
    """Show a matrix using :meth:`~matplotlib.axes.Axes.imshow`. Curves on x- and y-axis can be added.

    :arg matrix: matrix to be displayed
    :type matrix: :class:`~numpy.ndarray`

    :arg x_array: data to be plotted above the matrix
    :type x_array: :class:`~numpy.ndarray`

    :arg y_array: data to be plotted on the left side of the matrix
    :type y_array: :class:`~numpy.ndarray`

    :arg percentile: a percentile threshold to remove outliers, i.e. only showing data within *p*-th 
                     to *100-p*-th percentile
    :type percentile: float

    :arg interactive: turn on or off the interactive options
    :type interactive: bool

    :arg xtickrotation: how much to rotate the xticklabels in degrees
                        default is 0
    :type xtickrotation: float
    """

    from matplotlib import ticker
    from matplotlib.gridspec import GridSpec
    from matplotlib.collections import LineCollection
    from matplotlib.pyplot import gca, sca, sci, colorbar, subplot

    from .drawtools import drawTree, IndexFormatter

    p = kwargs.pop('percentile', None)
    vmin = vmax = None
    if p is not None:
        vmin = np.percentile(matrix, p)
        vmax = np.percentile(matrix, 100-p)
    
    vmin = kwargs.pop('vmin', vmin)
    vmax = kwargs.pop('vmax', vmax)
    vcenter = kwargs.pop('vcenter', None)
    norm = kwargs.pop('norm', None)

    if vcenter is not None and norm is None:
        if PY3K:
            try:
                from matplotlib.colors import DivergingNorm
            except ImportError:
                from matplotlib.colors import TwoSlopeNorm as DivergingNorm

            norm = DivergingNorm(vmin=vmin, vcenter=0., vmax=vmax)
        else:
            LOGGER.warn('vcenter cannot be used in Python 2 so norm remains None')

    lw = kwargs.pop('linewidth', 1)
    
    W = H = kwargs.pop('ratio', 6)

    ticklabels = kwargs.pop('ticklabels', None)
    xticklabels = kwargs.pop('xticklabels', ticklabels)
    yticklabels = kwargs.pop('yticklabels', ticklabels)

    xtickrotation = kwargs.pop('xtickrotation', 0.)

    show_colorbar = kwargs.pop('colorbar', True)
    cb_extend = kwargs.pop('cb_extend', 'neither')
    allticks = kwargs.pop('allticks', False) # this argument is temporary and will be replaced by better implementation
    interactive = kwargs.pop('interactive', True)

    cmap = kwargs.pop('cmap', 'jet')
    origin = kwargs.pop('origin', 'lower')

    try: 
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')
    tree_mode_y = isinstance(y_array, Phylo.BaseTree.Tree)
    tree_mode_x = isinstance(x_array, Phylo.BaseTree.Tree)

    if x_array is not None and y_array is not None:
        nrow = 2; ncol = 2
        i = 1; j = 1
        width_ratios = [1, W]
        height_ratios = [1, H]
        aspect = 'auto'
    elif x_array is not None and y_array is None:
        nrow = 2; ncol = 1
        i = 1; j = 0
        width_ratios = [W]
        height_ratios = [1, H]
        aspect = 'auto'
    elif x_array is None and y_array is not None:
        nrow = 1; ncol = 2
        i = 0; j = 1
        width_ratios = [1, W]
        height_ratios = [H]
        aspect = 'auto'
    else:
        nrow = 1; ncol = 1
        i = 0; j = 0
        width_ratios = [W]
        height_ratios = [H]
        aspect = kwargs.pop('aspect', None)

    main_index = (i, j)
    upper_index = (i-1, j)
    left_index = (i, j-1)

    complex_layout = nrow > 1 or ncol > 1

    ax1 = ax2 = ax3 = None

    if complex_layout:
        gs = GridSpec(nrow, ncol, width_ratios=width_ratios, 
                      height_ratios=height_ratios, hspace=0., wspace=0.)

    ## draw matrix
    if complex_layout:
        ax3 = subplot(gs[main_index])
    else:
        ax3 = gca()
    
    im = ax3.imshow(matrix, aspect=aspect, vmin=vmin, vmax=vmax, 
                    norm=norm, cmap=cmap, origin=origin, **kwargs)
                    
    #ax3.set_xlim([-0.5, matrix.shape[0]+0.5])
    #ax3.set_ylim([-0.5, matrix.shape[1]+0.5])

    if xticklabels is not None:
        ax3.xaxis.set_major_formatter(IndexFormatter(xticklabels))
    if yticklabels is not None and ncol == 1:
        ax3.yaxis.set_major_formatter(IndexFormatter(yticklabels))

    if allticks:
        ax3.xaxis.set_major_locator(ticker.IndexLocator(offset=0.5, base=1.))
        ax3.yaxis.set_major_locator(ticker.IndexLocator(offset=0.5, base=1.))
    else:
        locator = ticker.AutoLocator()
        locator.set_params(integer=True)
        minor_locator = ticker.AutoMinorLocator()

        ax3.xaxis.set_major_locator(locator)
        ax3.xaxis.set_minor_locator(minor_locator)

        locator = ticker.AutoLocator()
        locator.set_params(integer=True)
        minor_locator = ticker.AutoMinorLocator()
        
        ax3.yaxis.set_major_locator(locator)
        ax3.yaxis.set_minor_locator(minor_locator)

    if ncol > 1:
        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
    
    ## draw x_ and y_array
    lines = []

    if nrow > 1:
        ax1 = subplot(gs[upper_index])

        if tree_mode_x:
            Y, X = drawTree(x_array, label_func=None, orientation='vertical', 
                            inverted=True)
            miny = min(Y.values())
            maxy = max(Y.values())

            minx = min(X.values())
            maxx = max(X.values())

            ax1.set_xlim(minx-.5, maxx+.5)
            ax1.set_ylim(miny, 1.05*maxy)
        else:
            ax1.set_xticklabels([])
            
            y = x_array
            xp, yp = interpY(y)
            points = np.array([xp, yp]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcy = LineCollection(segments, array=yp, linewidths=lw, cmap=cmap)
            lines.append(lcy)
            ax1.add_collection(lcy)

            ax1.set_xlim(xp.min()-.5, xp.max()+.5)
            ax1.set_ylim(yp.min(), yp.max())

        if ax3.xaxis_inverted():
            ax2.invert_xaxis()

        ax1.axis('off')

    if ncol > 1:
        ax2 = subplot(gs[left_index])
        
        if tree_mode_y:
            X, Y = drawTree(y_array, label_func=None, inverted=True)
            miny = min(Y.values())
            maxy = max(Y.values())

            minx = min(X.values())
            maxx = max(X.values())

            ax2.set_ylim(miny-.5, maxy+.5)
            ax2.set_xlim(minx, 1.05*maxx)
        else:
            ax2.set_xticklabels([])
            
            y = y_array
            xp, yp = interpY(y)
            points = np.array([yp, xp]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lcx = LineCollection(segments, array=yp, linewidths=lw, cmap=cmap)
            lines.append(lcx)
            ax2.add_collection(lcx)
            ax2.set_xlim(yp.min(), yp.max())
            ax2.set_ylim(xp.min()-.5, xp.max()+.5)
        
        ax2.invert_xaxis()

        if ax3.yaxis_inverted():
            ax2.invert_yaxis()

        ax2.axis('off')

    ## draw colorbar
    sca(ax3)
    cb = None
    if show_colorbar:
        if nrow > 1:
            axes = [ax1, ax2, ax3]
            while None in axes:
                axes.remove(None)
            s = H / (H + 1.)
            cb = colorbar(mappable=im, ax=axes, anchor=(0, 0), shrink=s, extend=cb_extend)
        else:
            cb = colorbar(mappable=im, extend=cb_extend)

    sca(ax3)
    sci(im)

    if interactive:
        from prody.utilities import ImageCursor
        from matplotlib.pyplot import connect
        cursor = ImageCursor(ax3, im)
        connect('button_press_event', cursor.onClick)

    ax3.tick_params(axis='x', rotation=xtickrotation)

    return im, lines, cb

def reorderMatrix(names, matrix, tree, axis=None):
    """
    Reorder a matrix based on a tree and return the reordered matrix 
    and indices for reordering other things.

    :arg names: a list of names associated with the rows of the matrix
        These names must match the ones used to generate the tree
    :type names: list

    :arg matrix: any square matrix
    :type matrix: :class:`~numpy.ndarray`

    :arg tree: any tree from :func:`calcTree`
    :type tree: :class:`~Bio.Phylo.BaseTree.Tree`

    :arg axis: along which axis the matrix should be reordered. 
               Default is **None** which reorder along all the axes
    :type axis: int
    """

    try:
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')

    try:
        if matrix.ndim != 2:
            raise ValueError('matrix should be a 2D matrix.')
    except AttributeError:
        raise TypeError('matrix should be a numpy array.')

    if np.shape(matrix)[0] != np.shape(matrix)[1]:
        raise ValueError('matrix should be a square matrix')
    
    names = np.asarray(names)

    if np.isscalar(names):
        raise TypeError('names should be list-like')
    
    if not len(names):
        raise TypeError('names is empty')

    if not isinstance(tree, Phylo.BaseTree.Tree):
        raise TypeError('tree should be a BioPython Tree')

    if len(names) != len(matrix):
        raise ValueError('names should have entries for each matrix row/column')
    
    terminals = tree.get_terminals()
    if len(names) != len(terminals):
        raise ValueError('names should have entries for each tree terminal')

    if len(terminals) != len(matrix):
        raise ValueError('matrix should have a row for each tree terminal')

    indices = []
    for terminal in terminals:
        name = terminal.name
        locs = np.where(names == name)[0]
        if not len(locs):
            raise ValueError('inconsistent names and tree: %s not in names'%name)

        if len(locs) > 1:
            raise ValueError('inconsistent names and tree: duplicate name %s in names'%name)
        indices.append(locs[0])

    if axis is not None:
        I = [np.arange(s) for s in matrix.shape] 
        axes = [axis] if np.isscalar(axis) else axis
        for ax in axes:
            I[ax] = indices
    else:
        I = [indices] * matrix.ndim
    
    rmatrix = matrix[np.ix_(*I)]
    
    return rmatrix, indices

def findSubgroups(tree, c, method='naive', **kwargs):
    """
    Divide **tree** into subgroups using a criterion **method** and a cutoff **c**.
    Returns a list of lists with labels divided into subgroups.
    """

    method = method.lower().strip()
    terminals = tree.get_terminals()
    names = [clade.name for clade in terminals]
    Z = None

    if method != 'naive':
        try:
            Z = getLinkage(names, tree)
        except LinkageError:
            print('Failed to build linkage; fall back to naive criterion')
            method = 'naive'
    
    if method == 'naive':
        subgroups = [[names[0]]]
        for i in range(len(terminals)-1):
            curr_clade = terminals[i]
            next_clade = terminals[i + 1]
            d = tree.distance(curr_clade, next_clade)
            if d > c:
                subgroups.append([])
            subgroups[-1].append(next_clade.name)
    else:
        from scipy.cluster.hierarchy import fcluster
        
        T = fcluster(Z, c, criterion=method, **kwargs)
        labels = np.unique(T)
        subgroups = [[] for _ in range(len(labels))]

        for i, t in enumerate(T):
            subgroups[t-1].append(names[i])

    return subgroups
