"""This module defines miscellaneous utility functions."""

import numpy as np

from numpy import unique, linalg, diag, sqrt, dot
import scipy.cluster.hierarchy as sch
from scipy import spatial
from .misctools import addBreaks

__all__ = ['calcTree', 'clusterMatrix', 'showData', 'reorderMatrix', 'findSubgroups']

def calcTree(names, distance_matrix, method='nj'):
    """ Given a distance matrix for an ensemble, it creates an returns a tree structure.
    :arg names: an list of names. 
    :type names: list-like
    :arg distance_matrix: a square matrix with length of ensemble. If numbers does not mismatch
    it will raise an error. 
    :type distance_matrix: numpy.ndarray 
    """
    try: 
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')
    
    if len(names) != distance_matrix.shape[0] or len(names) != distance_matrix.shape[1]:
        raise ValueError("Mismatch between the sizes of matrix and names.")

    matrix = []
    k = 1
    for row in distance_matrix:
        matrix.append(list(row[:k]))
        k = k + 1
    from Bio.Phylo.TreeConstruction import _DistanceMatrix
    dm = _DistanceMatrix(names, matrix)
    constructor = Phylo.TreeConstruction.DistanceTreeConstructor()

    method = method.strip().lower()
    if method == 'nj':
        tree = constructor.nj(dm)
    elif method == 'upgma':
        tree = constructor.upgma(dm)
    else:
        raise ValueError('Method can be only either "nj" or "upgma".')

    for node in tree.get_nonterminals():
        node.name = None
    return tree

def clusterMatrix(similarity_matrix=None, distance_matrix=None, labels=None, no_plot=True, **kwargs):
    """
    Cluster a similarity matrix or a distance matrix using scipy.cluster.hierarchy and 
    return the linkage matrix, indices, sorted matrix, sorted labels, and 
    dendrogram dict.
    
    :arg similarity_matrix: an N-by-N matrix containing some measure of similarity 
        such as sequence identity, mode-mode overlap, or spectral overlap
    :type similarity_matrix: array
    
    :arg distance_matrix: an N-by-N matrix containing some measure of distance 
        such as 1. - seqid_matrix, rmsds, or distances in PCA space
    :type similarity_matrix: array
    
    :arg labels: labels for each matrix row that can be returned sorted
    :type labels: list
    
    Other arguments for scipy.hierarchy.linkage and scipy.hierarchy.dendrogram
        can also be provided and will be taken as kwargs.
        
    :arg no_plot: If True, don't plot the dendrogram.
        default is True
    :type no_plot: bool
    """
    if similarity_matrix is None and distance_matrix is None:
        raise ValueError('Please provide a similarity matrix or a distance matrix')
    elif distance_matrix is None:
        distance_matrix = 1. - similarity_matrix
    
    orientation = kwargs.pop('orientiation','right')
    
    formatted_distance_matrix = spatial.distance.squareform(distance_matrix)
    linkage_matrix = sch.linkage(formatted_distance_matrix, **kwargs)
    sorting_dendrogram = sch.dendrogram(linkage_matrix, orientation=orientation, labels=labels, no_plot=no_plot)

    indices = sorting_dendrogram['leaves']
    sorted_labels = sorting_dendrogram['ivl']
    
    if similarity_matrix is None:
        sorted_matrix = distance_matrix[indices,:]
    else:
        sorted_matrix = similarity_matrix[indices,:]
    sorted_matrix = sorted_matrix[:,indices]
    
    return linkage_matrix, indices, sorted_matrix, sorted_labels, sorting_dendrogram

def showData(*args, **kwargs):
    """
    Show data using :func:`~matplotlib.axes.Axes.plot`. 
    
    :arg x: (optional) x coordinates. *x* can be an 1-D array or a 2-D matrix of 
    column vectors.
    :type x: `~numpy.ndarray`

    :arg y: data array. *y* can be an 1-D array or a 2-D matrix of 
    column vectors.
    :type y: `~numpy.ndarray`

    :arg dy: an array of variances of *y* which will be plotted as a 
    band along *y*. It should have the same shape with *y*.
    :type dy: `~numpy.ndarray`

    :arg alpha: the transparency of the band(s).
    :type alpha: float

    :arg ticklabels: user-defined tick labels for x-axis.
    :type ticklabels: list
    """
    
    # note for developers: this function serves as a low-level 
    # plotting function which provides basic utilities for other 
    # plotting functions. Therefore showFigure and new_fig are 
    # not handled in this function as it should be already handled in 
    # the caller.

    ticklabels = kwargs.pop('ticklabels', None)
    dy = kwargs.pop('dy', None)
    alpha = kwargs.pop('alpha', 0.5)
    gap = kwargs.pop('gap', False)

    from matplotlib import cm, ticker
    from matplotlib.pyplot import figure, gca, xlim

    ax = gca()
    lines = ax.plot(*args, **kwargs)

    if dy is not None:
        dy = np.array(dy)
        if dy.ndim == 1:
            n, = dy.shape; m = 1
        elif dy.ndim == 2:
            n, m = dy.shape
        else:
            raise ValueError('dy should be either 1-D or 2-D.')
        
        for i, line in enumerate(lines):
            color = line.get_color()
            x, y = line.get_data()
            if m != 1 and m != len(lines) or n != len(y):
                raise ValueError('The shapes of dy and y do not match.')

            if dy.ndim == 1:
                _dy = dy
            else:
                _dy = dy[:, i]
            
            if gap:
                x_new, y_new = addBreaks(x, y)
                line.set_data(x_new, y_new)
                _, _dy = addBreaks(x, _dy)
            else:
                x_new, y_new = x, y

            ax.fill_between(x_new, y_new-_dy, y_new+_dy,
                    alpha=alpha, facecolor=color,
                    linewidth=1, antialiased=True)

    if ticklabels is not None:
        ax.get_xaxis().set_major_formatter(ticker.IndexFormatter(ticklabels))
    return ax

def reorderMatrix(names, matrix, tree):
    """
    Reorder a matrix based on a tree and return the reordered matrix 
    and indices for reordering other things.

    :arg names: a list of names associated with the rows of the matrix
        These names must match the ones used to generate the tree.
    :type names: a list of strings

    :arg matrix: any square matrix
    :type matrix: 2D array

    :arg tree: any tree from calcTree
    :type tree: Bio.Phylo.BaseTree.Tree
    """
    try:
        from Bio import Phylo
    except ImportError:
        raise ImportError('Phylo module could not be imported. '
            'Reinstall ProDy or install Biopython '
            'to solve the problem.')

    if type(names) is not list:
        raise TypeError('names should be a list.')

    if type(names[0]) is not str:
        raise TypeError('names should be a list of strings.')    

    if type(matrix) is not np.ndarray:
        raise TypeError('matrix should be a numpy array.')

    if matrix.ndim != 2:
        raise ValueError('matrix should be a 2D matrix.')

    if np.shape(matrix)[0] != np.shape(matrix)[1]:
        raise ValueError('matrix should be a square matrix')

    if type(tree) is not Phylo.BaseTree.Tree:
        raise TypeError('tree should be a BioPython Tree')

    if len(names) != len(matrix):
        raise ValueError('names should have entries for each matrix row/column')

    if len(names) != len(tree.get_terminals()):
        raise ValueError('names should have entries for each tree terminal')

    if len(tree.get_terminals()) != len(matrix):
        raise ValueError('matrix should have a row for each tree terminal')

    indices = []
    for terminal in tree.get_terminals():
        indices.append(np.where(np.array(names) == str(terminal))[0][0])

    reordered_matrix = matrix[:,indices]
    reordered_matrix = reordered_matrix[indices,:]
    
    return reordered_matrix, indices

def findSubgroups(tree, cutoff=0.8):
    """
    Divide a tree into subgroups using a distance cutoff.
    Returns a list of lists with labels divided into subgroups.
    """

    subgroups = [[]]

    for i, target_i in enumerate(tree.get_terminals()):
        subgroups[-1].append(str(target_i))
        for j, target_j in enumerate(tree.get_terminals()[:-1]):
            if i == j+1:
                neighbour_distance = tree.distance(target_i,target_j)
                if neighbour_distance > cutoff:
                    subgroups.append([])

    subgroups[-1].append(str(target_i))

    return subgroups
