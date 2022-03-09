# -*- coding: utf-8 -*-
"""This module defines functions for comparing normal modes from different
models."""

import numpy as np
from numbers import Integral
from prody import LOGGER, SETTINGS
from prody.utilities import openFile, isListLike

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .gnm import ZERO
from .analysis import calcFractVariance, calcSqFlucts

__all__ = ['calcOverlap', 'calcCumulOverlap', 'calcSubspaceOverlap', 'calcSpectralOverlap', 
           'calcCovOverlap', 'printOverlapTable', 'writeOverlapTable', 
           'calcSquareInnerProduct','pairModes', 'matchModes', 
           'calcRMSIP', 'calcSIP', 'calcRWSIP']

SO_CACHE = {}
WO_CACHE = {}

def calcOverlap(rows, cols, diag=False):
    """Returns overlap (or correlation) between two sets of modes (*rows* and
    *cols*).  Returns a matrix whose rows correspond to modes passed as *rows*
    argument, and columns correspond to those passed as *cols* argument.
    Both rows and columns are normalized prior to calculating overlap.
    
    This function can now return the diagonal of the overlap matrix if *diag*
    is set to **True**."""

    if not isinstance(rows, (NMA, ModeSet, Mode, Vector, np.ndarray)):
        raise TypeError('rows must be NMA, ModeSet, Mode, Vector, or array, not {0}'
                        .format(type(rows)))
    if not isinstance(cols, (NMA, ModeSet, Mode, Vector, np.ndarray)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, or array, not {0}'
                        .format(type(cols)))

    if isinstance(rows, np.ndarray):
        num_rows = rows.shape[0]
    else:
        num_rows = rows.numEntries()

    if isinstance(cols, np.ndarray):
        num_cols = cols.shape[0]
    else:
        num_cols = cols.numEntries()

    if num_rows != num_cols:
        raise ValueError('the length of vectors in rows and '
                         'cols must be the same')
    
    if not isinstance(rows, np.ndarray):
        rows = rows.getArray()
    rows *= 1 / (rows ** 2).sum(0) ** 0.5

    if not isinstance(cols, np.ndarray):
        cols = cols.getArray()
    cols *= 1 / (cols ** 2).sum(0) ** 0.5

    if diag:
        overlaps = np.einsum('ij,ji->i', rows.T, cols)
    else:
        overlaps = np.dot(rows.T, cols)
    
    return overlaps


def printOverlapTable(rows, cols):
    """Print table of overlaps (correlations) between two sets of modes.
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the printed table.  This function may be used to take
    a quick look into mode correspondences between two models.

    >>> # Compare top 3 PCs and slowest 3 ANM modes
    >>> printOverlapTable(p38_pca[:3], p38_anm[:3]) # doctest: +SKIP
    Overlap Table
                            ANM 1p38
                        #1     #2     #3
    PCA p38 xray #1   -0.39  +0.04  -0.71
    PCA p38 xray #2   -0.78  -0.20  +0.22
    PCA p38 xray #3   +0.05  -0.57  +0.06"""

    print(getOverlapTable(rows, cols))


def writeOverlapTable(filename, rows, cols):
    """Write table of overlaps (correlations) between two sets of modes to a
    file.  *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the overlap table.  See also :func:`.printOverlapTable`."""

    assert isinstance(filename, str), 'filename must be a string'
    out = openFile(filename, 'w')
    out.write(getOverlapTable(rows, cols))
    out.close()
    return filename


def getOverlapTable(rows, cols):
    """Make a formatted string of overlaps between modes in *rows* and *cols*.
    """

    overlap = calcOverlap(rows, cols)
    if isinstance(rows, Mode):
        rids = [rows.getIndex()]
        rname = str(rows.getModel())
    elif isinstance(rows, NMA):
        rids = np.arange(len(rows))
        rname = str(rows)
    elif isinstance(rows, ModeSet):
        rids = rows.getIndices()
        rname = str(rows.getModel())
    else:
        rids = [0]
        rname = str(rows)
    rlen = len(rids)
    if isinstance(cols, Mode):
        cids = [cols.getIndex()]
        cname = str(cols.getModel())
    elif isinstance(cols, NMA):
        cids = np.arange(len(cols))
        cname = str(cols)
    elif isinstance(cols, ModeSet):
        cids = cols.getIndices()
        cname = str(cols.getModel())
    else:
        cids = [0]
        cname = str(cols)
    clen = len(cids)
    overlap = overlap.reshape((rlen, clen))
    table = 'Overlap Table\n'
    table += (' '*(len(rname)+5) + cname.center(clen*7)).rstrip() + '\n'
    line = ' '*(len(rname)+5)
    for j in range(clen):
        line += ('#{0}'.format(cids[j]+1)).center(7)
    table += line.rstrip() + '\n'
    for i in range(rlen):
        line = rname + (' #{0}'.format(rids[i]+1)).ljust(5)
        for j in range(clen):
            if abs(overlap[i, j]).round(2) == 0.00:
                minplus = ' '
            elif overlap[i, j] < 0:
                minplus = '-'
            else:
                minplus = '+'
            line += (minplus+'{0:-.2f}').format(abs(overlap[i, j])).center(7)
        table += line.rstrip() + '\n'
    return table


def calcCumulOverlap(modes1, modes2, array=False):
    """Returns cumulative overlap of modes in *modes2* with those in *modes1*.
    Returns a number if *modes1* contains a single :class:`.Mode` or a
    :class:`.Vector` instance. If *modes1* contains multiple modes, returns an
    array. Elements of the array correspond to cumulative overlaps for modes
    in *modes1* with those in *modes2*.  If *array* is **True**, returns an array
    of cumulative overlaps. Returned array has the shape ``(len(modes1),
    len(modes2))``.  Each row corresponds to cumulative overlaps calculated for
    modes in *modes1* with those in *modes2*.  Each value in a row corresponds
    to cumulative overlap calculated using up to that many number of modes from
    *modes2*."""

    overlap = calcOverlap(modes1, modes2)
    if array:
        return np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    else:
        return np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))


def calcSubspaceOverlap(modes1, modes2):
    """Returns subspace overlap between two sets of modes (*modes1* and
    *modes2*).  Also known as the root mean square inner product (RMSIP)
    of essential subspaces [AA99]_.  This function returns a single number.

    .. [AA99] Amadei A, Ceruso MA, Di Nola A. On the convergence of the
       conformational coordinates basis set obtained by the essential
       dynamics analysis of proteins' molecular dynamics simulations.
       *Proteins* **1999** 36(4):419-424."""

    overlap = calcOverlap(modes1, modes2)
    if isinstance(modes1, Mode):
        length = 1
    else:
        length = len(modes1)
    rmsip = np.sqrt(np.power(overlap, 2).sum() / length)
    return rmsip

calcRMSIP = calcSubspaceOverlap


def calcRWSIP(modes1, modes2):
    """Returns root weighted square inner product (RWSIP)
    of essential subspaces [VC07]_.  This function returns a single number.

    .. [VC07] Carnevale V, Pontiggia F, Micheletti C. Structural and dynamical 
       alignment of enzymes with partial structural similarity.
       *J Phys Condens Matter.* **2007** 19:285206."""

    overlap = calcOverlap(modes1, modes2)
    if not isinstance(overlap, np.ndarray):
        overlap = np.array(overlap)
    
    if isinstance(modes1, Mode):
        length1 = 1
    else:
        length1 = len(modes1)

    if isinstance(modes2, Mode):
        length2 = 1
    else:
        length2 = len(modes2)

    vars1 = modes1.getVariances()
    if length1 == 1:
        vars1 = [vars1]
        overlap = overlap.reshape(-1, 1)

    vars2 = modes2.getVariances()
    if length2 == 1:
        vars2 = [vars2]
        overlap = overlap.reshape(1, -1)

    numerator = np.sum([vars1[l] * vars2[m] * overlap[l, m]**2
                        for l in range(length1)
                        for m in range(length2)])
    
    denominator = np.dot(vars1[:length1], vars2[:length1])
    
    rwsip = np.sqrt(numerator/denominator)
    return rwsip


def calcSquareInnerProduct(modes1, modes2):
    """Returns the square inner product (SIP) of fluctuations [SK02]_.  
    This function returns a single number.

    .. [SK02] Kundu S, Melton JS, Sorensen DC, Phillips GN: Dynamics of 
        proteins in crystals: comparison of experiment with simple models. 
        *Biophys J.* **2002**, 83:723-732.
        
    """
    if isinstance(modes1, (NMA, ModeSet)):
        w1 = calcSqFlucts(modes1)
    elif isListLike(modes1):
        w1 = modes1
    else:
        raise TypeError('modes1 should be a profile or an NMA or ModeSet object')

    if isinstance(modes2, (NMA, ModeSet)):
        w2 = calcSqFlucts(modes2)
    elif isListLike(modes2):
        w2 = modes2
    else:
        raise TypeError('modes2 should be a profile or an NMA or ModeSet object')

    return np.dot(w1, w2)**2 / (np.dot(w1, w1) * np.dot(w2, w2))

calcSIP = calcSquareInnerProduct


def calcSpectralOverlap(modes1, modes2, weighted=False, turbo=False):
    """Returns overlap between covariances of *modes1* and *modes2*.  Overlap
    between covariances are calculated using normal modes (eigenvectors),
    hence modes in both models must have been calculated.  This function
    implements equation 11 in [BH02]_.

    .. [BH02] Hess B. Convergence of sampling in protein simulations.
        *Phys Rev E* **2002** 65(3):031910.
    
    :arg weighted: if **True** then covariances are weighted by the trace.
    :type weighted: bool
    """

    if modes1.is3d() ^ modes2.is3d():
        raise TypeError('models must be either both 1-dimensional or 3-dimensional')
    if modes1.numAtoms() != modes2.numAtoms():
        raise ValueError('modes1 and modes2 must have same number of atoms')

    if isinstance(modes1, Mode):
        if weighted:
            varA = np.array([calcFractVariance(modes1)])
        else:
            varA = np.array([modes1.getVariance()])
        I = np.array([modes1.getIndex()])
    else:
        if weighted:
            varA = calcFractVariance(modes1)
        else:
            varA = modes1.getVariances()

        try:
            I = modes1.getIndices()
        except:
            try:
                modes1 = modes1[:]
                I = modes1.getIndices()
            except:
                raise TypeError('modes1 should be ModeSet or an object from which a ModeSet can be obtained')

    if isinstance(modes2, Mode):
        if weighted:
            varB = np.array([calcFractVariance(modes2)])
        else:
            varB = np.array([modes2.getVariance()])
        J = np.array([modes2.getIndex()])
    else:
        if weighted:
            varB = calcFractVariance(modes2)
        else:
            varB = modes2.getVariances()

        try:
            J = modes2.getIndices()
        except:
            try:
                modes2 = modes2[:]
                J = modes2.getIndices()
            except:
                raise TypeError('modes2 should be ModeSet or an object from which a ModeSet can be obtained')

    if turbo:
        model1 = modes1.getModel()
        model2 = modes2.getModel()
        if weighted:
            CACHE = WO_CACHE
        else:
            CACHE = SO_CACHE
        if (model1, model2) in CACHE:
            weights = CACHE[(model1, model2)]
        elif (model2, model1) in CACHE:
            weights = CACHE[(model2, model1)]
        else:
            farrayA = model1._getArray()
            farrayB = model2._getArray()

            fvarA = model1.getVariances()
            fvarB = model2.getVariances()

            dotAB = np.dot(farrayA.T, farrayB)**2
            outerAB = np.outer(fvarA**0.5, fvarB**0.5)
            CACHE[(model1, model2)] = weights = outerAB * dotAB
        
        weights = weights[I, :][:, J]
    else:
        arrayA = modes1._getArray()
        arrayB = modes2._getArray()

        dotAB = np.dot(arrayA.T, arrayB)**2
        outerAB = np.outer(varA**0.5, varB**0.5)
        weights = outerAB * dotAB

    diff = (np.sum(varA.sum() + varB.sum()) - 2 * np.sum(weights))

    if diff < ZERO:
        diff = 0
    else:
        diff = diff ** 0.5
    return 1 - diff / np.sqrt(varA.sum() + varB.sum())

def calcCovOverlap(modes1, modes2, turbo=False):
    """Returns overlap between covariances of *modes1* and *modes2*.  Overlap
    between covariances are calculated using normal modes (eigenvectors),
    hence modes in both models must have been calculated.  This function
    implements equation 11 in [BH02]_."""
    return calcSpectralOverlap(modes1, modes2, turbo=turbo)

def pairModes(modes1, modes2, **kwargs):
    """Returns the optimal matches between *modes1* and *modes2*. *modes1* 
    and *modes2* should have equal number of modes, and the function will 
    return a nested list where each item is a list containing a pair of modes.

    :arg index: if **True** then indices of modes will be returned instead of 
        :class:`Mode` instances.
    :type index: bool
    """

    index = kwargs.pop('index', False)
    method = kwargs.pop('method', None)

    if method is None:
        from scipy.optimize import linear_sum_assignment
        method = linear_sum_assignment

    if not (isinstance(modes1, (ModeSet, NMA)) \
        and isinstance(modes2, (ModeSet, NMA))):
        raise TypeError('modes1 and modes2 should be ModeSet or NMA instances')

    if len(modes1) != len(modes2):
        raise ValueError('the same number of modes should be provided')
    overlaps = calcOverlap(modes1, modes2)

    costs = 1 - abs(overlaps)
    row_ind, col_ind = method(costs)

    if index:
        return row_ind, col_ind

    if isinstance(modes1, ModeSet):
        row_ind = modes1._indices[row_ind]

    if isinstance(modes2, ModeSet):
        col_ind = modes2._indices[col_ind]

    outmodes1 = ModeSet(modes1.getModel(), row_ind)
    outmodes2 = ModeSet(modes2.getModel(), col_ind)

    return outmodes1, outmodes2

def _pairModes_wrapper(args):
    modeset0, modesets, index = args

    ret = []
    for modeset in modesets:
        _, reordered_modeset = pairModes(modeset0, modeset, index=index)
        ret.append(reordered_modeset)
    return ret

def matchModes(*modesets, **kwargs):
    """Returns the matches of modes among *modesets*. Note that the first 
    modeset will be treated as the reference so that only the matching 
    of each modeset to the first modeset is guaranteed to be optimal.
    
    :arg index: if **True** then indices of modes will be returned instead of 
                :class:`Mode` instances
    :type index: bool

    :arg turbo: if **True** then the computation will be performed in parallel. 
                The number of threads is set to be the same as the number of 
                CPUs. Assigning a number will specify the number of threads to be 
                used. Note that if writing a script, ``if __name__ == '__main__'`` 
                is necessary to protect your code when multi-tasking. 
                See https://docs.python.org/2/library/multiprocessing.html for details.
                Default is **False**
    :type turbo: bool, int
    """

    index = kwargs.pop('index', False)
    turbo = kwargs.pop('turbo', False)

    n_worker = None
    if not isinstance(turbo, bool):
        try:
            n_worker = int(turbo)
        except TypeError:
            raise TypeError('turbo should be Boolean or a number')

    modeset0 = modesets[0]
    if index:
        ret = [modeset0.getIndices()]
    else:
        ret = [modeset0]

    n_modes = len(modeset0)
    n_sets = len(modesets)
    if n_sets == 1:
        return ret
    elif n_sets == 0:
        raise ValueError('at least one modeset should be given')

    if turbo:
        from multiprocessing import Pool, cpu_count
        from math import ceil
        
        if not n_worker:
            n_worker = cpu_count()

        LOGGER.info('Matching {0} modes across {1} modesets with {2} threads...'
                        .format(n_modes, n_sets, n_worker))

        pool = Pool(n_worker)
        n_sets_per_worker = ceil((n_sets - 1) / n_worker)
        args = []
        for i in range(n_worker):
            start = i*n_sets_per_worker + 1
            end = (i+1)*n_sets_per_worker + 1
            subset = modesets[start:end]
            args.append((modeset0, subset, index))
        nested_ret = pool.map(_pairModes_wrapper, args)
        for entry in nested_ret:
            ret.extend(entry)

        pool.close()
        pool.join()
    else:
        LOGGER.progress('Matching {0} modes across {1} modesets...'
                        .format(n_modes, n_sets), n_sets, '_prody_matchModes')
        for i, modeset in enumerate(modesets):
            LOGGER.update(i, label='_prody_matchModes')
            if i > 0:
                _, reordered_modeset = pairModes(modeset0, modeset, index=index, **kwargs)
                ret.append(reordered_modeset)
        LOGGER.finish()
    
    return ret
