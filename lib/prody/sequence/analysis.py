# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2012 Ahmet Bakan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines MSA analysis functions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import dtype, zeros
from numpy import indices, tril_indices
from prody import LOGGER

__all__ = ['calcShannonEntropy', 'buildMutinfoMatrix', 'calcMSAOccupancy',
           'applyMutinfoCorr', 'applyMutinfoNorm', 'calcRankorder',
           'buildSeqidMatrix', 'uniqueSequences', 'buildOMESMatrix',
           'buildSCAMatrix' ]


doc_turbo = """

    By default, *turbo* mode, which uses memory as large as the MSA array
    itself but runs four to five times faster, will be used.  If memory
    allocation fails, the implementation will fall back to slower and
    memory efficient mode."""


def getMSA(msa):
    """Return MSA character array."""

    try:
        msa = msa._getArray()
    except AttributeError:
        pass

    try:
        dtype_, ndim, shape = msa.dtype, msa.ndim, msa.shape
    except AttributeError:
        raise TypeError('msa must be an MSA instance or a 2D character array')

    if dtype_ != dtype('|S1') or ndim != 2:
        raise TypeError('msa must be an MSA instance or a 2D character array')

    return msa


def calcShannonEntropy(msa, ambiguity=True, omitgaps=True, **kwargs):
    """Return Shannon entropy array calculated for *msa*, which may be
    an :class:`.MSA` instance or a 2D Numpy character array.  Implementation
    is case insensitive and handles ambiguous amino acids as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.

    All non-alphabet characters are considered as gaps, and they are handled
    in two ways:

      * non-existent, the probability of observing amino acids in a given
        column is adjusted, by default
      * as a distinct character with its own probability, when *omitgaps* is
        **False**"""

    msa = getMSA(msa)
    length = msa.shape[1]
    entropy = zeros(length, float)
    from .msatools import msaentropy
    return msaentropy(msa, entropy,
                      ambiguity=bool(ambiguity), omitgaps=bool(omitgaps))


def buildMutinfoMatrix(msa, ambiguity=True, turbo=True, **kwargs):
    """Return mutual information matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.  Implementation
    is case insensitive and handles ambiguous amino acids as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps.

    Mutual information matrix can be normalized or corrected using
    :func:`applyMINormalization` and :func:`applyMICorrection` methods,
    respectively.  Normalization by joint entropy can performed using this
    function with *norm* option set **True**."""

    msa = getMSA(msa)

    from .msatools import msamutinfo
    LOGGER.timeit('_mutinfo')
    length = msa.shape[1]
    mutinfo = zeros((length, length), float)
    mutinfo = msamutinfo(msa, mutinfo,
                         ambiguity=bool(ambiguity), turbo=bool(turbo),
                         norm=bool(kwargs.get('norm', False)),
                         debug=bool(kwargs.get('debug', False)))
    LOGGER.report('Mutual information matrix was calculated in %.2fs.',
                  '_mutinfo')

    return mutinfo

buildMutinfoMatrix.__doc__ += doc_turbo


def calcMSAOccupancy(msa, occ='res', count=False):
    """Return occupancy array calculated for residue positions (default,
    ``'res'`` or ``'col'`` for *occ*) or sequences (``'seq'`` or ``'row'``
    for *occ*) of *msa*, which may be an :class:`.MSA` instance or a 2D
    Numpy character array.  By default, occupancy [0-1] will be calculated.
    If *count* is **True**, count of non-gap characters will be returned.
    Implementation is case insensitive."""

    from .msatools import msaocc

    msa = getMSA(msa)

    try:
        dim = occ.startswith('res') or occ.startswith('col')
    except AttributeError:
        raise TypeError('occ must be a string')
    occ = zeros(msa.shape[int(dim)], float)
    return msaocc(msa, occ, dim, count=bool(count))


def applyMutinfoNorm(mutinfo, entropy, norm='sument'):
    """Apply one of the normalizations discussed in [MLC05]_ to *mutinfo*
    matrix.  *norm* can be one of the following:

      * ``'sument'``: :math:`H(X) + H(Y)`, sum of entropy of columns
      * ``'minent'``: :math:`min\{H(X), H(Y)\}`, minimum entropy
      * ``'maxent'``: :math:`max\{H(X), H(Y)\}`, maximum entropy
      * ``'mincon'``: :math:`min\{H(X|Y), H(Y|X)\}`, minimum conditional
         entropy
      * ``'maxcon'``: :math:`max\{H(X|Y), H(Y|X)\}`, maximum conditional
         entropy

    where :math:`H(X)` is the entropy of a column, and
    :math:`H(X|Y) = H(X) - MI(X, Y)`.  Normalization with joint entropy, i.e.
    :math:`H(X, Y)`, can be done using :func:`.buildMutinfoMatrix` *norm*
    argument.

    .. [MLC05] Martin LC, Gloor GB, Dunn SD, Wahl LM. Using information theory
       to search for co-evolving residues in proteins. *Bioinformatics*
       **2005** 21(22):4116-4124."""

    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        raise TypeError('mutinfo must be a 2D square array')

    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('mutinfo must be a 2D square array')

    try:
        ndim, shapent = entropy.ndim, entropy.shape
    except AttributeError:
        raise TypeError('entropy must be a numpy array')

    if ndim != 1:
        raise ValueError('entropy must be a 1D array')

    if shapent[0] != shape[0]:
        raise ValueError('shape of mutinfo and entropy does not match')

    try:
        sw = norm.startswith
    except AttributeError:
        raise TypeError('norm must be a string')

    if sw('sument'):
        norm = lambda i_val, j_val, val: i_val + j_val

    elif sw('minent'):
        norm = lambda i_val, j_val, val: min(i_val, j_val)

    elif sw('maxent'):
        norm = lambda i_val, j_val, val: max(i_val, j_val)

    elif sw('mincon'):
        norm = lambda i_val, j_val, val: min(i_val - val, j_val - val)

    elif sw('maxcon'):
        norm = lambda i_val, j_val, val: max(i_val - val, j_val - val)

    elif sw('joint'):
        raise ValueError('for joint entropy normalization, use '
                         'buildMutinfoMatrix function')
    else:
        raise ValueError('norm={0} is not a valid normalization type'
                         .format(norm))

    mi = mutinfo.copy()
    for i, i_val in enumerate(entropy):
        for j, j_val in enumerate(entropy):
            val = mi[i, j]
            div = norm(i_val, j_val, val)
            if div == 0:
                mi[i, j] = 0
            else:
                mi[i, j] /= div

    return mi


def applyMutinfoCorr(mutinfo, corr='prod'):
    """Return a copy of *mutinfo* array after average product correction
    (default) or average sum correction is applied.  See [DSD08]_ for details.

    .. [DSD08] Dunn SD, Wahl LM, Gloor GB. Mutual information without the
       influence of phylogeny or entropy dramatically improves residue
       contact prediction. *Bioinformatics* **2008** 24(3):333-340."""

    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        raise TypeError('mutinfo must be a 2D square array')

    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('mutinfo must be a 2D square array')

    try:
        sw = corr.startswith
    except AttributeError:
        raise TypeError('correction must be a string')

    avg_mipos = mutinfo.sum(1) / (shape[0] - 1)
    avg_mi = avg_mipos.mean()

    mi = mutinfo.copy()
    if sw('prod') or sw('apc'):
        for i, i_avg in enumerate(avg_mipos):
            for j, j_avg in enumerate(avg_mipos):
                mi[i, j] -= (i_avg * j_avg) / avg_mi
    elif sw('sum') or sw('asc'):
        for i, i_avg in enumerate(avg_mipos):
            for j, j_avg in enumerate(avg_mipos):
                mi[i, j] -= i_avg + j_avg - avg_mi
    else:
        raise ValueError('correction must be prod or sum, not ' + corr)

    return mi


def buildSeqidMatrix(msa, turbo=True):
    """Return sequence identity matrix for *msa*."""

    msa = getMSA(msa)

    LOGGER.timeit('_seqid')
    from .seqtools import msaeye

    seqid = msaeye(msa, turbo=bool(turbo))

    LOGGER.report('Sequence identity matrix was calculated in %.2fs.',
                  '_seqid')
    return seqid

buildSeqidMatrix.__doc__ += doc_turbo


def uniqueSequences(msa, seqid=0.98, turbo=True):
    """Return a boolean array marking unique sequences in *msa*.  A sequence
    sharing sequence identity of *sqid* or more with another sequence coming
    before itself in *msa* will have a **False** value in the array."""

    msa = getMSA(msa)

    from .seqtools import msaeye

    if not (0 < seqid <= 1):
        raise ValueError('seqid must satisfy 0 < seqid <= 1')

    return msaeye(msa, unique=seqid, turbo=bool(turbo))

uniqueSequences.__doc__ += doc_turbo


def calcRankorder(matrix, zscore=False, **kwargs):
    """Returns indices of elements and corresponding values sorted in
    descending order, if *descend* is **True** (default). Can apply a zscore
    normalization; by default along *axis* - 0 such that each column has
    mean=0 and std=1.  If *zcore* analysis is used, return value contains the
    zscores. If matrix is smymetric only lower triangle indices will be
    returned, with diagonal elements if *diag* is **True** (default)."""

    try:
        ndim, shape = matrix.ndim, matrix.shape
    except AttributeError:
        raise TypeError('matrix must be a 2D array')

    if ndim != 2:
        raise ValueError('matrix must be a 2D array')

    threshold = kwargs.get('thredhold', 0.0001)
    try:
        symm = abs((matrix.transpose() - matrix).max()) < threshold
    except:
        symm = False

    if zscore:
        axis = int(bool(kwargs.get('axis', 0)))
        matrix = (matrix - matrix.mean(axis)) / matrix.std(axis)
        LOGGER.info('Zscore normalization has been applied.')

    descend = kwargs.get('descend', True)
    if not symm:
        if descend:
            sorted_index = matrix.argsort(axis=None)[::-1]
        else:
            sorted_index = matrix.argsort(axis=None)
        row = indices(shape)[0].flatten()[sorted_index]
        column = indices(shape)[1].flatten()[sorted_index]
    else:
        LOGGER.info('Matrix is symmetric, only lower triangle indices '
                    'will be returned.')
        if kwargs.get('diag', True):
            k = 0
        else:
            k = -1
        ind_row, ind_column = tril_indices(shape[0], k=k)
        matrix_lt = matrix[ind_row, ind_column]
        if descend:
            sorted_index = matrix_lt.argsort(axis=None)[::-1]
        else:
            sorted_index = matrix_lt.argsort(axis=None)
        row = ind_row[sorted_index]
        column = ind_column[sorted_index]

    return (row, column, matrix[row, column])


def buildOMESMatrix(msa, ambiguity=True, turbo=True, **kwargs):
    """Return OMES (Observed Minus Expected Squared) covariance matrix
    calculated for *msa*, which may be an :class:`.MSA` instance or a 2D
    NumPy character array. OMES is defined as::

                        (N_OBS - N_EX)^2              (f_i,j - f_i * f_j)^2
      OMES_(i,j) = sum(------------------) = N * sum(-----------------------)
                             N_EX                           f_i * f_j

    Implementation is case insensitive and handles ambiguous amino acids
    as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps."""

    msa = getMSA(msa)

    from .msatools import msaomes
    LOGGER.timeit('_omes')
    length = msa.shape[1]
    omes = zeros((length, length), float)
    omes = msaomes(msa, omes, ambiguity=bool(ambiguity), turbo=bool(turbo),
                   debug=bool(kwargs.get('debug', False)))
    LOGGER.report('OMES matrix was calculated in %.2fs.',
                  '_omes')

    return omes


def buildSCAMatrix(msa, turbo=True, **kwargs):
    """Return SCA matrix calculated for *msa*, which may be an :class:`.MSA`
    instance or a 2D Numpy character array.

    Implementation is case insensitive and handles ambiguous amino acids
    as follows:

      * **B** (Asx) count is allocated to *D* (Asp) and *N* (Asn)
      * **Z** (Glx) count is allocated to *E* (Glu) and *Q* (Gln)
      * **J** (Xle) count is allocated to *I* (Ile) and *L* (Leu)
      * **X** (Xaa) count is allocated to the twenty standard amino acids
      * Joint probability of observing a pair of ambiguous amino acids is
        allocated to all potential combinations, e.g. probability of **XX**
        is allocated to 400 combinations of standard amino acids, similarly
        probability of **XB** is allocated to 40 combinations of *D* and *N*
        with the standard amino acids.

    Selenocysteine (**U**, Sec) and pyrrolysine (**O**, Pyl) are considered
    as distinct amino acids.  When *ambiguity* is set **False**, all alphabet
    characters as considered as distinct types.  All non-alphabet characters
    are considered as gaps."""

    msa = getMSA(msa)
    from .msatools import msasca
    LOGGER.timeit('_sca')
    length = msa.shape[1]
    sca = zeros((length, length), float)
    sca = msasca(msa, sca, turbo=bool(turbo))
    LOGGER.report('SCA matrix was calculated in %.2fs.', '_sca')
    return sca