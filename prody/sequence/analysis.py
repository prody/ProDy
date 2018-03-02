# -*- coding: utf-8 -*-
"""This module defines MSA analysis functions."""

__author__ = 'Anindita Dutta, Ahmet Bakan, Wenzhi Mao'

import os
from numpy import dtype, zeros, empty, ones, where, ceil, shape
from numpy import indices, tril_indices, array, ndarray
from prody import LOGGER
from prody.utilities import which
from prody.sequence.msa import MSA, refineMSA
from prody.sequence.msafile import parseMSA, writeMSA
from prody.sequence.sequence import Sequence
from prody.atomic import Atomic
from Bio import pairwise2
import sys

__all__ = ['calcShannonEntropy', 'buildMutinfoMatrix', 'calcMSAOccupancy',
           'applyMutinfoCorr', 'applyMutinfoNorm', 'calcRankorder', 'filterRankedPairs',
           'buildSeqidMatrix', 'uniqueSequences', 'buildOMESMatrix',
           'buildSCAMatrix', 'buildDirectInfoMatrix', 'calcMeff', 
           'buildPCMatrix', 'alignMultipleSequences', 'showAlignment', 
           'alignSequenceToMSA', 'calcPercentIdentities', 'alignSequencesByChain',]


doc_turbo = """

    By default, *turbo* mode, which uses memory as large as the MSA array
    itself but runs four to five times faster, will be used.  If memory
    allocation fails, the implementation will fall back to slower and
    memory efficient mode."""

def calcPercentIdentities(msa):
    percent_ids = []
    aas = ['A','C','D','E','F','G','H','I','J','K','L', \
           'M','N','P','Q','R','S','T','V','W','Y','-']

    for i in range(len(msa)):
        col_list = list(msa.getArray()[:,i])
        max_count = 0
        for aa in aas:
            if col_list.count(aa) > max_count:
                max_count = col_list.count(aa)
        percent_ids.append(float(max_count)/float(len(col_list))*100)

    return percent_ids

def getMSA(msa):
    """Returns MSA character array."""

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
    """Returns Shannon entropy array calculated for *msa*, which may be
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
    if msa.shape[0]<100:
        LOGGER.warning('SCA performs the best with higher number of sequences, and '
                       'minimal number of sequences is recommended as 100.')
    entropy = empty(length, float)
    from .msatools import msaentropy
    return msaentropy(msa, entropy,
                      ambiguity=bool(ambiguity), omitgaps=bool(omitgaps))


def buildMutinfoMatrix(msa, ambiguity=True, turbo=True, **kwargs):
    """Returns mutual information matrix calculated for *msa*, which may be an
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
    mutinfo = empty((length, length), float)
    mutinfo = msamutinfo(msa, mutinfo,
                         ambiguity=bool(ambiguity), turbo=bool(turbo),
                         norm=bool(kwargs.get('norm', False)),
                         debug=bool(kwargs.get('debug', False)))
    LOGGER.report('Mutual information matrix was calculated in %.2fs.',
                  '_mutinfo')

    return mutinfo

buildMutinfoMatrix.__doc__ += doc_turbo


def calcMSAOccupancy(msa, occ='res', count=False):
    """Returns occupancy array calculated for residue positions (default,
    ``'res'`` or ``'col'`` for *occ*) or sequences (``'seq'`` or ``'row'``
    for *occ*) of *msa*, which may be an :class:`.MSA` instance or a 2D
    NumPy character array.  By default, occupancy [0-1] will be calculated.
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
    """Returns a copy of *mutinfo* array after average product correction
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

def filterRankedPairs(indices, msa_indices, rank_row, rank_col, zscore_sort, \
                      num_of_pairs=20, seqDistance=5, resi_range=None, \
                      pdbDistance=8, chain1='A', chain2='A'):
    '''
    indices and msa_indices are lists output from alignSequenceToPDB
    
    rank_row, rank_col and zscore_sort are the outputs from calcRankorder
    
    :arg num_of_pairs: The number of pairs to be output, if no value is given
        then all pairs are output. Default is 20
    :type num_of_pairs: int
    
    :arg seqDistance: Remove pairs that are closer than this in the reference sequence
        Default is 5
    :type seqDistance: int
    
    :arg pdbDistance: Remove pairs with Calpha atoms further apart than this in the PDB
        Default is 8
    :type pdbDistance: int
    
    :arg chain1: The chain used for the residue specified by rank_row when measuring distances
    :type chain1: str
    
    :arg chain2: The chain used for the residue specified by rank_col when measuring distances
    :type chain2: str
    '''
    
    if indices is None:
        raiseValueError('Please provide indices output from alignSequenceToPDB')
    elif type(indices) != list:
        raiseTypeError('Please provide a valid indices list')
        
    if msa_indices is None:
        raiseValueError('Please provide msa_indices output from alignSequenceToPDB')
    elif type(indices) != list:
        raiseTypeError('Please provide valid msa_indices, which should be a list')
        
    if rank_row is None:
        raiseValueError('Please provide ranked row from calcRankorder')
        
    if rank_col is None:
        raiseValueError('Please provide ranked col from calcRankorder')
    
    if zscore_sort is None:
        raiseValueError('Please provide sorted Z scores from calcRankorder')
    
    if num_of_pairs is None:
        num_of_pairs = len(rank_row)
        
    pairList = []
    i = -1
    j = 0
    while j < num_of_pairs:
        
        i += 1
        
        if type(indices[where(msa_indices == rank_row[i])[0][0]]) != np.int64 or \
        type(indices[where(msa_indices == rank_col[i])[0][0]]) != np.int64:
            continue
        
        if indices[where(msa_indices == rank_row[i])[0][0]] - \
        indices[where(msa_indices == rank_col[i])[0][0]] < seqDistance:
            continue        
        
        distance = calcDistance(pdb.select('chain %s and resid %s' % (chain1, \
                                           indices[where(msa_indices == \
                                           rank_row[i])[0][0]])).copy(), \
                                pdb.select('chain %s and resid %s' % (chain2, \
                                           indices[where(msa_indices == \
                                           rank_col[i])[0][0]])).copy())
        if distance > pdbDistance:
            continue
            
        if resi_range is not None:
            if not indices[where(msa_indices == rank_row[i])[0][0]] in resi_range and \
            not indices[where(msa_indices == rank_col[i])[0][0]] in resi_range:
                continue
            
        pairList.append('%3d' % i + ':\t%3d' % indices[where(msa_indices == \
        rank_row[i])[0][0]] + '\t' + '%3d' % indices[where(msa_indices == \
        rank_col[i])[0][0]] + '\t' + '%5.1f' % zscore_sort[i] + '\t' + \
        '%5.1f' % distance + '\n')
        
        j += 1
    
    return pairList

def buildSeqidMatrix(msa, turbo=True):
    """Returns sequence identity matrix for *msa*."""

    msa = getMSA(msa)

    LOGGER.timeit('_seqid')
    from .seqtools import msaeye

    dim = msa.shape[0]
    seqid = msaeye(msa, ones((dim, dim), float), turbo=bool(turbo))

    LOGGER.report('Sequence identity matrix was calculated in %.2fs.',
                  '_seqid')
    return seqid

buildSeqidMatrix.__doc__ += doc_turbo


def uniqueSequences(msa, seqid=0.98, turbo=True):
    """Returns a boolean array marking unique sequences in ``msa``.  A sequence
    sharing sequence identity of ``seqid`` or more with another sequence coming
    before itself in ``msa`` will have a **False** value in the array."""

    msa = getMSA(msa)

    from .seqtools import msaeye

    if not (0 < seqid <= 1):
        raise ValueError('seqid must satisfy 0 < seqid <= 1')

    return msaeye(msa, zeros(msa.shape[0], bool),
                  unique=seqid, turbo=bool(turbo))

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
    """Returns OMES (Observed Minus Expected Squared) covariance matrix
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
    omes = empty((length, length), float)
    omes = msaomes(msa, omes, ambiguity=bool(ambiguity), turbo=bool(turbo),
                   debug=bool(kwargs.get('debug', False)))
    LOGGER.report('OMES matrix was calculated in %.2fs.',
                  '_omes')

    return omes

buildOMESMatrix.__doc__ += doc_turbo


def buildSCAMatrix(msa, turbo=True, **kwargs):
    """Returns SCA matrix calculated for *msa*, which may be an :class:`.MSA`
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

buildSCAMatrix.__doc__ += doc_turbo

def buildPCMatrix(msa, turbo=False, **kwargs):
    """Returns PC matrix calculated for *msa*, which may be an :class:`.MSA`
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
    are considered as gaps.
    """

    msa = getMSA(msa)
    from .msatools import msapsicov
    LOGGER.timeit('_psicov')
    length = msa.shape[1]
    pc = zeros((length, length), float)
    pc = msapsicov(msa, pc, turbo=bool(turbo))
    LOGGER.report('PC matrix was calculated in %.2fs.', '_psicov')
    return pc

def buildDirectInfoMatrix(msa, seqid=.8, pseudo_weight=.5, refine=False,
                          **kwargs):
    """Returns direct information matrix calculated for *msa*, which may be an
    :class:`.MSA` instance or a 2D Numpy character array.

    Sequences sharing sequence identity of *seqid* or more with another
    sequence are regarded as similar sequences for calculating their weights
    using :func:`.calcMeff`.

    *pseudo_weight* are the weight for pseudo count probability.

    Sequences are not refined by default. When *refine* is set **True**,
    the MSA will be refined by the first sequence and the shape of direct
    information matrix will be smaller.
    """

    msa = getMSA(msa)
    from .msatools import msadipretest, msadirectinfo1, msadirectinfo2
    from numpy import matrix

    LOGGER.timeit('_di')
    if msa.shape[0]<250:
        LOGGER.warning('DI performs the best with higher number of sequences, and '
                       'minimal number of sequences is recommended as 250.')
    refine = 1 if refine else 0
    # msadipretest get some parameter from msa to set matrix size
    length, q = msadipretest(msa, refine=refine)
    c = matrix.dot(matrix(zeros((length*q, 1), float)),
                   matrix(zeros((1, length*q), float)))
    prob = zeros((length, q+1), float)
    # msadirectinfo1 return c to be inversed and prob to be used
    meff, n, length, c, prob = msadirectinfo1(msa, c, prob, theta=1.-seqid,
                                              pseudocount_weight=pseudo_weight,
                                              refine=refine, q=q+1)

    c = c.I

    di = zeros((length, length), float)
    # get final DI
    di = msadirectinfo2(n, length, c, prob, di, q+1)
    del prob, c
    LOGGER.report('DI matrix was calculated in %.2fs.', '_di')
    return di


def calcMeff(msa, seqid=.8, refine=False, weight=False, **kwargs):
    """Returns the Meff for *msa*, which may be an :class:`.MSA`
    instance or a 2D Numpy character array.

    Since similar sequences in an *msa* decreases the diversity of *msa*,
    *Meff* gives a weight for sequences in the *msa*.

    For example: One sequence in MSA has 5 other similar sequences in this
    MSA(itself included). The weight of this sequence is defined as 1/5=0.2.
    Meff is the sum of all sequence weights. In another word, Meff can be
    understood as the effective number of independent sequences.

    Sequences sharing sequence identity of *seqid* or more with another
    sequence are regarded as similar sequences to calculate Meff.

    Sequences are not refined by default. When *refine* is set **True**, the
    MSA will be refined by the first sequence.

    The weight for each sequence are returned when *weight* is **True**."""

    msa = getMSA(msa)
    from .msatools import msameff
    LOGGER.timeit('_meff')
    refine = 1 if refine else 0
    weight = 0 if weight else 1  # A Mark for return weighted array.
    if (not weight):
        w = zeros((msa.shape[0]), float)
        meff = msameff(msa, theta=1.-seqid, meff_only=weight,
                       refine=refine, w=w)
    else:
        meff = msameff(msa, theta=1.-seqid, meff_only=weight, refine=refine)
    LOGGER.report('Meff was calculated in %.2fs.', '_meff')
    return meff

def msaeye(msa, unique, turbo):
    tic1 = timeit.default_timer()
    length = msa.shape[1]
    number = msa.shape[0]
    # number = 5
    array = eye(int(number))

    seqs = []
    for i in xrange(number):
        seqs.append(msa[i,:])
    iseq = zeros((number, length), dtype=int)

    for i in xrange(0,number-1):
        if i == 0:
            for k in xrange(length):
                if ord(seqs[i][k])>90:
                    iseq[i,k]=ord(seqs[i][k])-96 if ord(seqs[i][k])-96 > 0 and ord(seqs[i][k])-96 < 26 else 0
                else:
                    iseq[i,k]=ord(seqs[i][k])-64 if ord(seqs[i][k])-64 > 0 and ord(seqs[i][k])-64 < 26 else 0
            for j in xrange(i+1,number):
                score=0.
                ncols=0.
                for k in xrange(length):
                    if ord(seqs[j][k])>90:
                        iseq[j,k]=ord(seqs[j][k])-96 if ord(seqs[j][k])-96 > 0 and ord(seqs[j][k])-96 < 26 else 0
                    else:
                        iseq[j,k]=ord(seqs[j][k])-64 if ord(seqs[j][k])-64 > 0 and ord(seqs[j][k])-64 < 26 else 0
                    if iseq[i,k] or iseq[j,k]:
                        ncols += 1
                        if iseq[i,k]==iseq[j,k]:
                            score+=1
                array[i,j]=float(score)/ncols
                array[j,i]=array[i,j]
            # print iseq[0]
            # print seqs[0]
            # raw_input()
        else:
            for j in xrange(i+1,number):
                score=0.
                ncols=0.
                for k in xrange(length):
                    if iseq[i,k] or iseq[j,k]:
                        ncols += 1
                        if iseq[i,k]==iseq[j,k]:
                            score+=1
                array[i,j]= float(score)/ncols#float(sum((iseq[i] == iseq[j])*(iseq[i]*iseq[j]!=0))) / sum(iseq[i]*iseq[j]!=0)
                array[j,i]=array[i,j]

    toc1 = timeit.default_timer()
    elapsed1 = toc1 - tic1
    print(elapsed1)

def alignSequencesByChain(PDBs, **kwargs):
    """
    Runs alignMultipleSequences for each chain and optionally joins the results.
    Returns either a single MSA or a dictionary containing an MSA for each chain.

    :arg PDBs: a list or array of :class:`AtomGroup` objects or PDB IDs
        a mixed list containing both is acceptable
    :type PDBs: list or :class:`~numpy.ndarray`

    :arg join_chains: whether to join chain alignments
        default is True
    :type join_chains: bool 

    :arg join_char: a character for joining chain alignments
        default is '/' as used by PIR format alignments
    :type join_char: str
    """
    if not (isinstance(PDBs, list) or isinstance(PDBs, ndarray)):
        raise TypeError('sequences should be a list or array')

    pdbs = []
    chains = []
    for i, pdb in enumerate(PDBs):
        if isinstance(pdb, str):
            if len(pdb) in [4, 5, 6]:
                pdbs.append(parsePDB(pdb))
            else:
                raise ValueError('string entries in PDBs should be of length 4, 5 or 6')
        elif isinstance(pdb, Atomic):
            pdbs.append(pdb)
        else:
            raise TypeError('each entry in PDBs must be a :class:`Atomic` instance or a PDB ID')

        chains.append([])
        for chain in list(pdbs[i].getHierView()):
            chains[i].append(chain)

        if i != 0 and len(chains[i]) != len(chains[0]):
            raise ValueError('all pdbs should have the same number of chains')

    chains = array(chains)
    chain_alignments = []
    alignments = {}
    for j in range(len(chains[0])):
        msa = alignMultipleSequences(chains[:,j], prefix=pdbs[0].getChids()[j])
        msa = refineMSA(msa, colocc=1e-9) # remove gap-only cols
        chain_alignments.append(msa)
        alignments[pdbs[0].getChids()[j]] = msa

    join_chains = kwargs.get('join_chains', True)
    join_char = kwargs.get('join_char','/')
    if join_chains:
        aligned_sequences = list(zeros(shape(chain_alignments)).T)
        for j in range(shape(chain_alignments)[1]):
            aligned_sequences[j] = list(aligned_sequences[j])
        
        orig_labels = []
        for i, chain_alignment in enumerate(chain_alignments):
            for j, sequence in enumerate(chain_alignment):
                aligned_sequences[j][i] = str(sequence)
                if i == 0: orig_labels.append(sequence.getLabel())

        joined_msaarr = []
        for j in range(shape(chain_alignments)[1]):
            joined_msaarr.append(array(list(join_char.join(aligned_sequences[j]))))
        joined_msaarr = array(joined_msaarr)
        
        result = MSA(joined_msaarr, title='joined_chains', labels=orig_labels)

    else:
        result = alignments
            
    return result

def alignMultipleSequences(sequences, **kwargs):
    """
    Aligns sequences with clustalw or clustalw2 and returns the resulting MSA.

    :arg sequences: a file, MSA object or a list or array containing sequences
       as Atomic objects with :func:`getSequence` or Sequence objects or strings. 
       If strings are used then labels must be provided using ``labels``
    :type sequences: :class:`Atomic`, :class:`.MSAFile`, :class:`.MSA`, 
        :class:`~numpy.ndarray`, str

    :arg labels: a list of labels to go with the sequences
    :type labels: list

    :arg prefix: a prefix for filenames
    :type prefix: str
    """
    # 1. check if sequences are in a fasta file and if not make one
    fetched_labels = []

    if not (isinstance(sequences, list) or isinstance(sequences, ndarray)):
        raise TypeError('sequences should be a list or array')
    else:
        if isinstance(sequences[0], Atomic):
            msa = []
            for sequence in sequences:
                msa.append(sequence.getSequence())
                fetched_labels.append(sequence.getTitle())
            sequences = msa

        if isinstance(sequences[0], Sequence):
            msa = []
            for sequence in sequences:
                msa.append(str(sequence))
                fetched_labels.append(sequence.getLabel())
            sequences = msa

        if isinstance(sequences[0], str):
            max_len = 0
            for sequence in sequences:
                if len(sequence) > max_len:
                    max_len = len(sequence)

            msa = []
            for sequence in sequences:
                sequence = sequence + '-'*(max_len - len(sequence))
                msa.append(array(list(sequence)))
            sequences = array(msa)

        labels = kwargs.get('labels',None)
        if labels is None or len(labels) != len(sequences):
            if fetched_labels is not []:
                labels = fetched_labels
            else:
                raise ValueError('sequences can only be provided as a list of '
                                 'strings if corresponding labels are provided')

        sequences = MSA(msa=sequences, labels=labels)

    prefix = kwargs.get('prefix',None)
    if isinstance(sequences, MSA):
        if prefix is not None:
            sequences = writeMSA(prefix + '.fasta', sequences)
        else:
            raise ValueError('please provide a prefix for the MSA file to be '
                             'fed to clustalw')

    try:
        msa = parseMSA(sequences)
    except:
        raise ValueError('sequences is not an MSA file and could not be made into one')

    if prefix is None: 
        pos = sequences.rfind('.')
        prefix = sequences[:pos]

    # 2. find and run alignment method
    clustalw = which('clustalw')
    if clustalw is None:
        if which('clustalw2') is not None:
            clustalw = which('clustalw2')
        else:
            raise EnvironmentError("The executable for clustalw was not found, \
                                    install clustalw or add it to the path.")

    os.system(clustalw + " " + sequences)

    # 3. parse and return the new MSA
    return parseMSA(prefix + '.aln')

def showAlignment(alignment, row_size=60, max_seqs=5, **kwargs):
    """
    Prints out an alignment as sets of short rows with labels.

    :arg alignment: any object with aligned sequence
    :type alignment: :class: `.MSA`, tuple or list

    :arg row_size: the size of each row
        default 60
    :type row_size: int

    :arg max_seqs: the maximum number of sequences to show
        default 5
    :type max_seqs: int

    :arg indices: a set of indices for some or all sequences
        that will be shown above the relevant sequences
    :type indices: array, list or tuple of arrays, lists or tuples

    :arg index_start: how far along the alignment to start putting indices
        default 0
    :type index_start: int

    :arg index_stop: how far along the alignment to stop putting indices
        default the point when the shortest sequence stops
    :type index_stop: int
    """
    indices = kwargs.get('indices',None)
    index_start = kwargs.get('index_start',0)
    index_stop = kwargs.get('index_start',0)

    if index_stop == 0 and indices is not None:
        locs = []
        maxes = []
        for index in indices:
            int_index = []
            for i in index:
                if i == '':
                    int_index.append(0)
                else:
                    int_index.append(int(i))
            int_index = array(int_index)
            maxes.append(max(int_index))
            locs.append(where(int_index == max(int_index))[0][0])
        index_stop = locs[where(maxes == min(maxes))[0][0]]

    if len(alignment) < max_seqs:
        max_seqs = len(alignment)

    for i in range(int(ceil(len(alignment[0])/float(row_size)))):
        for j in range(max_seqs):
            if indices is not None:
                sys.stdout.write('\n' + ' '*15 + '\t')
                for k in range(row_size*i+10,row_size*(i+1)+10,10):
                    try:
                        if k > index_start + 10 and k < index_stop + 10:
                            sys.stdout.write('{:10d}'.format(int(indices[j][k-1])))
                        elif k > index_start:
                            sys.stdout.write(' '*(k-index_start))
                        else:
                            sys.stdout.write(' '*10)
                    except:
                            sys.stdout.write(' '*10)
                sys.stdout.write('\n')

            sys.stdout.write(alignment[j].getLabel()[:15] + \
                             ' ' * (15-len(alignment[j].getLabel()[:15])) + \
                             '\t' + str(alignment[j])[60*i:60*(i+1)] + '\n')

        sys.stdout.write('\n')
        
    return

def alignSequenceToMSA(seq, msa, label, match=5, mismatch=-1, gap_opening=-10, gap_extension=-1):
    """
    Align a sequence from a PDB or Sequence to a sequence from an MSA
    and create two sets of indices. 

    The sequence from the MSA (refSeq), the alignment and 
    the two sets of indices are returned. 
    
    The first set (indices) maps the residue numbers in the PDB to 
    the reference sequence. The second set (msa_indices) indexes the 
    reference sequence in the msa and is used for retrieving values 
    from the first indices.

    :arg seq: an object with an associated sequence string 
         or a sequence string itself
    :type seq: :class:`Atomic`, :class:`Sequence`, or str
    
    :arg msa: MSA object
    :type msa: :class:`.MSA`
    
    :arg label: a label for a sequence in msa or a PDB ID
        ``msa.getIndex(label)`` must return a sequence index
    :type label: str
    
    :arg chain: which chain from pdb to use for alignment
    :type chain: str
    
    :arg match: a positive integer, used to reward finding a match
        The default is 5, which we found to work in a test case.
    :type match: int
    
    :arg mismatch: a negative integer, used to penalise finding a mismatch
        The default is -1, which we found to work in a test case
    :type mismatch: int
    
    :arg gap_opening: a negative integer, used to penalise opening a gap
        The default is -10, which we found to work in a test case
    :type gap_opening: int
    
    :arg gap_extension: a negative integer, used to penalise extending a gap
        The default is -1, which we found to work in a test case
    :type gap_extension: int
    """
    if isinstance(seq, Atomic):
        ag = seq
        sequence = ag.calpha.getSequence()

    elif isinstance(seq, Sequence):
         sequence = str(seq)
         ag = None

    elif isinstance(seq, str):
        if len(seq) == 4 or len(seq) == 5:
            ag = parsePDB(seq)
            sequence = ag.calpha.getSequence()
        else:
            sequence = seq
            ag = None

    else:
        raise TypeError('seq must be an atomic class, sequence class, or str not {0}'
                        .format(type(seq)))

    if not isinstance(msa, MSA):
        raise TypeError('msa must be an MSA instance')

    if label is None:
        if ag:
            label = ag.getTitle().split('_')[0]
        elif isinstance(seq, Sequence):
            label = seq.getLabel()
        else:
            raise ValueError('A label cannot be extracted from seq so please provide one.')

    try:
        seqIndex = msa.getIndex(label)
    except:
        raise ValueError('Please provide a label that can be found in msa.')

    if isinstance(seqIndex, int):
        refMsaSeq = str(msa[seqIndex]).upper().replace('-','.')

    else:
        raise TypeError('The output from querying that label against msa is not a single sequence.')

    alignment = pairwise2.align.globalms(sequence, str(refMsaSeq), \
                                         match, mismatch, gap_opening, gap_extension)

    seq_indices = [0]
    msa_indices = [0]

    for i in range(len(alignment[0][0])):
        if alignment[0][0][i] != '-':
            seq_indices.append(seq_indices[i]+1)
        else:
            seq_indices.append(seq_indices[i])

        if alignment[0][1][i] != '-':
            msa_indices.append(msa_indices[i]+1)
        else:
            msa_indices.append(msa_indices[i])

    seq_indices = array(seq_indices)
    msa_indices = array(msa_indices)

    alignment = MSA(msa=array([array(list(alignment[0][0])), \
                               array(list(alignment[0][1]))]), \
                    labels=[ag.getTitle(), label])

    return alignment, seq_indices, msa_indices

