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

from numpy import all, zeros, dtype, array, char, fromstring

from prody import LOGGER

__all__ = ['calcShannonEntropy', 'buildMutinfoMatrix', 'calcMSAOccupancy', 
           'applyMICorrection', 'applyMINormalization']


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
        
    from .msatools import calcShannonEntropy
    return calcShannonEntropy(msa, ambiguity=bool(ambiguity), 
                              omitgaps=bool(omitgaps))
    
    
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
    
    By default, the will try to run in the *turbo* mode, which uses memory
    as large as the MSA array itself but runs four to five times faster.  If
    memory allocation fails, the implementation will switch to slower and
    memory efficient mode.
    
    Mutual information matrix can be normalized or corrected using
    :func:`applyMINormalization` and :func:`applyMICorrection` methods,
    respectively.  Normalization by joint entropy can performed using this
    function with *norm* option set **True**."""
    
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
        
    from .msatools import buildMutinfoMatrix
    LOGGER.timeit('_mutinfo')
    mutinfo = buildMutinfoMatrix(msa, ambiguity=bool(ambiguity), 
                                 turbo=bool(turbo), 
                                 norm=bool(kwargs.get('norm', False)), 
                                 debug=bool(kwargs.get('debug', False)))
    LOGGER.report('Mutual information matrix was calculated in %.2fs.', 
                  '_mutinfo')
        
    return mutinfo


def calcMSAOccupancy(msa, occ='res', count=False):
    """Return occupancy array calculated for residue positions (default, 
    ``'res'`` or ``'col'`` for *occ*) or sequences (``'seq'`` or ``'row'`` 
    for *occ*) of *msa*, which may be an :class:`.MSA` instance or a 2D 
    Numpy character array.  By default, occupancy [0-1] will be calculated.  
    If *count* is **True**, count of non-gap characters will be returned. 
    Implementation is case insensitive."""
    
    from .msatools import calcMSAOccupancy
    
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

    try:
        occ = occ.startswith('res') or occ.startswith('col')
    except AttributeError:
        raise TypeError('occ must be a string')
    return calcMSAOccupancy(msa, occ, count=bool(count))


def applyMINormalization(mutinfo, entropy, norm='sument'):
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
    argument."""
    
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
        
    mi = zeros(shape)
    if sw('sument'):
        for i, i_val in enumerate(entropy):
            for j, j_val in enumerate(entropy):
                mi[i, j] /= (i_val + j_val)
    
    elif sw('minent'):
        for i, i_val in enumerate(entropy):
            for j, j_val in enumerate(entropy):
                mi[i, j] /= min(i_val, j_val)
    
    elif sw('maxent'):
        for i, i_val in enumerate(entropy):
            for j, j_val in enumerate(entropy):
                mi[i, j] /= max(i_val, j_val)
    
    elif sw('mincon'):
        for i, i_val in enumerate(entropy):
            for j, j_val in enumerate(entropy):
                val = mi[i, j]
                mi[i, j] /= min(i_val - val, j_val - val)
    
    elif sw('maxcon'):
        for i, i_val in enumerate(entropy):
            for j, j_val in enumerate(entropy):
                val = mi[i, j]
                mi[i, j] /= max(i_val - val, j_val - val)
    
    elif sw('joint'):
        raise ValueError('for joint entropy normalization, use '
                         'buildMutinfoMatrix function')
    else:    
        raise ValueError('norm={0:s} is not a valid normalization type'
                         .format(norm))


def applyMICorrection(mutinfo, correction='prod'):
    """Return a copy of *mutinfo* array after average product correction 
    (default) or average sum correction is applied.  See [DSD08]_ for details.
    """
    
    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        raise TypeError('mutinfo must be a 2D square array')
    
    if ndim != 2 or shape[0] != shape[1]:
        raise ValueError('mutinfo must be a 2D square array')
    
    try:
        sw = norm.startswith
    except AttributeError:
        raise TypeError('correction must be a string')
    
    avg_mipos = mutinfo.sum(1) / (shape[0] - 1)
    avg_mi = avg_mipos.mean()
    
    mi = zeros(shape)
    if sw('prod'):
        for i, i_avg in enumerate(avg_mipos):
            for j, j_avg in enumerate(avg_mipos):
                mi[i, j] -= (i_avg * j_avg) / avg_mi
    elif sw('sum'):
        for i, i_avg in enumerate(avg_mipos):
            for j, j_avg in enumerate(avg_mipos):
                mi[i, j] -= i_avg + j_avg - avg_mi
    else:
        raise ValueError('correction must be prod or sum, not ' + correction)
    
    return micorrect
            
            
    
    
    
    
    
    
    
