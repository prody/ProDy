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

__all__ = ['calcShannonEntropy', 'buildMutinfoMatrix', 'calcMSAOccupancy']


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
    memory efficient mode."""
    
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
                        turbo=bool(turbo), debug=kwargs.get('debug', False))
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
