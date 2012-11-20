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

from numpy import arange

from .analysis import *

__all__ = ['showShannonEntropy']


def showShannonEntropy(entropy, indices=None, **kwargs):
    """Show a bar plot of Shannon *entropy* array.  :class:`MSA` instances 
    or Numpy character arrays storing sequence alignments are also accepted 
    as *entropy* argument, in which case :func:`.calcShannonEntropy` will 
    be used for calculations.  *indices* may be residue numbers, if **None**
    is given numbers starting from 1 will be used.
    
    Entropy is plotted using :func:`~matplotlib.pyplot.bar` function."""
    
    msa = None
    try:
        ndim = entropy.ndim
    except AttributeError:
        msa = entropy
        entropy = calcShannonEntropy(msa)
        ndim = entropy.ndim
    
    if ndim != 1:
        raise ValueError('entropy must be a 1D array')

    msa = kwargs.pop('msa', msa)
    xlabel = kwargs.pop('xlabel', None)
    if indices is None:
        length = len(entropy)    
        if msa is not None and False:
            try:
                occ = calcMSAOccupancy(msa, 'row')
            except TypeError:
                pass
            else:
                split, msa.split = msa.split, True
                rows = (occ == 1.0).nonzero()[0]
                for row in rows: 
                    #try:
                    label, seq, start, end = msa[row]
                    #except:
                    #    break
                    #else:
                    if end - start - 1 == length:
                        indices = arange(start, end + 1)
                        xlabel = 'Residue number ({0:s})'.format(label)
                        break
                msa.split = split
            
        if indices is None:  
            indices = arange(1, length + 1)
        xlabel = xlabel or 'MSA column index'

    ylabel = kwargs.pop('ylabel', 'Shannon entropy')  
    import matplotlib.pyplot as plt
    show = plt.bar(indices, entropy, **kwargs)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)  
    return show
