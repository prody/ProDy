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
from prody import LOGGER

__all__ = ['showMSAOccupancy', 'showShannonEntropy', 'showMutinfoMatrix']


def pickSequence(msa):
    """Pick a sequence without gaps and deletions and return its residue 
    numbers and labels to be used as indices and X-axis label, or a pair
    of **None** at failure."""
    
    try:
        counts = calcMSAOccupancy(msa, 'row', count=True)
    except TypeError:
        return None, None
    else:
        length = msa.numResidues()
        split, msa.split = msa.split, True
        rows = (counts == length).nonzero()[0]
        for row in rows:
            try:
                label, seq, start, end = msa[row]
            except:
                break
            if end - start + 1 == length:
                msa.split = split
                return (arange(start, end + 1),
                        'Residue number ({0})'.format(label)) 
        msa.split = split
        return None, None


def showMSAOccupancy(msa, occ='res', indices=None, count=False, **kwargs):
    """Show a bar plot of occupancy calculated for :class:`.MSA` instance *msa*
    using :func:`.calcMSAOccupancy`.  *occ* may be ``'res'`` or ``'col'``, or a
    a pre-calculated occupancy array.  If x-axis *indices* are not specified,
    they will be inferred from *msa*.
    
    Occupancy is plotted using :func:`~matplotlib.pyplot.bar` function."""
    
    try:
        numseq, lenseq = msa.numSequences(), msa.numResidues()
    except AttributeError:
        raise TypeError('msa must be an MSA instance')
    
    try:
        length = len(occ)
    except TypeError:
        raise TypeError("occ must be 'res', 'col', or an occupancy array")
    
    xlabel = None
    try:
        sw = occ.startswith
    except TypeError:
        try:
            ndim = occ.ndim
        except AttributeError:
            raise TypeError("occ must be 'res', 'col', or an occupancy array")
        else:
            if ndim != 1:
                raise ValueError('occ must be a 1-dimensional array')
        if length == numseq and indices is None:
            indices = arange(1, length + 1)  
            xlabel = 'MSA sequence index'       
    else:
        occ = calcMSAOccupancy(msa, occ, count)
        length = len(occ)
        if indices is None and (sw('row') or sw('seq')):
            indices = arange(1, length + 1)
            xlabel = 'MSA sequence index'
    

    xlabel = kwargs.pop('xlabel', xlabel)
    if indices is None:
        indices, xlabel = pickSequence(msa)
        if indices is None:
            indices = arange(1, length + 1)
        xlabel = xlabel or 'MSA column index'
    
    ylabel = kwargs.pop('ylabel', 'Count' if count else 'Occupancy')  
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)
    import matplotlib.pyplot as plt
    show = plt.bar(indices, occ, **kwargs)
    if format:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title is None:
            title = 'Occupancy: ' + str(msa)
        plt.title(title)
    return show


def showShannonEntropy(entropy, indices=None, **kwargs):
    """Show a bar plot of Shannon *entropy* array.  :class:`.MSA` instances 
    or Numpy character arrays storing sequence alignments are also accepted 
    as *entropy* argument, in which case :func:`.calcShannonEntropy` will 
    be used for calculations.  *indices* may be residue numbers, when not 
    specified they will be inferred from *msa* or indices starting from 1 
    will be used.
    
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
        if msa is not None:
            indices, xlabel = pickSequence(msa)
        if indices is None:
            indices = arange(1, length + 1)
        xlabel = xlabel or 'MSA column index'
    
    ylabel = kwargs.pop('ylabel', 'Shannon entropy')  
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)
    import matplotlib.pyplot as plt
    show = plt.bar(indices, entropy, **kwargs)
    if format:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title is None:
            if msa is None:  
                title = 'Entropy'
            else:
                title = 'Entropy: ' + str(msa)
        plt.title(title)
    return show


def showMutinfoMatrix(mutinfo, *args, **kwargs):
    """Show a heatmap of mutual information array.  :class:`.MSA` instances 
    or Numpy character arrays storing sequence alignment are also accepted 
    as *mutinfo* argument, in which case :func:`.buildMutinfoMatrix` will 
    be used for calculations.  Note that x, y axes contain indices of the
    matrix starting from 1.
    
    Mutual information is plotted using :func:`~matplotlib.pyplot.imshow`
    function. vmin and vmax values can be set by user to achieve better
    signals using this function."""
    
    msa = None
    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        msa = mutinfo
        mutinfo = buildMutinfoMatrix(mutinfo)
        ndim, shape = mutinfo.ndim, mutinfo.shape

    msa = kwargs.pop('msa', msa)
    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')
    y, x = shape
    if x != y:
        raise ValueError('mutinfo matrix must be a square matrix')
    
    kwargs.setdefault('interpolation', 'nearest') 
    kwargs.setdefault('origin', 'lower')
    
    if msa is not None:
        indices, msalabel = pickSequence(msa)
        if indices is not None:
            start = indices[0] + 0.5
            end = start + x
            extent = [start, end, start, end]
        else:
            extent = [0.5, x + 0.5, 0.5, y + 0.5]
    else:
        msalabel = None
        extent = [0.5, x + 0.5, 0.5, y + 0.5]
    
    xlabel = kwargs.pop('xlabel', None)
    if xlabel is None:
        xlabel = msalabel or 'MSA column index'
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)
    
    import matplotlib.pyplot as plt
    show = plt.imshow(mutinfo, extent=extent, *args, **kwargs), plt.colorbar()
    
    if format:
        plt.xlabel(xlabel)
        plt.ylabel(xlabel)
        if title is None:
            if msa is None:  
                title = 'Mutual Information'
            else:
                title = 'Mutual Information: ' + str(msa)
        plt.title(title)
    return show





if __name__ == '__main__':
    from prody import *
    msa = parseMSA('piwi_seed.sth')
    print(repr(msa))
    msa = refineMSA(msa, label=msa[0][0])
    print(repr(msa))
    print(calcMSAOccupancy(msa, 'row', count=True))
    print(pickSequence(msa))
