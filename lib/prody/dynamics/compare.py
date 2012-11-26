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

"""This module defines functions for comparing normal modes from different 
models.""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody.utilities import openFile

from .nma import NMA
from .modeset import ModeSet
from .mode import VectorBase, Mode, Vector
from .gnm import ZERO

__all__ = ['calcOverlap', 'calcCumulOverlap', 'calcSubspaceOverlap', 
           'calcCovOverlap', 'printOverlapTable', 'writeOverlapTable',]
         
           
def calcOverlap(rows, cols):
    """Return overlap (or correlation) between two sets of modes (*rows* and 
    *cols*).  Returns a matrix whose rows correspond to modes passed as *rows* 
    argument, and columns correspond to those passed as *cols* argument.
    Both rows and columns are normalized prior to calculating overlap."""
    
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(rows)))
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(cols)))
    
    if rows.numDOF() != cols.numDOF(): 
        raise ValueError('number of degrees of freedom of rows and '
                         'cols must be the same')
    rows = rows.getArray()
    rows *= 1 / (rows ** 2).sum(0) ** 0.5
    cols = cols.getArray()
    cols *= 1 / (cols ** 2).sum(0) ** 0.5
    return np.dot(rows.T, cols)


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
    """Return cumulative overlap of modes in *modes2* with those in *modes1*.
    Returns a number of *modes1* contains a single :class:`.Mode` or a 
    :class:`.Vector` instance. If *modes1* contains multiple modes, returns an
    array. Elements of the array correspond to cumulative overlaps for modes 
    in *modes1* with those in *modes2*.  If *array* is **True**, Return array 
    of cumulative overlaps. Returned array has the shape ``(len(modes1), 
    len(modes2))``.  Each row corresponds to cumulative overlaps calculated for
    modes in *modes1* with those in *modes2*.  Each value in a row corresponds
    to cumulative overlap calculated using upto that many number of modes from 
    *modes2*."""
    
    overlap = calcOverlap(modes1, modes2)
    if array:
        return np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    else:
        return np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))


def calcSubspaceOverlap(modes1, modes2):
    """Return subspace overlap between two sets of modes (*modes1* and 
    *modes2*).  Also known as the root mean square inner product (RMSIP) 
    of essential subspaces [AA99]_.  This function returns a single number."""
    
    overlap = calcOverlap(modes1, modes2)
    if isinstance(modes1, Mode):
        length = 1
    else:
        length = len(modes1)
    rmsip = np.sqrt(np.power(overlap, 2).sum() / length)
    return rmsip
    
    
def calcCovOverlap(modelA, modelB):
    """Return overlap between covariances of *modelA* and *modelB*.  Overlap 
    between covariances are calculated using normal modes (eigenvectors), 
    hence modes in both models must have been calculated.  This function 
    implements equation 11 in [BH02]_."""
    
    if not modelA.is3d() or not modelB.is3d(): 
        raise TypeError('both models must be 3-dimensional') 
    if len(modelA) == 0 or len(modelB) == 0:  
        raise TypeError('modes must be calculated for both models, '
                        'try calcModes method')
    if modelA.numAtoms() != modelB.numAtoms(): 
        raise ValueError('modelA and modelB must have same number of atoms')
    arrayA = modelA._getArray()
    varA = modelA.getVariances()
    arrayB = modelB._getArray()
    varB = modelB.getVariances()
    
    dotAB = np.dot(arrayA.T, arrayB)**2
    outerAB = np.outer(varA**0.5, varB**0.5)
    diff = (np.sum(varA.sum() + varB.sum()) - 2 * np.sum(outerAB * dotAB))
    if diff < ZERO:
        diff = 0
    else:
        diff = diff ** 0.5
    return 1 - diff / np.sqrt(varA.sum() + varB.sum())
