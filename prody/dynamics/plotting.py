# ProDy: A Python Package for Protein Structural Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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

""":mod:`plotting` module defines plotting functions for assisting interactive
dynamics analysis."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np
pl = None

from prody import ProDyLogger as LOGGER
import prody
from .functions import getUniformModel
from .functions import *
from .nma import *

__all__ = ['showFractOfVariances', 'showProjection',
           'showSumOfWeights', 'showOverlapMatrix',
           'showCrossCorrelations', 'showMode', 'showSqFlucts',
           'showContactMap', 'showOverlap', 'showCumulativeOverlap',
           'showCumFractOfVariances'
           ]

def showFractOfVariances(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`matplotlib.pyplot.bar`.
    
    Note that mode indices are increased by 1.
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA, or ModeSet, not {0:s}'.format(type(modes)))
    
    fracts = [(mode.getIndex(), mode.getFractOfVariance()) for mode in modes]
    fracts = np.array(fracts)
    show = pl.bar(fracts[:,0]+0.5, fracts[:,1], *args, **kwargs)
    axis = list(pl.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    pl.axis(axis)
    pl.xlabel('Mode index')
    pl.ylabel('Fraction of variance')
    return show

def showCumFractOfVariances(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`matplotlib.pyplot.plot`.
    
    Note that mode indices are increased by 1.
    
    """
    if pl is None: prody.importPyPlot()
    modes = getUniformModel(modes)
    if modes.simple:
        LOGGER.warning('Vector instances are not accepted')
        return None
    
    fracts = np.array([mode.getFractOfVariance() for mode in modes.modes]).cumsum()
    show = pl.plot(modes.indices+0.5, fracts, *args, **kwargs)
    axis = list(pl.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    pl.axis(axis)
    pl.xlabel('Mode index')
    pl.ylabel('Fraction of variance')
    return show

def showProjection(ensemble, modes, *args, **kwargs):
    """Show projection of conformational deviations onto given modes.
    
    :arg ensemble: an :class:`Ensemble` instance
    
    Matplotlib function used for plotting depends on the number of modes:
        
      * 1 mode: :func:`matplotlib.pyplot.hist`
      * 2 modes: :func:`matplotlib.pyplot.scatter`
      * 3 modes: not implemented yet
     
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    if isinstance(modes, Mode):
        n_modes = 1
    else:
        n_modes = len(modes)
    if not 0 < n_modes < 4: 
        raise Exception('1, 2, or 3 modes are accepted')
    projection = getProjection(ensemble, modes)
    if n_modes == 1:
        show = pl.hist(projection.flatten(), *args, **kwargs)
        pl.xlabel(str(modes.modes[0])+ ' coordinate')
        pl.ylabel('Number of conformations')
    elif n_modes == 2:
        show = pl.scatter(projection[:,0], projection[:,1], *args, **kwargs)
        modes = [m for m in modes]
        pl.xlabel(str(modes[0])+ ' coordinate')
        pl.ylabel(str(modes[1])+ ' coordinate')
    else:
        return
    return show

#def showSumOfWeights(ensemble, indices=None, *args, **kwargs):
def showSumOfWeights(ensemble, *args, **kwargs):
    """Show sum of weights from an ensemble using :func:`matplotlib.pyplot.plot`.
    
    Weights are summed for each atom over conformations in the ensemble.
    Size of the plotted array will be equal to the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    """
    """
    *indices*, if given, will be used as X values. Otherwise, X axis will
    start from 0 and increase by 1 for each atom. 
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble instance')
    
    weights = getSumOfWeights(ensemble)
    
    if weights is None:
        return None
    
    show = pl.plot(weights, *args, **kwargs)
    
    axis = list(pl.axis())
    axis[2] = 0
    axis[3] += 1
    pl.axis(axis)
    pl.xlabel('Atom index')
    pl.ylabel('Sum of weights')
    return show
    
    
def showOverlapMatrix(rows, cols, *args, **kwargs):
    """Show overlap matrix using :func:`matplotlib.pyplot.pcolor`.
    
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the displayed matrix.
    
    Note that mode indices are increased by 1. List of modes should contain
    a set of contiguous modes from the same model. 
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(rows, (NMA, ModeSet)):
        raise TypeError('rows must be an NMA model or a ModeSet, not {0:s}'.format(type(rows)))
    if not isinstance(rows, (NMA, ModeSet)):
        raise TypeError('cols must be an NMA model or a ModeSet, not {0:s}'.format(type(cols)))
        
    overlap = abs(getOverlap(rows, cols))
    
    if isinstance(rows, NMA):
        rows = rows[:]
    if isinstance(cols, NMA):
        cols = cols[:]
    
    X = rows.getIndices()+0.5#np.arange(modes1.indices[0]+0.5, len(modes1)+1.5)
    Y = cols.getIndices()+0.5#np.arange(modes2.indices[0]+0.5, len(modes2)+1.5)
    axis = (X[0], X[-1], Y[0], Y[-1])
    X, Y = np.meshgrid(X, Y)

    show = pl.pcolor(X, Y, overlap, cmap=pl.cm.jet, *args, **kwargs), pl.colorbar()
    pl.xlabel(str(cols))
    pl.ylabel(str(rows))
    pl.axis(axis)
    return show

def showCrossCorrelations(modes, *args, **kwargs):
    """Show cross-correlations for given modes using :func:`matplotlib.pyplot.imshow`.
    
    See also :func:`getCrossCorrelations`. 
    
    By default, *origin=lower* and *interpolation=bilinear* keyword
    arguments are passed to imshow function. User can overwrite these
    parameters.
    
    """
    if pl is None: prody.importPyPlot()
    modes = getUniformModel(modes)
    arange = np.arange(modes.n_atoms)
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:, 
                       arange[0]+1:] = getCrossCorrelations(modes)
    if not kwargs.has_key('interpolation'):
        kwargs['interpolation'] = 'bilinear'
    if not kwargs.has_key('origin'):
        kwargs['origin'] = 'lower'
    show = pl.imshow(cross_correlations, *args, **kwargs), pl.colorbar()
    pl.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    pl.title('{0:s} cross-correlations ({1:d} modes)'
             .format(str(modes), len(modes))) 
    pl.xlabel('Atom indices')
    pl.ylabel('Atom indices')
    return show

def showMode(mode, *args, **kwargs):
    """Show mode array using :func:`matplotlib.pyplot.plot`."""
    if pl is None: prody.importPyPlot()
    modes = getUniformModel(mode)

    for mode in modes:
        if modes.is3d:
            a3d = mode.getArrayNx3()
            show = pl.plot(a3d[:, 0], *args, label='x-component', **kwargs)
            pl.plot(a3d[:, 1], *args, label='y-component', **kwargs)
            pl.plot(a3d[:, 2], *args, label='z-component', **kwargs)
        else:
            show = pl.plot(mode.getArray(), *args, **kwargs)
    return show

def showSqFlucts(modes, *args, **kwargs):
    """Show square fluctuations using :func:`matplotlib.pyplot.imshow`."""
    if pl is None: prody.importPyPlot()
    modes = getUniformModel(modes)
    sqf = getSqFlucts(modes)
    show = pl.plot(sqf, *args, **kwargs)
    pl.xlabel('Atom indices')
    pl.ylabel('Square fluctuations (A^2)')
    return show

def showContactMap(enm, *args, **kwargs):
    """Show Kirchhoff matrix using :func:`matplotlib.pyplot.spy`."""
    if pl is None: prody.importPyPlot()
    if not isinstance(enm, GNMBase):
        raise TypeError('model argument must be an ENM instance')
    kirchhoff = enm.getKirchhoff()
    if kirchhoff is None:
        LOGGER.warning('kirchhoff matrix is not set')
        return None
    enm = getUniformModel(enm)
    show = pl.spy(kirchhoff, *args, **kwargs)
    pl.title('{0:s} connectivity'.format(str(enm))) 
    pl.xlabel('Residue index')
    pl.ylabel('Residue index')
    return show

def showOverlap(mode, modes, *args, **kwargs):
    """Show overlap :func:`matplotlib.pyplot.bar`.
    
    :arg mode: a single mode/vector
    :type mode: :class:`Mode`, :class:`Vector` 
    :arg modes: multiple modes
    :type modes: :class:`ModeSet`, :class:`ANM`, :class:`GNM`, or :class:`PCA` 
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0:s}'.format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    overlap = abs(getOverlap(mode, modes))
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+1.5)
    else:
        arange = modes.getIndices() + 0.5
    show = pl.bar(arange, overlap, *args, **kwargs)
            
    pl.title('Overlap: {0:s} & {1:s}'.format(str(mode), str(modes)))
    pl.xlabel('Mode index')
    pl.ylabel('Overlap')
    return show
    
def showCumulativeOverlap(mode, modes, *args, **kwargs):
    """Show cumulative overlap using :func:`matplotlib.pyplot.plot`.
   :type mode: :class:`Mode`, :class:`Vector` 
    :arg modes: multiple modes
    :type modes: :class:`ModeSet`, :class:`ANM`, :class:`GNM`, or :class:`PCA` 
    
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0:s}'.format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    cumov = (getOverlap(mode, modes) ** 2).cumsum(1) ** 0.5
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+1.5)
    else:
        arange = modes.getIndices() + 0.5
    show = pl.plot(indices, cumov, *args, **kwargs)
    pl.title('Cumulative overlap: {0:s} & {1:s}'.format(str(mode), str(modes)))
    pl.xlabel('Mode index')
    pl.ylabel('Cumulative overlap')
    pl.axis((indices[0]-0.5, indices[-1]+0.5, 0, 1))
    return show
    
