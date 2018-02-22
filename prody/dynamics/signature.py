# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import numpy as np

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure, showMatrix, copy
from prody.ensemble import Ensemble, Conformation

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes

from .analysis import calcSqFlucts, calcCrossCorr
from .plotting import showAtomicData, showAtomicMatrix
from .anm import ANM
from .gnm import GNM

__all__ = ['Signature', 'calcEnsembleENMs', 'getSignatureProfile', 'calcEnsembleSpectralOverlaps',
           'showSignatureProfile', 'calcAverageCrossCorr', 'showAverageCrossCorr', 'showMatrixAverageCrossCorr']

class Signature(object):
    """
    A class for signature dynamics calculated from an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. The class is a collection of column vectors, 
    and each is associated with a value.
    """

    __slots__ = ['_array', '_vars', '_title', '_is3d']

    def __init__(self, vecs, vals=None, title=None, is3d=False):
        vecs = np.array(vecs)
        self._vars = vals

        ndim = vecs.ndim
        if ndim == 1:
            vecs = vecs[:, np.newaxis]
        elif ndim > 2:
            raise ValueError('array can be only 1-D or 2-D')
        self._array = vecs
        shape = vecs.shape

        if vals is None:
            vals = np.zeros(shape[1])
        self._vars = vals

        if title is None:
            self._title = '%d vectors of length %d'%shape
        else:
            self._title = title

        self._is3d = is3d

    def __len__(self):
        """Returns the number of vectors."""
        return self._array.shape[1]

    def __iter__(self):
        for mode in zip(self._array.T, self._vars):
            yield mode

    def __repr__(self):
        return '<Signature: %d vectors of length %d>'%self._array.shape

    def __str__(self):
        return self.getTitle()
    
    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        vecs = self._array[:, index]
        vals = self._vars[index]
        return Signature(vecs, vals, is3d=self.is3d)

    def is3d(self):
        """Returns **True** is model is 3-dimensional."""
        
        return self._is3d

    def numAtoms(self):
        """Returns number of atoms."""

        return self._array.shape[0]

    def numVectors(self):
        """Returns number of modes in the instance (not necessarily maximum
        number of possible modes)."""

        return len(self)

    numModes = numVectors

    def getTitle(self):
        """Returns title of the signature."""

        return self._title

    def getValues(self):
        """Returns variances of vectors. """

        return self._vars.copy()

    def getArray(self):
        """Returns a copy of vectors."""

        return self._array.copy()

    getVectors = getArray

    def _getArray(self):
        """Returns vectors."""

        return self._array

    def getMean(self):
        return self._array.mean(axis=1)
    mean = getMean

    def getVariance(self):
        return self._array.var(axis=1)
    var = getVariance
    
    def getStd(self):
        return self._array.std(axis=1)
    std = getStd

    def getMin(self):
        return self._array.min(axis=1)
    min = getMin

    def getMax(self):
        return self._array.max(axis=1)
    max = getMax

def calcEnsembleENMs(ensemble, model='gnm', trim='trim', n_modes=20):
    """Description"""

    if model is GNM:
        model_type = 'GNM'
    elif model is ANM:
        model_type = 'ANM'
    else:
        model_type = str(model).strip().upper()

    atoms = ensemble.getAtoms()
    select = None
    if ensemble._indices is not None:
        select = atoms
        atoms = atoms.getAtomGroup()
        
    labels = ensemble.getLabels()

    verb = LOGGER.verbosity
    LOGGER.verbosity = 'info'
    ### ENMs ###
    ## ENM for every conf
    enms = []
    n_confs = ensemble.numConfs()
    LOGGER.progress('Calculating {0} {1} modes for {2} conformations...'
                    .format(n_modes, model_type, n_confs), n_confs)

    for i in range(n_confs):
        coords = ensemble.getCoordsets(i, selected=False)
        if atoms is not None:
            atoms.setCoords(coords)
        else:
            atoms = coords
        enm, _ = calcENM(atoms, select, model=model, trim=trim, 
                            n_modes=n_modes, title=labels[i])
        enms.append(enm)

        lbl = labels[i] if labels[i] != '' else '%d-th conformation'%(i+1)
        LOGGER.update(i)
    
    LOGGER.update(n_confs, 'Finished.')
    LOGGER.verbosity = verb

    LOGGER.info('{0} {1} modes were calculated for each of the {2} conformations.'
                        .format(n_modes, model_type, n_confs))
    return enms

def _getEnsembleENMs(ensemble, **kwargs):
    if isinstance(ensemble, Ensemble):
        enms = calcEnsembleENMs(ensemble, **kwargs)
    if isinstance(ensemble, Conformation):
        enms = calcEnsembleENMs([ensemble], **kwargs)
    else:
        try:
            enms = []
            for enm in ensemble:
                if not isinstance(enm, (Mode, NMA, ModeSet)):
                    raise TypeError('ensemble can be a list of Mode, '
                                    'NMA, or ModeSet instances, '
                                    'not {0}'.format(type(enm)))
                enms.append(enm)
        except TypeError:
            raise TypeError('ensemble must be an Ensemble instance, '
                            'or a list of NMA, Mode, or ModeSet instances.')
    return enms

def calcEnsembleSpectralOverlaps(ensemble, distance=False, **kwargs):
    """Description"""

    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    overlaps = np.zeros((len(enms), len(enms)))
    for i, enmi in enumerate(enms):
        for j, enmj in enumerate(enms):
            covlap = calcSpectralOverlap(enmi, enmj)
            overlaps[i, j] = covlap

    if distance:
        overlaps = np.arccos(overlaps)

    return overlaps

def getSignatureProfile(ensemble, index, **kwargs):
    """
    Get the signature profile of *ensemble*. If *ensemble* is an instance of 
    :class:`Ensemble` then the ENMs will be first calculated using 
    :func:`calcEnsembleENMs`. 
    
    :arg ensemble: an ensemble of structures or ENMs 
    :type ensemble: :class: `Ensemble` or list

    :arg index: mode index for displaying the mode shape or a list 
                of mode indices for displaying the mean square fluctuations. 
                The list can contain only one index.
    :type index: int or list
    """

    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    matches = matchModes(*enms)

    V = []; W = []; is3d = None
    if np.isscalar(index):
        modes = matches[index]
        v0 = modes[0].getEigvec()
        for mode in modes:
            v = mode.getEigvec()
            c = np.dot(v, v0)
            if c < 0:
                v *= -1
            w = mode.getVariance()
            V.append(v); W.append(w)
    else:
        for j in range(len(enms)):
            model = None
            indices = []
            for i in index:
                mode = matches[i][j]
                indices.append(mode.getIndex())
                if model is None: 
                    model = mode.getModel()

            modes = ModeSet(model, indices)
            sqfs = calcSqFlucts(modes)
            vars = modes.getVariances()
            V.append(sqfs); W.append(np.sum(vars))
    V = np.vstack(V); V = V.T

    try:
        title = ensemble.getTitle()
    except AttributeError:
        title = None
    sig = Signature(V, W, title=title, is3d=mode.is3d())

    return sig
    
def showSignatureProfile(ensemble, index, linespec='-', **kwargs):
    """
    Show the signature profile of *ensemble* using :func:`showAtomicData`. 
    
    :arg ensemble: an ensemble of structures or ENMs 
    :type ensemble: :class: `Ensemble` or list

    :arg index: mode index for displaying the mode shape or a list 
                of mode indices for displaying the mean square fluctuations. 
                The list can contain only one index.
    :type index: int or list

    :arg atoms: an object with method :func:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 

    :arg alpha: the transparency of the band(s).
    :type alpha: float
    """

    from matplotlib.pyplot import figure, plot, fill_between, gca
    from .signature import getSignatureProfile

    V = getSignatureProfile(ensemble, index, **kwargs)
    meanV, stdV, minV, maxV = V.mean(), V.std(), V.min(), V.max()
    x = range(meanV.shape[0])
    
    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = ensemble.getAtoms()
        except:
            pass

    zero_line = kwargs.pop('show_zero', None)
    if zero_line is None:
        zero_line = np.isscalar(index)
    lines, _, bars, _ = showAtomicData(meanV, atoms=atoms, linespec=linespec, show_zero=zero_line, **kwargs)
    line = lines[-1]
    color = line.get_color()
    x, _ = line.get_data()
    polys = []
    poly = fill_between(x, minV, maxV,
                        alpha=0.3, facecolor=color,
                        linewidth=1, antialiased=True)
    polys.append(poly)
    poly = fill_between(x, meanV-stdV, meanV+stdV,
                        alpha=0.5, facecolor=color,
                        linewidth=1, antialiased=True)
    polys.append(poly)
    return lines, polys, bars
    
def calcAverageCrossCorr(modeEnsemble, modeIndex, *args, **kwargs):
    """Calculate average cross-correlations for a modeEnsemble (a list of modes)."""
    
    matches = matchModes(*modeEnsemble, index=True)
    n_sets = len(modeEnsemble)
    CCs = []
    for i in range(n_sets):
        CC = calcCrossCorr(modeEnsemble[i][modeIndex])
        CCs.append(CC)
    C = np.vstack(CCs)
    n_atoms = modeEnsemble[0].numAtoms()
    C = C.reshape(len(CCs), n_atoms, n_atoms)
    mean = C.mean(axis=0)
    std = C.std(axis=0)
        
    return C, mean, std

def showAverageCrossCorr(modeEnsemble, modeIndex, plotStd=False, *args, **kwargs):
    """Show average cross-correlations using :func:`~matplotlib.pyplot.imshow`.  By
    default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcAverageCrossCorr`."""

    import matplotlib.pyplot as plt
    if SETTINGS['auto_show']:
        plt.figure()
    arange = np.arange(modeEnsemble[0].numAtoms())
    C, mean, std = calcAverageCrossCorr(modeEnsemble, modeIndex)
    if plotStd:
        matrixData = std
    else:
        matrixData = mean
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    cmap = kwargs.pop('cmap', 'jet')
    show = showAtomicMatrix(matrixData, cmap=cmap, *args, **kwargs), plt.colorbar()
    if np.isscalar(modeIndex):
        title_str = ', mode '+str(modeIndex+1)
    else:
        modeIndexStr = ','.join([str(x+1) for x in modeIndex])
        if len(modeIndexStr) > 8:
            title_str = ', '+str(len(modeIndex))+' modes '+modeIndexStr[:5]+'...'
        else:
            title_str = ', modes '+modeIndexStr
        # title_str = ', '+str(len(modeIndex))+' modes'
    if plotStd:
        plt.title('Std - Cross-correlations'+title_str, size=14)
    else:
        plt.title('Avg - Cross-correlations'+title_str, size=14)
    plt.xlabel('Indices', size=14)
    plt.ylabel('Indices', size=14)
    if SETTINGS['auto_show']:
        showFigure()
    return show

def showMatrixAverageCrossCorr(modeEnsemble, modeIndex, plotStd=False, *args, **kwargs):
    """Show average cross-correlations using :func:`showMatrix`.  By
    default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcAverageCrossCorr`."""

    import matplotlib.pyplot as plt
    C, mean, std = calcAverageCrossCorr(modeEnsemble, modeIndex)
    if plotStd:
        matrixData = std
    else:
        matrixData = mean
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    cmap = kwargs.pop('cmap', 'jet')
    show = showMatrix(matrixData, cmap=cmap, *args, **kwargs)
    if np.isscalar(modeIndex):
        title_str = ', mode '+str(modeIndex+1)
    else:
        # modeIndexStr = ','.join([str(x+1) for x in modeIndex])
        # if len(modeIndexStr) > 8:
            # title_str = ', '+str(len(modeIndex))+' modes '+modeIndexStr[:5]+'...'
        # else:
            # title_str = ', modes '+modeIndexStr
        title_str = ', '+str(len(modeIndex))+' modes'
    if plotStd:
        plt.title('Std - Cross-correlations'+title_str, size=12)
    else:
        plt.title('Avg - Cross-correlations'+title_str, size=12)
    plt.xlabel('Indices', size=14)
    plt.ylabel('Indices', size=14)
    if SETTINGS['auto_show']:
        showFigure()
    return show
