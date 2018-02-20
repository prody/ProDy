# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import numpy as np

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure, showMatrix
from prody.ensemble import Ensemble

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes

from .analysis import calcSqFlucts, calcCrossCorr
from .plotting import showAtomicData, showAtomicMatrix
from .anm import ANM
from .gnm import GNM

__all__ = ['calcEnsembleENMs', 'getSignatureProfile', 'calcEnsembleSpectralOverlaps',
           'showSignatureProfile', 'calcAverageCrossCorr', 'showAverageCrossCorr', 'showMatrixAverageCrossCorr']

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
        atoms.setCoords(ensemble.getCoordsets(i, selected=False))
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

    if np.isscalar(index):
        modes = matches[index]
        V = []
        v0 = modes[0].getEigvec()
        for mode in modes:
            v = mode.getEigvec()
            c = np.dot(v, v0)
            if c < 0:
                v *= -1
            V.append(v)
    else:
        V = []
        for j in range(len(matches)):
            modes = []
            for i in index:
                mode = matches[i][j]
                modes.append(mode)
            sqfs = calcSqFlucts(modes)
            V.append(sqfs)
    V = np.vstack(V); V = V.T

    meanV = V.mean(axis=1)
    stdV = V.std(axis=1)

    return V, (meanV, stdV)
    
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

    V, (meanV, stdV) = getSignatureProfile(ensemble, index, **kwargs)
    minV = V.min(axis=1)
    maxV = V.max(axis=1)

    x = range(meanV.shape[0])
    
    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = ensemble.getAtoms()
        except:
            pass

    lines, _, bars, _ = showAtomicData(meanV, atoms=atoms, linespec=linespec, **kwargs)
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
