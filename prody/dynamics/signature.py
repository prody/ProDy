# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import numpy as np

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure
from prody.ensemble import Ensemble

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes

from .analysis import calcSqFlucts, calcCrossCorr
from .plotting import showAtomicData
from .anm import ANM
from .gnm import GNM
from .plotting import showMatrix

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

    new_fig = kwargs.pop('new_fig', True)

    if new_fig:
        figure()
    
    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = ensemble.getAtoms()
        except:
            pass

    ax = showAtomicData(meanV, atoms=atoms, linespec=linespec, new_fig=False, **kwargs)
    line = ax.lines[0]
    color = line.get_color()
    x, _ = line.get_data()
    fill_between(x, minV, maxV,
                alpha=0.3, facecolor=color,
                linewidth=1, antialiased=True)
    fill_between(x, meanV-stdV, meanV+stdV,
                alpha=0.5, facecolor=color,
                linewidth=1, antialiased=True)
    
    if SETTINGS['auto_show']:
        showFigure()
    return gca()
    
def calcAverageCrossCorr(modesEnsemble, modeIndex, *args, **kwargs):
    """Calculate average cross-correlations for a modesEnsemble (a list of modes)."""
    
    matches = matchModes(*modesEnsemble)
    CCs = []
    for mode_i in matches[modeIndex]:
        CC = calcCrossCorr(mode_i)
        CCs.append(CC)
    C = np.vstack(CCs)

    n_atoms = modesEnsemble[0].numAtoms()
    C = C.reshape(len(CCs), n_atoms, n_atoms)
    mean = C.mean(axis=0)
    std = C.std(axis=0)
    return C, mean, std

def showAverageCrossCorr(modesEnsemble, modeIndex, plotStd=False, *args, **kwargs):
    """Show average cross-correlations using :func:`~matplotlib.pyplot.imshow`.  By
    default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcAverageCrossCorr`."""

    import matplotlib.pyplot as plt
    if kwargs.pop('new_fig', True):
        plt.figure()
    arange = np.arange(modesEnsemble[0].numAtoms())
    C, mean, std = calcAverageCrossCorr(modesEnsemble, modeIndex)
    if plotStd:
        matrixData = std
    else:
        matrixData = mean
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    show = plt.imshow(matrixData, *args, **kwargs), plt.colorbar()
    plt.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    if plotStd:
        plt.title('Std - Average Cross-correlations')
    else:
        plt.title('Average Cross-correlations')
    plt.xlabel('Indices')
    plt.ylabel('Indices')
    if SETTINGS['auto_show']:
        showFigure()
    return show

def showMatrixAverageCrossCorr(modesEnsemble, modeIndex, plotStd=False, *args, **kwargs):
    """Show average cross-correlations using :func:`showMatrix`.  By
    default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcAverageCrossCorr`."""

    C, mean, std = calcAverageCrossCorr(modesEnsemble, modeIndex)
    if plotStd:
        matrixData = std
    else:
        matrixData = mean
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'
    if not 'origin' in kwargs:
        kwargs['origin'] = 'lower'
    show = showMatrix(matrixData, *args, **kwargs)
    # if plotStd:
        # title('Std - Average Cross-correlations')
    # else:
        # title('Average Cross-correlations')
    # xlabel('Indices')
    # ylabel('Indices')
    if SETTINGS['auto_show']:
        showFigure()
    return show
