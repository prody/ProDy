# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import numpy as np

from prody import LOGGER
from prody.ensemble import Ensemble

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes
from .analysis import calcSqFlucts

__all__ = ['calcEnsembleENMs', 'getSignatureProfile', 'calcSpectralDistances',
           'showSignatureProfile']

def calcEnsembleENMs(ensemble, model='gnm', trim='trim', n_modes=20):
    """Description"""

    atoms = ensemble.getAtoms()
    select = None
    if ensemble._indices is not None:
        select = atoms
        atoms = atoms.getAtomGroup()
        
    labels = ensemble.getLabels()

    ### ENMs ###
    ## ENM for every conf
    enms = []
    for i in range(ensemble.numConfs()):
        atoms.setCoords(ensemble.getCoordsets(i, selected=False))
        enm, _ = calcENM(atoms, select, model=model, trim=trim, 
                            n_modes=n_modes, title=labels[i])
        enms.append(enm)
    
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

def calcSpectralDistances(ensemble, **kwargs):
    """Description"""

    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    overlaps = np.zeros((len(enms), len(enms)))
    for i, enmi in enumerate(enms):
        for j, enmj in enumerate(enms):
            covlap = calcSpectralOverlap(enmi, enmj)
            overlaps[i, j] = covlap

    ### build tree based on similarity matrix ###
    dist_mat = np.arccos(overlaps)

    return dist_mat

def getSignatureProfile(ensemble, index, **kwargs):
    """Description"""

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
    """Description"""

    from matplotlib.pyplot import figure, plot, fill_between, gca
    from .signature import getSignatureProfile

    V, (meanV, stdV) = getSignatureProfile(ensemble, index, **kwargs)
    minV = V.min(axis=1)
    maxV = V.max(axis=1)

    x = range(meanV.shape[0])

    new_fig = kwargs.pop('new_fig', True)

    if new_fig:
        figure()
    line = plot(x, meanV, linespec)[0]
    color = line.get_color()
    fill_between(x, minV, maxV,
                alpha=0.3, facecolor=color,
                linewidth=1, antialiased=True)
    fill_between(x, meanV-stdV, meanV+stdV,
                alpha=0.5, facecolor=color,
                linewidth=1, antialiased=True)
    
    if SETTINGS['auto_show']:
        showFigure()
    return gca()