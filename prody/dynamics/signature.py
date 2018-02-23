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

from .analysis import calcSqFlucts, calcCrossCorr, calcFractVariance
from .plotting import showAtomicData, showAtomicMatrix
from .anm import ANM
from .gnm import GNM

__all__ = ['Signature', 'calcEnsembleENMs', 'calcSignatureMobility', 'calcEnsembleSpectralOverlaps',
           'showSignatureMobility', 'calcSignatureCrossCorr', 'showSignatureCrossCorr',
           'showSignatureVariances']

class Signature(object):
    """
    A class for signature dynamics calculated from an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. The class is a collection of row vectors, and 
    each is associated with a value. The vectors can have more than one 
    dimension.
    """

    __slots__ = ['_array', '_vars', '_title', '_is3d']

    def __init__(self, vecs, vars=None, title=None, is3d=False):
        vecs = np.array(vecs)

        ndim = vecs.ndim
        if ndim == 1:
            vecs = np.reshape(1, len(vecs))

        self._array = vecs
        shape = vecs.shape

        if vars is None:
            vars = np.zeros(shape[1])
        if np.isscalar(vars):
            vars = [vars]
        if len(vars) != shape[0]:
            raise ValueError('the number of vars does not match the number of vecs')
        self._vars = np.array(vars)

        if title is None:
            self._title = '{0} vectors of size {1}'.format(shape[0], shape[1:])
        else:
            self._title = title

        self._is3d = is3d

    def __len__(self):
        """Returns the number of vectors."""
        return self._array.shape[0]

    def __iter__(self):
        for mode in zip(self._array, self._vars):
            yield mode

    def __repr__(self):
        shape = self._array.shape
        return '<Signature: {0} vectors of size {1}>'.format(shape[0], shape[1:])

    def __str__(self):
        return self.getTitle()
    
    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        vecs = self._array[index]
        vars = self._vars[index]
        if np.isscalar(index):
            vecs = np.array([vecs])
        return Signature(vecs, vars, is3d=self.is3d)

    def is3d(self):
        """Returns **True** is model is 3-dimensional."""
        
        return self._is3d

    def numAtoms(self):
        """Returns number of atoms."""

        return self._array.shape[1]

    def numVectors(self):
        """Returns number of modes in the instance (not necessarily maximum
        number of possible modes)."""

        return len(self)

    numModes = numVectors

    def getShape(self):
        return self._array.shape[1:]
    shape = getShape

    def getTitle(self):
        """Returns title of the signature."""

        return self._title

    def getVariances(self, index=None):
        """Returns variances of vectors. """
        if index is not None:
            return self._vars[index]
        return self._vars.copy()

    def getArray(self, index=None):
        """Returns a copy of row vectors."""

        if index is not None:
            return self._array[index].copy()

        return self._array.copy()

    getVectors = getArray

    def _getArray(self, index=None):
        """Returns row vectors."""

        if index is not None:
            return self._array[index]

        return self._array

    def getMean(self):
        return self._array.mean(axis=0)
    mean = getMean
    
    def getStd(self):
        return self._array.std(axis=0)
    std = getStd

    def getMin(self):
        return self._array.min(axis=0)
    min = getMin

    def getMax(self):
        return self._array.max(axis=0)
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

        #lbl = labels[i] if labels[i] != '' else '%d-th conformation'%(i+1)
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

def calcSignatureMobility(ensemble, index, **kwargs):
    """
    Get the signature mobility of *ensemble*. If *ensemble* is an instance of 
    :class:`Ensemble` then the ENMs will be first calculated using 
    :func:`calcEnsembleENMs`. 
    
    :arg ensemble: an ensemble of structures or ENMs 
    :type ensemble: :class: `Ensemble` or list

    :arg index: mode index for displaying the mode shape or a list 
                of mode indices for displaying the mean square fluctuations. 
                The list can contain only one index.
    :type index: int or list

    :arg fractional: if set to ``True``, mode variances are normalized, default
                    is ``False``
    :type fractional: bool
    """

    fract = kwargs.pop('fractional', False)
    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    modesets = matchModes(*enms)

    V = []; W = []
    if np.isscalar(index):
        v0 = modesets[0][index].getEigvec()
        for modeset in modesets:
            mode = modeset[index]
            v = mode.getEigvec()
            c = np.dot(v, v0)
            if c < 0:
                v *= -1
            if not fract:
                w = mode.getVariance()
            else:
                w = calcFractVariance(mode)
            V.append(v); W.append(w)
        is3d = mode.is3d()
    else:
        for modeset in modesets:
            modes = modeset[index]
            sqfs = calcSqFlucts(modes)
            if not fract:
                w = modes.getVariances().sum()
            else:
                w = calcFractVariance(modes).cumsum()[-1]
            V.append(sqfs); W.append(w)
        is3d = modeset.is3d()
    V = np.vstack(V)

    try:
        title_str = ensemble.getTitle()
    except AttributeError:
        if np.isscalar(index):
            title_str = 'mode %d'%(index+1)
        else:
            title_str = '%d modes'%len(index)
    sig = Signature(V, W, title=title_str, is3d=is3d)

    return sig
    
def showSignatureMobility(ensemble, index, linespec='-', **kwargs):
    """
    Show the signature mobility of *ensemble* using :func:`showAtomicData`. 
    
    :arg ensemble: an ensemble of structures or ENMs, or a signature profile 
    :type ensemble: :class: `Ensemble`, list, :class:`Signature`

    :arg index: mode index for displaying the mode shape or a list 
                of mode indices for displaying the mean square fluctuations. 
                The list can contain only one index.
    :type index: int or list

    :arg linespec: line specifications that will be passed to :func:`showAtomicData`
    :type linespec: str

    :arg atoms: an object with method :func:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 

    :arg alpha: the transparency of the band(s).
    :type alpha: float
    """

    from matplotlib.pyplot import figure, plot, fill_between, \
                                  gca, xlabel, ylabel, title

    if isinstance(ensemble, Signature):
        V = ensemble
    else:
        V = calcSignatureMobility(ensemble, index, **kwargs)

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

    xlabel('Residues')

    if isinstance(ensemble, Ensemble):
        title_str = ensemble.getTitle()
    else:
        if np.isscalar(index):
            title_str = 'mode %d'%(index+1)
        else:
            title_str = '%d modes'%len(index)
    title('Signature profile of ' + title_str)

    return lines, polys, bars
    
def calcSignatureCrossCorr(ensemble, index, *args, **kwargs):
    """Calculate average cross-correlations for a modeEnsemble (a list of modes)."""
    
    enms = _getEnsembleENMs(ensemble, **kwargs)
    matches = matchModes(*enms)
    n_atoms = enms[0].numAtoms()
    n_sets = len(enms)

    C = np.zeros((n_sets, n_atoms, n_atoms))
    W = []; is3d = None
    for i in range(n_sets):
        m = matches[i][index]
        c = calcCrossCorr(m)
        C[i, :, :] = c
        if np.isscalar(index):
            var = m.getVariance()
        else:
            var = np.sum(m.getVariances())
        W.append(var)
        if is3d is None:
            is3d = m.is3d()
    
    try:
        title = ensemble.getTitle()
    except AttributeError:
        title = None

    sig = Signature(C, W, title=title, is3d=is3d)
        
    return sig

def showSignatureCrossCorr(ensemble, index, show_std=False, **kwargs):
    """Show average cross-correlations using :func:`showAtomicMatrix`. 
    By default, *origin=lower* and *interpolation=bilinear* keyword  arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcSignatureCrossCorr`.
    
    :arg ensemble: an ensemble of structures or ENMs, or a signature profile 
    :type ensemble: :class: `Ensemble`, list, :class:`Signature`

    :arg index: mode index for displaying the mode shape or a list 
                of mode indices for displaying the mean square fluctuations. 
                The list can contain only one index.
    :type index: int or list

    :arg atoms: an object with method :func:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 
    """

    import matplotlib.pyplot as plt
    
    if isinstance(ensemble, Signature):
        C = ensemble
    else:
        C = calcSignatureCrossCorr(ensemble, index, **kwargs)

    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = ensemble.getAtoms()
        except:
            pass

    if show_std:
        matrixData = C.std()
    else:
        matrixData = C.mean()
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'

    show = showAtomicMatrix(matrixData, atoms=atoms, **kwargs)
    if np.isscalar(index):
        title_str = ', mode '+str(index+1)
    else:
        modeIndexStr = ','.join([str(x+1) for x in index])
        if len(modeIndexStr) > 8:
            title_str = ', '+str(len(index))+' modes '+modeIndexStr[:5]+'...'
        else:
            title_str = ', modes '+modeIndexStr
        # title_str = ', '+str(len(modeIndex))+' modes'
    if show_std:
        plt.title('Cross-correlations (standard deviation)'+title_str)
    else:
        plt.title('Cross-correlations (average)'+title_str)
    plt.xlabel('Residues')
    plt.ylabel('Residues')
    
    return show

def showSignatureVariances(*signatures, **kwargs):
    """
    Show the distribution of signature variances using 
    :func:`~matplotlib.pyplot.hist`.
    """
    
    from matplotlib.pyplot import figure, hist, annotate, legend, xlabel, ylabel
    from matplotlib.figure import Figure

    fig = kwargs.pop('figure', None)

    if isinstance(fig, Figure):
        fig_num = fig.number
    elif fig is None or isinstance(fig, (int, str)):
        fig_num = fig
    else:
        raise TypeError('figure can be either an instance of matplotlib.figure.Figure '
                        'or a figure number.')
    if SETTINGS['auto_show']:
        figure(fig_num)
    elif fig_num is not None:
        figure(fig_num)

    show_legend = kwargs.pop('legend', True)

    W = []; legends = []; weights = []
    for signature in signatures:
        vars = signature.getVariances()
        W.append(vars)
        legends.append(signature.getTitle())
        weight = np.ones_like(vars)/float(len(vars))
        weights.append(weight)

    W = np.vstack(W[::-1])  # reversed to accommodate with matplotlib.pyplot.hist
    weights = np.vstack(weights[::-1]) 
    legends = legends[::-1]

    bins = kwargs.pop('bins', 'auto')
    if bins == 'auto':
        _, bins = np.histogram(W.flatten(), bins='auto')
    elif np.isscalar(bins) and isinstance(bins, (int, np.integer)):
        step = (W.max() - W.min())/bins
        bins = np.arange(W.min(), W.max(), step)

    histtype = kwargs.pop('histtype', 'stepfilled')
    label = kwargs.pop('label', legends)
    weights = kwargs.pop('weights', weights)
    n, bins, patches = hist(W.T, bins=bins, weights=weights.T, 
                            histtype=histtype, label=label, **kwargs)
    if show_legend:
        legend()

    xlabel('Variance')
    ylabel('Probability')

    if SETTINGS['auto_show']:
        showFigure()

    return n, bins, patches
