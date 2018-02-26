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

__all__ = ['ModeEnsemble', 'Signature', 'calcEnsembleENMs', 'showSignature', 'showSignatureMode', 
           'showSignatureSqFlucts', 'calcEnsembleSpectralOverlaps', 'calcSignatureSqFlucts', 
           'calcSignatureCrossCorr', 'showSignatureCrossCorr', 'showVarianceBar',
           'showSignatureVariances']

class ModeEnsemble(object):
    """
    A class for ENMs calculated for conformations in an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. 
    """

    __slots__ = ['_modesets', '_title', '_labels', '_atoms']

    def __init__(self, title=None):
        self._modesets = []
        self._title = 'Unknown' if title is None else title
        self._labels = None
        self._atoms = None

    def __len__(self):
        """Returns the number of vectors."""

        return len(self._modesets)

    def __iter__(self):
        for modeset in self._modesets:
            yield modeset

    def __repr__(self):
        return '<ModeEnsemble: {0} modesets ({1} modes, {2} atoms)>'\
                .format(len(self), self.numModes(), self.numAtoms())

    def __str__(self):
        return self.getTitle()
    
    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        if isinstance(index, slice):
            modesets = self._modesets[index]
            labels = self._labels[index] if self._labels else None
        elif isinstance(index, (list, tuple)):
            modesets = []; labels = []
            for i in index:
                assert isinstance(i, int), 'all indices must be integers'
                modesets.append(self._modesets[i])
                if self._labels is not None:
                    labels.append(self._labels[i])
        else:
            try:
                index = int(index)
            except Exception:
                raise IndexError('indices must be int, slice, list, or tuple')
            else:
                return self._modesets[index]
        
        ens = ModeEnsemble(title=self.getTitle())
        ens.addModeSet(modesets, label=labels)
        return ens

    def __add__(self, other):
        """Concatenate two mode ensembles. """

        if not isinstance(other, ModeEnsemble):
            raise TypeError('an ModeEnsemble instance cannot be added to an {0} '
                            'instance'.format(type(other)))
        if self.numAtoms() != other.numAtoms():
            raise ValueError('Ensembles must have same number of atoms.')
        if self.numModes() != other.numModes():
            raise ValueError('Ensembles must have same number of modes.')

        ensemble = ModeEnsemble('{0} + {1}'.format(self.getTitle(),
                                                  other.getTitle()))
        ensemble.addModeSet(self._modesets, self._labels)
        ensemble.setAtoms(self.getAtoms())
        
        ensemble.addModeSet(other.getModeSets(), label=other.getLabels())
        return ensemble

    def is3d(self):
        """Returns **True** is model is 3-dimensional."""
        
        if self._modesets:
            return self._modesets[0].is3d()
        return False

    def numAtoms(self):
        """Returns number of atoms."""

        if self._modesets:
            return self._modesets[0].numAtoms()
        return 0

    def numModes(self):
        """Returns number of modes in the instance (not necessarily maximum
        number of possible modes)."""

        if self._modesets:
            return self._modesets[0].numModes()
        return 0

    def numModeSets(self):
        """Returns number of modesets in the instance."""

        return len(self)

    def getTitle(self):
        """Returns title of the signature."""

        return self._title

    def getModeSets(self, index=None):
        if index is None:
            return self._modesets
        return self[index]._modesets

    def getArray(self, mode_index=0):
        """Returns a copy of row vectors."""

        if np.isscalar(mode_index):
            sig = self.getEigvec(mode_index)
        else:
            sig = self.getEigvecs(mode_index)
        return sig
    
    _getArray = getArray
    
    def _getModeData(self, name, mode_index=0, sign_correction=False):
        modesets = self._modesets

        V = []
        for i, modeset in enumerate(modesets):
            mode = modeset[mode_index]
            func = getattr(mode, name)
            v = func()
            if sign_correction:
                if i == 0:
                    v0 = v
                else:
                    c = np.dot(v, v0)
                    if c < 0:
                        v *= -1
            V.append(v)
        is3d = mode.is3d()
        V = np.vstack(V)

        title = 'mode %d'%(mode_index+1)
        sig = Signature(V, title=title, labels=self.getLabels(), is3d=is3d)
        return sig

    def _getData(self, name, mode_indices=None):
        modesets = self._modesets
        V = []
        for modeset in modesets:
            modes = modeset if mode_indices is None else modeset[mode_indices]
            func = getattr(modes, name)
            vecs = func()
            V.append(vecs)
        is3d = modeset.is3d()
        V = np.array(V)

        title = '%d modes'%len(V)
        sig = Signature(V, title=title, labels=self.getLabels(), is3d=is3d)
        return sig

    def getEigvec(self, mode_index=0):
        """Returns a copy of eigenvector."""

        return self._getModeData('getEigvec', mode_index, sign_correction=True)

    def getEigvecs(self, mode_indices=None):
        """Returns a copy of eigenvectors."""

        return self._getData('getEigvecs', mode_indices)

    def getVariance(self, mode_index=0):
        
        return self._getModeData('getVariance', mode_index)

    def getVariances(self, mode_indices=None):
        
        return self._getData('getVariances', mode_indices)

    def getIndex(self, mode_index=0):
        """Returns indices of modes matched to the reference modeset."""

        return self._getModeData('getIndex', mode_index)

    def getIndices(self, mode_indices=None):
        """Returns indices of modes in the mode ensemble."""
        
        return self._getData('getIndices', mode_indices)

    def getAtoms(self):
        return self._atoms

    def setAtoms(self, atoms):
        if atoms is not None and self._atoms is not None:
            if len(atoms) != self.numAtoms():
                raise ValueError('atoms should have %d atoms'%self.numAtoms())
        self._atoms = atoms

    def getLabels(self):
        return self._labels

    def match(self):
        if self._modesets:
            self._modesets = matchModes(*self._modesets)
        return

    def addModeSet(self, modeset, label=None):
        if isinstance(modeset, (NMA, ModeSet)):
            modesets = [modeset]
        else:
            modesets = modeset

        if label is not None:
            if np.isscalar(label):
                labels = [label]
            else:
                labels = label
            if len(labels) != len(modesets):
                raise ValueError('labels should have the same length as modesets')
        
        if not self._labels and labels:
            self._labels = ['']*len(self._modesets)
            self._labels.extend(labels)

        if not labels and self._labels:
            labels = ['']*len(modesets)
            self._labels.extend(labels)

        for i in range(len(modesets)):
            modeset = modesets[i]
            if isinstance(modeset, NMA):
                modeset = modeset[:]
            if not isinstance(modeset, ModeSet):
                raise TypeError('modesets should be a list of ModeSet instances')
            if self._modesets:
                if modeset.numAtoms() != self.numAtoms():
                    raise ValueError('to be added, modesets should contain exactly %d atoms'%self.numAtoms())
                if modeset.numModes() < self.numModes():
                    raise ValueError('to be added, modesets should contain at least %d modes'%self.numModes())
                self._modesets.append(modeset[:self.numModes()])
            else:
                self._modesets = [modeset]


class Signature(object):
    """
    A class for signature dynamics calculated from an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. The class is a collection of row vectors, and 
    each is associated with a value. The vectors can have more than one 
    dimension.
    """

    __slots__ = ['_array', '_title', '_labels', '_is3d']

    def __init__(self, vecs, labels=None, title=None, is3d=False):
        vecs = np.array(vecs)

        self._array = vecs
        shape = vecs.shape

        if title is None:
            self._title = '{0} vectors of size {1}'.format(shape[0], shape[1:])
        else:
            self._title = title
        
        self._labels = labels
        self._is3d = is3d

    def __len__(self):
        """Returns the number of vectors."""
        return self._array.shape[0]

    def __iter__(self):
        for mode in self._array:
            yield mode

    def __repr__(self):
        shape = self._array.shape
        return '<Signature: {0} vectors of size {1}>'.format(shape[0], shape[1:])

    def __str__(self):
        return self.getTitle()
    
    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        vecs = self._array[index]
        labels = self._labels
        if np.isscalar(index):
            vecs = np.array([vecs])
        if labels is not None:
            labels = list(np.array(labels)[index])
        return Signature(vecs, labels=labels, is3d=self.is3d)

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

    def getLabels(self):
        return self._labels

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

    str_modes = 'all' if n_modes is None else str(n_modes)
    LOGGER.progress('Calculating {0} {1} modes for {2} conformations...'
                    .format(str_modes, model_type, n_confs), n_confs)

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
                        .format(str_modes, model_type, n_confs))

    modeens = ModeEnsemble(title=ensemble.getTitle())
    modeens.addModeSet(enms, label=ensemble.getLabels())
    modeens.setAtoms(ensemble.getAtoms())
    return modeens

def _getEnsembleENMs(ensemble, **kwargs):
    if isinstance(ensemble, Ensemble):
        enms = calcEnsembleENMs(ensemble, **kwargs)
    elif isinstance(ensemble, Conformation):
        enms = calcEnsembleENMs([ensemble], **kwargs)
    elif isinstance(ensemble, ModeEnsemble):
        enms = ensemble
    else:
        try:
            enms = ModeEnsemble()
            enms.addModeSet(ensemble)
        except TypeError:
            raise TypeError('ensemble must be an Ensemble or a ModeEnsemble instance,'
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

def calcSignatureSqFlucts(mode_ensemble, indices, **kwargs):
    """
    Get the signature square fluctuations of *mode_ensemble*. 
    
    :arg mode_ensemble: an ensemble of structures or ENMs 
    :type mode_ensemble: :class: `ModeEnsemble`

    :arg indices: mode indices for displaying the mean square fluctuations. 
    :type indices: int or list
    """

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')
    modesets = mode_ensemble
    modesets.match()

    V = []
    for modeset in modesets:
        modes = modeset[indices]
        sqfs = calcSqFlucts(modes)
        V.append(sqfs)
    is3d = modeset.is3d()
    V = np.vstack(V)

    title_str = '%d modes'%len(indices)
    sig = Signature(V, title=title_str, is3d=is3d)

    return sig
    
def showSignature(signature, linespec='-', **kwargs):
    """
    Show the signature dynamics using :func:`showAtomicData`. 
    
    :arg signature: the signature dynamics to be plotted 
    :type signature: :class:`Signature`

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

    V = signature
        
    meanV, stdV, minV, maxV = V.mean(), V.std(), V.min(), V.max()
    x = range(meanV.shape[0])
    
    atoms = kwargs.pop('atoms', None)

    zero_line = kwargs.pop('show_zero', False)
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
    title('Signature profile of ' + V.getTitle())

    return lines, polys, bars

def showSignatureMode(mode_ensemble, index):
    mode = mode_ensemble.getEigvec(index)
    return showSignature(mode, atoms=mode_ensemble.getAtoms(), show_zero=True)

def showSignatureSqFlucts(mode_ensemble, indices):
    sqf = calcSignatureSqFlucts(mode_ensemble, indices)
    return showSignature(sqf, atoms=mode_ensemble.getAtoms(), show_zero=False)

def calcSignatureCrossCorr(mode_ensemble, index):
    """Calculate average cross-correlations for a modeEnsemble (a list of modes)."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    mode_ensemble.match()
    matches = mode_ensemble
    n_atoms = matches.numAtoms()
    n_sets = len(matches)

    C = np.zeros((n_sets, n_atoms, n_atoms))
    is3d = None
    for i in range(n_sets):
        m = matches[i][index]
        c = calcCrossCorr(m)
        C[i, :, :] = c
        
        if is3d is None:
            is3d = m.is3d()

    sig = Signature(C, title=mode_ensemble.getTitle(), is3d=is3d)
        
    return sig

def calcSignatureFractVariance(mode_ensemble, index):
    """Calculate average cross-correlations for a modeEnsemble (a list of modes)."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    mode_ensemble.match()
    matches = mode_ensemble
    n_sets = len(matches)

    W = []; is3d = None
    for i in range(n_sets):
        m = matches[i][index]
        var = calcFractVariance(m)
        if np.isscalar(index):
            var = var.sum()
        W.append(var)
        if is3d is None:
            is3d = m.is3d()

    sig = Signature(W, title=mode_ensemble.getTitle(), is3d=is3d)
        
    return sig

def showSignatureCrossCorr(mode_ensemble, index, show_std=False, **kwargs):
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
    
    C = calcSignatureCrossCorr(mode_ensemble, index, **kwargs)

    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = mode_ensemble.getAtoms()
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

def showSignatureVariances(mode_ensemble, index, **kwargs):
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

    fract = kwargs.pop('fraction', True)
    show_legend = kwargs.pop('legend', True)

    if fract:
        sig = calcSignatureFractVariance(mode_ensemble, index)
    else:
        sig = mode_ensemble.getVariances(index) 
    W = sig.getArray()[:, ::-1] # reversed to accommodate with matplotlib.pyplot.hist
    weights = np.ones_like(W)/float(len(W))
    if np.isscalar(index):
        index = [index]
    legends = ['mode %d'%(i+1) for i in index][::-1]

    bins = kwargs.pop('bins', 'auto')
    if bins == 'auto':
        _, bins = np.histogram(W.flatten(), bins='auto')
    elif np.isscalar(bins) and isinstance(bins, (int, np.integer)):
        step = (W.max() - W.min())/bins
        bins = np.arange(W.min(), W.max(), step)

    histtype = kwargs.pop('histtype', 'stepfilled')
    label = kwargs.pop('label', legends)
    weights = kwargs.pop('weights', weights)
    n, bins, patches = hist(W, bins=bins, weights=weights, 
                            histtype=histtype, label=label, **kwargs)
    if show_legend:
        legend()

    xlabel('Variance')
    ylabel('Probability')

    if SETTINGS['auto_show']:
        showFigure()

    return n, bins, patches

def showVarianceBar(variances, labels=None, **kwargs):

    from matplotlib.pyplot import figure, gca
    from matplotlib.figure import Figure
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize
    from matplotlib import cm, colors
    
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

    meanVar = variances.mean()
    stdVar = variances.std()
    
    variances = (variances - meanVar)/stdVar

    maxVar = variances.max()
    minVar = variances.min()

    # Set the colormap and norm to correspond to the data for which
    # the colorbar will be used.
    cmap = kwargs.pop('cmap', 'jet')
    norm = Normalize(vmin=minVar, vmax=maxVar)

    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    cb = ColorbarBase(gca(), cmap=cmap, norm=norm,
                      orientation='horizontal')
    cb.set_label('Variances')

    if SETTINGS['auto_show']:
        showFigure()
    return cb
