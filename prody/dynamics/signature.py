# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import time
from numbers import Integral
from numpy import ndarray
import numpy as np

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure, showMatrix, copy, checkWeights, openFile, getValue
from prody.ensemble import Ensemble, Conformation

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes, calcOverlap

from .analysis import calcSqFlucts, calcCrossCorr, calcFractVariance, calcCollectivity
from .plotting import showAtomicLines, showAtomicMatrix
from .anm import ANM
from .gnm import GNM

__all__ = ['ModeEnsemble', 'sdarray', 'calcEnsembleENMs', 'showSignature1D', 'showSignatureAtomicLines', 
           'showSignatureMode', 'showSignatureDistribution', 'showSignatureCollectivity',
           'showSignatureSqFlucts', 'calcEnsembleSpectralOverlaps', 'calcSignatureSqFlucts', 
           'calcSignatureCollectivity', 'calcSignatureFractVariance',
           'calcSignatureCrossCorr', 'showSignatureCrossCorr', 'showVarianceBar',
           'showSignatureVariances', 'calcSignatureOverlaps', 'showSignatureOverlaps',
           'saveModeEnsemble', 'loadModeEnsemble']

class ModeEnsemble(object):
    """
    A collection of ENMs calculated for conformations in an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. 
    """

    __slots__ = ['_modesets', '_title', '_labels', '_atoms', '_weights', '_matched']

    def __init__(self, title=None):
        self._modesets = []
        self._title = 'Unknown' if title is None else title
        self._labels = None
        self._atoms = None
        self._weights = None
        self._matched = False

    def __len__(self):
        """Returns the number of modesets."""

        return len(self._modesets)

    def __iter__(self):
        for modeset in self._modesets:
            yield modeset

    iterModeSets = __iter__

    def iterModes(self):
        for i in range(self.numModes()):
            modesets = []
            for modeset in self._modesets:
                modesets.append(modeset[i])
            ens = ModeEnsemble(self.getTitle())
            ens._modesets = modesets
            yield ens

    def __repr__(self):
        return '<ModeEnsemble: {0} modesets ({1} modes, {2} atoms)>'\
                .format(len(self), self.numModes(), self.numAtoms())

    def __str__(self):
        return self.getTitle()
    
    def _getitembase(self, modeset_index, mode_index):

        if isinstance(modeset_index, slice):
            modesets = self._modesets[modeset_index]
            labels = None if self._labels is None else self._labels[modeset_index]
        elif not np.isscalar(modeset_index):
            modesets = []; labels = []
            for i in modeset_index:
                assert isinstance(i, Integral), 'all indices must be integers'
                modesets.append(self._modesets[i])
                if self._labels is not None:
                    labels.append(self._labels[i])
        else:
            try:
                modeset_index = int(modeset_index)
            except Exception:
                raise IndexError('indices must be int, slice, or array-like objects')
            else:
                return self._modesets[modeset_index][mode_index]
        
        if np.isscalar(mode_index):
            mode_index = [mode_index]
        for i in range(len(modesets)):
            modesets[i] = modesets[i][mode_index]

        if self._weights is None:
            weights = None
        else:
            weights = self._weights[modeset_index, :, :]

        ens = ModeEnsemble(title=self.getTitle())
        ens.addModeSet(modesets, weights=weights, label=labels)
        ens.setAtoms(self.getAtoms())
        ens._matched = self._matched
        return ens

    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        modeset_index = slice(None, None, None)
        mode_index = slice(None, None, None)

        if isinstance(index, tuple):
            if len(index) >= 1:
                modeset_index = index[0]
            if len(index) >= 2:
                mode_index = index[1]
            if len(index) >= 3:
                raise ValueError('ModeEnsemble indexing supports up to two arguments')
        else:
            modeset_index = index
        ens = self._getitembase(modeset_index, mode_index)
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
        
        ensemble.addModeSet(other.getModeSets(), weights=other.getWeights(),
                            label=other.getLabels())
        ensemble._matched = self._matched and other._matched
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

    def getWeights(self):
        """Returns a copy of weights."""
        return copy(self._weights)

    def setWeights(self, weights):
        """Set atomic weights."""

        n_atoms = self.numAtoms()
        n_modesets = self.numModeSets()

        if n_atoms > 0 and n_modesets > 0:
            self._weights = checkWeights(weights, n_atoms, n_modesets)
        else:
            raise RuntimeError('modesets must be set before weights can be assigned')

    def getModeSets(self, index=None):
        """
        Returns the modeset of the given index. If index is **None** then all 
        modesets are returned.
        """

        if index is None:
            return self._modesets
        return self[index]._modesets

    def getArray(self, mode_index=0):
        """Returns a sdarray of row arrays."""

        if np.isscalar(mode_index):
            sig = self.getEigvec(mode_index)
        else:
            sig = self.getEigvecs(mode_index)
        return sig
    
    _getArray = getArray
    
    def _getModeData(self, name, mode_index=0, weights=None, sign_correction=False):
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

        mode_num = modesets[0][mode_index].getIndex()
        title = 'mode %d'%(mode_num+1)
        sig = sdarray(V, title=title, weights=weights, labels=self.getLabels(), is3d=is3d)
        return sig

    def _getData(self, name, mode_indices=None, weights=None, sign_correction=False):
        modesets = self._modesets
        V = []
        for i, modeset in enumerate(modesets):
            modes = modeset if mode_indices is None else modeset[mode_indices]
            func = getattr(modes, name)
            vecs = func()
            if sign_correction:
                if i == 0:
                    vecs0 = vecs
                else:
                    c = np.sign(np.dot(vecs.T, vecs0))
                    c = np.diag(np.diag(c))
                    vecs = np.dot(vecs, c)
            V.append(vecs)
        is3d = modeset.is3d()
        V = np.array(V)

        title = '%d modes'%len(V)
        sig = sdarray(V, title=title, weights=weights, labels=self.getLabels(), is3d=is3d)
        return sig

    def getEigvec(self, mode_index=0, sign_correction=True):
        """Returns a sdarray of eigenvector across modesets."""

        weights = self.getWeights()
        if weights is not None:
            weights = weights[:, :, 0]
            if self.is3d():
                weights = np.repeat(weights, 3, axis=1)
        return self._getModeData('getEigvec', mode_index, weights=weights, sign_correction=True)

    def getEigvecs(self, mode_indices=None, sign_correction=True):
        """Returns a sdarray of eigenvectors across modesets."""

        weights = self.getWeights()
        if weights is not None:
            weights = weights[:, :, 0]
            if mode_indices is not None:
                n_modes = len(mode_indices)
            else:
                n_modes = self.numModes()
            W = [weights.copy() for _ in range(n_modes)]  # almost equivalent to W = [weights]*n_modes
            W = np.stack(W, axis=-1)
            if self.is3d():
                W = np.repeat(W, 3, axis=1)

        return self._getData('getEigvecs', mode_indices, weights=W, sign_correction=True)

    def getVariance(self, mode_index=0):
        """Returns variances of a given mode index with respect to the reference."""
        
        return self._getModeData('getVariance', mode_index)

    def getVariances(self, mode_indices=None):
        """Returns a sdarray of variances across modesets."""

        return self._getData('getVariances', mode_indices)

    def getEigval(self, mode_index=0):
        """Returns eigenvalue of a given mode index with respect to the reference."""
        
        return self._getModeData('getEigval', mode_index)

    def getEigvals(self, mode_indices=None):
        """Returns a sdarray of eigenvalues across modesets."""

        return self._getData('getEigvals', mode_indices)

    def getIndex(self, mode_index=0):
        """Returns indices of modes matched to the reference modeset."""

        return self._getModeData('getIndex', mode_index)

    def getIndices(self, mode_indices=None):
        """Returns indices of modes in the mode ensemble."""
        
        return self._getData('getIndices', mode_indices)

    def getAtoms(self):
        """Returns associated atoms of the mode ensemble."""

        return self._atoms

    def setAtoms(self, atoms):
        """Sets the atoms of the mode ensemble."""

        if atoms is not None and self._atoms is not None:
            if len(atoms) != self.numAtoms():
                raise ValueError('atoms should have %d atoms'%self.numAtoms())
        self._atoms = atoms

    def getLabels(self):
        """Returns the labels of the mode ensemble."""

        return self._labels

    def match(self):
        """Matches the modes across mode sets according the mode overlaps."""

        if self._modesets:
            #LOGGER.debug('Matching {0} modes across {1} modesets...'
            #                .format(self.numModes(), self.numModeSets()))
            start = time.time()
            self._modesets = matchModes(*self._modesets)
            LOGGER.debug('{0} modes across {1} modesets were matched in {2:.2f}s.'
                            .format(self.numModes(), self.numModeSets(), time.time()-start))
        else:
            LOGGER.warn('Mode ensemble has no modesets')
        self._matched = True
        return

    def reorder(self):
        """Reorders the modes across mode sets according to their collectivity"""
        if not self._matched:
            LOGGER.warn('Mode ensemble has not been matched')
        else:
            vals = self.getEigvals()
            
            mean_val = vals.mean()
            order = np.argsort(mean_val)

            ret = []
            for modeset in self._modesets:
                ret.append(ModeSet(modeset.getModel(), order))

            self._modesets = ret

    def addModeSet(self, modeset, weights=None, label=None):
        """Adds a modeset or modesets to the mode ensemble."""

        if isinstance(modeset, (NMA, ModeSet, Mode)):
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
            elif isinstance(modeset, Mode):
                modeset = ModeSet(modeset.getModel(), [modeset.getIndex()])
            if not isinstance(modeset, ModeSet):
                raise TypeError('modesets should be a list of Mode or ModeSet instances')
            if self._modesets:
                if modeset.numAtoms() != self.numAtoms():
                    raise ValueError('to be added, modesets should contain exactly %d atoms'%self.numAtoms())
                if modeset.numModes() < self.numModes():
                    raise ValueError('to be added, modesets should contain at least %d modes'%self.numModes())
                if modeset.numModes() > self.numModes():
                    modeset = modeset[:self.numModes()]
                self._modesets.append(modeset)
            else:
                self._modesets = [modeset]

        if weights is None:
            weights = np.ones((len(modesets), self.numAtoms(), 1))
            if self._weights is not None:
                self._weights = np.concatenate((self._weights, weights), axis=0)
            else:
                self._weights = weights
        else:
            weights = checkWeights(weights, self.numAtoms(), len(modesets))
            if self._weights is None:
                self._weights = np.ones((self.numModeSets()-len(modesets), self.numAtoms(), 1))
            self._weights = np.concatenate((self._weights, weights), axis=0)

    def delModeSet(self, index):
        """Removes a modeset or modesets from the mode ensemble."""

        if np.isscalar(index):
            index = [index]
        else:
            index = list(index)
        
        n_modesets = self.numModeSets()

        for i in reversed(index):
            self._modesets.pop(i)
            if self._labels:
                self._labels.pop(i)

        if self._weights is not None:
            torf = np.ones(n_modesets, dtype=bool)
            torf[index] = False
            self._weights = self._weights[torf, :, :]

    def isMatched(self):
        """Returns whether the modes are matched across modesets in the mode ensemble"""

        return self._matched

class sdarray(ndarray):
    """
    A class for representing a collection of arrays. It is derived from 
    :class:`~numpy.ndarray`, and the first axis is reserved for indexing 
    the collection. 
    
    :class:`sdarray` functions exactly the same as :class:`~numpy.ndarray`, 
    except that :meth:`sdarray.mean`, :meth:`sdarray.std`, 
    :meth:`sdarray.max`, :meth:`sdarray.min` are overriden. 
    Average, standard deviation, minimum, maximum, etc. are weighted and 
    calculated over the first axis by default. "sdarray" stands for 
    "signature dynamics array".

    Note for developers: please read the following article about subclassing 
    :class:`~numpy.ndarray` before modifying this class:

    https://docs.scipy.org/doc/numpy-1.14.0/user/basics.subclassing.html
    """

    __slots__ = ['_title', '_labels', '_is3d', '_weights', '_oneset']

    def __new__(self, array, weights=None, labels=None, title=None, is3d=False, 
                oneset=False):

        array = np.asarray(array)
        obj = array.view(self)
        shape = array.shape

        if title is None:
            if oneset:
                obj._title = '1 array of size {0}'.format(shape)
            else:
                obj._title = '{0} array{2} of size {1}'.format(
                             shape[0], shape[1:], 's' if shape[0]>1 else '')
        else:
            obj._title = title
        
        if labels is not None and np.isscalar(labels):
            labels = [labels]
        
        n_modesets = 1 if oneset else len(array)
        if len(labels) != n_modesets:
            raise ValueError('the number of labels does not match the size of the array')
        obj._labels = labels
        obj._is3d = is3d

        if weights is not None:
            weights = np.asarray(weights)
        obj._weights = weights
        obj._oneset = oneset
        return obj

    # def __getitem__(self, index):
    #     if isinstance(index, tuple):
    #         index0 = index[0]
    #         index1 = index[1:]
    #     else:
    #         index0 = index
    #         index1 = ()

    #     arr = np.asarray(self)[index0]
    #     w = self._weights
    #     if w is not None:
    #         w = w[index0]
    #     if arr.ndim != self.ndim:
    #         arr = np.expand_dims(arr, axis=0)
    #         if w is not None:
    #             w = np.expand_dims(w, axis=0)
    #     new_index = [slice(None, None, None)]
    #     new_index.extend(index1)

    #     arr = arr[new_index]
    #     if w is not None:
    #         w = w[new_index]
        
    #     labels = self._labels
    #     if labels is not None:
    #         labels = np.array(labels)[index0].tolist()
    #     return sdarray(arr, weights=w, labels=labels, title=self._title, is3d=self.is3d)

    def __getitem__(self, index):

        if self._oneset:
            index0 = 0
        else:
            if isinstance(index, tuple):
                index0 = index[0]
            else:
                index0 = index

        arr = np.asarray(self)[index]
        if np.ndim(arr) == 0:
            return arr

        oneset = np.isscalar(index0)

        w = self._weights
        if w is not None:
            w = w[index]
        
        labels = self._labels
        if labels is not None:
            labels = np.array(labels)[index0].tolist()
        return sdarray(arr, weights=w, labels=labels, title=self._title, 
                       is3d=self.is3d, oneset=oneset)

    def __array_finalize__(self, obj):
        if obj is None:
            return

        self._title = getattr(obj, '_title', None)
        self._weights = getattr(obj, '_weights', None)
        self._is3d = getattr(obj, '_is3d', None)
        self._labels = getattr(obj, '_labels', None)
        self._oneset = getattr(obj, '_oneset', None)

    def __str__(self):
        return self.getTitle()

    def is3d(self):
        """Returns **True** is model is 3-dimensional."""
        
        return self._is3d

    def __repr__(self):
        arr_repr = ndarray.__repr__(self)
        weights = self._weights
        if weights is not None:
            weights_repr = repr(weights)
            arr_repr += '\nweights=\n{0}'.format(weights_repr)
        return arr_repr

    def numAtoms(self):
        """Returns the number of atoms assuming it is represented by the second axis."""

        try:
            if self._oneset:
                n_atoms = self.shape[0]
            else:
                n_atoms = self.shape[1]
            if self.is3d():
                n_atoms /= 3
            return n_atoms
        except IndexError:
            LOGGER.warn('{0} is not related to the number of atoms'.format(self.getTitle()))
            return 0

    def numModeSets(self):
        """Returns the number of modesets in the instance """

        if self._oneset:
            return 1
        return self.shape[0]

    def getTitle(self):
        """Returns the title of the signature."""

        return self._title

    def getLabels(self):
        """Returns the labels of the signature."""

        return self._labels

    def mean(self, axis=0, **kwargs):
        """Calculates the weighted average of the sdarray over modesets (`axis=0`)."""

        arr = np.asarray(self)
        return np.average(arr, axis=axis, weights=self._weights)
    
    def std(self, axis=0, **kwargs):
        """Calculates the weighted standard deviations of the sdarray over modesets (`axis=0`)."""

        arr = np.asarray(self)
        mean = np.average(arr, weights=self._weights, axis=axis)
        variance = np.average((arr - mean)**2, 
                              weights=self._weights, axis=axis)
        return np.sqrt(variance)

    def min(self, axis=0, **kwargs):
        """Calculates the minimum values of the sdarray over modesets (`axis=0`)."""

        arr = np.asarray(self)
        return arr.min(axis=axis)

    def max(self, axis=0, **kwargs):
        """Calculates the maximum values of the sdarray over modesets (`axis=0`)."""

        arr = np.asarray(self)
        return arr.max(axis=axis)

    def getWeights(self):
        """Returns the weights of the signature."""

        return self._weights

    def setWeights(self, weights):
        """Sets the weights of the signature."""

        self._weights = weights

    weights = property(getWeights, setWeights)

    def getArray(self):
        """Returns the signature as an numpy array."""

        return np.asarray(self)

def calcEnsembleENMs(ensemble, model='gnm', trim='reduce', n_modes=20, **kwargs):
    """Description"""

    match = kwargs.pop('match', True)
    if isinstance(ensemble, Conformation):
        conformation = ensemble
        ensemble = conformation.getEnsemble()
        index = conformation.getIndex()
        ensemble = ensemble[index:index+1]
    if model is GNM:
        model_type = 'GNM'
    elif model is ANM:
        model_type = 'ANM'
    else:
        model_type = str(model).strip().upper()

    start = time.time()

    atoms = ensemble.getAtoms()
    select = None
    if ensemble.isSelected():
        select = atoms
        atoms = ensemble.getAtoms(selected=False)
        
    labels = ensemble.getLabels()

    ### ENMs ###
    ## ENM for every conf
    enms = []
    n_confs = ensemble.numConfs()

    str_modes = 'all' if n_modes is None else str(n_modes)
    LOGGER.progress('Calculating {0} {1} modes for {2} conformations...'
                    .format(str_modes, model_type, n_confs), n_confs, '_prody_calcEnsembleENMs')

    for i in range(n_confs):
        LOGGER.update(i, label='_prody_calcEnsembleENMs')
        coords = ensemble.getCoordsets(i, selected=False)
        nodes = coords[0, :, :]
        if atoms is not None:
            atoms.setCoords(nodes)
            nodes = atoms
        enm, _ = calcENM(nodes, select, model=model, trim=trim, 
                            n_modes=n_modes, title=labels[i], **kwargs)
        enms.append(enm)

        #lbl = labels[i] if labels[i] != '' else '%d-th conformation'%(i+1)
    LOGGER.finish()

    min_n_modes = ensemble.numAtoms() * 3
    for enm in enms:
        n_modes = enm.numModes()
        if n_modes < min_n_modes:
            min_n_modes = n_modes

    for i in range(len(enms)):
        n_modes = enms[i].numModes()
        if n_modes > min_n_modes:
            enms[i] = enms[i][:min_n_modes]
            LOGGER.warn('last {0} modes for {1} has been discarded because at least one '
                        'conformation has only {2} modes'.format(n_modes-min_n_modes, 
                        enms[i].getTitle(), min_n_modes))

    LOGGER.info('{0} {1} modes were calculated for each of the {2} conformations in {3:.2f}s.'
                        .format(str_modes, model_type, n_confs, time.time()-start))

    modeens = ModeEnsemble(title=ensemble.getTitle())
    modeens.addModeSet(enms, weights=ensemble.getWeights(), 
                             label=ensemble.getLabels())
    modeens.setAtoms(ensemble.getAtoms())

    if match:
        modeens.match()
    return modeens

def _getEnsembleENMs(ensemble, **kwargs):
    if isinstance(ensemble, (Ensemble, Conformation)):
        enms = calcEnsembleENMs(ensemble, **kwargs)
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
    """Calculate the spectral overlaps between each pair of conformations in the 
    *ensemble*.
    
    :arg ensemble: an ensemble of structures or ENMs 
    :type ensemble: :class: `Ensemble`, :class: `ModeEnsemble`

    :arg distance: if set to **True**, spectral overlap will be converted to spectral 
                   distance via arccos.
    :type distance: bool
    """

    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    overlaps = np.zeros((len(enms), len(enms)))
    for i, enmi in enumerate(enms):
        for j, enmj in enumerate(enms):
            covlap = calcSpectralOverlap(enmi, enmj)
            overlaps[i, j] = covlap

    if distance:
        overlaps = np.arccos(overlaps)

    return overlaps

def calcSignatureSqFlucts(mode_ensemble, **kwargs):
    """
    Get the signature square fluctuations of *mode_ensemble*. 
    
    :arg mode_ensemble: an ensemble of ENMs 
    :type mode_ensemble: :class: `ModeEnsemble`
    """

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')
    
    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    modesets = mode_ensemble
    V = []
    for modes in modesets:
        sqfs = calcSqFlucts(modes)
        V.append(sqfs)
    V = np.vstack(V)

    title_str = '%d modes'%mode_ensemble.numModes()
    weights = mode_ensemble.getWeights()
    if weights is not None:
        weights = weights[:, :, 0]
    labels = mode_ensemble.getLabels()

    # even the original model is 3d, sqfs are still 1d
    sig = sdarray(V, title=title_str, weights=weights, labels=labels, is3d=False)

    return sig

def showSignatureAtomicLines(y, std=None, min=None, max=None, atoms=None, **kwargs):
    """
    Show the signature dynamics data using :func:`showAtomicLines`. 
    
    :arg y: the mean values of signature dynamics to be plotted 
    :type y: :class:`~numpy.ndarray`

    :arg std: the standard deviations of signature dynamics to be plotted 
    :type std: :class:`~numpy.ndarray`

    :arg min: the minimum values of signature dynamics to be plotted 
    :type min: :class:`~numpy.ndarray`

    :arg max: the maximum values of signature dynamics to be plotted 
    :type max: :class:`~numpy.ndarray`

    :arg linespec: line specifications that will be passed to :func:`showAtomicLines`
    :type linespec: str

    :arg atoms: an object with method :meth:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 
    """

    from matplotlib.pyplot import figure, plot, fill_between, \
                                  gca, xlabel, ylabel, title, ylim

    linespec = kwargs.pop('linespec', '-')
    zero_line = kwargs.pop('zero_line', False)

    x = range(y.shape[0])
    lines, polys, bars, texts = showAtomicLines(y, atoms=atoms, dy=std, lower=max, upper=min, 
                                        linespec=linespec, show_zero=zero_line, **kwargs)
        
    return lines, polys, bars, texts

def showSignature1D(signature, linespec='-', **kwargs):
    """
    Show *signature* using :func:`showAtomicLines`. 
    
    :arg signature: the signature dynamics to be plotted 
    :type signature: :class:`sdarray`

    :arg linespec: line specifications that will be passed to :func:`showAtomicLines`
    :type linespec: str

    :arg atoms: an object with method :func:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 

    :arg alpha: the transparency of the band(s).
    :type alpha: float

    :arg range: whether shows the minimum and maximum values. 
                Default is **True**
    :type range: bool
    """

    from matplotlib.pyplot import figure, plot, fill_between, \
                                  gca, xlabel, ylabel, title, ylim

    V = signature
        
    meanV, stdV, minV, maxV = V.mean(), V.std(), V.min(), V.max()

    atoms = kwargs.pop('atoms', None)
    zero_line = kwargs.pop('show_zero', False)
    zero_line = kwargs.pop('zero', zero_line)
    show_range = kwargs.pop('range', True)

    bars = []; polys = []; lines = []; texts = []

    if V.is3d():
        meanV = np.reshape(meanV, (V.numAtoms(), 3)).T
        stdV = np.reshape(stdV, (V.numAtoms(), 3)).T
        minV = np.reshape(minV, (V.numAtoms(), 3)).T
        maxV = np.reshape(maxV, (V.numAtoms(), 3)).T

        atoms_ = None; zero_line_ = False
        for i in range(3):
            if i == 2:
                atoms_ = atoms
                zero_line_ = zero_line
            if not show_range:
                minV[i] = maxV[i] = None
            _lines, _polys, _bars, _texts = showSignatureAtomicLines(meanV[i], stdV[i], minV[i], maxV[i], 
                                                   atoms=atoms_, zero_line=zero_line_,
                                                   linespec=linespec, **kwargs)
            lines.extend(_lines)
            bars.extend(_bars)
            polys.extend(_polys)
            texts.extend(_texts)

    else:
        if not show_range:
            minV = maxV = None
        _lines, _polys, _bars, _texts = showSignatureAtomicLines(meanV, stdV, minV, maxV, 
                                               atoms=atoms, zero_line=zero_line,
                                               linespec=linespec, **kwargs)
        lines.extend(_lines)
        bars.extend(_bars)
        polys.extend(_polys)
        texts.extend(_texts)

    xlabel('Residues')
    title('Signature profile of ' + V.getTitle())

    return lines, polys, bars, texts

def showSignatureMode(mode_ensemble, **kwargs):

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    scale = kwargs.pop('scale', 1.0)
    mode = mode_ensemble.getEigvec() * scale
    show_zero = kwargs.pop('show_zero', True)
    return showSignature1D(mode, atoms=mode_ensemble.getAtoms(), show_zero=show_zero, **kwargs)

def showSignatureSqFlucts(mode_ensemble, **kwargs):

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    scale = kwargs.pop('scale', 1.0)
    sqf = calcSignatureSqFlucts(mode_ensemble) * scale
    show_zero = kwargs.pop('show_zero', False)
    return showSignature1D(sqf, atoms=mode_ensemble.getAtoms(), show_zero=show_zero, **kwargs)

def calcSignatureCrossCorr(mode_ensemble, norm=True):
    """Calculate average cross-correlations for a ModeEnsemble."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')
    matches = mode_ensemble
    n_atoms = matches.numAtoms()
    n_sets = len(matches)

    C = np.zeros((n_sets, n_atoms, n_atoms))
    for i in range(n_sets):
        m = matches[i]
        c = calcCrossCorr(m, norm=norm)
        C[i, :, :] = c

    title_str = '%d modes'%mode_ensemble.numModes()
    weights = mode_ensemble.getWeights()
    if weights is not None:
        W = np.zeros((mode_ensemble.numModeSets(), 
                      mode_ensemble.numAtoms(), 
                      mode_ensemble.numAtoms()))
        for i, w in enumerate(weights):
            w2 = np.outer(w, w)
            W[i, :, :] = w2
    labels = mode_ensemble.getLabels()

    # even the original model is 3d, cross-correlations are still 1d
    sig = sdarray(C, title=title_str, weights=W, labels=labels, is3d=False)
        
    return sig

def calcSignatureCollectivity(mode_ensemble, masses=None):
    """Calculate average collectivities for a ModeEnsemble."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')
    
    n_modes = mode_ensemble.numModes()
    n_sets = len(mode_ensemble)

    C = np.zeros((n_sets, n_modes))
    for i in range(n_sets):
        m = mode_ensemble[i]
        c = calcCollectivity(m, masses=masses)
        C[i, :] = c

    title_str = 'collectivities of %d modes'%mode_ensemble.numModes()
    labels = mode_ensemble.getLabels()

    # even the original model is 3d, cross-correlations are still 1d
    sig = sdarray(C, title=title_str, weights=None, labels=labels, is3d=False)
        
    return sig

def calcSignatureOverlaps(mode_ensemble, diag=True):
    """Calculate average mode-mode overlaps for a ModeEnsemble."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    n_sets = mode_ensemble.numModeSets()
    n_modes = mode_ensemble.numModes()

    if diag:
        overlaps = np.zeros((n_modes, n_sets, n_sets))
        for m in range(mode_ensemble.numModes()):
            for i, modeset_i in enumerate(mode_ensemble):
                mode_i = modeset_i[m]
                for j, modeset_j in enumerate(mode_ensemble):
                    mode_j = modeset_j[m]
                    if j >= i:
                        overlaps[m, i, j] = overlaps[m, j, i] = abs(calcOverlap(mode_i, mode_j))
                    
    else:
        overlaps = np.zeros((n_modes, n_modes, n_sets, n_sets))

        for i, modeset_i in enumerate(mode_ensemble):
            for j, modeset_j in enumerate(mode_ensemble):
                if j >= i:                
                    overlaps[:,:,i,j] = overlaps[:,:,j,i] = abs(calcOverlap(modeset_i, modeset_j))

    return overlaps

def showSignatureOverlaps(mode_ensemble):

    from matplotlib.pyplot import xlabel, ylabel

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    overlaps = calcSignatureOverlaps(mode_ensemble, diag=True)
    r, c = np.triu_indices(overlaps.shape[1], k=1)
    overlap_triu = overlaps[:, r, c]

    meanV = overlap_triu.mean(axis=1)
    stdV = overlap_triu.std(axis=1)

    show = showSignatureAtomicLines(meanV, stdV)
    xlabel('Mode index')
    ylabel('Overlap')
    
    return show

def calcSignatureFractVariance(mode_ensemble):
    """Calculate signature fractional variance for a ModeEnsemble."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    matches = mode_ensemble
    n_sets = len(matches)

    W = []; is3d = None
    for i in range(n_sets):
        m = matches[i]
        var = calcFractVariance(m)
        W.append(var)
        if is3d is None:
            is3d = m.is3d()

    title_str = '%d modes'%mode_ensemble.numModes()
    labels = mode_ensemble.getLabels()
    sig = sdarray(W, title=title_str, weights=None, labels=labels, is3d=is3d)
        
    return sig

def showSignatureCrossCorr(mode_ensemble, std=False, **kwargs):
    """Show average cross-correlations using :func:`showAtomicMatrix`. 
    By default, *origin=lower* and *interpolation=bilinear* keyword arguments
    are passed to this function, but user can overwrite these parameters.
    See also :func:`.calcSignatureCrossCorr`.
    
    :arg ensemble: an ensemble of structures or ENMs, or a signature profile 
    :type ensemble: :class: `Ensemble`, list, :class:`sdarray`

    :arg atoms: an object with method :func:`getResnums` for use 
                on the x-axis.
    :type atoms: :class:`Atomic` 
    """

    import matplotlib.pyplot as plt
    
    norm = kwargs.pop('norm', True)
    C = calcSignatureCrossCorr(mode_ensemble, norm=norm)

    atoms = kwargs.pop('atoms', None)
    if atoms is None:
        try:
            atoms = mode_ensemble.getAtoms()
        except:
            pass

    if std:
        matrixData = C.std()
    else:
        matrixData = C.mean()
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'bilinear'

    show = showAtomicMatrix(matrixData, atoms=atoms, **kwargs)

    indices = mode_ensemble.getIndices()[0]
    if len(indices) == 1:
        title_str = ', mode '+str(indices[0]+1)
    else:
        modeIndexStr = ','.join([str(x+1) for x in indices])
        if len(modeIndexStr) > 8:
            title_str = ', '+str(len(indices))+' modes '+modeIndexStr[:5]+'...'
        else:
            title_str = ', modes '+modeIndexStr
        # title_str = ', '+str(len(modeIndex))+' modes'
    if std:
        plt.title('Cross-correlations (standard deviation)'+title_str)
    else:
        plt.title('Cross-correlations (average)'+title_str)
    plt.xlabel('Residues')
    plt.ylabel('Residues')
    
    return show

def showSignatureDistribution(signature, **kwargs):
    """
    Show the distribution of signature values using 
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
        
    W = signature.getArray()[:, ::-1] # reversed to accommodate with matplotlib.pyplot.hist
    weights = np.ones_like(W)/float(len(W))
    
    show_legend = kwargs.pop('legend', True)
    bins = kwargs.pop('bins', 'auto')
    if bins == 'auto':
        _, bins = np.histogram(W.flatten(), bins='auto')
    elif np.isscalar(bins) and isinstance(bins, (int, np.integer)):
        step = (W.max() - W.min())/bins
        bins = np.arange(W.min(), W.max(), step)

    histtype = kwargs.pop('histtype', 'stepfilled')
    label = kwargs.pop('label', None)
    labels = kwargs.pop('labels', label)
    weights = kwargs.pop('weights', weights)
    n, bins, patches = hist(W, bins=bins, weights=weights, 
                            histtype=histtype, label=labels, **kwargs)

    colors = []
    for patch_i in patches:
        colors.append(patch_i[0].get_facecolor())

    for i, patch_i in enumerate(patches):
        for j, patch_j in enumerate(patch_i):
            patch_j.set_color(colors[-(i+1)])

    if show_legend:
        legend()

    xlabel('Signature value')
    ylabel('Probability')

    if SETTINGS['auto_show']:
        showFigure()

    return n, bins, patches

def showSignatureVariances(mode_ensemble, **kwargs):
    """
    Show the distribution of signature variances using 
    :func:`showSignatureDistribution`.
    """

    from matplotlib.pyplot import xlabel

    fract = kwargs.pop('fraction', True)

    if fract:
        sig = calcSignatureFractVariance(mode_ensemble)
    else:
        sig = mode_ensemble.getVariances() 

    indices = np.asarray(mode_ensemble.getIndices())[0]
    labels = ['mode %d'%(i+1) for i in indices][::-1]

    show = showSignatureDistribution(sig, label=labels, **kwargs)

    xlabel('Variance')
    return show

def showSignatureCollectivity(mode_ensemble, **kwargs):
    """
    Show the distribution of signature variances using 
    :func:`showSignatureDistribution`.
    """

    from matplotlib.pyplot import xlabel

    sig = calcSignatureCollectivity(mode_ensemble)

    indices = np.asarray(mode_ensemble.getIndices())[0]
    labels = ['mode %d'%(i+1) for i in indices][::-1]

    show = showSignatureDistribution(sig, label=labels, **kwargs)

    xlabel('Collectivity')
    return show

def showVarianceBar(mode_ensemble, highlights=None, **kwargs):

    from matplotlib.pyplot import figure, gca, annotate, subplots_adjust, plot
    from matplotlib.figure import Figure
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize, NoNorm
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
        if fig_num is None:
            figure(figsize=(6, 2))
        else:
            figure(fig_num)
    elif fig_num is not None:
        figure(fig_num)
    ax = gca()

    # adjust layouts
    box = ax.get_position()
    _, _, _, height = box.bounds
    ratio = 2.5
    box.y1 = box.y0 + height/ratio
    #box.y0 += height/7.
    ax.set_position(box)

    fract = kwargs.pop('fraction', True)

    #defarrow = {'width':1, 'headwidth':2, 
    #            'facecolor':'black',
    #            'headlength': 4}
    defarrow = {'arrowstyle': '->'}
    arrowprops = kwargs.pop('arrowprops', defarrow)

    if fract:
        sig = calcSignatureFractVariance(mode_ensemble)
    else:
        sig = mode_ensemble.getVariances() 

    variances = sig.getArray().sum(axis=1)
    #meanVar = variances.mean()
    #stdVar = variances.std()
    
    #variances = (variances - meanVar)/stdVar

    maxVar = variances.max()
    minVar = variances.min()

    cmap = kwargs.pop('cmap', 'jet')
    norm = Normalize(vmin=minVar, vmax=maxVar)
    cb = ColorbarBase(ax, cmap=cmap, norm=norm,
                      orientation='horizontal')

    if not highlights:
        highlights = []

    indices = []; labels = []
    ens_labels = mode_ensemble.getLabels()
    for hl in highlights:
        if isinstance(hl, str):
            if not ens_labels:
                raise TypeError('highlights should be a list of integers because '
                                    'mode_ensemble has no label')
            indices.append(ens_labels.index(hl))
            labels.append(hl)
        else:
            try:
                index = int(hl)
            except:
                raise TypeError('highlights should be a list of integers or strings') 
            indices.append(index)
            if ens_labels:
                labels.append(ens_labels[index])
            else:
                labels.append(str(index))

    annotations = []
    for i, label in zip(indices, labels):
        x = norm(variances[i])
        an = annotate(label, xy=(x, 1), xytext=(x, ratio), arrowprops=arrowprops)
        annotations.append(an)

    for i in range(len(variances)):
        x = norm(variances[i])
        plot([x, x], [0, 1], 'w')

    cb.set_label('Variances')

    if SETTINGS['auto_show']:
        showFigure()
    return cb, annotations

def saveModeEnsemble(mode_ensemble, filename=None, atoms=False, **kwargs):
    """Save *mode_ensemble* as :file:`filename.modeens.npz`.  If *filename* 
    is **None**, title of the ModeEnsemble instance will be used as the 
    filename, after ``" "`` (white spaces) in the title are replaced with 
    ``"_"`` (underscores).  Upon successful completion of saving, filename 
    is returned. This function makes use of :func:`~numpy.savez_compressed` 
    function."""

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('invalid type for mode_ensemble, {0}'
                        .format(type(mode_ensemble)))
    if len(mode_ensemble) == 0:
        raise ValueError('mode_ensemble instance does not contain data')

    attr_list = ['_modesets', '_title', '_labels', '_weights', '_matched']
    attr_dict = {}

    if atoms:
        attr_list.append('_atoms')
    
    for attr in attr_list:
        value = getattr(mode_ensemble, attr)
        if value is not None:
            if attr == '_atoms':
                value = [value, None]
            if attr == '_modesets':
                value = list(value)
                value.append(None)
            attr_dict[attr] = value

    if filename is None:
        filename = mode_ensemble.getTitle().replace(' ', '_')
    
    suffix = '.modeens'
    if not filename.lower().endswith('.npz'):
        if not filename.lower().endswith(suffix):
            filename += suffix + '.npz'
        else:
            filename += '.npz'
            
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez_compressed(ostream, **attr_dict)
    ostream.close()

    return filename

def loadModeEnsemble(filename, **kwargs):
    """Returns ModeEnsemble instance after loading it from file (*filename*).
    This function makes use of :func:`numpy.load` function.  See
    also :func:`saveModeEnsemble`."""

    if not 'encoding' in kwargs:
        kwargs['encoding'] = 'latin1'
    data = np.load(filename, **kwargs)
    
    weights = getValue(data, '_weights', None)
    labels = getValue(data, '_labels', None)
    matched = getValue(data, '_matched', False)
    title = getValue(data, '_title', None)
    modesets = getValue(data, '_modesets', [])
    atoms = getValue(data, '_atoms', [None])[0]

    if isinstance(title, np.ndarray):
        title = np.asarray(title, dtype=str)
    title = str(title)

    if isinstance(modesets, np.ndarray):
        modesets = modesets.tolist()
    while (None in modesets):
        modesets.remove(None)

    if labels is not None:
        labels = labels.tolist()

    modeens = ModeEnsemble(title=title)
    modeens._weights = weights
    modeens._labels = labels
    modeens._matched = matched
    modeens._modesets = modesets
    modeens._atoms = atoms

    return modeens
