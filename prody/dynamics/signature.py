# -*- coding: utf-8 -*-
"""This module defines functions for analyzing normal modes obtained 
for conformations in an ensemble."""

import time
from numbers import Integral
from numpy import ndarray
import numpy as np
import warnings

from prody import LOGGER, SETTINGS
from prody.utilities import showFigure, showMatrix, copy, checkWeights, openFile, DTYPE
from prody.utilities import getValue, importLA, wmean, div0, isListLike
from prody.ensemble import Ensemble, Conformation
from prody.atomic import AtomGroup

from .nma import NMA
from .modeset import ModeSet
from .mode import Mode, Vector
from .functions import calcENM
from .compare import calcSpectralOverlap, matchModes, calcOverlap

from .analysis import calcSqFlucts, calcCrossCorr, calcFractVariance, calcCollectivity
from .plotting import showAtomicLines, showAtomicMatrix, showDomainBar
from .perturb import calcPerturbResponse
from .anm import ANM
from .gnm import GNM

__all__ = ['ModeEnsemble', 'sdarray', 'calcEnsembleENMs', 
           'showSignature1D', 'psplot', 'showSignatureAtomicLines', 
           'showSignatureMode', 'showSignatureDistribution', 'showSignatureCollectivity',
           'showSignatureSqFlucts', 'calcEnsembleSpectralOverlaps', 'calcSignatureSqFlucts', 
           'calcSignatureCollectivity', 'calcSignatureFractVariance', 'calcSignatureModes', 
           'calcSignatureCrossCorr', 'showSignatureCrossCorr', 'showVarianceBar',
           'showSignatureVariances', 'calcSignatureOverlaps', 'showSignatureOverlaps',
           'saveModeEnsemble', 'loadModeEnsemble', 'saveSignature', 'loadSignature',
           'calcSubfamilySpectralOverlaps','showSubfamilySpectralOverlaps', 'calcSignaturePerturbResponse']

class ModeEnsemble(object):
    """
    A collection of ENMs calculated for conformations in an :class:`Ensemble`. 
    or :class:`PDBEnsemble`. 
    """

    __slots__ = ['_modesets', '_title', '_labels', '_atoms', '_weights', '_matched', '_reweighted']

    def __init__(self, title=None):
        self._modesets = []
        self._title = 'Unknown' if title is None else title
        self._labels = None
        self._atoms = None
        self._weights = None
        self._matched = False
        self._reweighted = False

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
        if self.numModeSets() == 1:
            modesets_str = '1 modeset'
        else:
            modesets_str = '{0} modesets'.format(self.numModeSets())

        if self.numModes() == 1:
            modes_str = '1 mode'
        else:
            modes_str = '{0} modes'.format(self.numModes())

        if self.numAtoms() == 1:
            atoms_str = '1 atom'
        else:
            atoms_str = '{0} atoms'.format(self.numAtoms())

        return '<ModeEnsemble: {0} ({1}, {2})>'\
                .format(modesets_str, modes_str, atoms_str)

    def __str__(self):
        return self.getTitle()
    
    def _getitembase(self, modeset_index, mode_index):

        if isinstance(modeset_index, slice):
            modesets = self._modesets[modeset_index]
            modeset_indices = modeset_index
            labels = None if self._labels is None else self._labels[modeset_index]
            matched = self.getMatchingStatus()[modeset_index]
            reweighted = self.getReweightingStatus()[modeset_index]
        elif not np.isscalar(modeset_index):
            modesets = []; modeset_indices = []; labels = []; matched = []; reweighted = []
            for i in modeset_index:
                if isinstance(i, Integral):
                    j = i
                elif isinstance(i, str):
                    try:
                        j = self._labels.index(i)
                    except:
                        raise IndexError('invalid label: %s'%i)
                else:
                    raise IndexError('all indices must be integers or strings (labels)')
                modesets.append(self._modesets[j])
                modeset_indices.append(j)
                if self._labels is not None:
                    labels.append(self._labels[j])
                
                allstatus = self.getMatchingStatus()
                matched.append(allstatus[j])

                allstatus = self.getReweightingStatus()
                reweighted.append(allstatus[j])
        else:
            if isinstance(modeset_index, Integral):
                pass
            elif isinstance(modeset_index, str):
                try:
                    modeset_index = self._labels.index(modeset_index)
                except:
                    raise IndexError('invalid label: %s'%modeset_index)
            else:
                raise IndexError('indices must be int, slice, or array-like objects')
            return self._modesets[modeset_index][mode_index]
        
        if np.isscalar(mode_index):
            mode_index = [mode_index]
        for i in range(len(modesets)):
            modesets[i] = modesets[i][mode_index]

        if self._weights is None:
            weights = None
        else:
            weights = self._weights[modeset_indices, :, :]

        ens = ModeEnsemble(title=self.getTitle())
        ens.addModeSet(modesets, weights=weights, label=labels, matched=matched, reweighted=reweighted)
        ens.setAtoms(self.getAtoms())
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
        ensemble.addModeSet(self._modesets, weights=self._weights, label=self._labels,
                            matched=self.getMatchingStatus(), reweighted=self.getReweightingStatus())
        ensemble.setAtoms(self.getAtoms())
        
        ensemble.addModeSet(other.getModeSets(), weights=other.getWeights(),
                            label=other.getLabels(), matched=other.getMatchingStatus(),
                            reweighted=other.getReweightingStatus())
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

        subset = self[index]
        if isinstance(subset, ModeSet):
            return subset
        else:
            return subset._modesets

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

    def setLabels(self, labels):
        """Returns the labels of the mode ensemble."""
        if len(labels) != self.numModeSets():
            raise ValueError('the number of labels and mode sets mismatch')
        self._labels = labels

    def getMatchingStatus(self):
        """Returns the matching status of each mode ensemble."""

        if isinstance(self._matched, bool):
            matched = [self._matched for _ in self]
        else:
            matched = self._matched
        return matched

    def setMatchingStatus(self, status):
        """Returns the matching status of each mode ensemble."""

        if not isinstance(status, bool):
            if len(status) != self.numModeSets():
                raise ValueError('the number of status and mode sets mismatch')
        self._matched = status

    def getReweightingStatus(self):
        """Returns the reweighting status of each mode ensemble."""

        if isinstance(self._reweighted, bool):
            reweighted = [self._reweighted for _ in self]
        else:
            reweighted = self._reweighted
        return reweighted

    def setReweightingStatus(self, status):
        """Returns the reweighting status of each mode ensemble."""

        if not isinstance(status, bool):
            if len(status) != self.numModeSets():
                raise ValueError('the number of status and mode sets mismatch')
        self._reweighted = status

    def match(self, turbo=False, method=None):
        """Matches the modes across mode sets according the mode overlaps.

        :arg turbo: if **True** then the computation will be performed in parallel. 
                The number of threads is set to be the same as the number of 
                CPUs. Assigning a number will specify the number of threads to be 
                used. Default is **False**
        :type turbo: bool, int
        """

        if self._modesets:
            #LOGGER.debug('Matching {0} modes across {1} modesets...'
            #                .format(self.numModes(), self.numModeSets()))
            start = time.time()

            matched = self.getMatchingStatus()
            if np.any(matched):
                matched[0] = False
                indices = [i for i in range(len(matched)) if not matched[i]]
                modesets = [self._modesets[i] for i in indices]

                modesets = matchModes(*modesets, turbo=turbo, method=method)

                for n, i in enumerate(indices):
                    if n > 0:
                        self._modesets[i] = modesets[n]

                n_modesets = len(modesets)
            else: # if all not matched, start from scratch
                self._modesets = matchModes(*self._modesets, turbo=turbo, method=method)
                n_modesets = len(self._modesets)

            LOGGER.debug('{0} modes across {1} modesets were matched in {2:.2f}s.'
                            .format(self.numModes(), n_modesets, time.time()-start))
        else:
            LOGGER.warn('Mode ensemble has no modesets')
        self._matched = True
        return

    def reweight(self):
        """Reweight the modes based on matched orders"""

        self._reweighted = []
        for i in range(self.numModeSets()):
            modes = self._modesets[i]
            model = modes._model

            I = modes.getIndices()   # matched order for the subset of modes
            J = np.sort(I) # original order for the subset
            V = model.getEigvecs()

            W = model.getEigvals()
            W[I] = W[J]
            model.setEigens(V, W)
            self._reweighted.append(True)

    def undoMatching(self):
        """Restores the original orders of modes"""

        n_modes = self.numModes()
        matched = self.getMatchingStatus()
        for i, enm in enumerate(self):
            if matched[i]:
                model = enm._model
                self._modesets[i] = model[:n_modes]

        self.setMatchingStatus(False)

    def undoReweighting(self):
        """Restores the original weighting of modes"""

        reweighted = self.getReweightingStatus()
        for i, enm in enumerate(self):
            model = enm._model
            if reweighted[i]:
                vars = model.getVariances()
                I = np.argsort(vars)[::-1]

                V = model.getEigvecs()

                W = model.getEigvals()
                W = W[I]
                model.setEigens(V, W)

        self.setReweightingStatus(False)

    def reorder(self):
        """Reorders the modes across mode sets according to their collectivity"""
        if not self.isMatched():
            LOGGER.warn('Mode ensemble has not been all matched')
        else:
            vals = self.getEigvals()
            
            mean_val = vals.mean()
            order = np.argsort(mean_val)

            ret = []
            for modeset in self._modesets:
                ret.append(ModeSet(modeset.getModel(), order))

            self._modesets = ret

    def addModeSet(self, modeset, weights=None, label=None, matched=False, reweighted=False):
        """Adds a modeset or modesets to the mode ensemble."""

        if isinstance(modeset, (NMA, ModeSet, Mode)):
            modesets = [modeset]
        else:
            modesets = modeset

        matchingstatus = self.getMatchingStatus()
        if isinstance(matched, bool):
            matched = [matched for _ in range(len(modesets))]
        matchingstatus.extend(matched)
        self._matched = matchingstatus

        reweightingstatus = self.getReweightingStatus()
        if isinstance(reweighted, bool):
            reweighted = [reweighted for _ in range(len(modesets))]
        reweightingstatus.extend(reweighted)
        self._reweighted = reweightingstatus

        if label is None:
            labels = ['']*len(modesets)
        elif np.isscalar(label):
            labels = [label]
        else:
            labels = label
            
        if len(labels) != len(modesets):
            raise ValueError('labels should have the same length as modesets')

        if self._labels is None:
            self._labels = ['']*len(self._modesets)

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

            if isinstance(self._matched, list):
                self._matched.pop(i)

            if isinstance(self._reweighted, list):
                self._reweighted.pop(i)

        if self._weights is not None:
            torf = np.ones(n_modesets, dtype=bool)
            torf[index] = False
            self._weights = self._weights[torf, :, :]

    def isMatched(self):
        """Returns whether the modes are matched across ALL modesets in the mode ensemble"""

        return np.all(self._matched)

    def isReweighted(self):
        """Returns whether the modes are matched across ALL modesets in the mode ensemble"""

        return np.all(self._reweighted)

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
            return int(n_atoms)
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
        weights = self._weights
        return wmean(arr, weights, axis)
    
    def std(self, axis=0, **kwargs):
        """Calculates the weighted standard deviations of the sdarray over modesets (`axis=0`)."""

        arr = np.asarray(self)
        mean = wmean(arr, self._weights, axis)

        if axis is not None:
            reps = list(arr.shape)
            axes = []
            if isinstance(axis, tuple):
                axes.extend(axis)
            else:
                axes = [axis]
            
            for ax in axes:
                reps[ax] = 1
            mean = np.reshape(mean, reps)
        variance = wmean((arr - mean)**2, self._weights, axis)
        return np.sqrt(variance)

    def min(self, axis=0, **kwargs):
        """Calculates the minimum values of the sdarray over modesets (`axis=0`)."""

        # np.array instead asarray is used to make sure a copy of the original data is created
        arr = np.array(self)
        weights = self._weights
        if weights is not None:
            weights = weights.astype(bool)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rmin = np.nanmin(arr, axis=axis)

        return rmin

    def max(self, axis=0, **kwargs):
        """Calculates the maximum values of the sdarray over modesets (`axis=0`)."""

        # np.array instead asarray is used to make sure a copy of the original data is created
        arr = np.array(self)
        weights = self._weights
        if weights is not None:
            weights = weights.astype(bool)
            arr[~weights] = np.nan

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rmax = np.nanmax(arr, axis=axis)

        return rmax

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

    def transpose(self, axes=None):
        a = np.asarray(self)
        return np.transpose(a, axes=axes)


def calcEnsembleENMs(ensemble, model='gnm', trim='reduce', n_modes=20, **kwargs):
    """Calculates normal modes for each member of *ensemble*.
    
    :arg ensemble: normal modes of whose members to be computed
    :type ensemble: :class:`.PDBEnsemble`
    
    :arg model: type of ENM that will be performed. It can be either 'anm' 
                or 'gnm'
    :type model: str

    :arg trim: type of method that will be used to trim the model. It can 
               be either 'trim' , 'slice', or 'reduce'. If set to 'trim', the parts 
               that is not in the selection will simply be removed
    :type trim: str

    :arg n_modes: number of modes to be computed
    :type trim: int

    :arg turbo: if **True** then the computation will be performed in parallel. 
                The number of threads is set to be the same as the number of 
                CPUs. Assigning a number to specify the number of threads to be 
                used. Default is **False**
    :type turbo: bool, int

    :arg match: whether the modes should be matched using :func:`.matchModes`. 
                Default is **True**
    :type match: bool

    :arg method: the alternative function that is used to match the modes. 
                Default is **None**
    :type method: function

    :arg turbo: whether use :class:`~multiprocessing.Pool` to accelerate the computation. 
                Note that if writing a script, ``if __name__ == '__main__'`` is necessary 
                to protect your code when multi-tasking. 
                See https://docs.python.org/2/library/multiprocessing.html for details.
                Default is **False**
    :type turbo: bool

    :returns: :class:`.ModeEnsemble`
    """

    match = kwargs.pop('match', True)
    method = kwargs.pop('method', None)
    turbo = kwargs.pop('turbo', False)

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
    select = ensemble.getIndices()
    labels = ensemble.getLabels()
    n_atoms = ensemble.numAtoms(selected=False)

    if select is None:
        torf_selected = np.ones(n_atoms, dtype=bool)
    else:
        torf_selected = np.zeros(n_atoms, dtype=bool)
        torf_selected[select] = True

    ### ENMs ###
    ## ENM for every conf
    enms = []
    n_confs = ensemble.numConfs()

    str_modes = 'all' if n_modes is None else str(n_modes)
    LOGGER.progress('Calculating {0} {1} modes for {2} conformations...'
                    .format(str_modes, model_type, n_confs), n_confs, '_prody_calcEnsembleENMs')

    coordsets = ensemble.getCoordsets(selected=False)
    weights = ensemble.getWeights(selected=False)
    for i in range(n_confs):
        LOGGER.update(i, label='_prody_calcEnsembleENMs')
        coords = coordsets[i]

        if weights.ndim == 3:
            weight = weights[i].flatten()
        else:
            weight = weights.flatten()
        torf_mapped = weight != 0

        coords = coords[torf_mapped, :]
        system = torf_selected[torf_mapped]
        mask = torf_mapped[torf_selected]

        enm, _ = calcENM(coords, system, model=model, mask=mask, trim=trim, 
                         n_modes=n_modes, title=labels[i], **kwargs)
        enm.masked = False
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

    if n_modes > 1:
        LOGGER.info('{0} {1} modes were calculated for each of the {2} conformations in {3:.2f}s.'
                            .format(str_modes, model_type, n_confs, time.time()-start))
    else:
        LOGGER.info('{0} {1} mode was calculated for each of the {2} conformations in {3:.2f}s.'
                            .format(str_modes, model_type, n_confs, time.time()-start))

    modeens = ModeEnsemble(title=ensemble.getTitle())
    modeens.addModeSet(enms, weights=ensemble.getWeights(), 
                             label=ensemble.getLabels())
    modeens.setAtoms(ensemble.getAtoms())
    
    if match:
        modeens.match(turbo=turbo, method=method)
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
            raise TypeError('ensemble must be an Ensemble or a ModeEnsemble instance, '
                            'or a list of NMA, Mode, or ModeSet instances.')
    return enms

def calcEnsembleSpectralOverlaps(ensemble, distance=False, turbo=False, **kwargs):
    """Calculate the spectral overlaps between each pair of conformations in the 
    *ensemble*.
    
    :arg ensemble: an ensemble of structures or ENMs 
    :type ensemble: :class: `Ensemble`, :class: `ModeEnsemble`

    :arg distance: if set to **True**, spectral overlap will be converted to spectral 
                   distance via arccos.
    :type distance: bool

    :arg turbo: if **True**, extra memory will be used to remember previous calculation 
                results to accelerate the next calculation, so this option is particularly 
                useful if spectral overlaps of the same ensemble are calculated repeatedly, 
                e.g. using different number of modes. Note that for single calculation, 
                *turbo* will compromise the speed.
                Default is **False**
    :type turbo: bool
    """

    enms = _getEnsembleENMs(ensemble, **kwargs)
    
    overlaps = np.ones((len(enms), len(enms)))
    for i in range(enms.numModeSets()):
        for j in range(i+1, enms.numModeSets()):
            covlap = calcSpectralOverlap(enms[i, :], enms[j, :], turbo=turbo)
            overlaps[i, j] = overlaps[j, i] = covlap

    if distance:
        overlaps = np.arccos(overlaps)

    return overlaps

def calcSignatureSqFlucts(mode_ensemble, **kwargs):
    """
    Get the signature square fluctuations of *mode_ensemble*. 
    
    :arg mode_ensemble: an ensemble of ENMs 
    :type mode_ensemble: :class: `ModeEnsemble`

    :keyword norm: whether to normalize the square fluctuations. Default is **True**
    :type norm: bool

    :keyword scale: whether to rescale the square fluctuations based on the reference. 
                    Default is **False**
    :type scale: bool
    """

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')
    
    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    ifnorm = kwargs.pop('norm', True)
    ifscale = kwargs.pop('scale', False)

    norm = importLA().norm

    V = []
    for i, modes in enumerate(mode_ensemble):
        sqfs = calcSqFlucts(modes)

        if ifnorm:
            sqfs = div0(sqfs, norm(sqfs))
        elif ifscale:
            if i == 0:
                norm0 = norm(sqfs)
            else:
                sqfs = div0(sqfs, norm(sqfs) * norm0)
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
    :type atoms: :class:`.Atomic` 
    """

    from matplotlib.pyplot import figure, plot, fill_between, \
                                  gca, xlabel, ylabel, title, ylim

    linespec = kwargs.pop('linespec', '-')
    zero_line = kwargs.pop('zero_line', False)

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
    :type atoms: :class:`.Atomic` 

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

psplot = showSignature1D

def showSignatureMode(mode_ensemble, **kwargs):
    """Show signature mode profile.

    :arg mode_ensemble: mode ensemble from which to extract an eigenvector
                        If this is not indexed already then index 0 is used by default
    :type mode_ensemble: :class:`ModeEnsemble`    

    :arg atoms: atoms for showing residues along the x-axis
                Default option is to use mode_ensemble.getAtoms()
    :type atoms: :class:`.Atomic`

    :arg scale: scaling factor. Default is 1.0
    :type scale: float    
    """

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    atoms = kwargs.pop('atoms', mode_ensemble.getAtoms())
    scale = kwargs.pop('scale', 1.0)
    mode = mode_ensemble.getEigvec() * scale
    show_zero = kwargs.pop('show_zero', True)
    return showSignature1D(mode, atoms=atoms, show_zero=show_zero, **kwargs)

def showSignatureSqFlucts(mode_ensemble, **kwargs):
    """Show signature profile of square fluctations.

    :arg mode_ensemble: mode ensemble from which to calculate square fluctutations
    :type mode_ensemble: :class:`ModeEnsemble`    

    :arg atoms: atoms for showing residues along the x-axis
                Default option is to use mode_ensemble.getAtoms()
    :type atoms: :class:`.Atomic`

    :arg scale: scaling factor. Default is 1.0
    :type scale: float  

    :arg show_zero: where to show a grey line at y=0
                    Default is False
    :type show_zero: bool      
    """

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    atoms = kwargs.pop('atoms', mode_ensemble.getAtoms())
    scale = kwargs.pop('scale', 1.0)
    sqf = calcSignatureSqFlucts(mode_ensemble) * scale
    show_zero = kwargs.pop('show_zero', False)
    return showSignature1D(sqf, atoms=atoms, show_zero=show_zero, **kwargs)

def calcSignatureCrossCorr(mode_ensemble, norm=True):
    """Calculate the signature cross-correlations based on a :class:`ModeEnsemble` instance.
    
    :arg mode_ensemble: an ensemble of ENMs 
    :type mode_ensemble: :class: `ModeEnsemble`

    :keyword norm: whether to normalize the cross-correlations. Default is **True**
    :type norm: bool
    """
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')
    n_atoms = mode_ensemble.numAtoms()
    n_sets = len(mode_ensemble)

    C = np.zeros((n_sets, n_atoms, n_atoms))
    for i in range(n_sets):
        modes = mode_ensemble[i]
        c = calcCrossCorr(modes, norm=norm)
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

def calcSignaturePerturbResponse(mode_ensemble, **kwargs):
    """Calculate the signature perturbation response scanning based on a :class:`ModeEnsemble` instance.
    
    :arg mode_ensemble: an ensemble of ENMs 
    :type mode_ensemble: :class: `ModeEnsemble`

    """
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')
    n_atoms = mode_ensemble.numAtoms()
    n_sets = len(mode_ensemble)

    P = np.zeros((n_sets, n_atoms, n_atoms))
    E = np.zeros((n_sets, n_atoms))
    S = np.zeros((n_sets, n_atoms))
    for i in range(n_sets):
        modes = mode_ensemble[i]
        prs_mat, eff, sen = calcPerturbResponse(modes, **kwargs)
        P[i, :, :] = prs_mat
        E[i, :] = eff
        S[i, :] = sen

    title_str = '%d modes'%mode_ensemble.numModes()
    weights = mode_ensemble.getWeights()
    if weights is not None:
        W2 = np.zeros((mode_ensemble.numModeSets(), 
                       mode_ensemble.numAtoms(), 
                       mode_ensemble.numAtoms()))
        for i, w in enumerate(weights):
            w2 = np.outer(w, w)
            W2[i, :, :] = w2

        W = weights[:, :, 0]
    labels = mode_ensemble.getLabels()

    # even the original model is 3d, cross-correlations are still 1d
    sig_prs_mat = sdarray(P, title=title_str, weights=W2, labels=labels, is3d=False)
    sig_eff = sdarray(E, title=title_str, weights=W, labels=labels, is3d=False)
    sig_sen = sdarray(S, title=title_str, weights=W, labels=labels, is3d=False)
        
    return sig_prs_mat, sig_eff, sig_sen

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

def calcSignatureOverlaps(mode_ensemble, diag=True, collapse=False):
    """Calculate average mode-mode overlaps for a ModeEnsemble.
    
    If *diag* is **True** (default) then only diagonal values will be calculated. 
    Otherwise, the whole overlap matrices will be calculated.
    
    By default (*collapse* is **False**), the whole overlap matrices are returned as 
    a 4-dimensional sdarray that is a matrix of overlap matrices. 
    
    If *collapse* is **True** then these will be collapsed together, giving a 2-dimensional 
    array for full matrices. This operation is not defined for diagonal values."""
    
    if isinstance(mode_ensemble, ModeEnsemble):
        if not mode_ensemble.isMatched():
            LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                        'Consider running mode_ensemble.match() prior to using this function')

        n_sets = mode_ensemble.numModeSets()
        n_modes = mode_ensemble.numModes()
    else:
        if not isListLike(mode_ensemble):
            raise TypeError('mode_ensemble should be list-like or an instance of ModeEnsemble')

        n_sets = len(mode_ensemble)

        n_modes = np.array([len(modeset) for modeset in mode_ensemble])
        if np.all(n_modes == n_modes[0]):
            n_modes = n_modes[0]
        else:
            raise ValueError('all mode sets in mode_ensemble should have the same number of modes')

    if diag:
        if collapse:
            LOGGER.warn('cannot collapse diagonal values')
        overlaps = np.zeros((n_modes, n_sets, n_sets))
    else:
        if collapse:
            overlaps = np.zeros((n_modes*n_sets, n_modes*n_sets))
        else:
            overlaps = np.zeros((n_modes, n_modes, n_sets, n_sets))

    for i, modeset_i in enumerate(mode_ensemble):
        for j, modeset_j in enumerate(mode_ensemble):
            if j >= i:
                if diag:
                    overlaps[:,i,j] = overlaps[:,j,i] = abs(calcOverlap(modeset_i, 
                                                                        modeset_j, 
                                                                        diag=True))
                else:
                    if collapse:
                        overlaps[i*n_modes:(i+1)*n_modes,
                                 j*n_modes:(j+1)*n_modes] = np.abs(calcOverlap(modeset_i,
                                                                               modeset_j))
                        overlaps[j*n_modes:(j+1)*n_modes,
                                 i*n_modes:(i+1)*n_modes] = np.abs(calcOverlap(modeset_j,
                                                                               modeset_i))
                    else:
                        overlaps[:, :, i, j] = abs(calcOverlap(modeset_i,
                                                               modeset_j,
                                                               diag=False))
                        overlaps[:, :, j, i] = abs(calcOverlap(modeset_j,
                                                               modeset_i,
                                                               diag=False))

    return overlaps


def calcSignatureModes(mode_ensemble):
    """Calculate mean eigenvalues and eigenvectors
    and return a new GNM or ANM object containing them."""

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be a ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    eigvecs = mode_ensemble.getEigvecs().mean()
    eigvals = mode_ensemble.getEigvals().mean()

    if isinstance(mode_ensemble[0].getModel(), GNM):
        ret = GNM('mean of ' + mode_ensemble.getTitle())
    elif isinstance(mode_ensemble[0].getModel(), ANM):
        ret = ANM('mean of ' + mode_ensemble.getTitle())
    else:
        ret = NMA('mean of ' + mode_ensemble.getTitle())

    ret.setEigens(eigvecs, eigvals)
    return ret


def showSignatureOverlaps(mode_ensemble, **kwargs):
    """Show a curve of mode-mode overlaps against mode number
    with shades for standard deviation and range

    :arg diag: Whether to calculate the diagonal values only.
               Default is **False** and :func:`showMatrix` is used.
               If set to **True**, :func:`showSignatureAtomicLines` is used.
    :type diag: bool

    :arg std: Whether to show the standard deviation matrix
              when **diag** is **False** (and whole matrix is shown).
              Default is **False**, meaning the mean matrix is shown.
    type: std: bool
    """
    diag = kwargs.pop('diag', False)
    std = kwargs.pop('std', False)
    from matplotlib.pyplot import xlabel, ylabel, Normalize

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    overlaps = calcSignatureOverlaps(mode_ensemble, diag=diag)

    if diag:
        r, c = np.triu_indices(overlaps.shape[1], k=1)
        overlap_triu = overlaps[:, r, c]

        meanV = overlap_triu.mean(axis=1)
        stdV = overlap_triu.std(axis=1)

        show = showSignatureAtomicLines(meanV, stdV, **kwargs)
        xlabel('Mode index (ref)')
        ylabel('Overlap')
    else:
        r, c = np.triu_indices(overlaps.shape[2], k=1)
        overlap_triu = overlaps[:, :, r, c]

        if std:
            stdV = overlap_triu.std(axis=-1)
            show = showMatrix(stdV, **kwargs)
        else:
            meanV = overlap_triu.mean(axis=-1)
            norm = kwargs.pop('norm', Normalize(0, 1))
            show = showMatrix(meanV, norm=norm, **kwargs)

        xlabel('Mode index (ref)')
        ylabel('Mode index (ref)')
    
    return show

def calcSignatureFractVariance(mode_ensemble):
    """Calculate signature fractional variance for a ModeEnsemble."""
    
    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('mode_ensemble should be an instance of ModeEnsemble')

    if not mode_ensemble.isMatched():
        LOGGER.warn('modes in mode_ensemble did not match cross modesets. '
                    'Consider running mode_ensemble.match() prior to using this function')

    n_sets = mode_ensemble.numModeSets()
    n_modes = mode_ensemble.numModes()

    W = []; is3d = None
    for i in range(n_sets):
        m = mode_ensemble[i]
        if n_modes == 1:
            var = np.array([calcFractVariance(m)])
        else:
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
    :type atoms: :class:`.Atomic` 
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
        if isinstance(patch_i, list):
            colors.append(patch_i[0].get_facecolor())
        else:
            colors.append(patch_i.get_facecolor())

    for i, patch_i in enumerate(patches):
        if isinstance(patch_i, list):
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
    """Show the distribution of variances (cumulative if multiple modes) using 
    :func:`~numpy.histogram`. 
    
    :arg mode_ensemble: an ensemble of modes whose variances are displayed
    :type mode_ensemble: :class:`.ModeEnsemble`

    :arg highlights: labels of conformations whose locations on the bar 
                     will be highlighted by arrows and texts
    :type highlights: list

    :arg fraction: whether the variances should be weighted or not. 
                   Default is **True**
    :type fraction: bool
    """

    from matplotlib.pyplot import figure, gca, annotate, subplots_adjust
    from matplotlib.pyplot import fill_between, xlabel, yticks, xlim
    from matplotlib.figure import Figure
    from matplotlib.colors import Normalize, NoNorm
    from matplotlib import cm, colors
    
    fig = kwargs.pop('figure', None)
    fract = kwargs.pop('fraction', True)
    bins = kwargs.pop('bins', 50)
    cmap = kwargs.pop('cmap', 'Reds')

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
    
    hist, edges = np.histogram(variances, bins=bins)
    color_norm  = colors.Normalize(vmin=hist.min(), vmax=hist.max())
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap=cmap)
    colors = scalar_map.to_rgba(hist)

    areas = []
    for i in range(len(hist)):
        x = [edges[i], edges[i+1]]
        area = fill_between(x, [0, 0], [1, 1], color=colors[i])
        areas.append(area)

    if not highlights:
        highlights = []

    indices = []; labels = []
    ens_labels = mode_ensemble.getLabels()
    for hl in highlights:
        if isinstance(hl, str):
            if not ens_labels:
                raise TypeError('highlights should be a list of integers because '
                                    'mode_ensemble has no label')
            index = ens_labels.index(hl)
            if isinstance(highlights, dict):
                label = highlights[hl]
            else:
                label = hl
        else:
            try:
                index = int(hl)
            except:
                raise TypeError('highlights should be a list of integers or strings') 
            if isinstance(highlights, dict):
                label = highlights[hl]
            else:
                label = ens_labels[index] if ens_labels else str(index)
        indices.append(index)
        labels.append(label)

    annotations = []
    for i, label in zip(indices, labels):
        x = variances[i]
        an = annotate(label, xy=(x, 1), xytext=(x, ratio), arrowprops=arrowprops)
        annotations.append(an)

    # for i in range(len(variances)):
    #     x = variances[i]
    #     plot([x, x], [0, 1], 'w')

    xlabel('Variances')
    yticks([])
    xlim([variances.min(), variances.max()])

    if SETTINGS['auto_show']:
        showFigure()
    return areas, annotations

def saveModeEnsemble(mode_ensemble, filename=None, atoms=False, **kwargs):
    """Save *mode_ensemble* as :file:`filename.modeens.npz`.  If *filename* 
    is **None**, title of the *mode_ensemble* will be used as the 
    filename, after ``" "`` (white spaces) in the title are replaced with 
    ``"_"`` (underscores).  Upon successful completion of saving, filename 
    is returned. This function makes use of :func:`~numpy.savez_compressed` 
    function."""

    if not isinstance(mode_ensemble, ModeEnsemble):
        raise TypeError('invalid type for mode_ensemble, {0}'
                        .format(type(mode_ensemble)))
    if len(mode_ensemble) == 0:
        raise ValueError('mode_ensemble instance does not contain data')

    attr_list = ['_modesets', '_title', '_labels', '_weights', '_matched', '_reweighted']
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

    if not 'allow_pickle' in kwargs:
        kwargs['allow_pickle'] = True

    data = np.load(filename, **kwargs)
    
    weights = getValue(data, '_weights', None)
    labels = getValue(data, '_labels', None)
    matched = getValue(data, '_matched', False)
    title = getValue(data, '_title', None)
    modesets = getValue(data, '_modesets', [])
    atoms = getValue(data, '_atoms', [None])[0]
    reweighted = getValue(data, '_reweighted', False)

    if isinstance(title, np.ndarray):
        title = np.asarray(title, dtype=str)
    title = str(title)

    if isinstance(modesets, np.ndarray):
        modesets = modesets.tolist()
    while (None in modesets):
        modesets.remove(None)

    if labels is not None:
        char = labels.dtype.char
        if char in 'SU' and char != DTYPE:
            labels = labels.astype(str)
        labels = labels.tolist()

    if isinstance(matched, np.ndarray):
        matched = matched.tolist()

    if isinstance(reweighted, np.ndarray):
        reweighted = reweighted.tolist()

    modeens = ModeEnsemble(title=title)
    modeens._weights = weights
    modeens._labels = labels
    modeens._matched = matched
    modeens._reweighted = reweighted
    modeens._modesets = modesets

    if atoms is not None:
        if isinstance(atoms, AtomGroup):
            data = atoms._data
        else:
            data = atoms._ag._data
            
        for key in data:
            arr = data[key]
            char = arr.dtype.char
            if char in 'SU' and char != DTYPE:
                arr = arr.astype(str)
                data[key] = arr

    modeens._atoms = atoms

    return modeens

def saveSignature(signature, filename=None, **kwargs):
    """Save *signature* as :file:`filename.sdarray.npz`.  If *filename* 
    is **None**, title of the *signature* will be used as the 
    filename, after ``" "`` (white spaces) in the title are replaced with 
    ``"_"`` (underscores).  Upon successful completion of saving, filename 
    is returned. This function makes use of :func:`~numpy.savez_compressed` 
    function."""

    if not isinstance(signature, sdarray):
        raise TypeError('invalid type for signature, {0}'
                        .format(type(signature)))

    attr_list = ['_title', '_labels', '_is3d', '_weights', '_oneset', '_array']
    attr_dict = {}
    
    for attr in attr_list:
        if attr == '_array':
            value = np.asarray(signature)
        else:
            value = getattr(signature, attr)
        if value is not None:
            attr_dict[attr] = value

    if filename is None:
        filename = signature.getTitle().replace(' ', '_')
    
    suffix = '.sdarray'
    if not filename.lower().endswith('.npz'):
        if not filename.lower().endswith(suffix):
            filename += suffix + '.npz'
        else:
            filename += '.npz'
            
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez_compressed(ostream, **attr_dict)
    ostream.close()

    return filename

def loadSignature(filename, **kwargs):
    """Returns :class:`sdarray` instance after loading it from file (*filename*).
    This function makes use of :func:`numpy.load` function.  See
    also :func:`saveSignature`."""

    if not 'encoding' in kwargs:
        kwargs['encoding'] = 'latin1'
    data = np.load(filename, **kwargs)
    
    weights = getValue(data, '_weights', None)
    labels = getValue(data, '_labels', None)
    title = getValue(data, '_title', None)
    is3d = getValue(data, '_is3d', False)
    oneset = getValue(data, '_oneset', False)
    array = getValue(data, '_array', None)

    if isinstance(title, np.ndarray):
        title = np.asarray(title, dtype=str)
    title = str(title)

    if isinstance(is3d, np.ndarray):
        is3d = bool(is3d)

    if isinstance(oneset, np.ndarray):
        oneset = bool(oneset)

    if labels is not None:
        labels = labels.tolist()

    signature = sdarray(array, weights=weights, labels=labels, title=title, 
                      is3d=is3d, oneset=oneset)

    return signature

def calcSubfamilySpectralOverlaps(mode_ens, subfamily_dict, **kwargs):
    """Calculate average spectral overlaps (or distances) within and between 
    subfamilies in a mode ensemble defined using a dictionary where each key is an 
    ensemble member and the associate value is a subfamily name.
    
    To use a range of modes, please index the mode ensemble e.g. 
    mode_ens=mode_ensemble[:,3:20] to use modes 4 to 20 inclusive. 
    Alternatively, there is the option to provide first and last 
    keyword arguments, which would be used as the 3 and 20 above.
    
    :arg mode_ensemble: an ensemble of modes corresponding to a set of modes 
        for each family member
    :type mode_ensemble: :class:`.ModeEnsemble`
    
    :arg subfamily_dict: a dictionary providing a subfamily label for 
        each family member
    :type subfamily_dict: dict
    
    :keyword first: the first index for a range of modes
    :type first: int
    
    :keyword last: the last index for a range of modes
    :type last: int
    
    :keyword remove_small: whether to remove small subfamilies with 
        fewer than 4 members. Default is True
    :type remove_small: bool
    
    :keyword return_reordered_subfamilies: whether to return the reordered 
        subfamilies in addition to the matrix. Default is False
    type return_reordered_subfamilies: bool 
    
    """

    if not isinstance(mode_ens, ModeEnsemble):
        raise TypeError('mode_ens should be a mode ensemble')

    if not isinstance(subfamily_dict, dict):
        raise TypeError('subfamily_dict should be a dictionary')

    if any([label not in list(subfamily_dict.keys()) for label in mode_ens.getLabels()]):
        raise ValueError('The are member labels in mode_ens with no associated entry in subfamily_dict')

    first = kwargs.get('first', 0)
    if first is not None:
        if not isinstance(first,int):
            raise TypeError('first should be an integer')

    last = kwargs.get('last', -1)
    if last is not None:
        if not isinstance(last,int):
            raise TypeError('last should be an integer')

    try:
        mode_ens = mode_ens[:,first:last]
    except:
        try:
            mode_ens = mode_ens[:,first:]
        except:
            raise ValueError('first is not a valid index for indexing mode_ens')
        
        try:
            mode_ens = mode_ens[:,:last]
        except:
            raise ValueError('last is not a valid index for indexing mode_ens')

    first_mode_index = mode_ens.getIndices()[0,0]
    last_mode_index = mode_ens.getIndices()[0,-1]
    LOGGER.info('The mode range used for this analysis is {0} to {1}'
                .format(first_mode_index+1,last_mode_index+1))

    tree_labels = mode_ens.getLabels()
    distance = kwargs.get('distance',True)
    distm = calcEnsembleSpectralOverlaps(mode_ens, distance=distance)

    subfamilies = np.unique(list(subfamily_dict.values()))
    reverse_dict = dict()
    for i in subfamilies:
        reverse_dict[i] = []
    for i in range(len(tree_labels)):
        subfamily_i = subfamily_dict[tree_labels[i]]
        if subfamily_i in subfamilies:
            reverse_dict[subfamily_i].append(i)

    remove_small = kwargs.get('remove_small',True)
    if remove_small:
        temp_dict = dict()
        for key, value in reverse_dict.items():
            if len(value) >= 4:
                temp_dict[key] = value
        reverse_dict = temp_dict

    N_group = len(reverse_dict)

    subfamily_overlap_matrix = []
    for subfamily_i in subfamilies:
        index_i = reverse_dict[subfamily_i]
        for subfamily_j in subfamilies:
            index_j = reverse_dict[subfamily_j]
            temp_sub_matrix = distm[np.ix_(index_i, index_j)]
            if subfamily_i == subfamily_j:
                subfamily_overlap_matrix.append(temp_sub_matrix[np.triu_indices(
                    np.shape(temp_sub_matrix)[0])].mean())
            else:
                subfamily_overlap_matrix.append(temp_sub_matrix.mean())
    subfamily_overlap_matrix = np.asarray(subfamily_overlap_matrix)
    subfamily_overlap_matrix = subfamily_overlap_matrix.reshape([N_group, N_group])

    return_reordered_subfamilies = kwargs.get('return_reordered_subfamilies',False)
    if return_reordered_subfamilies:
        return subfamily_overlap_matrix, subfamilies

    return subfamily_overlap_matrix

def showSubfamilySpectralOverlaps(mode_ens, subfamily_dict, **kwargs):
    """Calculate and show the matrix of spectral overlaps or distances averaged 
    over subfamilies. Inputs are the same as calcSubfamilySpectralOverlaps
    plus the following and those of showDomainBar if you wish.

    :keyword show_subfamily_bar: whether to show the subfamilies as colored bars 
        using showDomainBar. Default is False
    :type show_subfamily_bar: bool
    """
    kwargs['return_reordered_subfamilies'] = True

    subfamily_overlap_matrix, subfamilies = calcSubfamilySpectralOverlaps(mode_ens, subfamily_dict, **kwargs)
    show = showMatrix(subfamily_overlap_matrix, origin='lower',
                      xticklabels=subfamilies, yticklabels=subfamilies,
                      allticks=True, vmin=0., vmax=1.6)

    show_subfamily_bar = kwargs.get('show_subfamily_bar',False)
    text = kwargs.pop('text',False)
    if show_subfamily_bar:
        showDomainBar(subfamilies, axis='x', text=text, **kwargs)
        showDomainBar(subfamilies, axis='y', text=text, **kwargs)

    if SETTINGS['auto_show']:
        showFigure()

    return show
