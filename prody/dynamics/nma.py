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
""":mod:`nma` module defines classes for normal mode analysis. 

Classes:

  * Normal Mode Analysis (:class:`NMA`) 
  * Base class for GNM/ANM (:class:`GNMBase`)
  * Gaussian Network Model (:class:`GNM`)
  * Anisotropic Network Model (:class:`ANM`)
  * Principal Component Analysis (:class:`PCA`)
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import time

import numpy as np

import prody
from prody import ProDyLogger as LOGGER
from .modes import Mode


__all__ = ['GNM', 'ANM', 'PCA', 'NMA', 'ModeSet']


ZERO = 1e-8

class NMAError(Exception):
    pass

class NMA(object):
    
    """Base class for ENM, GNM, and PCA."""


    def __init__(self, name):
        """Initialize a Normal Mode analysis with a given name."""
        self._name = str(name)
        self._modes = []
        self._n_modes = 0
        self._cov = None

        self._n_atoms = 0
        self._dof = 0
        
        self._array = None                # modes/eigenvectors
        self._eigvals = None
        self._vars = None                # evs for PCA, inverse evs for ENM
        self._trace = None
        
        self._is3d = True               # is set to false for GNM


    def __len__(self):
        return self._n_modes
        
    def __getitem__(self, index):
        if isinstance(index, int):
            return self.getMode(index)
        elif isinstance(index, slice):
            modes = self._modes[index]
            return ModeSet(self, np.arange(*index.indices(len(self))))
        else:        
            raise IndexError('indices may be an integer or a slice')
        
    def __iter__(self):
        for i in xrange(self._n_modes):
            yield self.getMode(i)
    
    def __repr__(self):
        return '<NMA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._name, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'NMA {0:s}'.format(self._name)


    def _reset(self):
        
        self._n_modes = 0        
        self._cov = None
        
        self._n_atoms = 0
        self._dof = 0
        
        self._array = None
        self._eigvals = None
        self._vars = None
        self._trace = None
        
        self._is3d = True
    
    def is3d(self):
        return self._is3d
    
    def getNumOfAtoms(self):
        """Return number of modes."""
        return self._n_atoms
    
    def getNumOfModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes)."""
        return self._n_modes
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom."""
        return self._dof
        
    def getName(self):
        """Return name."""
        return self._name
    
    def setName(self, name):
        """Set name."""
        self._name = str(name)
    
    def getMode(self, index):
        """Return mode at given index."""
        mode = self._modes[index]
        if mode is None:
            mode = Mode(self, index)
            self._modes[index] = mode
        return mode

    def getModes(self):
        """Return all modes in a list."""
        getMode = self.getMode
        return [getMode(i) for i in range(len(self))]

    
    def getEigenvalues(self):
        """Return eigenvalues."""
        return self._eigvals.copy()

    def getEigenvectors(self):
        """Return eigenvectors."""
        return self._array.copy()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues)."""
        return self._vars.copy()

    def getArray(self):
        """Return eigenvectors."""
        return self._array.copy()
        
    def getCovariance(self):
        """Return covariance matrix, if needed after calculating it using available modes."""
        if self._cov is None:
            if self._array is None:
                return None
            self._cov = np.dot(self._array,np.dot(np.diag(self._vars), self._array.T))
        return self._cov
        
    def calcModes(self):
        pass

class ModeSet(object):
    """A class for providing access to data for a subset of modes.
    
    Instances are obtained by slicing an NMA model (:class:`ANM`, :class:`GNM`, 
    or :class:`PCA`).
    
    ModeSet's contain a reference to the model and a list of mode indices.
    Methods common to NMA models are also defined for mode sets.
    
    """
    
    __slots__ = ['_model', '_indices', '_slice']
    
    def __init__(self, model, indices):
        if not isinstance(model, NMA):
            raise TypeError('model must be an NMA, not {0:s}'.format(type(model)))
        self._model = model
        self._indices = np.array(indices, np.int64)
        self._slice = prody.getIntAsStr(indices+1, sep=' to ')
        
    def __len__(self):
        return len(self._indices)
        
    def __iter__(self):
        for i in self._indices:
            yield self._model.getMode(i)
    
    def __repr__(self):
        return '<ModeSet: {0:s} from {1:s} ({2:d} modes)>'.format(self._slice,
                self._model._name, len(self))

    def __str__(self):
        return 'ModeSet {0:s} from {1:s}'.format(self._slice, self._model._name)
    
    def is3d(self):
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Return number of nodes."""
        return self._model._n_atoms
    
    def getNumOfModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes)."""
        return len(self._indices)
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom."""
        return self._model._dof
        
    def getModes(self):
        """Return all modes in the subset in a list."""
        getMode = self._model.getMode
        return [getMode(i) for i in self._indices]
    
    def getName(self):
        """Return name."""
        return str(self)
    
    def getModel(self):
        """The NMA instance that the mode belongs to"""
        return self._model
    
    def getIndices(self):
        """Return indices of modes in the mode set."""
        return self._indices
    
    def getEigenvalues(self):
        """Return eigenvalues."""
        return self._model._eigvals[self._indices].copy()

    def getEigenvectors(self):
        """Return eigenvectors."""
        return self._model._array[:, self._indices].copy()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues)."""
        return self._model._vars[self._indices].copy()

    def getArray(self):
        """Return eigenvectors."""
        return self._model._array[:, self._indices].copy()
        
    def getCovariance(self):
        """Return covariance matrix calculated for modes in the set."""
        array = self.getArray()
        return np.dot(array, np.dot(np.diag(self.getVariances()), array.T))
   
    

class GNMBase(NMA):
    
    """Class for Gaussian Network Model analysis of proteins."""
    

    def __init__(self, name):
        NMA.__init__(self, name)
        self._is3d = False
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        
    def __repr__(self):
        return '<GNM: {0:s} ({1:d} modes, {2:d} nodes)>'.format(
                                    self._name, self.__len__(), self._n_atoms)
    
    def __str__(self):
        return 'GNM {0:s}'.format(self._name)
    
    def _reset(self):
        NMA._reset(self)
        self._cutoff = None
        self._gamma = None
        
    
    def getCutoff(self):
        """Return cutoff distance."""
        return self._cutoff
    
    def getGamma(self):
        """Return spring constant."""
        return self._gamma

    def getKirchhoff(self):
        """Return Kirchhoff matrix."""
        return self._kirchhoff.copy()    

class GNM(GNMBase):
    
    """A class for Gaussian Network Model (GNM) analysis of proteins ([IB97]_, [TH97]_)."""
    
    def setKirchhoff(self, kirchhoff):
        """Set Kirchhoff matrix."""
        if not isinstance(kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be an ndarray')
        elif not (kirchhoff.ndim == 2 and kirchhoff.shape[0] == kirchhoff.shape[1]):
            raise TypeError('kirchhoff must be square matrix')
        self._reset()
        self._kirchhoff = kirchhoff
        self._n_atoms = kirchhoff.shape[0]
        self._dof = kirchhoff.shape[0]
    
    def buildKirchhoff(self, coords, cutoff=10., gamma=1., masses=None):
        """Build Kirchhoff matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        
        :arg cutoff: cutoff distance (A) for pairwise interactions.
        :type cutoff: float, default is 10.0
        
        :arg gamma: spring constant
        :type gamma: float, default is 1.0
        
        *masses* is not used yet.

        """
        if prody.KDTree is None:
            prody.importBioKDTree()
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoordinates()
            except AttributeError:
                raise TypeError('coords must be an ndarray instance or '
                                'must contain getCoordinates as an attribute')

        if coords.ndim != 2:
            raise ValueError('coords must be a 2d array')
        elif coords.shape[1] != 3:
            raise ValueError('shape of coords must be (n_atoms,3)')
        elif coords.dtype != np.float64:
            try:
                coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
                                 
        cutoff = float(cutoff)
        gamma = float(gamma)
        n_atoms = coords.shape[0]
        kdtree = prody.KDTree(3)
        start = time.time()
        kdtree.set_coords(coords) 
        kdtree.all_search(cutoff)
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        for i, j in kdtree.all_get_indices():
            kirchhoff[i, j] = -gamma
            kirchhoff[j, i] = -gamma
            kirchhoff[i, i] += gamma
            kirchhoff[j, j] += gamma
        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'
                         .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._cutoff = cutoff
        self._gamma = gamma
        self._n_atoms = n_atoms
        self._dof = n_atoms
        
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if prody.la is None:
            prody.importScipyLinalg()
        
        if self._kirchhoff is None:
            raise RuntimeError('Kirchhoff matrix is not set')
            
        start = time.time()
        shift = 0
        if n_modes is None:
            eigvals = None
            n_modes = self._dof 
        else: 
            n_modes = int(n_modes)
            eigvals = (0, n_modes + shift)
        if eigvals: 
            turbo = False
        values, vectors = prody.la.eigh(self._kirchhoff, turbo=turbo, eigvals=eigvals)
        if zeros:
            shift = -1
        self._eigvals = values[1+shift:]
        self._vars = 1 / self._eigvals
        self._trace = self._vars.sum()
        self._array = vectors[:, 1+shift:]
        self._n_modes = len(self._eigvals)
        self._modes = [None] * self._n_modes
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                          ''.format(self._n_modes, time.time()-start))


class ANM(GNMBase):

    """Class for Anisotropic Network Model (ANM) analysis of proteins ([PD00]_, [ARA01]_)"""

    def __init__(self, name):
        GNMBase.__init__(self, name)
        self._is3d = True
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        self._hessian = None


    def __repr__(self):
        return '<ANM: {0:s} ({1:d} modes, {2:d} nodes)>'.format(
                                    self._name, self.__len__(), self._n_atoms)
                                    
    def __str__(self):
        return 'ANM {0:s}'.format(self._name)

    def _reset(self):
        GNMBase._reset(self)
        self._kirchhoff = None
        
    def getHessian(self):
        """Return a copy of Hessian matrix."""
        return self._hessian.copy()
    
    def setHessian(self, hessian):
        """Set Hessian matrix."""
        if not isinstance(hessian, np.ndarray):
            raise TypeError('hessian must be an ndarray')
        elif not (hessian.ndim == 2 and hessian.shape[0] == hessian.shape[1]):
            raise TypeError('hessian must be square matrix')
        self._reset()
        self._hessian = hessian
        self._dof = hessian.shape[0]
        self._n_atoms = self._dof / 3 

    def buildHessian(self, coords, cutoff=15., gamma=1., masses=None):
        """Build Hessian matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        
        :arg cutoff: cutoff distance (A) for pairwise interactions
        :type cutoff: float, default is 15.0
        
        :arg gamma: spring constant
        :type gamma: float, default is 1.0
        
        *masses* is not used yet.
                        
        """
        if prody.KDTree is None:
            prody.importBioKDTree()
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoordinates()
            except AttributeError:
                raise TypeError('coords must be an ndarray instance or '
                                'must contain getCoordinates as an attribute')

        if coords.ndim != 2:
            raise ValueError('coords must be a 2d array')
        elif coords.shape[1] != 3:
            raise ValueError('shape of coords must be (n_atoms,3)')
        elif coords.dtype != np.float64:
            try:
                coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        
        cutoff = float(cutoff)
        gamma = float(gamma)
        n_atoms = coords.shape[0]
        dof = n_atoms * 3
        kdtree = prody.KDTree(3)
        start = time.time()
        kdtree.set_coords(coords) 
        kdtree.all_search(cutoff)
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        hessian = np.zeros((dof, dof), 'd')

        for i, j in kdtree.all_get_indices():
            i2j = coords[j] - coords[i]
            dist2 = np.dot(i2j, i2j)
            super_element = np.outer(i2j, i2j) / dist2 * gamma 
            res_i3 = i*3
            res_i33 = res_i3+3
            res_j3 = j*3
            res_j33 = res_j3+3
            hessian[res_i3:res_i33, res_j3:res_j33] = -super_element
            hessian[res_j3:res_j33, res_i3:res_i33] = -super_element
            hessian[res_i3:res_i33, res_i3:res_i33] += super_element
            hessian[res_j3:res_j33, res_j3:res_j33] += super_element
            kirchhoff[i, j] = -gamma
            kirchhoff[j, i] = -gamma
            kirchhoff[i, i] += gamma
            kirchhoff[j, j] += gamma

        LOGGER.info('Hessian was built in {0:.2f}s.'
                         .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._gamma = gamma
        self._cutoff = cutoff
        self._n_atoms = n_atoms
        self._dof = dof

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if prody.la is None:
            prody.importScipyLinalg()

        if self._hessian is None:
            raise RuntimeError('Hessian matrix is not set')
            
        start = time.time()
        shift = 5
        if n_modes is None:
            eigvals = None
            n_modes = self._dof 
        else: 
            n_modes = int(n_modes)
            eigvals = (0, n_modes + shift)
        if eigvals: 
            turbo = False
        values, vectors = prody.la.eigh(self._hessian, turbo=turbo, eigvals=eigvals)
        n_zeros = sum(values < ZERO)
        if n_zeros < 6: 
            LOGGER.warning('Less than 6 zero eigenvalues are calculated.')
            shif = n_zeros - 1
        elif n_zeros > 6: 
            LOGGER.warning('More than 6 zero eigenvalues are calculated.')
            shif = n_zeros - 1
        if zeros:
            shift = -1
        self._eigvals = values[1+shift:]
        self._vars = 1 / self._eigvals
        self._trace = self._vars.sum()
        self._array = vectors[:, 1+shift:]
        self._n_modes = len(self._eigvals)
        self._modes = [None] * self._n_modes
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                          ''.format(self._n_modes, time.time()-start))

class PCAError(Exception):
    pass
    

class PCA(NMA):
    """A class for Principal Component Analysis (PCA) of conformational ensembles 
    (also known as Essential Dynamics Analysis (EDA) in [AA93]_)."""

    def __init__(self, name):
        NMA.__init__(self, name)
    
    def __repr__(self):
        return '<PCA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._name, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'PCA {0:s}'.format(self._name)
    
    def setCovariance(self, covariance):
        if not isinstance(covariance, np.ndarray):
            raise TypeError('covariance must be an ndarray')
        elif not (covariance.ndim == 2 and covariance.shape[0] == covariance.shape[1]):
            raise TypeError('covariance must be square matrix')
        self._reset()
        self._cov = covariance
        self._dof = covariance.shape[0]
        self._n_atoms = self._dof / 3
        self._trace = self._cov.trace()
    

    def buildCovariance(self, coordsets, weights=None):
        """Build a weighted covariance matrix for coodsets.
        
        *coordsets* must have the following methods: 
            * ``getCoordinates``
            * ``getCoordsets``
            * ``getNumOfCoordsets``
            * ``getNumOfAtoms``
        
        :class:`Ensemble`, :class:`AtomGroup`, :class:`Selection`,
        :class:`Chain`, and :class:`Residue` instances are acceptable.
        
        If *weights* is ``None``, but *coordsets* has getWeights method,
        weights from that method will be used. 
        
        """
        pass
        start = time.time()
        n_atoms = coordsets.getNumOfAtoms()
        dof = n_atoms * 3
        coordinates = coordsets.getCoordinates()
        try:
            acsi = coordsets.getActiveCoordsetIndex()
            indices = range(coordsets.getNumOfCoordsets())
            indices.pop(acsi)
            conformations = coordsets.getCoordsets(indices)
        except:
            conformations = coordsets.getCoordsets()
        n_confs = conformations.shape[0]
        if weights is None:
            try:
                weights = coordsets.getWeights()
            except:
                pass
        if weights is None:
            d_xyz = (conformations - coordinates)
            d_xyz = d_xyz.reshape((n_confs, dof))
            self._cov = np.dot(d_xyz.T, d_xyz) / n_confs
        else:
            d_xyz = ((conformations - coordinates) * weights)
            d_xyz = d_xyz / (weights + (weights == 0)) 
            d_xyz = d_xyz.reshape((n_confs, dof))
            which_axis = weights.ndim-1
            divide_by = weights.repeat(3, axis=which_axis).reshape((n_confs, dof))
            self._cov = np.dot(d_xyz.T, d_xyz) / np.dot(divide_by.T, divide_by)
        self._trace = self._cov.trace()
        self._dof = dof
        self._n_atoms = n_atoms
        LOGGER.info('Covariance matrix was calculated in {0:.2f}s.'
                    .format(time.time()-start))
        
    def calcModes(self, n_modes=20, turbo=True):
        """Calculate essential modes.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if prody.la is None:
            prody.importScipyLinalg()

        start = time.time()
        dof = self._dof
        if n_modes: 
            eigvals = (dof - n_modes, dof - 1)
        else: 
            eigvals = None
            n_modes = dof
        values, vectors = prody.la.eigh(self._cov, turbo=turbo, 
                                  eigvals=eigvals)
        # Order by descending SV
        revert = range(len(values)-1, -1, -1)
        values = values[revert]
        vectors = vectors[:, revert]
        which = values > 1e-8
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._n_modes = len(self._eigvals)
        self._modes = [None] * self._n_modes
        LOGGER.debug('{0:d} essential modes were calculated in {1:.2f}s.'
                         .format(self._n_modes, time.time()-start))
