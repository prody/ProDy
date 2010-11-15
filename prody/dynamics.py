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
"""This module defines classes and functions for dynamics analysis. 

Classes
-------

  * :class:`ANM`
  * :class:`GNM`
  * :class:`PCA`
  * :class:`Mode`
  * :class:`ModeSet`
  * :class:`Vector`
  
Base Classes
------------

  * :class:`NMA`
  * :class:`GNMBase`
  * :class:`VectorBase`

Inheritance Diagram
-------------------

.. inheritance-diagram:: prody.dynamics
   :parts: 1

Functions
---------

Many of the functions documented in this page accepts a *modes* argument 
(may also have different names).
This argument may be:

  * an NMA model, which may be an instance of one of :class:`ANM`, :class:`GNM`, :class:`PCA`.
  * a :class:`Mode` instance obtained by indexing an NMA model, e.g. ``nma[0]``
  * a :class:`ModeSet` instance obtained by slicing an NMA model, e.g. ``nma[:10]``

Some of these functions may also accept :class:`Vector` instances as *mode* argument.

**Short-hand functions**:

  * :func:`getANM`
  * :func:`getGNM`

**Analysis**:

  * :func:`getCovariance`
  * :func:`getCrossCorrelations`
  * :func:`getSqFlucts`  

**Write data**:

  * :func:`writeArray`
  * :func:`writeModes`
  * :func:`writeNMD`
  * :func:`writeOverlapTable`
        
**Visualization**:

  * :func:`getVMDpath`
  * :func:`setVMDpath`
  * :func:`viewNMDinVMD`
    
**Compare NMA models**:

  * :func:`getOverlap`
  * :func:`getCumulativeOverlap`
  * :func:`getCumulativeOverlapArray`
  * :func:`getSubspaceOverlap`
  * :func:`printOverlapTable`
  
**Other**:

  * :func:`getProjection`
  * :func:`reduceModel`
  * :func:`sampleModes`

**Plotting Functions**:

Plotting functions are called by the name of the plotted data/property and are
prefixed with "show". Function documentations include the :mod:`matplotlib.pyplot` 
function utilized for actual plotting. Arguments and keyword arguments are passed to the
Matplotlib functions.  


  * :func:`showContactMap`
  * :func:`showCrossCorrelations`
  * :func:`showCumulativeOverlap`
  * :func:`showCumFractOfVariances`
  * :func:`showFractOfVariances`
  * :func:`showMode`
  * :func:`showOverlap`
  * :func:`showOverlapTable`
  * :func:`showProjection`
  * :func:`showSqFlucts`
    
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import os.path
import time
import os
import gzip

import numpy as np
pl = None

import prody
from .atomic import *
from .ensemble import *
from prody import ProDyLogger as LOGGER

__all__ = ['ANM', 'GNM', 'PCA', 'EDA', 'Mode', 'ModeSet', 'Vector', 
           
           'NMA', 'GNMBase', 'VectorBase',
           
           'getANM', 'getGNM', 
           
           'getCovariance', 'getCrossCorrelations', 'getSqFlucts', 
           
           'writeArray', 'writeModes', 'writeNMD', 'writeOverlapTable',
           
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD', 
           
           'getOverlap', 'getCumulativeOverlap', 'getCumulativeOverlapArray', 
           'getSubspaceOverlap', 'printOverlapTable',
           
           'getProjection', 'reduceModel', 'sampleModes',
            
           'showContactMap', 'showCrossCorrelations', 'showCumulativeOverlap', 
           'showCumFractOfVariances', 'showFractOfVariances', 'showMode', 
           'showOverlap', 'showOverlapTable', 'showProjection', 'showSqFlucts',
           ]

VMDPATH = '/usr/local/bin/vmd'
ZERO = 1e-8

class NMAError(Exception):
    pass

class VectorBase(object):
    """A base class for :class:`Mode` and :class:`Vector`.
    
    This base class defines some shared methods, such as scalar multiplication
    or addition of mode instances.
    
    Defined operations are:
        
        * Absolute value (abs(mode)) returns mode length
        * Additive inverse (-mode) 
        * Mode addition (mode1 + mode2)
        * Mode substraction (mode1 - mode2)
        * Scalar mulitplication (x*mode or mode*x)
        * Division by a scalar (mode/x)
        * Dot product (mode1*mode2)
        * Power (mode**x)
    
    """
    
    def __abs__(self):
        return np.sqrt((self.getArray()**2).sum())
    
    def __neg__(self):
        return Vector('-({0:s})'.format(str(self)), -self.getArray(), self.is3d())
    
    def __div__(self, other):
        if isinstance(other, (int, float, long)):
            return Vector('({1:s})/{0}'.format(other, str(self)), 
                              self.getArray() / other, self.is3d())
        else:
            raise TypeError('{0} is not a scalar'.format(other))
    
    def __idiv__(self, other):
        return self.__div__(other)
    
    def __mul__(self, other):
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector('{0}*({1:s})'.format(other, str(self)),
                              other * self.getArray(), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self.getArray(), other.getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
    
    def __rmul__(self, other):   
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector('{0}*({1:s})'.format(other, str(self)),
                              other * self.getArray(), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self.getArray(), other.getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
            
    def __imul__(self, other):
        return self.__mul__(other)    

    def __add__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector('{0:s} + {1:s}'.format(str(self), str(other)),
                              self.getArray() + other.getArray(), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __radd__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector('{0:s} + {1:s}'.format(str(other), str(self)),
                               self.getArray() + other.getArray(), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))
        
    def __iadd__(self, other):
        return self.__add__(other)   

    def __sub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector('{0:s} - {1:s}'.format(str(self), str(other)), 
                                self.getArray() - other.getArray(), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __rsub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return  Vector('{0:s} - {1:s}'.format(str(other), str(self)),
                               other.getArray() - self.getArray(), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __isub__(self, other):
        return self.__sub__(other)   

    def __pow__(self, other):
        if isinstance(other, (int, float, long)): 
            return Vector('({0:s})**{1}'.format(str(self), other),
                              self.getArray() ** other, self.is3d())
        else:
            raise TypeError('{0} is not a scalar'.format(other))

    def getArray(self):
        """Return array."""
        pass
        
    def getArrayNx3(self):
        """Return array with shape (N, 3)."""
        if self.is3d():
            return self.getArray().reshape((len(self)/3, 3))
        else:
            return self.getArray()

class Mode(VectorBase):
    
    """A class to provide access to and operations on mode data.
        
    """
    __slots__ = ['_model', '_index']
    
    def __init__(self, model, index):
        """Initialize mode object as part of an NMA model.
        
        :arg model: a normal mode analysis instance
        :type model: :class:`ANM`, :class:`GNM`, or :class:`PCA` 
        :arg index: index of the mode 
        :type index: int
            
        """
        self._model = model
        self._index = index
        
    def __len__(self):
        return self._model._dof
    
    def __repr__(self):
        return '<Mode: {0:d} from {1:s}>'.format(self._index+1, self._model._name)

    def __str__(self):
        return 'Mode {0:d} from {1:s}'.format(self._index+1, self._model._name)

    def is3d(self):
        return self._model._is3d

    
    def getNumOfAtoms(self):
        """Return number of nodes."""
        return self._model._n_atoms
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom."""
        return self._model._dof
    
    def getName(self):
        """A descriptive name for the mode instance."""
        return str(self)
    
    def getIndex(self):
        """Return index of the mode."""
        return self._index
    
    def getModel(self):
        """The NMA instance that the mode belongs to"""
        return self._model
    
    def getArray(self):
        """Normal mode array (eigenvector)"""
        return self._model._array[:, self._index].copy()
    getEigenvector = getArray
    
    
    def getEigenvalue(self):
        """Normal mode eigenvalue"""
        return self._model._eigvals[self._index]
    
    def getVariance(self):
        """Variance along the mode. 
        
        If the model is an ENM, inverse of the eigenvalue is returned.
        
        """
        return self._model._vars[self._index]

    def getFractOfVariance(self):
        """Return fraction of variance explained by the mode.
        
        Fraction of variance is the ratio of the variance along this mode to 
        the trace of the covariance matrix.
        
        See :meth:`getVariance`
        
        """
        return self.getVariance() / self._model._trace
    
    def getCollectivity(self, masses=None):
        """Return collectivity of the mode.
        
        :arg masses: atomic masses
        :type masses: :class:`numpy.ndarray`
        
        This function implements collectivity as defined in equation 5 of [BR95]_.
        
        If *masses* are provided, they will be incorporated in the calculation.
        Otherwise, atoms are assumed to have uniform masses.
        
        """
        if masses is not None:
            if len(masses) != self._model._n_atoms: 
                raise ValueError('length of massesmust be equal to number of atoms')
            u2in = (self.getArrayNx3() ** 2 / masses).sum(1)
        else:
            u2in = (self.getArrayNx3() ** 2 ).sum(1)
        u2in = u2in / u2in.sum()
        coll = np.exp(-(u2in * np.log(u2in)).sum()) / self._model._n_atoms
        return coll
        
    
    def getCovariance(self):
        """Return covariance matrix calculated for this mode instance."""
        array = self.getArray()
        return np.outer(array, array) * self.getVariance()
    
    def getSqFlucts(self):
        """Return square fluctuations.
        
        Square fluctuations are obtained by multiplying the squared the mode 
        array with the variance (:meth:`getVariance`) along the mode.
        
        """
        if self.is3d():
            return (self.getArrayNx3()**2).sum(axis=1) * self.getVariance()
        else:
            return (self.getArray() ** 2)  * self.getVariance()


class Vector(VectorBase):
    
    """A class to provide operations on a modified mode array.
    
    This class holds only mode array (i.e. eigenvector) data, and has no
    associations with an NMA instance.
    
    Scalar multiplication of :class:`Mode` instance or addition of two
    :class:`Mode` instances results in a :class:`Vector` instance. 
    
    """
    
    __slots__ = ['_name', '_array', '_is3d']
    
    def __init__(self, name, array, is_3d=True):
        """Instantiate with a name, an array, and a 3d flag."""
        if not isinstance(name, str):
            raise TypeError('name must be a string')
        if not isinstance(array, np.ndarray):
            raise TypeError('array must be an ndarray')
        if not isinstance(is_3d, bool):
            raise TypeError('is_3d must be a boolean')
        self._name = name
        self._array = array
        self._is3d = is_3d
        
    def __len__(self):
        return len(self._array)
    
    def __repr__(self):
        return '<Vector: {0:s}>'.format(self._name)
    
    def __str__(self):
        return self._name 
    
    def is3d(self):
        return self._is3d
    
    def getName(self):
        """A descriptive name for the vector instance."""
        return self._name
    
    def getArray(self):
        """Normal mode array"""
        return self._array.copy()
    
    def getNormed(self):
        """Return mode after normalizing it."""
        return Vector('({0:s})/||{0:s}||'.format(self._name), 
                          self._array/(self._array**2).sum()**0.5, self._is3d)

    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom."""
        return len(self._array)

    def getNumOfAtoms(self):
        """Return number of atoms."""
        if self._is3d: 
            return len(self._array)/3
        else:
            return len(self._array)

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
        return 'Modes {0:s} from {1:s}'.format(self._slice, str(self._model))
    
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

EDA = PCA

def getVMDpath():
    """Return path to the VMD executable."""
    return VMDPATH

def setVMDpath(path):
    """Set the path to VMD executable."""
    if not path.startswith('vmd') and not os.path.isfile(path):
        LOGGER.warning('{0:s} may not be a valid path for VMD executable.')
    VMDPATH = path

def writeNMD(filename, modes, atoms):
    """Writes an NMD file for given *modes* and includes applicable data from 
    *atoms*.
    
    Returns *filename*, if file is successfully written.
    
    This function skips modes with zero eigenvalues.    
    """
    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(modes)))
    if modes.getNumOfAtoms() != atoms.getNumOfAtoms():
        raise Exception('number of atoms do not match')
    out = open(filename, 'w')
    
    #out.write('#!{0:s} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0:s}\n'.format(os.path.abspath(filename)))
    out.write('name {0:s}\n'.format(modes.getName()))
    try:
        coords = atoms.getCoordinates()
    except:
        raise RuntimeError('coordinates could not be retrived from atoms instance')
    if coords is None:
        raise RuntimeError('coordinates could not be retrived from atoms instance')
    
    try:
        data = atoms.getAtomNames()
        if data is not None:
            out.write('atomnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResidueNames()
        if data is not None:
            out.write('resnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResidueNumbers()
        if data is not None:
            out.write('resids {0:s}\n'.format(' '.join(data.astype('|S5'))))
    except:
        pass
    try:
        data = atoms.getChainIdentifiers()
        if data is not None:
            out.write('chainids {0:s}\n'.format(' '.join(data)))
    except:
        pass
    
    out.write('coordinates {0:s}\n'.format(' '.join(['{0:.3f}'.format(x) for x in coords.flatten()])))
    
    count = 0
    if isinstance(modes, Vector):
        out.write('mode {0:s}\n'.format(' '.join(['{0:.3f}'.format(x) for x in modes.getArray()])))
        count += 1
    else:
        for mode in modes:
            if mode.getEigenvalue() < ZERO:
                continue
            out.write('mode {0:d} {1:.2f} {2:s}\n'.format(mode.getIndex()+1, mode.getVariance()**0.5, ' '.join(['{0:.3f}'.format(x) for x in mode.getArray()])))
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. Given modes might have 0 eigenvalues.')
    out.close() 
    return filename  

def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""
    os.system('{0:s} -e {1:s}'.format(VMDPATH, os.path.abspath(filename)))
    
def getANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Return an ANM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`prody.atomic.AtomGroup`, :class:`prody.atomic.Selection`,  
    or :class:`prody.atomic.Chain` instance.  
    
    """
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, (AtomGroup, Selection, Chain, AtomMap)):
        ag = pdb
        if isinstance(pdb, prody.AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    anm = ANM(name + ' ANM')
    anm.buildHessian(ag.select(selstr), cutoff, gamma)
    anm.calcModes(n_modes)
    return anm

def getGNM(pdb, selstr='all', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Return an GNM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`prody.atomic.AtomGroup`, :class:`prody.atomic.Selection`,  
    or :class:`prody.atomic.Chain` instance.  
    
    """
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, (prody.AtomGroup, prody.Selection, prody.Chain, prody.AtomMap)):
        ag = pdb
        if isinstance(pdb, prody.AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    gnm = GNM(name + ' GNM')
    gnm.buildKirchhoff(ag.select(selstr), cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm

def getProjection(ensemble, modes):
    """Return projection of conformational deviations onto given modes.

    For K conformations and M modes, a (K,M) matrix is returned.     
                   
    """
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    if not modes.is3d(): 
        raise ValueError('modes must be 3-dimensional')
    if ensemble.getNumOfAtoms() != modes.getNumOfAtoms():
        raise ValueError('number of atoms are not the same')
    deviations = ensemble.getDeviations()
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0], 
                                         deviations.shape[1] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    return np.dot(deviations, modes.getArray())


def getOverlap(rows, cols):
    """Return overlap (or correlation) between two sets of modes (*rows* and *cols*).
    
    Returns a matrix whose rows correspond to modes passed as 
    *rows* argument, and columns correspond to those passed as *cols* 
    argument.
    
    """
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(rows)))
    if not isinstance(rows, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(cols)))
    
    if rows.getNumOfDegOfFreedom() != cols.getNumOfDegOfFreedom(): 
        raise ValueError('number of defrees of freedom of rows and cols must be the same')
        
    return np.dot(rows.getArray().T, cols.getArray())

def printOverlapTable(rows, cols):
    """Print table of overlaps (correlations) between two sets of modes.
    
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the printed table.

    This function may be used to take a quick look into mode correspondences 
    between two models.

    """
    print getOverlapTable(rows, cols)

def writeOverlapTable(filename, rows, cols):
    """Write table of overlaps (correlations) between two sets of modes to a file.

    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the overlap table.
    
    See also :func:`printOverlapTable`.

    """
    out = open(filename, 'w')
    out.write(getOverlapTable(rows, cols))
    out.close()
    return filename
    
def getOverlapTable(rows, cols):
    overlap = getOverlap(rows, cols)
    if isinstance(rows, Mode):
        rids = [rows.getIndex()]
        rname = str(rows.getModel())
    elif isinstance(rows, NMA): 
        rids = np.arange(len(rows))
        rname = str(rows)
    elif isinstance(rows, ModeSet): 
        rids = rows.getIndices()
        rname = str(rows.getModel())
    else:
        rids = [0]
        rname = str(rows)
    rlen = len(rids)
    if isinstance(cols, Mode):
        cids = [cols.getIndex()]
        cname = str(cols.getModel())
    elif isinstance(cols, NMA): 
        cids = np.arange(len(cols))    
        cname = str(cols)        
    elif isinstance(cols, ModeSet):
        cids = cols.getIndices()
        cname = str(cols.getModel())
    else:
        cids = [0]
        cname = str(cols)        
    clen = len(cids)
    overlap = overlap.reshape((rlen, clen)) 
    table = 'Overlap Table\n'
    table += ' '*(len(rname)+5) + cname.center(clen*7)+'\n'  
    table += ' '*(len(rname)+5)
    for j in range(clen):
        table += ('#{0}'.format(cids[j]+1)).center(7)
    table += '\n'
    for i in range(rlen):
        table += rname + (' #{0}'.format(rids[i]+1)).ljust(5)
        for j in range(clen):
            if overlap[i, j] < 0: 
                minplus = '-'
            else: 
                minplus = '+'
            table += (minplus+'{0:-.2f}').format(abs(overlap[i, j])).center(7)
        table += '\n'
    return table

def reduceModel(model, atoms, selstr):
    """Reduce dynamics model to a subset of *atoms* matching a selection *selstr*.

    Returns a tuple containing reduced model and atom selection.
    
    This function behaves depending on the type of the model.
    
    :arg model: dynamics model
    :type model: :class:`ANM`, :class:`GNM`, or :class:`PCA`
    :arg atoms: atoms that were used to build the model
    :arg selstr: a selection string specifying subset of atoms  

    For ANM and GNM:    
       This function implements [KH00]_. Selected atoms constitute the system 
       and the rest is the environment.
    
    For PCA:
       This function simply takes the sub-covariance matrix for the selected
       atoms.
       
    
    """
    if prody.la is None:
        prody.importScipyLinalg()

    #LOGGER.warning('Implementation of this function is not finalized. Use it with caution.')
    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance, not {0:s}'.format(type(model)))
    if not isinstance(atoms, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap)):
        raise TypeError('atoms type is not valid')
    if len(atoms) <= 1:
        raise TypeError('atoms must contain more than 1 atoms')

    if isinstance(model, GNM):
        matrix = model._kirchhoff
    elif isinstance(model, ANM):
        matrix = model._hessian
    elif isinstance(model, PCA):
        matrix = model._cov
    else:
        raise TypeError('model does not have a valid type derived from NMA')

    if matrix is None:
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not built')

    if isinstance(atoms, prody.AtomGroup):
        indices = np.arange(atoms.getNumOfAtoms())
    else:        
        indices = atoms.getIndices()
    
    selection = atoms.select(selstr)
    if len(selection) == 0:
        LOGGER.warning('selection has 0 atoms')
        return None
    sel = selection.getIndices()
    if len(indices) == len(sel):
        LOGGER.warning('selection results in same number of atoms, model is not reduced')
        return None
    ndim = 1
    if model._is3d:
        ndim = 3
    system = [] 
    other = []
    index = 0
    for i in indices:
        if i in sel:
            system.extend(range(index*ndim, (index+1)*ndim))
        else:
            other.extend(range(index*ndim, (index+1)*ndim))
        index += 1
    ss = matrix[system,:][:,system]
    if isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(ss)
        return eda, selection
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(prody.la.inv(oo), os))
    
    if isinstance(model, GNM):
        gnm = GNM(model.getName()+' reduced')
        gnm.setKirchhoff(matrix)
        return gnm, selection
    elif isinstance(model, ANM):
        anm = ANM(model.getName()+' reduced')
        anm.setHessian(matrix)
        return anm, selection
    elif isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(matrix)
        return eda, selection

def writeModes(filename, modes, format='g', sep=' ', compressed=False):
    """Write *modes* (eigenvectors) into a plain text file with name *filename*.
    
    See also :func:`writeArray`.
    
    """
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    return writeArray(filename, modes.getArray(), format=format, sep=sep)
    
def writeArray(filename, array, format='g', sep=' ', compressed=False):
    """Write 1-d or 2-d array data into a delimited text file.
    
    If *filename* is ``None``, output will be returned as a formatted string. 
    Default *format* argument is "g", which writes enough many significant 
    digits automatically.  Default seperator (*sep*) argument is white 
    space " ". Optionally, a gzip *compressed* file may be outputted.
    
    If array is written into a file, *filename* will be returned upon
    successful writing. 
     
    """
    if not array.ndim in (1, 2):
        raise ValueError('array must be 1- or 2-dimensional')
    if array.ndim == 2:
        length = array.shape[1]
        line = '{' + '0[{0:d}]:{1:s}'.format(0, format) + '}'
    else:
        length = 1
        line = '{' + '0:{0:s}'.format(format) + '}'
    for j in range(1, length):
        line += sep + '{' + '0[{0:d}]:{1:s}'.format(j, format) + '}' 
    line += '\n'
    output = ''
    for row in array:
        output += line.format(row)
    if compressed:
        out = gzip.open(filename + '.gz', 'w')
    else:
        out = open(filename, 'w')
    out.write(output)
    out.close()
    return filename

def sampleModes(modes, atoms=None, **kwargs):
    """Samples conformations along given *modes* and returns an ensemble.
    
    .. versionadded:: 0.2
    
    If *atoms* are provided, sampling is around the active coordinate set of 
    *atoms*. Otherwise, sampling is around the 0 coordinates, Conformations
    are deformations.
    
    :arg modes: Modes along which sampling will be performed.
    :type modes: :class:`Mode`, :class:`ModeSet`, :class:`PCA`, or :class:`ANM`  
    
    :arg atoms: Atoms whose active coordinate will be used as the initial 
        conformation.
    :type atoms: :class:`prody.atomic.AtomGroup`, :class:`prody.atomic.Chain`, 
        :class:`prody.atomic.Selection`  
    
    :keyword n_steps: Number of steps to take along each direction of each
        mode. For example, for ``n_steps=10``, 20 conformations will be 
        generated along the first mode. Default is 10.
    :type n_steps: int 
    
    :keyword rmsd_step: RMSD change along the slowest mode at each step.
        Default is 0.2 A. 
    :type rmsd_step: float    
    
    :keyword rmsd_max: The maximum RMSD that the conformations will have with 
        respect to the initial conformation.  
    :type rmsd_max: float 

    :returns: :class:`prody.ensemble.Ensemble`
    
    If *rmsd_max* and *n_steps* are provided, *rmsd_step* is 
    calculated based on them. If *n_steps* and *rmsd_step* are given,
    *rmsd_max* is determined based on them.
    
    Note that *rmsd_step* is for the slowest mode. Effective *rmsd_step* is 
    smaller for faster modes. It is scaled down by the ratio of 
    the inverse eigenvalue of the faster mode to that of the slowest mode.
    """

    if not isinstance(modes, (Mode, PCA, ANM, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be a from a 3-dimentional model.')
    n_atoms = modes.getNumOfAtoms()
    if atoms is not None:
        if not isinstance(atoms, (AtomGroup, Selection, Chain)):
            raise TypeError('{0:s} is not correct type for atoms'.format(type(modes)))

        if atoms.getNumOfAtoms() != n_atoms:
            raise ValueError('Number of atoms must match.')
        confs = [atoms.getCoordinates()]
    else:    
        confs = [np.zeros((n_atoms, 3))]
    if isinstance(modes, Mode):
        modes = [modes]
    elif isinstance(modes, ModeSet):
        modes = list(modes)
    if 'rmsd_max' in kwargs:
        rmsd_max = float(kwargs['rmsd_max']) + 0.000004
        LOGGER.info('Parameter: rmsd_max={0:.2f} A'.format(rmsd_max))
        if 'rmsd_step' in kwargs:
            rmsd_step = float(kwargs['rmsd_step'])
            n_steps = int(round(rmsd_max / rmsd_step))
        else:
            n_steps = int(kwargs.get('n_steps', 10))
            rmsd_step = rmsd_max / n_steps
    else:
        rmsd_max = None
        n_steps = int(kwargs.get('n_steps', 10))
        rmsd_step = float(kwargs.get('rmsd_step', 0.2))
    LOGGER.info('Parameter: n_steps={0:d}'.format(n_steps))
    LOGGER.info('Parameter: rmsd_step={0:.2f} A'.format(rmsd_step))
    
    m = modes[0]
    arr = m.getArrayNx3()
    var = m.getVariance()
    scale = ((n_atoms * rmsd_step**2) / var) **0.5
    LOGGER.info('Modes are scaled by s={0:g}'.format(scale))
    
    arrays = []
    drmsds = []
    for m in modes:
        arr = m.getArrayNx3()
        var = m.getVariance()
        arrays.append( arr * scale * var**0.5 )
        if rmsd_max is not None:         
            drmsds.append(((var * scale**2) / n_atoms)**0.5)
    extend = confs.extend
    rmsds = [0]
    for mdarr, drmsd in zip(arrays, drmsds):
        newconfs = []
        append = newconfs.append 
        for r, c in enumerate(confs):
            for s in range(-n_steps,0)+range(1, n_steps+1):
                if rmsd_max is not None:
                    rmsd = rmsds[r] + abs(s)*drmsd  
                    if rmsd > rmsd_max:
                        continue
                append(c +  mdarr * s)
                rmsds.append(rmsd)
        extend(newconfs)
    ensemble = Ensemble('Conformations along {0:s}'.format(modes))
    ensemble.setCoordinates(confs[0])    
    ensemble.addCoordset(np.array(confs[1:]))
    return ensemble

def getSqFlucts(modes):
    """Return sum of square-fluctuations for given set of normal *modes*."""
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Mode):
        modes = [modes]
    square_fluctuations = np.zeros(modes.getNumOfAtoms()) 
    for mode in modes.modes:
        square_fluctuations += mode.getSqFlucts()
    return square_fluctuations
 
def getCrossCorrelations(modes, n_cpu=1):
    """Return cross-correlations matrix.
    
    For a 3-d model, cross-correlations matrix is an NxN matrix, where N is the 
    number of atoms. Each element of this matrix is the trace of the 
    submatrix corresponding to a pair of atoms.
    
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.

    For large systems, calculation of cross-correlations matrix may be time 
    consuming. Optionally, multiple processors may be employed to perform
    calculations by passing ``n_cpu=2`` or more. 

    :arg n_cpu: number of CPUs to use 
    :type n_cpu: int, default is 1
    
    """
    if not isinstance(n_cpu, int):
        raise NMAError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise NMAError('n_cpu must be equal to or greater than 1')
        
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
        
    if modes.is3d():
        model = modes
        if isinstance(modes, (Mode, ModeSet)):
            model = modes._model
            if isinstance(modes, (Mode)):
                indices = [modes.getIndex()]
                n_modes = 1
            else:
                indices = modes.getIndices()
                n_modes = len(modes)
        else:
            n_modes = len(modes)
            indices = np.arange(n_modes)
        array = model._array
        n_atoms = model._n_atoms
        variances = model._vars
        if n_cpu == 1:
            arvar = (array[:, indices]*variances[indices]).T.reshape((n_modes,
                                                                   n_atoms, 3))
            array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
            covariance = np.tensordot(array.transpose(2, 0, 1),
                                      arvar.transpose(0, 2, 1),
                                      axes=([0, 1], [1, 0]))
        else:
            import multiprocessing
            n_cpu = min(multiprocessing.cpu_count(), n_cpu)
            queue = multiprocessing.Queue()
            size = n_modes / n_cpu
            for i in range(n_cpu):
                if n_cpu - i == 1:
                    indices = modes.indices[i*size:]
                else:
                    indices = modes.indices[i*size:(i+1)*size]
                process = multiprocessing.Process(target=_cross_correlations, 
                              args=(queue, n_atoms, array, variances, indices))
                process.start()
            while queue.qsize() < n_cpu:
                time.sleep(0.05)
            covariance = queue.get()
            while queue.qsize() > 0:
                covariance += queue.get()
    else:
        covariance = modes.getCovariance()
    diag = np.power(covariance.diagonal(), 0.5)
    return covariance / np.outer(diag, diag)

def _crossCorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""
    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)


def getCumulativeOverlap(modes1, modes2):
    """Return cumulative overlap of modes in *modes2* with those in *modes1*.
    
    Returns a number of *modes1* contains a single :class:`Mode` or a 
    :class:`Vector` instance. If *modes1* contains multiple modes, returns an
    array. Elements of the array correspond to cumulative overlaps for modes 
    in *modes1* with those in *modes2*."""
    overlap = getOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    return cumov

def getCumulativeOverlapArray(modes1, modes2):
    """Return array of cumulative overlaps.
   
    .. versionadded: 0.2
    
    Returned array has the shape ``(len(modes1), len(modes2))``. Each row
    corresponds to cumulative overlaps calculated for modes in *modes1* with
    those in *modes2*. Each value in a row corresponds to cumulative overlap
    calculated using upto that many number of modes from *modes2*.
    """
    overlap = getOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))
    return cumov


def getSubspaceOverlap(modes1, modes2):
    """Return subspace overlap between two sets of modes (*modes1* and *modes2*).
    
    Also known as the root mean square inner product (RMSIP) of essential 
    subspaces [AA99]_.
    
    This function returns a single number.
        
    """
    modes1 = get_dict(modes1)
    modes2 = get_dict(modes2)
    overlap = getOverlap(modes1, modes2)
    rmsip = np.sqrt(np.power(overlap, 2).sum() /
                               len(adict['modes']))
    return rmsip

def getCovariance(modes):
    """Calculate covariance matrix from given modes and return it."""
    return modes.getCovariance()

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
    if not isinstance(modes, (Mode, NMA, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Mode):
        modes = [modes]
    fracts = np.array([mode.getFractOfVariance() for mode in modes]).cumsum()
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
    
    :arg ensemble: a :class:`prody.ensemble.Ensemble` instance
    
    Matplotlib function used for plotting depends on the number of modes:
        
      * 1 mode: :func:`matplotlib.pyplot.hist`
      * 2 modes: :func:`matplotlib.pyplot.scatter`
      * 3 modes: :meth:`mpl_toolkits.mplot3d.Axes3D.scatter`
    """
    if pl is None: prody.importPyPlot()
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    if isinstance(modes, Mode):
        projection = getProjection(ensemble, modes)
        show = pl.hist(projection.flatten(), *args, **kwargs)
        pl.xlabel('Mode {0:d} coordinate'.format(modes.getIndex()+1))
        pl.ylabel('Number of conformations')
    elif len(modes) == 2:
        projection = getProjection(ensemble, modes)
        show = pl.scatter(projection[:,0], projection[:,1], *args, **kwargs)
        modes = [m for m in modes]
        pl.xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        pl.ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
    elif len(modes) == 3:
        projection = getProjection(ensemble, modes)
        modes = [m for m in modes]
        show = Axes3D(pl.gcf())
        show.scatter(projection[:, 0], projection[:, 1], projection[:, 2], *args, **kwargs)
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    else:
        raise ValueError('Projection onto upto 3 modes can be shown. You have given {0:d} mode.'.format(len(modes)))
    
    return show
    
def showOverlapTable(rows, cols, *args, **kwargs):
    """Show overlap table using :func:`matplotlib.pyplot.pcolor`.
    
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

showOverlapMatrix = showOverlapTable

def showCrossCorrelations(modes, *args, **kwargs):
    """Show cross-correlations for given modes using :func:`matplotlib.pyplot.imshow`.
    
    See also :func:`getCrossCorrelations`. 
    
    By default, *origin=lower* and *interpolation=bilinear* keyword
    arguments are passed to imshow function. User can overwrite these
    parameters.
    
    """
    if pl is None: prody.importPyPlot()
    arange = np.arange(modes.getNumOfAtoms())
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:, 
                       arange[0]+1:] = getCrossCorrelations(modes)
    if not kwargs.has_key('interpolation'):
        kwargs['interpolation'] = 'bilinear'
    if not kwargs.has_key('origin'):
        kwargs['origin'] = 'lower'
    show = pl.imshow(cross_correlations, *args, **kwargs), pl.colorbar()
    pl.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    pl.title('Cross-correlations for {0:s}'
             .format(str(modes))) 
    pl.xlabel('Indices')
    pl.ylabel('Indices')
    return show

def showMode(mode, *args, **kwargs):
    """Show mode array using :func:`matplotlib.pyplot.plot`."""
    if pl is None: prody.importPyPlot()
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode, not {0:s}'.format(type(modes)))
    if mode.is3d():
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
    sqf = getSqFlucts(modes)
    show = pl.plot(sqf, *args, **kwargs)
    pl.xlabel('Indices')
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
    show = pl.spy(kirchhoff, *args, **kwargs)
    pl.title('{0:s} contact map'.format(enm.getName())) 
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
        raise TypeError('mode must be Mode or Vector, not {0:s}'.format(type(mode)))
    if not isinstance(modes, (NMA, ModeSet)):
        raise TypeError('modes must be NMA or ModeSet, not {0:s}'.format(type(modes)))
    overlap = abs(getOverlap(mode, modes))
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+0.5)
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
    cumov = (getOverlap(mode, modes) ** 2).cumsum() ** 0.5
    if isinstance(modes, NMA):
        arange = np.arange(0.5, len(modes)+1.5)
    else:
        arange = modes.getIndices() + 0.5
    show = pl.plot(arange, cumov, *args, **kwargs)
    pl.title('Cumulative overlap: {0:s} & {1:s}'.format(str(mode), str(modes)))
    pl.xlabel('Mode index')
    pl.ylabel('Cumulative overlap')
    pl.axis((arange[0]-0.5, arange[-1]+0.5, 0, 1))
    return show
    
