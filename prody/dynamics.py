# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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
  * :class:`EDA`
  * :class:`NMA`
  * :class:`Mode`
  * :class:`ModeSet`
  * :class:`Vector`
  * :class:`GammaStructureBased`
  * :class:`GammaVariableCutoff`
  
Base Classes
------------

  * :class:`NMABase`
  * :class:`GNMBase`
  * :class:`VectorBase`
  * :class:`Gamma` 

Inheritance Diagram
-------------------

.. inheritance-diagram:: prody.dynamics
   :parts: 1

Functions
---------

Many of the functions documented in this page accepts a *modes* argument (may 
also appear in different names). One of the following may be accepted as this 
argument:

  * an NMA model, which may be an instance of one of :class:`ANM`, 
    :class:`GNM`, :class:`NMA`, :class:`PCA`.
  * a :class:`Mode` instance obtained by indexing an NMA model, e.g. ``nma[0]``
  * a :class:`ModeSet` instance obtained by slicing an NMA model, e.g. ``nma[:10]``

Some of these functions may also accept :class:`Vector` instances as *mode* 
argument. These are noted in function documentations. 

**Short-hand functions**:

  * :func:`calcANM`
  * :func:`calcGNM`

**Analysis**:

  * :func:`calcCollectivity`
  * :func:`calcCovariance`
  * :func:`calcCrossCorrelations`
  * :func:`calcSqFlucts`  
  * :func:`calcProjection`
  * :func:`calcTempFactors`

**Parse/write data**:

  * :func:`parseArray`
  * :func:`parseModes`
  * :func:`parseSparseMatrix`
  * :func:`parseNMD`
  * :func:`writeArray`
  * :func:`writeModes`
  * :func:`writeNMD`
  * :func:`writeOverlapTable`
        
**Save/load models**:
    
  * :func:`saveModel`
  * :func:`loadModel`
  * :func:`saveVector`
  * :func:`loadVector`

**Visualization**:

  * :func:`getVMDpath`
  * :func:`setVMDpath`
  * :func:`viewNMDinVMD`
    
**Comparative analysis**:

  * :func:`calcOverlap`
  * :func:`calcCumulativeOverlap`
  * :func:`calcCumulativeOverlapArray`
  * :func:`calcSubspaceOverlap`
  * :func:`calcCovarianceOverlap`
  * :func:`printOverlapTable`
  
**Sampling**:

  * :func:`deform`
  * :func:`sampleModes`
  * :func:`traverseMode`

**Model extrapolation/reduction**:
  
  * :func:`extrapolateModel`  
  * :func:`reduceModel`
  * :func:`sliceVector`
  * :func:`sliceMode`
  * :func:`sliceModel`

**Plotting Functions**:

Plotting functions are called by the name of the plotted data/property and are
prefixed with "show". Function documentations include the :mod:`matplotlib.pyplot` 
function utilized for actual plotting. Arguments and keyword arguments are passed 
to the Matplotlib functions.  


  * :func:`showContactMap`
  * :func:`showCrossCorrelations`
  * :func:`showCumulativeOverlap`
  * :func:`showCumFractOfVariances`
  * :func:`showFractOfVariances`
  * :func:`showMode`
  * :func:`showOverlap`
  * :func:`showOverlapTable`
  * :func:`showProjection`
  * :func:`showCrossProjection`
  * :func:`showEllipsoid`
  * :func:`showSqFlucts`
  * :func:`showScaledSqFlucts`
  * :func:`showNormedSqFlucts`
  * :func:`resetTicks`
    

Examples
--------

Results from the example :ref:`pca-xray-calculations` will be used to
illustrate class methods and functions in the module.
  

>>> from prody import *
>>> import matplotlib.pyplot as plt
>>> import numpy as np

>>> p38_pca = loadModel('p38_xray.pca.npz')
>>> p38_anm = loadModel('1p38.anm.npz') 
>>> p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
>>> p38_structure = parsePDB('p38_ref_chain.pdb')


.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np
   
   plt.close('all')
   
   p38_pca = loadModel('p38_xray.pca.npz')
   p38_anm = loadModel('1p38.anm.npz') 
   p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
   p38_structure = parsePDB('p38_ref_chain.pdb')
   
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import os.path
import time
import os
from types import FunctionType

import numpy as np
linalg = None
scipyla = None
plt = None
KDTree = None
scipy_sparse = None
scipy_sparse_la = None

import prody
from .atomic import *
from .ensemble import *
from prody import ProDyLogger as LOGGER
from prody import ProDyAtomSelect as SELECT
from prody import ProDyException

__all__ = ['ANM', 'GNM', 'NMA', 'PCA', 'EDA', 'Mode', 'ModeSet', 'Vector', 
           
           'NMABase', 'GNMBase', 'VectorBase', 
           
           'Gamma', 'GammaStructureBased', 'GammaVariableCutoff',
           
           'calcANM', 'calcGNM', 
           
           'calcCollectivity', 'calcCovariance', 'calcCrossCorrelations', 
           
           'calcSqFlucts', 'calcTempFactors',
           
           'calcProjection',  
           
           'parseArray', 'parseModes', 'parseNMD',
           
           'writeArray', 'writeModes', 'parseSparseMatrix',
           'writeNMD', 'writeOverlapTable',
           
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD', 
           
           'calcOverlap', 'calcCumulativeOverlap', 
           'calcCumulativeOverlapArray', 'calcSubspaceOverlap', 
           'calcCovarianceOverlap', 'printOverlapTable',
           
           'deform', 'sampleModes', 'traverseMode',
            
           'extrapolateModel', 'reduceModel', 'sliceVector', 
           'sliceMode', 'sliceModel',
            
           'showContactMap', 'showCrossCorrelations', 'showCumulativeOverlap', 
           'showCumFractOfVariances', 'showFractOfVariances', 'showMode', 
           'showOverlap', 'showOverlapTable', 'showProjection', 
           'showCrossProjection', 'showEllipsoid', 'showSqFlucts', 
           'showScaledSqFlucts', 'showNormedSqFlucts', 'resetTicks'
           ]

ZERO = 1e-8

class VectorBase(object):
    """A base class for :class:`Mode` and :class:`Vector`.
    
    This base class defines some shared methods, such as scalar multiplication
    or addition of mode instances.
    
    Defined operations are:
        
        * Absolute value (abs(mode)) returns mode length
        * Additive inverse (-mode) 
        * Mode addition (mode1 + mode2)
        * Mode subtraction (mode1 - mode2)
        * Scalar multiplication (x*mode or mode*x)
        * Division by a scalar (mode/x)
        * Dot product (mode1*mode2)
        * Power (mode**x)
    
    """
    
    def __abs__(self):
        return np.sqrt((self.getArray()**2).sum())
    
    def __neg__(self):
        return Vector(-self.getArray(), '-({0:s})'.format(str(self)), self.is3d())
    
    def __div__(self, other):
        if isinstance(other, (int, float, long)):
            return Vector(self.getArray() / other, 
                          '({1:s})/{0}'.format(other, str(self)), self.is3d())
        else:
            raise TypeError('{0} is not a scalar'.format(other))
    
    def __idiv__(self, other):
        return self.__div__(other)
    
    def __mul__(self, other):
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector(other * self.getArray(), 
                          '{0}*({1:s})'.format(other, str(self)), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self.getArray(), other.getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
    
    def __rmul__(self, other):   
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector(other * self.getArray(), 
                          '{0}*({1:s})'.format(other, str(self)), self.is3d())
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
            return Vector(self.getArray() + other.getArray(), 
                          '{0:s} + {1:s}'.format(str(self), str(other)), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __radd__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self.getArray() + other.getArray(), 
                          '{0:s} + {1:s}'.format(str(other), str(self)), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))
        
    def __iadd__(self, other):
        return self.__add__(other)   

    def __sub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self.getArray() - other.getArray(), 
                          '{0:s} - {1:s}'.format(str(self), str(other)), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __rsub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return  Vector(other.getArray() - self.getArray(), 
                           '{0:s} - {1:s}'.format(str(other), str(self)), self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __isub__(self, other):
        return self.__sub__(other)   

    def __pow__(self, other):
        if isinstance(other, (int, float, long)): 
            return Vector(self.getArray() ** other, 
                          '({0:s})**{1}'.format(str(self), other), self.is3d())
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
        
        >>> mode = p38_anm[0]
        >>> mode
        <Mode: 1 from ANM 1p38>
        
        """
        
        self._model = model
        self._index = int(index)
        
    def __len__(self):
        return self._model._dof
    
    def __repr__(self):
        return '<Mode: {0:d} from {1:s}>'.format(self._index+1, str(self._model))

    def __str__(self):
        return 'Mode {0:d} from {1:s}'.format(self._index+1, str(self._model))

    def is3d(self):
        """Returns True if mode instance is from a 3-dimensional model.
        
        >>> print mode.is3d()
        True
        
        """
        
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Return number of atoms.
        
        >>> mode.getNumOfAtoms()
        321
        
        """
        
        return self._model._n_atoms
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom.
        
        >>> mode.getNumOfDegOfFreedom()
        963
        
        """
        
        return self._model._dof
    
    def getName(self):
        """A descriptive name for the mode instance.
        
        >>> mode.getName()
        'Mode 1 from ANM 1p38'
        
        """
        
        return str(self)
    
    def getIndex(self):
        """Return index of the mode.
        
        >>> mode.getIndex() # Note that this returns 0, i.e. anm[0]
        0
        
        """
        
        return self._index
    
    def getModel(self):
        """Return the model that the mode belongs to.
        
        >>> mode.getModel()
        <ANM: 1p38 (20 modes, 321 nodes)>
        
        """
        
        return self._model
    
    def getArray(self):
        """Return the normal mode array (eigenvector).
        
        >>> mode.getArray() # doctest: +SKIP
        array([ -2.17154375e-02,   5.19952720e-02,  -5.45241061e-02,
                -2.10540040e-02,   5.74771315e-02,  -5.74002951e-02,
                -3.27903774e-02,   7.57414975e-02,  -5.86200994e-02,
                 ...,
                 9.90919829e-03,  -8.97774222e-03,  -3.05413370e-02,
                 9.46978952e-03,  -8.34108970e-03,  -2.94414502e-02,
                 1.14453479e-02,  -1.91394533e-02,  -2.43839818e-02])        
        """
        
        return self._model._array[:, self._index].copy()
    
    getEigenvector = getArray
    
    def getEigenvalue(self):
        """Return normal mode eigenvalue.
        
        >>> mode.getEigenvalue() # doctest: +SKIP
        0.10652815070679021
        
        """
        
        return self._model._eigvals[self._index]
    
    def getVariance(self):
        """Variance along the mode. 
        
        If the model is not a PCA, inverse of the eigenvalue is returned.
        
        >>> mode.getVariance() # doctest: +SKIP
        9.3871900841723619
        >>> mode.getEigenvalue()**-1 # doctest: +SKIP
        9.3871900841723619
        
        """
        
        return self._model._vars[self._index]

    def getFractOfVariance(self):
        """Return fraction of variance explained by the mode.
        
        Fraction of variance is the ratio of the variance along this mode to 
        the trace of the covariance matrix.
        
        See :meth:`getVariance`
        
        >>> mode.getFractOfVariance() # doctest: +SKIP
        0.32090162902353198
        
        """
        
        return self.getVariance() / self._model._trace
    
    def getCollectivity(self, masses=None):
        """Return collectivity of the mode.
        
        This function implements collectivity as defined in equation 5 of 
        [BR95]_.
        
        To incorporate atomic *masses* in calculations, use 
        :func:`calcCollectivity` function.
        
        >>> mode.getCollectivity() # doctest: +SKIP
        0.64762987283788642
        
        """

        return calcCollectivity(self)

    def getCovariance(self):
        """Return covariance matrix calculated for this mode instance.
        
        >>> mode.getCovariance() # doctest: +SKIP
        array([[ 0.00442663, -0.01059908,  0.01111457, ..., -0.0023331 ,
                 0.00390152,  0.0049706 ],
               [-0.01059908,  0.02537835, -0.02661264, ...,  0.00558635,
                -0.00934177, -0.01190157],
               [ 0.01111457, -0.02661264,  0.02790697, ..., -0.00585805,
                 0.00979611,  0.01248041],
               ..., 
               [-0.0023331 ,  0.00558635, -0.00585805, ...,  0.00122968,
                -0.00205634, -0.00261981],
               [ 0.00390152, -0.00934177,  0.00979611, ..., -0.00205634,
                 0.0034387 ,  0.00438096],
               [ 0.0049706 , -0.01190157,  0.01248041, ..., -0.00261981,
                 0.00438096,  0.00558142]])
                 
        """
        
        array = self.getArray()
        return np.outer(array, array) * self.getVariance()
    
    def getSqFlucts(self):
        """Return square fluctuations.
        
        Square fluctuations are obtained by multiplying the squared the mode 
        array with the variance (:meth:`getVariance`) along the mode.
        
        >>> mode.getSqFlucts() # doctest: +SKIP
        array([  5.77119441e-02,   6.61016413e-02,   9.62027339e-02,
                 8.11873159e-02,   1.02183864e-01,   1.12300137e-01,
                 ...
                 6.55654114e-03,   6.14363088e-03,   3.52230202e-03,
                 1.04344751e-02,   9.63172338e-03,   1.02498093e-02])                 

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
    
    def __init__(self, array, name='Unknown', is_3d=True):
        """Instantiate with a name, an array, and a 3d flag."""
        
        if not isinstance(array, np.ndarray) or array.ndim != 1:
            raise TypeError('array must be a 1-dimensional numpy.ndarray')
        if not isinstance(is_3d, bool):
            raise TypeError('is_3d must be a boolean')
        self._name = str(name)
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
        """Get the descriptive name for the vector instance."""
        
        return self._name
    
    def setName(self, name):
        """Set the descriptive name for the vector instance."""
        
        self._name = str(name) 
    
    def getArray(self):
        """Normal mode array"""
        
        return self._array.copy()
    
    def getNormed(self):
        """Return mode after normalizing it."""
        
        return Vector(self._array/(self._array**2).sum()**0.5, 
                      '({0:s})/||{0:s}||'.format(self._name), self._is3d)

    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom."""
        return len(self._array)

    def getNumOfAtoms(self):
        """Return number of atoms.
        
        For a 3-dimensional vector, returns length of the vector divided by 3.
        
        """
        
        if self._is3d: 
            return len(self._array)/3
        else:
            return len(self._array)
    

class NMABase(object):
    
    """Base class for Normal Mode Analysis calculations.
    
    Derived classes are:
        
        * :class:`GNMBase`
        * :class:`NMA`
        * :class:`PCA`
        
    """
    
    def __init__(self, name):
        """Initialize a Normal Mode analysis with a given name.
        
        .. versionchanged:: 0.7
           When an empty string is passed as *name* argument, NMA instance 
           is called "Unnamed".
        
        """
        
        name = str(name)
        if name == '':
            name = 'Unnamed'
        self._name = name 
        
        self._modes = []
        self._n_modes = 0
        self._cov = None
        self._n_atoms = 0
        self._dof = 0
        self._array = None      # modes/eigenvectors
        self._eigvals = None
        self._vars = None       # evs for PCA, inverse evs for ENM
        self._trace = None
        self._is3d = True       # is set to false for GNM

    def __len__(self):
        return self._n_modes
        
    def __getitem__(self, index):
        """
        .. versionchanged:: 0.6
           A list or tuple of integers can be used for indexing. 
        
        """
        
        if self._n_modes == 0:
            raise ProDyException('{0:s} modes are not calculated, try '
                                 'calcModes() method'.format(str(self)))
        if isinstance(index, int):
            return self.getMode(index)
        elif isinstance(index, slice):
            #modes = self._modes[index]
            indices = np.arange(*index.indices(len(self)))
            if len(indices) > 0:
                return ModeSet(self, indices)
        elif isinstance(index, (list, tuple)):
            for i in index:
                assert isinstance(i, int), 'all indices must be integers'
            if len(index) == 1:
                return self.getMode(index[0])
            return ModeSet(self, index)
        else:        
            raise IndexError('indices may be an integer, slice, list, or tuple')
        
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
    
    def getModel(self):
        """Return self."""
        
        return self
    
    def is3d(self):
        """Return ``True`` is model is 3-dimensional."""
        
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
        if self._n_modes == 0:
            raise ProDyException('{0:s} modes are not calculated, try '
                                 'calcModes() method'.format(str(self)))
        if index >= self._n_modes or index < -self._n_modes:
            raise IndexError('{0:s} contains {1:d} modes, try 0 <= index < '
                             '{1:d}'.format(str(self), self._n_modes))
        if index < 0:
            index += self._n_modes
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
        
        if self._eigvals is None: return None
        return self._eigvals.copy()

    def getEigenvectors(self):
        """Return eigenvectors."""
        
        if self._array is None: return None
        return self._array.copy()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues)."""
        
        if self._vars is None: return None
        return self._vars.copy()

    def getArray(self):
        """Return eigenvectors."""
        
        if self._array is None: return None
        return self._array.copy()
        
    def getCovariance(self):
        """Return covariance matrix.
        
        If covariance matrix is not set or calculated yet, it will be 
        calculated using available modes and then returned.
        
        """
        
        if self._cov is None:
            if self._array is None:
                return None
            self._cov = np.dot(self._array, np.dot(np.diag(self._vars), 
                                                   self._array.T))
        return self._cov
        
    def calcModes(self):
        pass
    
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to the :class:`NMA` instance.
        
        .. versionadded:: 0.5.3
        
        If *eigenvalue* is not given, it will be set to 1. 
        Variances are set as the inverse eigenvalues.
        
        """
        vector = eigenvector
        value = eigenvalue
        if not isinstance(vector, np.ndarray):
            raise TypeError('eigenvector must be a Numpy ndarray, not {0:s}'
                            .format(type(vector)))
        if vector.ndim > 2:
            raise ValueError('eigenvector must be 1- or 2-dimensional array, '
                            'not {0:d}-d'.format(type(vector.ndim)))
        elif vector.ndim == 2 and vector.shape[0] >= vector.shape[1]:
            raise ValueError('eigenvectors must correspond to columns')
        else:
            vector = vector.reshape((vector.shape[0], 1))
        
        if value is None:
            if vector.ndim == 1:
                value = np.ones(1)
            else:
                value = np.ones(vector.shape[2])
        elif isinstance(value, (int, float)):
            value = np.array([value])
        elif isinstance(value, np.ndarray):
            if value.ndim > 1:
                raise ValueError('eigenvalue must be a 1-d array')
            elif value.shape[0] != vector.shape[0]:
                raise ValueError('number of eigenvectors and eigenvalues '
                                 'must match')
        
        if self._array is None:
            self._array = vector
            self._eigvals = value
            self._dof = vector.shape[0]
            if self._is3d:
                self._n_atoms = self._dof / 3
            else:
                self._n_atoms = self._dof
            self._n_modes = vector.shape[1]
            self._modes = [None] * self._n_modes
        else:
            if vector.shape[0] != self._array.shape[0]: 
                raise ValueError('shape of vector do not match shape of ' 
                                 'existing vectors')
            self._array = np.concatenate((self._array, vector), 1)
            self._eigvals = np.concatenate((self._eigvals, value))
            self._n_modes += vector.shape[1]            
            self._modes += [None] * vector.shape[1]
        
        self._vars = 1 / self._eigvals
    
    def setEigens(self, vectors, values=None):
        """Set eigenvectors and eigenvalues.
        
        .. versionadded:: 0.5.3
        
        :arg vectors: eigenvectors
        :type vectors: numpy.ndarray
        
        :arg values: Eigenvalues. When ``None`` is passed (default value), 
            all eigenvalues will be set to ``1``.
        :type values: numpy.ndarray
        
        For M modes and N atoms, *vectors* must have shape ``(3*N, M)``
        and values must have shape ``(M,)``.
        
        Variances are set as the inverse eigenvalues.
        
        """
        
        if not isinstance(vectors, np.ndarray):
            raise TypeError('vectors must be a numpy.ndarray, not {0:s}'
                            .format(type(vectors)))
        elif vectors.ndim != 2:
            raise ValueError('vectors must be a 2-dimensional array')
        else:
            dof = vectors.shape[0]
            if self._is3d:
                n_atoms = dof / 3
            else: 
                n_atoms = dof
            if self._n_atoms > 0 and n_atoms != self._n_atoms:
                    raise ValueError('vectors do not have the right shape, '
                                 'which is (M,{0:d})'.format(n_atoms*3))
            n_modes = vectors.shape[1]
        if values is not None:
            if not isinstance(vectors, np.ndarray):
                raise TypeError('values must be a numpy.ndarray, not {0:s}'
                                .format(type(vectors)))
            elif values.ndim != 1:
                raise ValueError('values must be a 1-dimensional array')
            else:
                if values.shape[0] != vectors.shape[1]:
                    raise ValueError('number of vectors and values do not match')
        else:
            values = np.ones(n_modes)
            
        self._array = vectors
        self._eigvals = values
        self._dof = dof
        self._n_atoms = n_atoms
        self._n_modes = n_modes
        self._modes = [None] * n_modes
        self._vars = 1 / values


class NMA(NMABase):
    
    """A class for analysis of externally calculated Hessian matrices and 
    normal modes.
    
    """
    
    def __init__(self, name):
        NMABase.__init__(self, name)

        
class ModeSet(object):
    """A class for providing access to data for a subset of modes.
    
    Instances are obtained by slicing an NMA model (:class:`ANM`, :class:`GNM`, 
    or :class:`PCA`).
    
    ModeSet's contain a reference to the model and a list of mode indices.
    Methods common to NMA models are also defined for mode sets.
    
    >>> modes = p38_anm[:3]
    >>> modes
    <ModeSet: 3 modes from ANM 1p38>
    
    """
    
    __slots__ = ['_model', '_indices']
    
    def __init__(self, model, indices):
        if not isinstance(model, NMABase):
            raise TypeError('model must be an NMA, not {0:s}'
                            .format(type(model)))
        self._model = model
        self._indices = np.array(indices, np.int64)
        
    def __len__(self):
        return len(self._indices)
        
    def __iter__(self):
        for i in self._indices:
            yield self._model.getMode(i)
    
    def __repr__(self):
        return '<ModeSet: {0:d} modes from {1:s}>'.format(len(self),
                                                       str(self._model))

    def __str__(self):
        return '{0:d} modes from {1:s}'.format(len(self._indices), 
                                               str(self._model))
    
    def is3d(self):
        """Return True if mode instance is from a 3-dimensional model.
        
        >>> modes.is3d()
        True
        
        """
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Return number of nodes.
        
        >>> modes.getNumOfAtoms()
        321
        
        """
        return self._model._n_atoms
    
    def getNumOfModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes).
        
        >>> modes.getNumOfModes()
        3
        
        """
        
        return len(self._indices)
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom.
        
        >>> modes.getNumOfDegOfFreedom()
        963
        
        """
        
        return self._model._dof
        
    def getModes(self):
        """Return all modes in the subset in a list.
        
        >>> modes.getModes() # doctest: +SKIP
        [<Mode: 1 from ANM 1p38>, <Mode: 2 from ANM 1p38>, 
        <Mode: 3 from ANM 1p38>]
        
        """
        
        getMode = self._model.getMode
        return [getMode(i) for i in self._indices]
    
    def getName(self):
        """Return name.
        
        >>> modes.getName()
        '3 modes from ANM 1p38'
        
        """
        
        return str(self)
    
    def getModel(self):
        """Return the model that the modes belongs to.
        
        >>> modes.getModel()
        <ANM: 1p38 (20 modes, 321 nodes)>
        
        """
        
        return self._model
    
    def getIndices(self):
        """Return indices of modes in the mode set.
        
        >>> modes.getIndices()
        array([0, 1, 2])
        
        """
        
        return self._indices
    
    def getEigenvalues(self):
        """Return eigenvalues.
        
        >>> modes.getEigenvalues() # doctest: +SKIP
        array([ 0.10652815,  0.1557703 ,  0.33017013])
        
        """
        
        return self._model._eigvals[self._indices].copy()

    def getEigenvectors(self):
        """Return eigenvectors.
        
        >>> modes.getEigenvectors() # doctest: +SKIP
        array([[-0.02171544, -0.04111048, -0.03200268],
               [ 0.05199527, -0.06430818, -0.01006391],
               [-0.05452411, -0.01514202, -0.05723971],
               ..., 
               [ 0.01144535, -0.04767575,  0.02667254],
               [-0.01913945, -0.0217476 , -0.00341831],
               [-0.02438398, -0.01942338,  0.00499727]])
               
        """
        
        return self._model._array[:, self._indices].copy()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues).
        
        >>> modes.getVariances() # doctest: +SKIP
        array([ 9.38719008,  6.41970913,  3.02874162])
        >>> modes.getEigenvalues()**-1 # doctest: +SKIP
        array([ 9.38719008,  6.41970913,  3.02874162])
        
        """
        
        return self._model._vars[self._indices].copy()

    def getArray(self):
        """Return eigenvectors.
        
        >>> modes.getArray() # doctest: +SKIP
        array([[-0.02171544, -0.04111048, -0.03200268],
               [ 0.05199527, -0.06430818, -0.01006391],
               [-0.05452411, -0.01514202, -0.05723971],
               ..., 
               [ 0.01144535, -0.04767575,  0.02667254],
               [-0.01913945, -0.0217476 , -0.00341831],
               [-0.02438398, -0.01942338,  0.00499727]])
               
        """
        
        return self._model._array[:, self._indices].copy()
        
    def getCovariance(self):
        """Return covariance matrix calculated for modes in the set.
        
        >>> modes.getCovariance() # doctest: +SKIP
        array([[ 0.01837834,  0.00734844,  0.02065894, ...,  0.00766404,
                 0.00997241,  0.00961239],
               [ 0.00734844,  0.05223408, -0.01861669, ...,  0.0244558 ,
                -0.0002593 , -0.00403514],
               [ 0.02065894, -0.01861669,  0.03930221, ..., -0.00584768,
                 0.01250275,  0.01350216],
               ..., 
               [ 0.00766404,  0.0244558 , -0.00584768, ...,  0.01797626,
                 0.00432368,  0.0037287 ],
               [ 0.00997241, -0.0002593 ,  0.01250275, ...,  0.00432368,
                 0.00651035,  0.00704099],
               [ 0.00961239, -0.00403514,  0.01350216, ...,  0.0037287 ,
                 0.00704099,  0.00807901]])
                
        """
        
        array = self.getArray()
        return np.dot(array, np.dot(np.diag(self.getVariances()), array.T))
   

class GNMBase(NMABase):
    """Class for Gaussian Network Model analysis of proteins."""

    def __init__(self, name):
        NMABase.__init__(self, name)
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
        NMABase._reset(self)
        self._cutoff = None
        self._gamma = None
        self._kirchhoff = None
    
    def getCutoff(self):
        """Return cutoff distance."""
        
        return self._cutoff
    
    def getGamma(self):
        """Return spring constant (or the gamma function or :class:`Gamma`
        instance).
        
        """
        return self._gamma

    def getKirchhoff(self):
        """Return Kirchhoff matrix."""
        
        if self._kirchhoff is None: return None
        return self._kirchhoff.copy()    

class GNM(GNMBase):
    
    """A class for Gaussian Network Model (GNM) analysis of proteins 
    ([IB97]_, [TH97]_).
    
    |example| See example :ref:`gnm`.
    
    """
    
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
    
    def buildKirchhoff(self, coords, cutoff=10., gamma=1., sparse=False):
        """Build Kirchhoff matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        :type coords: :class:`~numpy.ndarray` or :class:`~prody.atomic.Atomic`
        
        :arg cutoff: Cutoff distance (Å) for pairwise interactions
            default is 10.0 Å
        :type cutoff: float
        
        :arg gamma: Spring constant, default is 1.0.
        :type gamma: float
        
        .. versionchanged:: 0.6
            Instances of :class:`Gamma` classes and custom functions are
            accepted as *gamma* argument.        

        """
        
        if KDTree is None: 
            prody.importBioKDTree()
            if not KDTree:
                LOGGER.debug('Using a slower method for building the '
                             'Kirchhoff matrix.')
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
        assert cutoff > 0, 'cutoff distance must be greater than 0'
        self._cutoff = cutoff
        if isinstance(gamma, Gamma):
            self._gamma = gamma
            gamma = gamma.gamma
        elif isinstance(gamma, FunctionType):
            self._gamma = gamma
        else:
            g = float(gamma)
            assert g > 0, 'force constant (gamma) must be greater than 0'
            self._gamma = g
            gamma = lambda dist2, i, j: g
        n_atoms = coords.shape[0]
        start = time.time()
        if sparse:
            prody.importScipySparse()
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        if KDTree:
            kdtree = KDTree(3)
            kdtree.set_coords(coords) 
            kdtree.all_search(cutoff)
            radii = kdtree.all_get_radii()
            r = 0
            for i, j in kdtree.all_get_indices():
                g = gamma(radii[r]**2, i, j)
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] += g
                kirchhoff[j, j] += g
                r += 1
        else:
            cutoff2 = cutoff * cutoff
            for i in range(n_atoms):
                xyz_i = coords[i, :]
                for j in range(i+1, n_atoms):
                    i2j = coords[j, :] - xyz_i
                    dist2 = np.dot(i2j, i2j)
                    if dist2 > cutoff2:
                        continue             
                    g = gamma(dist2, i, j)
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] += g
                    kirchhoff[j, j] += g
            
        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'
                     .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._n_atoms = n_atoms
        self._dof = n_atoms
        
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize 
        Kirchhoff matrix. When Scipy is not found, :func:`numpy.linalg.eigh` 
        is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if self._kirchhoff is None:
            raise ProDyException('Kirchhoff matrix is not built or set')
        if linalg is None:
            prody.importLA()
        start = time.time()
        shift = 0
        if scipyla:
            if n_modes is None:
                eigvals = None
                n_modes = self._dof 
            else:
                n_modes = int(n_modes)
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = self._dof
                else: 
                    eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            if isinstance(self._kirchhoff, np.ndarray):            
                values, vectors = linalg.eigh(self._kirchhoff, turbo=turbo, 
                                              eigvals=eigvals)
            else:
                prody.importScipySparseLA()
                values, vectors = scipy_sparse_la.eigsh(self._kirchhoff, 
                                                       k=n_modes + 1,
                                                       which='SA')                
        else:
            values, vectors = linalg.eigh(self._kirchhoff)
        n_zeros = sum(values < ZERO)
        if n_zeros < 1: 
            LOGGER.warning('Less than 1 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        elif n_zeros > 1: 
            LOGGER.warning('More than 1 zero eigenvalues are calculated.')
            shift = n_zeros - 1
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

    """Class for Anisotropic Network Model (ANM) analysis of proteins 
    ([PD00]_, [ARA01]_)
    
    |example| See example :ref:`anm`.
    
    """

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
        
    def getHessian(self):
        """Return a copy of Hessian matrix."""
        
        if self._hessian is None: return None
        return self._hessian.copy()
    
    def setHessian(self, hessian):
        """Set Hessian matrix.
        
        A symmetric matrix is expected, i.e. not a lower- or upper-triangular
        matrix.
        
        """
        
        if not isinstance(hessian, np.ndarray):
            raise TypeError('hessian must be an ndarray')
        elif not (hessian.ndim == 2 and hessian.shape[0] == hessian.shape[1]):
            raise TypeError('hessian must be square matrix')
        self._reset()
        self._hessian = hessian
        self._dof = hessian.shape[0]
        self._n_atoms = self._dof / 3 

    def buildHessian(self, coords, cutoff=15., gamma=1.):
        """Build Hessian matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        :type coords: :class:`~numpy.ndarray` or :class:`~prody.atomic.Atomic`
        
        :arg cutoff: Cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å. 
        :type cutoff: float
        
        :arg gamma: Spring constant, default is 1.0. 
        :type gamma: float, :class:`Gamma`
        
        When available, this method makes use of Bio.KDTree.
        
        .. versionchanged:: 0.6
            Instances of :class:`Gamma` classes and custom functions are
            accepted as *gamma* argument.        
                           
        """
        
        if KDTree is None: 
            prody.importBioKDTree()
            if not KDTree:
                LOGGER.debug('Using a slower method for building the '
                             'Hessian matrix.')
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
        assert cutoff > 0, 'cutoff distance must be greater than 0'
        self._cutoff = cutoff
        if isinstance(gamma, Gamma):
            self._gamma = gamma
            gamma = gamma.gamma
        elif isinstance(gamma, FunctionType):
            self._gamma = gamma
        else:
            g = float(gamma)
            assert g > 0, 'force constant (gamma) must be greater than 0'
            self._gamma = g
            gamma = lambda dist2, i, j: g 
         
        n_atoms = coords.shape[0]
        dof = n_atoms * 3
        start = time.time()
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        hessian = np.zeros((dof, dof), 'd')
        if KDTree:
            kdtree = KDTree(3)
            kdtree.set_coords(coords) 
            kdtree.all_search(cutoff)
            for i, j in kdtree.all_get_indices():
                #if k < i:
                #    j = i
                #    i = k
                #else:
                #    j = k
                i2j = coords[j] - coords[i]
                dist2 = np.dot(i2j, i2j)
                g = gamma(dist2, i, j)
                super_element = np.outer(i2j, i2j) * (- g / dist2)  
                res_i3 = i*3
                res_i33 = res_i3+3
                res_j3 = j*3
                res_j33 = res_j3+3
                hessian[res_i3:res_i33, res_j3:res_j33] = super_element
                hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                hessian[res_i3:res_i33, res_i3:res_i33] -= super_element
                hessian[res_j3:res_j33, res_j3:res_j33] -= super_element
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] += g
                kirchhoff[j, j] += g
        else:
            cutoff2 = cutoff * cutoff 
            for i in range(n_atoms):
                res_i3 = i*3
                res_i33 = res_i3+3
                xyz_i = coords[i, :]
                for j in range(i+1, n_atoms):
                    i2j = coords[j, :] - xyz_i
                    dist2 = np.dot(i2j, i2j)
                    if dist2 > cutoff2:
                        continue             
                    g = gamma(dist2, i, j)
                    res_j3 = j*3
                    res_j33 = res_j3+3
                    super_element = np.outer(i2j, i2j) * (- g / dist2) 
                    hessian[res_i3:res_i33, res_j3:res_j33] = super_element 
                    hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                    hessian[res_i3:res_i33, res_i3:res_i33] -= super_element
                    hessian[res_j3:res_j33, res_j3:res_j33] -= super_element
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] += g
                    kirchhoff[j, j] += g
        LOGGER.info('Hessian was built in {0:.2f}s.'.format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._n_atoms = n_atoms
        self._dof = dof

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize 
        Hessian matrix. When Scipy is not found, :func:`numpy.linalg.eigh` 
        is used.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        
        if self._hessian is None:
            raise ProDyException('Hessian matrix is not built or set')
        if linalg is None:
            prody.importLA()
            
        start = time.time()
        shift = 5
        if scipyla:
            if n_modes is None:
                eigvals = None
                n_modes = self._dof 
            else: 
                n_modes = int(n_modes)
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = self._dof 
                else:
                    eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            values, vectors = linalg.eigh(self._hessian, turbo=turbo, 
                                          eigvals=eigvals)
        else:
            values, vectors = linalg.eigh(self._hessian)
        n_zeros = sum(values < ZERO)
        if n_zeros < 6: 
            LOGGER.warning('Less than 6 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        elif n_zeros > 6: 
            LOGGER.warning('More than 6 zero eigenvalues are calculated.')
            shift = n_zeros - 1
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


class PCA(NMABase):
    
    """A class for Principal Component Analysis (PCA) of conformational ensembles 
    (also known as Essential Dynamics Analysis (EDA) in [AA93]_).
    
    
    |example| See examples in :ref:`pca`.
    
    """

    def __init__(self, name):
        NMABase.__init__(self, name)
    
    def __repr__(self):
        return '<PCA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._name, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'PCA {0:s}'.format(self._name)
    
    def setCovariance(self, covariance):
        """Set covariance matrix."""
        
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
        
        *coordsets* argument must be a :class:`~prody.atomic.Atomic` or a 
        :class:`~prody.ensemble.Ensemble` instance. 
        
        If *weights* is ``None``, but *coordsets* has a :meth:`getWeights` 
        method, weights from that method will be used. 
        
        """
        
        start = time.time()
        
        if not isinstance(coordsets, (prody.Ensemble, prody.Atomic)):
            raise TypeError('coordsets must be an Ensemble or Atomic instance')
        n_atoms = coordsets.getNumOfAtoms()
        n_confs = coordsets.getNumOfCoordsets()
        if n_confs < 3:
            raise ValueError('coordsets must have more than 3 coordinate sets')
        if n_atoms < 3:
            raise ValueError('coordsets must have more than 3 atoms')

        if isinstance(coordsets, prody.Atomic):
            coordsets = prody.Ensemble(coordsets)
        dof = n_atoms * 3
        if weights is None:
            try:
                weights = coordsets.getWeights()
            except AttributeError:
                pass
        if weights is None:
            #d_xyz = (conformations - coordinates)
            #d_xyz = d_xyz.reshape((n_confs, dof))
            #self._cov = np.dot(d_xyz.T, d_xyz) / n_confs
            self._cov = np.cov(
                            coordsets.getCoordsets().reshape((n_confs, dof)).T,
                            bias=1)
        else:
            conformations = coordsets.getCoordsets()
            coordinates = coordsets.getCoordinates()
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
        """Calculate principal (or essential) modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize
        covariance matrix. When Scipy is not found, :func:`numpy.linalg.eigh` 
        is used.
        
        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        
        if linalg is None:
            prody.importLA()
        if self._cov is None:
            raise ProDyException('covariance matrix is not built or set')
        start = time.time()
        dof = self._dof
        if scipyla:        
            if n_modes is None:
                eigvals = None
                n_modes = dof
            else:
                n_modes = int(n_modes)
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = dof
                else:
                    eigvals = (dof - n_modes, dof - 1)
            values, vectors = linalg.eigh(self._cov, turbo=turbo, 
                                          eigvals=eigvals)
        else:
            values, vectors = linalg.eigh(self._cov)
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
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def performSVD(self, coordsets):
        """Calculate principal modes using singular value decomposition (SVD).
        
        *coordsets* argument must be a :class:`~prody.atomic.Atomic` or a 
        :class:`~prody.ensemble.Ensemble` instance. 
        
        .. versionadded:: 0.6.2
        
        This is a considerably faster way of performing PCA calculations 
        compared to eigenvalue decomposition of covariance matrix, but is
        an approximate method when heterogeneous datasets are analyzed. 
        Covariance method should be preferred over this one for analysis
        of ensembles with missing atomic data. See :ref:`pca-xray-calculations`
        examples for comparison of results from SVD and covariance methods.
        
        
        """

        if linalg is None:
            prody.importLA()

        start = time.time()

        if not isinstance(coordsets, (prody.Ensemble, prody.Atomic)):
            raise TypeError('coordsets must be an Ensemble or Atomic instance')
        n_atoms = coordsets.getNumOfAtoms()
        n_confs = coordsets.getNumOfCoordsets()
        if n_confs < 3:
            raise ValueError('coordsets must have more than 3 coordinate sets')
        if n_atoms < 3:
            raise ValueError('coordsets must have more than 3 atoms')
        if isinstance(coordsets, prody.Ensemble):
            deviations = coordsets.getDeviations()
        elif isinstance(coordsets, prody.Atomic):
            deviations = coordsets.getCoordsets() - coordsets.getCoordinates()

        if isinstance(coordsets, prody.Atomic):
            coordsets = prody.Ensemble(coordsets)
        dof = n_atoms * 3        
        deviations = deviations.reshape((n_confs, dof)).T

        vectors, values, self._temp = linalg.svd(deviations, full_matrices=False)
        values = (values ** 2) / n_confs
        self._dof = dof
        self._n_atoms = n_atoms
        which = values > 1e-18
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._trace = self._vars.sum()
        self._n_modes = len(self._eigvals)
        self._modes = [None] * self._n_modes
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                         .format(self._n_modes, time.time()-start))
        
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to the :class:`NMA` instance.
        
        .. versionadded:: 0.7
        
        If *eigenvalue* is not given, it will be set to 1. 
        Eigenvalue is also set as the variance.
        
        """

        NMABase.addEigenpair(self, eigenvector, eigenvalue)
        self._vars = self._eigvals.copy()


    def setEigens(self, vectors, values=None):
        """Set eigenvectors and eigenvalues.
        
        .. versionadded:: 0.7
        
        :arg vectors: eigenvectors
        :type vectors: numpy.ndarray
        
        :arg values: Eigenvalues. When ``None`` is passed (default value), 
            all eigenvalues will be set to ``1``.
        :type values: numpy.ndarray
        
        For M modes and N atoms, *vectors* must have shape ``(3*N, M)``
        and values must have shape ``(M,)``.
        
        Eigenvalues are also set as the variances.
        
        """
        
        NMABase.setEigens(self, vectors, values)
        self._vars = self._eigvals.copy()

EDA = PCA

class Gamma(object):
    
    """Base class for facilitating use of atom type, residue type, or residue
    property dependent force constants (γ).
    
    .. versionadded:: 0.6    
        
    Derived classes:
        
    * :class:`GammaStructureBased`
    * :class:`GammaVariableCutoff`
    
    """
    
    
    def __init__(self):
        pass
    
    def gamma(self, dist2, i, j):
        """Return force constant.
        
        For efficiency purposes square of the distance between interacting
        atom/residue (node) pairs is passed to this function. In addition, 
        node indices are passed.
        
        """
        
        pass
    
    
class GammaStructureBased(Gamma):
    
    """Facilitate setting the spring constant based on the secondary structure 
    and connectivity of the residues.
    
    .. versionadded:: 0.6    
        
    A recent systematic study [LT10]_ of a large set of NMR-structures analyzed 
    using a method based on entropy maximization showed that taking into 
    consideration properties such as sequential separation between 
    contacting residues and the secondary structure types of the interacting 
    residues provides refinement in the ENM description of proteins.
    
    This class determines pairs of connected residues or pairs of proximal 
    residues in a helix or a sheet, and assigns them a larger user defined 
    spring constant value.
    
     DSSP single letter abbreviations are recognized: 
       * **H**: α-helix
       * **G**: 3-10-helix
       * **I**: π-helix
       * **E**: extended part of a sheet
    
    *helix*: 
        Applies to residue (or Cα atom) pairs that are in the same helical 
        segment, at most 7 Å apart, and separated by at most 
        3 (3-10-helix), 4 (α-helix), or 5 (π-helix) residues.
        
    *sheet*:  
        Applies to Cα atom pairs that are in different β-strands and at most 
        6 Å apart.
        
    *connected*:
        Applies to Cα atoms that are at most 4 Å apart.
        
    Note that this class does not take into account insertion codes.        
    
    **Example**:

    Let's parse coordinates and header data from a PDB file, and then
    assign secondary structure to the atoms. 
        
    >>> from prody import *
    >>> ubi, header = parsePDB('1aar', chain='A', subset='calpha', header=True)
    >>> assignSecondaryStructure(header, ubi)
    <AtomGroup: 1aar (76 atoms; 1 coordinate sets, active set index: 0)>

    In the above we parsed only the atoms needed for this calculation, i.e.
    Cα atoms from chain A. 
    
    We build the Hessian matrix using structure based force constants as 
    follows;
    
    >>> gamma = GammaStructureBased(ubi)
    >>> anm = ANM('')
    >>> anm.buildHessian(ubi, gamma=gamma)
    
    We can obtain the force constants assigned to residue pairs from the 
    Kirchhoff matrix as follows: 
    
    >>> k = anm.getKirchhoff()
    >>> k[0,1] # a pair of connected residues
    -10.0
    >>> k[0,16] # a pair of residues from a sheet
    -6.0
    
    """
    
    def __init__(self, atoms, gamma=1.0, helix=6.0, sheet=6.0, connected=10.0):
        """Setup the parameters.
        
        
        :arg atoms: A set of atoms with chain identifiers, residue numbers,
            and secondary structure assignments are set.
        :type atoms: :class:`~prody.atomic.Atomic`

        :arg gamma: Force constant in arbitrary units. Default is 1.0.
        :type gamma: float
            
        :arg helix: Force constant factor for residues hydrogen bonded in 
            α-helices, 3,10-helices, and π-helices. Default is 6.0, i.e.
            ``6.0`*gamma``.
        :type helix: float

        :arg sheet: Force constant factor for residue pairs forming a hydrogen 
            bond in a β-sheet. Default is 6.0, i.e. ``6.0`*gamma``.
        :type sheet: float
            
        :arg connected: Force constant factor for residue pairs that are
            connected. Default is 10.0, i.e. ``10.0`*gamma``.
        :type connected: float
        
        """
        
        if not isinstance(atoms, prody.Atomic):
            raise TypeError('atoms must be an Atomic instance')
        n_atoms = atoms.getNumOfAtoms()
        sstr = atoms.getSecondaryStrs()
        assert sstr is not None, 'secondary structure assignments must be set'
        chid = atoms.getChainIdentifiers()
        assert chid is not None, 'chain identifiers must be set'
        rnum = atoms.getResidueNumbers()
        assert rnum is not None, 'residue numbers must be set'
        gamma = float(gamma)
        assert gamma > 0, 'gamma must be greater than 0'
        helix = float(helix)
        assert helix > 0, 'helix must be greater than 0'
        sheet = float(sheet)
        assert sheet > 0, 'sheet must be greater than 0'
        connected = float(connected)
        assert connected > 0, 'connected must be greater than 0'
        
        ssid = np.zeros(n_atoms)
        for i in range(1, n_atoms):
            if (sstr[i-1] == sstr[i] and chid[i-1] == chid[i] and
               rnum[i]-rnum[i-1] == 1): 
                ssid[i] = ssid[i-1]
            else:
                ssid[i] = ssid[i-1] + 1
        self._sstr = sstr
        self._chid = chid
        self._rnum = rnum
        self._ssid = ssid
        self._gamma = gamma
        self._helix = gamma * helix
        self._sheet = gamma * sheet
        self._connected = gamma * connected
    
    def getSecondaryStr():
        """Return a copy of secondary structure assignments."""
        
        return self._sstr.copy()    
    
    def getChainIdentifiers():
        """Return a copy of chain identifiers."""
        
        return self._chid.socopypy()    

    def getResidueNumbers():
        """Return a copy of residue numbers."""
        
        return self._rnum.copy()    


    def gamma(self, dist2, i, j):
        """Return force constant."""
        
        if dist2 <= 16:
            return self._connected
        sstr = self._sstr
        ssid = self._ssid
        rnum = self._rnum
        if ssid[i] == ssid[j]:
            i_j = abs(rnum[j] - rnum[i])
            if ((i_j <= 4 and sstr[i] == 'H') or 
                (i_j <= 3 and sstr[i] == 'G') or 
                (i_j <= 5 and sstr[i] == 'I')) and dist2 <= 49: 
                return self._helix
        elif sstr[i] == sstr[j] == 'E' and dist2 <= 36:
            return self._sheet
        
        return self._gamma
    
    
class GammaVariableCutoff(Gamma):
    
    """Facilitate setting the cutoff distance based on user defined 
    atom/residue (node) radii.
    
    .. versionadded:: 0.6    
        
    Half of the cutoff distance can be thought of as the radius of a node. 
    This class enables setting different radii for different node types.
    
    **Example**:
    
    Let's think of a protein-DNA complex for which we want to use different
    radius for different residue types. Let's say, for protein Cα atoms we
    want to set the radius to 7.5 Å, and for nucleic acid phosphate atoms to 
    10 Å. We use the HhaI-DNA complex structure :file:`1mht`.

    >>> hhai = parsePDB('1mht')
    >>> ca_p = hhai.select('(protein and name CA) or (nucleic and name P)')
    >>> print ca_p.getAtomNames() # doctest: +SKIP
    ['P' 'P' 'P' ..., 'CA' 'CA' 'CA']
    
    We set the radii of atoms: 
     
    >>> variableCutoff = GammaVariableCutoff(ca_p.getAtomNames(), gamma=1, default_radius=7.5, debug=True, P=10)
    >>> print variableCutoff.getRadii() # doctest: +SKIP
    [ 10.   10.   10.   ...,   7.5   7.5   7.5]
    
    The above shows that for phosphate atoms radii is set to 10 Å, because
    we passed the ``P=10`` argument. As for Cα atoms, the default 7.5 Å
    is set as the radius (``default_radius=7.5``). Note we also passed
    ``debug=True`` argument for demonstration purposes. This argument 
    allows printing debugging information on the screen.
    
    We build :class:`ANM` Hessian matrix as follows:  
        
    >>> anm = ANM('HhaI-DNA')
    >>> anm.buildHessian(ca_p, gamma=variableCutoff, cutoff=20) # doctest: +SKIP
    CA_275 -- P_7 effective cutoff: 17.5 distance: 19.5948930081 gamma: 0
    CA_275 -- CA_110 effective cutoff: 15.0 distance: 13.5699586587 gamma: 1.0
    P_20 -- P_6 effective cutoff: 20.0 distance: 18.150633763 gamma: 1.0
    CA_275 -- P_6 effective cutoff: 17.5 distance: 18.7748343268 gamma: 0
    CA_275 -- CA_109 effective cutoff: 15.0 distance: 10.5352334573 gamma: 1.0
    ...
    
    Note that we set passed ``cutoff=20.0`` to the :meth:`ANM.buildHessian` 
    method. This is equal to the largest possible cutoff distance (between two
    phosphate atoms) for this system, and ensures that all of the potential 
    interactions are evaluated. 
    
    For pairs of atoms for which the actual distance is larger than the 
    effective cutoff, the :meth:`GammaVariableCutoff.gamma` method returns 
    ``0``. This annuls the interaction between those atom pairs.
    
    """
    
    def __init__(self, identifiers, gamma=1., default_radius=7.5, debug=False, 
                 **kwargs):
        """Set the radii of atoms.
        
        :arg identifiers: List of atom names or types, or residue names.
        :type identifiers: list or :class:`numpy.ndarray`
        
        :arg gamma: Uniform force constant value. Default is 1.0.
        :type gamma: float
        
        :arg default_radius: Default radius for atoms whose radii is not set
            as a keyword argument. Default is 7.5
        :value default_radius: float
        
        :arg debug: Print debugging information. Default is ``False``.
        :type debug: bool
        
        Keywords in keyword arguments must match those in *atom_identifiers*.
        Values of keyword arguments must be :class:`float`.  
        
        """
        
        self._identifiers = identifiers
        radii = np.ones(len(identifiers)) * default_radius
        
        for i, identifier in enumerate(identifiers): 
            radii[i] = kwargs.get(identifier, default_radius)
        self._radii = radii
        self._gamma = float(gamma)
        self._debug = bool(debug)
    
    def getRadii(self):
        """Return a copy of radii array."""
        
        return self._radii.copy()

    def getGamma(self):
        """Return the uniform force constant value."""
        
        return self._gamma_constant

    def gamma(self, dist2, i, j):
        """Return force constant."""
        
        cutoff = (self._radii[i] + self._radii[j])
        cutoff2 = cutoff ** 2
        
        if dist2 < cutoff2: 
            gamma = self._gamma
        else:
            gamma = 0
        if self._debug:
            print self._identifiers[i]+'_'+str(i), '--', \
                  self._identifiers[j]+'_'+str(j), \
                  'effective cutoff:', cutoff, 'distance:', dist2**0.5, \
                  'gamma:', gamma    
        return gamma


def saveModel(nma, filename=None, matrices=False):
    """Save *nma* model data as :file:`filename.nma.npz`. 
    
    .. versionadded:: 0.5
    
    By default, eigenvalues, eigenvectors, variances, trace of covariance 
    matrix, and name of the model will be saved. If *matrices* is ``True``,
    covariance, Hessian or Kirchhoff matrices are saved too, whichever are 
    available.
    
    If *filename* is ``None``, name of the NMA instance will be used as 
    the filename, after " " (blank spaces) in the name are replaced with "_" 
    (underscores) 
    
    Extension may differ based on the type of the NMA model. For ANM models,
    it is :file:`.anm.npz`.
    
    Upon successful completion of saving, filename is returned.
    
    This function makes use of :func:`numpy.savez` function.
    
    """
    if not isinstance(nma, NMABase):
        raise TypeError('invalid type for nma, {0:s}'.format(type(nma)))
    if len(nma) == 0:
        raise ValueError('nma instance does not contain data')
    
    dict_ = nma.__dict__
    attr_list = ['_name', '_trace', '_array', '_eigvals', '_vars', '_n_atoms',
                 '_dof', '_n_modes']
    if filename is None:
        filename = nma.getName().replace(' ', '_')
    if isinstance(nma, GNMBase):
        attr_list.append('_cutoff')
        attr_list.append('_gamma')
        if matrices:
            attr_list.append('_kirchhoff')
            if isinstance(nma, ANM):
                attr_list.append('_hessian')
        if isinstance(nma, ANM):
            type_ = 'ANM'
        else:
            type_ = 'GNM'
    elif isinstance(nma, PCA):
        type_ = 'PCA'
    else:
        type_ = 'NMA'  
    
    if matrices:
        attr_list.append('_cov')
    attr_dict = {'type': type_}
    for attr in attr_list:
        value = dict_[attr]
        if value is not None:
            attr_dict[attr] = value
    filename += '.' + type_.lower() + '.npz'
    np.savez(filename, **attr_dict)
    return filename



def loadModel(filename):
    """Return NMA instance after loading it from file (*filename*).
    
    .. versionadded:: 0.5
    
    .. seealso: :func:`saveModel`
    
    This function makes use of :func:`numpy.load` function.
    
    """
    
    attr_dict = np.load(filename)
    type_ = attr_dict['type']
    if type_ == 'ANM':
        nma = ANM(str(attr_dict['_name']))
    elif type_ == 'PCA':
        nma = PCA(str(attr_dict['_name']))
    elif type_ == 'GNM':
        nma = GNM(str(attr_dict['_name']))
    else:
        nma = NMA(str(attr_dict['_name']))
    dict_ = nma.__dict__ 
    for attr in attr_dict.files:
        if attr in ('type', '_name'): 
            continue
        elif attr in ('_trace', '_cutoff', '_gamma'):
            dict_[attr] = float(attr_dict[attr])
        elif attr in ('_dof', '_n_atoms', '_n_modes'):
            dict_[attr] = int(attr_dict[attr])
        else:
            dict_[attr] = attr_dict[attr]
    nma._modes = [None] * len(nma)
    return nma

def saveVector(vector, filename):
    """Save *vector* data as :file:`filename.vec.npz`. 
    
    Upon successful completion of saving, filename is returned.
    
    This function makes use of :func:`numpy.savez` function.
    
    """
    
    if not isinstance(vector, Vector):
        raise TypeError('invalid type for vector, {0:s}'.format(type(vector)))
    attr_dict = {}
    attr_dict['name'] = vector.getName()
    attr_dict['array'] = vector.getArray()
    attr_dict['is3d'] = vector.is3d()
    filename += '.vec.npz'
    np.savez(filename, **attr_dict)
    return filename

def loadVector(filename):
    """Return :class:`Vector` instance after loading it from file (*filename*).
    
    .. seealso: :func:`saveVector`
    
    This function makes use of :func:`numpy.load` function.
    
    """
    
    attr_dict = np.load(filename)
    return Vector(attr_dict['array'], str(attr_dict['name']), 
                  bool(attr_dict['is3d']))

def getVMDpath():
    """Return path to the VMD executable."""
    path = prody._ProDySettings.get('vmd')
    if path is None:
        LOGGER.warning('VMD path is not set.')
    return path

def setVMDpath(path):
    """Set the path to VMD executable."""
    
    if not os.path.isfile(path):
        LOGGER.warning('{0:s} is not a file.'.format(path))
        return
    prody._ProDySettings['vmd'] = path
    prody._saveProDySettings()
    
    

def parseNMD(filename, type=NMA):
    """Returns normal mode and atomic data parsed from an NMD file.
    
    .. versionadded:: 0.5.3
    
    .. versionchanged:: 0.7
       User can pass NMA type for the data, eg. :class:`ANM` or :class:`PCA`.
    
    Normal mode data is returned in an :class:`NMA` instance. Atomic
    data is returned in an :class:`~prody.atomic.AtomGroup` instance. 
    
    """
    assert not isinstance(type, NMABase), 'type must be NMA, ANM, GNM, or PCA'
    
    atomic = dict()
    modes = []
    nmd = open(filename)
    for line in nmd:
        split = line.find(' ')
        if line[:split] == 'mode':
            modes.append(line[split:].strip())
        elif line[:split] in ('coordinates', 'atomnames', 'resnames', 
                              'resnums', 'resids', 'chainids', 'bfactors',
                              'name'):
            atomic[line[:split]] = line[split:].strip()
    nmd.close()
    
    name = atomic.pop('name', os.path.splitext(os.path.split(filename)[1])[0])
    coords = atomic.pop('coordinates', None)
    dof = None
    if coords is not None:
        coords = np.fromstring( coords, dtype=np.float64, sep=' ')
        dof = coords.shape[0]
        ag = None
        n_atoms = dof / 3
        coords = coords.reshape((n_atoms, 3))
        ag = AtomGroup(name)
        ag.setCoordinates(coords)
        data = atomic.pop('atomnames', None)
        if data is not None:
            ag.setAtomNames(data.split())
        data = atomic.pop('resnames', None)
        if data is not None:
            ag.setResidueNames(data.split())
        data = atomic.pop('chainids', None)
        if data is not None:
            ag.setChainIdentifiers(data.split())
        data = atomic.pop('resnums', None)
        if data is not None:
            ag.setResidueNumbers(np.fromstring(data, np.int64, sep=' '))
        data = atomic.pop('resids', None)
        if data is not None:
            ag.setResidueNumbers(np.fromstring(data, np.int64, sep=' '))
        data = atomic.pop('bfactors', None)
        if data is not None:
            ag.setTempFactors(np.fromstring(data, np.float64, sep=' '))
    nma = type(name)
    for mode in modes:
        
        items = mode.split()
        diff = len(items) - dof
        mode = np.array(items[diff:]).astype(np.float64)
        if len(mode) != dof:
            pass
        if diff == 1 and not items[0].isdigit():
            value = float(items[0])
        else:
            if not items[0].isdigit():
                value = float(items[0])
            elif not items[1].isdigit():
                value = float(items[1])
            else:
                value = 1.0
        nma.addEigenpair(mode, value)
    return nma, ag
    

def writeNMD(filename, modes, atoms):
    """Writes an NMD file for given *modes* and includes applicable data from 
    *atoms*.
    
    Returns *filename*, if file is successfully written.
    
    NMD file format is described at 
    http://www.csb.pitt.edu/People/abakan/software/NMWiz/nmdformat.html.
    
    .. note:: 
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`Vector` instance is given, it will be normalized before
          it is written. It's length before normalization will be written
          as the scaling factor of the vector.
        
    """
    
    if not isinstance(modes, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0:s}'.format(type(modes)))
    if modes.getNumOfAtoms() != atoms.getNumOfAtoms():
        raise Exception('number of atoms do not match')
    out = open(filename, 'w')
    
    #out.write('#!{0:s} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0:s}\n'.format(os.path.abspath(filename)))
    name = modes.getName()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = os.path.splitext(os.path.split(filename)[1])[0]
    out.write('name {0:s}\n'.format(name))
    try:
        coords = atoms.getCoordinates()
    except:
        raise ProDyException('coordinates could not be retrived '
                             'from atoms instance')
    if coords is None:
        raise ProDyException('coordinates could not be retrived '
                             'from atoms instance')
    
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
    
    try:
        data = atoms.getTempFactors()
        if data is not None:
            out.write('bfactors {0:s}\n'.format(' '.join(['{0:.3f}'.format(x) for x in data.flatten()])))
    except:
        pass
    
    out.write('coordinates {0:s}\n'.format(
                    ' '.join(['{0:.3f}'.format(x) for x in coords.flatten()])))
    
    count = 0
    if isinstance(modes, Vector):
        out.write('mode 1 {0:.2f} {1:s}\n'.format(abs(modes), ' '.join(
                ['{0:.3f}'.format(x) for x in modes.getNormed().getArray()])))
        count += 1
    else:
        for mode in modes:
            if mode.getEigenvalue() < ZERO:
                continue
            out.write('mode {0:d} {1:.2f} {2:s}\n'.format(
                       mode.getIndex()+1, mode.getVariance()**0.5, 
                       ' '.join(
                            ['{0:.3f}'.format(x) for x in mode.getArray()])))
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. '
                       'Given modes might have 0 eigenvalues.')
    out.close() 
    return filename  

def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""
    
    os.system('{0:s} -e {1:s}'.format(prody._ProDySettings['vmd'], 
                                      os.path.abspath(filename)))
    
def calcANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, 
            zeros=False):
    """Return an :class:`ANM` instance and atoms used for the calculations.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`~prody.atomic.Atomic` instance.
    
    .. versionchanged:: 0.6
       Returns also the :class:`~prody.atomic.Selection` instance.
    
    """
    
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if isinstance(pdb, AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'
                        .format(type(pdb)))
    anm = ANM(name)
    sel = ag.select(selstr)
    anm.buildHessian(sel, cutoff, gamma)
    anm.calcModes(n_modes)
    return anm, sel 

def calcGNM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Return a :class:`GNM` instance and atoms used for the calculations.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`~prody.atomic.Atomic` instance.  
    
    .. versionchanged:: 0.6
       Returns also the :class:`~prody.atomic.Selection` instance.
    
    """
    
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        name = ag.getName()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if isinstance(pdb, AtomGroup):
            name = ag.getName()
        else: 
            name = ag.getAtomGroup().getName()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'
                        .format(type(pdb)))
    gnm = GNM(name)
    sel = ag.select(selstr)
    gnm.buildKirchhoff(sel, cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm, sel

def calcCollectivity(mode, masses=None):
    """Return collectivity of the mode.
    
    :arg mode: mode or vector
    :type mode: :class:`Mode` or :class:`Vector`
    
    :arg masses: atomic masses
    :type masses: :class:`numpy.ndarray`
    
    This function implements collectivity as defined in equation 5 of [BR95]_.
    
    If *masses* are provided, they will be incorporated in the calculation.
    Otherwise, atoms are assumed to have uniform masses.
    
    """
    
    is3d = mode.is3d()
    if masses is not None:
        if len(masses) != mode.getNumOfAtoms(): 
            raise ValueError('length of massesmust be equal to number of atoms')
        if is3d:
            u2in = (mode.getArrayNx3() ** 2).sum(1) / masses
    else:
        if is3d:
            u2in = (mode.getArrayNx3() ** 2 ).sum(1)
        else:
            u2in = (mode.getArrayNx3() ** 2 )
    u2in = u2in * (1 / u2in.sum() ** 0.5)
    coll = np.exp(-(u2in * np.log(u2in)).sum()) / mode.getNumOfAtoms()
    return coll
    

def calcProjection(ensemble, modes):
    """Return projection of conformational deviations onto given modes.

    For K conformations and M modes, a (K,M) matrix is returned.     
    
    >>> calcProjection(p38_ensemble, p38_pca[:3]) # doctest: +SKIP
    array([[ 11.41 ,   2.651,   1.116],
           [  7.846,  -1.856,   1.591],
           ...
           [  5.516,  -5.108,  -5.983],
           [  5.875,  -7.408,   2.208]])
    
    """
    
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'
                        .format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
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


def calcOverlap(rows, cols):
    """Return overlap (or correlation) between two sets of modes 
    (*rows* and *cols*).
    
    Returns a matrix whose rows correspond to modes passed as 
    *rows* argument, and columns correspond to those passed as *cols* 
    argument.
    
    >>> calcOverlap(p38_pca[0], p38_anm[2]) # doctest: +SKIP
    -0.71366564906422636
    
    .. versionchanged:: 0.7
       Both rows and columns are normalized prior to calculating overlap.       
    
    """
    
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(rows)))
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(cols)))
    
    if rows.getNumOfDegOfFreedom() != cols.getNumOfDegOfFreedom(): 
        raise ValueError('number of degrees of freedom of rows and '
                         'cols must be the same')
    rows = rows.getArray()
    rows *= 1 / (rows ** 2).sum(0) ** 0.5
    cols = cols.getArray()
    cols *= 1 / (cols ** 2).sum(0) ** 0.5
    return np.dot(rows.T, cols)

def printOverlapTable(rows, cols):
    """Print table of overlaps (correlations) between two sets of modes.
    
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the printed table.

    This function may be used to take a quick look into mode correspondences 
    between two models.
    
    >>> # Compare top 3 PCs and slowest 3 ANM modes
    >>> printOverlapTable(p38_pca[:3], p38_anm[:3]) # doctest: +SKIP   
    Overlap Table
                            ANM 1p38
                        #1     #2     #3
    PCA p38 xray #1   -0.39  +0.04  -0.71
    PCA p38 xray #2   -0.78  -0.20  +0.22
    PCA p38 xray #3   +0.05  -0.57  +0.06
    """
    print calcOverlapTable(rows, cols)

def writeOverlapTable(filename, rows, cols):
    """Write table of overlaps (correlations) between two sets of modes to a 
    file.

    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the overlap table.
    
    See also :func:`printOverlapTable`.
    """
    out = open(filename, 'w')
    out.write(calcOverlapTable(rows, cols))
    out.close()
    return filename
    
def calcOverlapTable(rows, cols):
    overlap = calcOverlap(rows, cols)
    if isinstance(rows, Mode):
        rids = [rows.getIndex()]
        rname = str(rows.getModel())
    elif isinstance(rows, NMABase): 
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
    elif isinstance(cols, NMABase): 
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
    table += (' '*(len(rname)+5) + cname.center(clen*7)).rstrip() + '\n'
    line = ' '*(len(rname)+5)
    for j in range(clen):
        line += ('#{0}'.format(cids[j]+1)).center(7)
    table += line.rstrip() + '\n'
    for i in range(rlen):
        line = rname + (' #{0}'.format(rids[i]+1)).ljust(5)
        for j in range(clen):
            if overlap[i, j] < 0: 
                minplus = '-'
            else: 
                minplus = '+'
            line += (minplus+'{0:-.2f}').format(abs(overlap[i, j])).center(7)
        table += line.rstrip() + '\n'
    return table

def extrapolateModel(enm, nodes, atoms):
    """Extrapolate *enm* built for *nodes* to *atoms*.
    
    .. version added:: 0.6.2
    
    This function is designed for extrapolating an NMA model built at coarse 
    grained level to all atom level. For each atom in *nodes* argument
    *atoms* argument must contain a corresponding residue.
    
    Note that modes in the extrapolated model will not be normalized.
    
    For a usage example see :ref:`extrapolate`.
    
    """
    
    if not isinstance(enm, NMABase):
        raise TypeError('enm must be an NMABase instance')
    if not isinstance(nodes, Atomic):
        raise TypeError('nodes must be an Atomic instance')
    if enm.getNumOfAtoms() != nodes.getNumOfAtoms():
        raise ValueError('enm and nodes must have same number of atoms')
    
    if isinstance(atoms, Atomic):
        is3d = enm.is3d()            
        atom_indices = []
        indices = []
        hierview = atoms.getHierView()
        for i, node in enumerate(nodes):
            res = hierview[node.getChainIdentifier(), node.getResidueNumber(), 
                           node.getInsertionCode()]
            if res is None:
                raise ValueError('hierview must contain a residue for all atoms')
            atom_indices.append(res.getIndices())
            if is3d:
                indices.append(range(i*3, (i+1)*3) * len(res))
            else:
                indices.append([i] * len(res))
        atom_indices = np.concatenate(atom_indices)
        indices = np.concatenate(indices)
        
        array = enm.getArray()[indices,:]
        extra = NMA('Extrapolated ' + str(enm))
        extra.setEigens(array, enm.getEigenvalues())
        if isinstance(atoms, AtomGroup):
            ag = atoms
        else: 
            ag = atoms.getAtomGroup()
        atommap = AtomMap(ag, atom_indices, np.arange(len(atom_indices)), 
                          np.array([]), str(atoms), 
                          atoms.getActiveCoordsetIndex())
        return extra, atommap
    else:
        raise TypeError('atoms must be an Atomic instance')
    

def sliceVector(vector, atoms, selstr):
    """Return a slice of *vector* matching *atoms* specified by *selstr*.
    
    .. versionadded:: 0.5
    
    .. versionchanged:: 0.7
       Returns the vector and the corresponding atom selection. 
    
    Note that returned :class:`Vector` instance is not normalized.
    
    :arg vector: vector instance to be sliced
    :type vector: :class:`VectorBase`
    
    :arg atoms: atoms for which *vector* describes a deformation, motion, etc.
    :type atoms: :class:`~prody.atomic.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`Vector`, :class:`~prody.atomic.Selection`)
    
    """
    
    if not isinstance(vector, VectorBase):
        raise TypeError('vector must be a VectorBase instance, not {0:s}'
                        .format(type(vector)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.getNumOfAtoms() != vector.getNumOfAtoms(): 
        raise ValueError('number of atoms in *vector* and *atoms* must be '
                         'equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getActiveCoordsetIndex())
    vec = Vector(vector.getArrayNx3()[
                 which, :].flatten(),
                 '{0:s} slice "{1:s}"'.format(str(vector), selstr), 
                 vector.is3d())
    return (vec, sel)

def sliceMode(mode, atoms, selstr):
    """Return a slice of *mode* matching *atoms* specified by *selstr*.
    
    .. versionadded:: 0.5
    
    .. versionchanged:: 0.7
       Returns the vector and the corresponding atom selection. 
    
    This works slightly difference from :func:`sliceVector`. Mode array 
    (eigenvector) is multiplied by square-root of the variance along the mode.
    If mode is from an elastic network model, variance is defined as the 
    inverse of the eigenvalue.
    
    Note that returned :class:`Vector` instance is not normalized.
    
    :arg mode: mode instance to be sliced
    :type mode: :class:`Mode`
    
    :arg atoms: atoms for which *mode* describes a deformation, motion, etc.
    :type atoms: :class:`~prody.atomic.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`Vector`, :class:`~prody.atomic.Selection`)
    
    """
    
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, not {0:s}'
                        .format(type(mode)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.getNumOfAtoms() != mode.getNumOfAtoms(): 
        raise ValueError('number of atoms in *mode* and *atoms* must be equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getActiveCoordsetIndex())
    vec = Vector(mode.getArrayNx3()[
                 which,:].flatten() * mode.getVariance()**0.5,
                 '{0:s} slice "{1:s}"'.format(str(mode), selstr), 
                 mode.is3d()) 
    return (vec, sel)

def sliceModel(model, atoms, selstr):
    """|new| Return a slice of *model* matching *atoms* specified by *selstr*.
    
    .. versionadded:: 0.7

    Note that sliced normal modes (eigenvectors) are not normalized.
    
    :arg mode: NMA model instance to be sliced
    :type mode: :class:`NMABase`
    
    :arg atoms: atoms for which the *model* was built
    :type atoms: :class:`~prody.atomic.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: (:class:`NMA`, :class:`~prody.atomic.Selection`)
    
    """
    
    if not isinstance(model, NMABase):
        raise TypeError('mode must be a NMABase instance, not {0:s}'
                        .format(type(model)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if atoms.getNumOfAtoms() != model.getNumOfAtoms(): 
        raise ValueError('number of atoms in *model* and *atoms* must be '
                         'equal')
    
    array = model.getArray()
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getActiveCoordsetIndex())

    nma = type(model)('{0:s} slice "{1:s}"'.format(model.getName(), selstr))
    if model.is3d():
        which = [which.reshape((len(which),1))*3]
        which.append(which[0]+1)
        which.append(which[0]+2)
        which = np.concatenate(which, 1).flatten()
    nma.setEigens( array[which, :], model.getEigenvalues() )
    return (nma, sel)
    
def reduceModel(model, atoms, selstr):
    """Return reduced NMA model.
    
    .. versionchanged:: 0.7
       Returns the reduced model and the corresponding atom selection. 
    
    Reduces a :class:`NMA` model to a subset of *atoms* matching a selection 
    *selstr*.
    
    This function behaves depending on the type of the model.
    
    For ANM and GNM or other NMA models, this functions derives the force 
    constant matrix for system of interest (specified by the *selstr*) from 
    the force constant matrix for the *model* by assuming that for any given 
    displacement of the system of interest, the other atoms move along in 
    such a way as to minimize the potential energy. This is based on the
    formulation in in [KH00]_.
       
    For PCA models, this function simply takes the sub-covariance matrix for 
    the selected atoms.

    :arg model: dynamics model
    :type model: :class:`ANM`, :class:`GNM`, or :class:`PCA`
    :arg atoms: atoms that were used to build the model
    :arg selstr: a selection string specifying subset of atoms  
       
    """
    
    if linalg is None:
        prody.importLA()

    if not isinstance(model, NMABase):
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

    system = SELECT.getBoolArray(atoms, selstr)
    other = np.invert(system)
    n_sel = sum(system) 
    if n_sel == 0:
        LOGGER.warning('selection has 0 atoms')
        return None
    if len(atoms) == n_sel:
        LOGGER.warning('selection results in same number of atoms, '
                       'model is not reduced')
        return None

    if model.is3d():
        system = np.tile(system, (3,1)).transpose().flatten()
        other = np.tile(other, (3,1)).transpose().flatten()
    ss = matrix[system,:][:,system]
    if isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(ss)
        return eda, system
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(linalg.inv(oo), os))
    
    if isinstance(model, GNM):
        gnm = GNM(model.getName()+' reduced')
        gnm.setKirchhoff(matrix)
        return gnm, system
    elif isinstance(model, ANM):
        anm = ANM(model.getName()+' reduced')
        anm.setHessian(matrix)
        return anm, system
    elif isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(matrix)
        return eda, system

def writeModes(filename, modes, format='%.18e', delimiter=' '):
    """Write *modes* (eigenvectors) into a plain text file with name *filename*.
    
    .. versionchanged:: 0.5.3
       A compressed file is not outputted.
    
    See also :func:`writeArray`.
        
    >>> writeModes('p38_pca_modes_1-3.txt', p38_pca[:3])
    'p38_pca_modes_1-3.txt'

    """
    
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    return writeArray(filename, modes.getArray(), format=format, delimiter=delimiter)

def parseModes(normalmodes, eigenvalues=None, nm_delimiter=None, nm_skiprows=0, 
               nm_usecols=None, ev_delimiter=None, ev_skiprows=0, ev_usecols=None, 
               ev_usevalues=None):
    """Return :class:`NMA` instance with normal modes parsed from *normalmodes*.
    
    .. versionadded:: 0.5.3
    
    In normal mode file *normalmodes*, columns must correspond to  
    modes (eigenvectors).

    Optionally, *eigenvalues* can be parsed from a separate file. If 
    eigenvalues are not provided, they will all be set to 1.
    
    
    :arg normalmodes: File or filename that contains normal modes. 
        If the filename extension is :file:`.gz` or :file:`.bz2`, the file is 
        first decompressed.
    :type normalmodes: str or file
    
    :arg eigenvalues: Optional, file or filename that contains eigenvalues. 
        If the filename extension is :file:`.gz` or :file:`.bz2`, 
        the file is first decompressed.
    :type eigenvalues: str or file

    :arg nm_delimiter: The string used to separate values in *normalmodes*. 
        By default, this is any whitespace.
    :type nm_delimiter: str

    :arg nm_skiprows: Skip the first *skiprows* lines in *normalmodes*. 
        Default is ``0``.
    :type nm_skiprows: 0

    :arg nm_usecols: Which columns to read from *normalmodes*, with 0 being the 
        first. For example, ``usecols = (1,4,5)`` will extract the 2nd, 5th and 
        6th columns. The default, ``None``, results in all columns being read.
    :type nm_usecols: list

    :arg ev_delimiter: The string used to separate values in *eigenvalues*. 
        By default, this is any whitespace.
    :type ev_delimiter: str

    :arg ev_skiprows: Skip the first *skiprows* lines in *eigenvalues*. 
        Default is ``0``.
    :type ev_skiprows: 0

    :arg ev_usecols: Which columns to read from *eigenvalues*, with 0 being the 
        first. For example, ``usecols = (1,4,5)`` will extract the 2nd, 5th and 
        6th columns. The default, ``None``, results in all columns being read.
    :type ev_usecols: list

    :arg ev_usevalues: Which columns to use after the eigenvalue column is
        parsed from *eigenvalues*, with 0 being the first. 
        This can be used if *eigenvalues* contains more values than the
        number of modes in *normalmodes*.
    :type ev_usevalues: list
    
    See :func:`parseArray` for details of parsing arrays from files.
    
    """
    
    modes = parseArray(normalmodes, delimiter=nm_delimiter, 
                       skiprows=nm_skiprows, usecols=nm_usecols)
    if eigenvalues is not None:
        values = parseArray(eigenvalues, delimiter=ev_delimiter, 
                            skiprows=ev_skiprows, usecols=ev_usecols)
        values = values.flatten()
        if ev_usevalues is not None:
            values = values[ev_usevalues]
    nma = NMA(os.path.splitext(os.path.split(normalmodes)[1])[0])
    nma.setEigens(modes, values)
    return nma
    
    
def writeArray(filename, array, format='%d', delimiter=' '):
    """Write 1-d or 2-d array data into a delimited text file.
    
    .. versionchanged:: 0.5.3
       A compressed file is not outputted.
       
    This function is using :func:`numpy.savetxt` to write the file, after 
    making some type and value checks.
    
    Default *format* argument is ``"%d"``.
    Default *delimiter* argument is white space, ``" "``.
    
    *filename* will be returned upon successful writing. 
     
    >>> writeArray('p38_cross-correlations.txt', calcCrossCorrelations(p38_pca))
    'p38_cross-correlations.txt'
    
    """
    if not isinstance(array, np.ndarray):
        raise TypeError('array must be a Numpy ndarray, not {0:s}'
                        .format(type(array)))
    elif not array.ndim in (1, 2):
        raise ValueError('array must be a 1 or 2-dimensional Numpy ndarray, '
                         'not {0:d}-d'.format(type(array.ndim)))
    np.savetxt(filename, array, format, delimiter)
    return filename

def parseArray(filename, delimiter=None, skiprows=0, usecols=None,
               dtype=np.float64):
    """Parse array data from a file.
    
    .. versionadded:: 0.5.3
    
    This function is using :func:`numpy.loadtxt` to parse the file.
    
    Each row in the text file must have the same number of values.
    
    :arg filename: File or filename to read. If the filename extension is 
        :file:`.gz` or :file:`.bz2`, the file is first decompressed.
    :type filename: str or file
    
    :arg delimiter: The string used to separate values. By default, 
        this is any whitespace.
    :type delimiter: str
    
    :arg skiprows: Skip the first *skiprows* lines, default is ``0``.
    :type skiprows: int
     
    :arg usecols: Which columns to read, with 0 being the first. For example, 
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns. 
        The default, ``None``, results in all columns being read.
    :type usecols: list
    
    :arg dtype: Data-type of the resulting array, default is 
        :class:`numpy.float64`. 
    :type dtype: :class:`numpy.dtype`.
    
    """

    array = np.loadtxt(filename, dtype=dtype, delimiter=delimiter, 
                       skiprows=skiprows, usecols=usecols)
    return array
        
def parseSparseMatrix(filename, symmetric=False, delimiter=None, skiprows=0,
                      irow=0, icol=1, first=1):
    """|new| Parse sparse matrix data from a file.
    
    .. versionadded:: 0.7
    
    This function is using :func:`parseArray` to parse the file.
    
    Input must have the following format::
        
       1       1    9.958948135375977e+00
       1       2   -3.788214445114136e+00
       1       3    6.236155629158020e-01
       1       4   -7.820609807968140e-01
    
    Each row in the text file must have the same number of values.
    
    :arg filename: File or filename to read. If the filename extension is 
        :file:`.gz` or :file:`.bz2`, the file is first decompressed.
    :type filename: str or file
    
    :arg symmetric: Set ``True`` if the file contains triangular part of a 
        symmetric matrix, default is ``False``.
    :type symmetric: bool
    
    :arg delimiter: The string used to separate values. By default, 
        this is any whitespace.
    :type delimiter: str
    
    :arg skiprows: Skip the first *skiprows* lines, default is ``0``.
    :type skiprows: int
    
    :arg irow: Index of the column in data file corresponding to row indices,
        default is ``0``. 
    :type irow: int 
        
    :arg icol: Index of the column in data file corresponding to row indices,
        default is ``0``. 
    :type icol: int
    
    :arg first: First index in the data file (0 or 1), default is ``1``. 
    :type first: int

    Data-type of the resulting array, default is :class:`numpy.float64`. 

    """
    irow = int(irow)
    icol = int(icol)
    first = int(first)
    assert 0 <= irow <= 2 and 0 <= icol <= 2, 'irow/icol may be 0, 1, or 2'
    assert icol != irow, 'irow and icol must not be equal' 
    idata = [0, 1, 2]
    idata.pop(idata.index(irow))
    idata.pop(idata.index(icol))
    idata = idata[0]
    sparse = parseArray(filename, delimiter, skiprows)
    dof = sparse[:,[irow, icol]].max() 
    matrix = np.zeros((dof,dof))
    irow = (sparse[:,irow] - first).astype(np.int64)
    icol = (sparse[:,icol] - first).astype(np.int64)
    matrix[irow, icol] = sparse[:,idata]
    if symmetric:
        matrix[icol, irow] = sparse[:,idata]
    return matrix

def sampleModes(modes, atoms=None, n_confs=1000, rmsd=1.0):
    """Return an ensemble of randomly sampled conformations along given *modes*.
    
    If *atoms* are provided, sampling will be around its active coordinate set. 
    Otherwise, sampling is around the 0 coordinate set.
    
    :arg modes: Modes along which sampling will be performed.
    :type modes: :class:`Mode`, :class:`ModeSet`, :class:`PCA`, :class:`ANM` 
                 or :class:`NMA`   
    
    :arg atoms: Atoms whose active coordinate set will be used as the initial 
        conformation.
    :type atoms: :class:`~prody.atomic.Atomic`  
    
    :arg n_confs: Number of conformations to generate. Default is 1000.
    :type n_steps: int 
    
    :arg rmsd: The average RMSD that the conformations will have with 
        respect to the initial conformation. Default is 1.0 A.
    :type rmsd: float 
    
    For given normal modes :math:`[u_1 u_2 ... u_m]` and their eigenvalues
    :math:`[\lambda_1 \lambda_2 ... \lambda_m]`, a new conformation 
    is sampled using the relation:
        
    .. math::
    
       R_k = R_0 + s \sum_{i=1}^{m} r_i^k \lambda^{-0.5}_i u_i 
    
    :math:`R_0` is the active coordinate set of *atoms*.
    :math:`[r_1^k r_2^k ... r_m^k]` are normally distributed random numbers 
    generated for conformation :math:`k` using :func:`numpy.random.randn`.
    
    RMSD of the new conformation from :math:`R_0` can be calculated as
     
    .. math::
        
      RMSD^k = \sqrt{ {\\left( s \sum_{i=1}^{m} r_i^k \lambda^{-0.5}_i u_i  \\right)}^{2} / N } = \\frac{s}{ \sqrt{N}} \sqrt{ \sum_{i=1}^{m} (r_i^k)^2 \lambda^{-1}_i  } 


    Average :math:`RMSD` of the generated conformations from the initial conformation is: 
        
    .. math::
        
      \\left< RMSD^k \\right> = \\frac{s}{ \sqrt{N}} \\left< \sqrt{ \sum_{i=1}^{m} (r_i^k)^2 \lambda^{-1}_i } \\right>

 
    From this relation :math:`s` scaling factor obtained using the relation 
    
    .. math::
       
       s =  \\left< RMSD^k \\right> \sqrt{N} {\\left< \sqrt{ \sum_{i=1}^{m} (r_i)^2 \lambda^{-1}_i} \\right>}^{-1}
       
     
    Note that random numbers are generated before conformations are 
    sampled, hence exact value of :math:`s` is known from this relation to
    ensure that the generated ensemble will have user given average *rmsd* value. 
     
    Note that if modes are from a :class:`PCA`, variances are used instead of 
    inverse eigenvalues, i.e. :math:`\sigma_i \sim \lambda^{-1}_i`.
    
    |more| See also :func:`showEllipsoid`.
    
    .. plot::
       :context:
       :include-source:
        
       # Generate 300 conformations using ANM modes 1-3
       ensemble = sampleModes( p38_anm[:3], n_confs=500 )
       # Project these conformations onto the space spanned by these modes
       plt.figure(figsize=(5,4))
       showProjection(ensemble, p38_anm[:3])
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    if not isinstance(modes, (Mode, NMABase, ModeSet)):
        raise TypeError('modes must be a NMA or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be from a 3-dimensional model')
    n_confs = int(n_confs)
    n_atoms = modes.getNumOfAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, (Atomic)):
            raise TypeError('{0:s} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.getNumOfAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoordinates()

    rmsd = float(rmsd)
    LOGGER.info('Parameter: rmsd = {0:.2f} A'.format(rmsd))
    n_confs = int(n_confs)
    LOGGER.info('Parameter: n_confs = {0:d}'.format(n_confs))
    
    
    if isinstance(modes, Mode):
        n_modes = 1
        variances = np.array([modes.getVariance()])
    else:
        n_modes = len(modes)
        variances = modes.getVariances()
    if np.any(variances == 0):
        raise ValueError('one or more modes has zero variance')
    randn = np.random.standard_normal((n_confs, n_modes))
    coef = ((randn ** 2 * variances).sum(1) ** 0.5).mean()
    scale = n_atoms**0.5 * rmsd / coef
    
    LOGGER.info('Modes are scaled by {0:g}.'.format(scale))
    
    confs = []
    append = confs.append
    scale = scale * variances ** 0.5
    array = modes.getArray()
    if array.ndim > 1:
        for i in range(n_confs):
            append( (array * scale * randn[i]).sum(1).reshape((n_atoms, 3)) )
    else:
        for i in range(n_confs):
            append( (array * scale * randn[i]).reshape((n_atoms, 3)) )

    ensemble = Ensemble('Conformations along {0:s}'.format(modes))
    if initial is None:
        ensemble.setCoordinates(np.zeros((n_atoms, 3)))
        ensemble.addCoordset(np.array(confs))
    else:    
        ensemble.setCoordinates(initial)
        ensemble.addCoordset(np.array(confs) + initial)
    return ensemble  

def showEllipsoid(modes, onto=None, n_std=2, scale=1., *args, **kwargs):
    """Show an ellipsoid using :meth:`~mpl_toolkits.mplot3d.Axes3D.plot_wireframe`.
    
    Ellipsoid volume gives an analytical view of the conformational space that
    given modes describe.
    
    :arg modes: 3 modes for which ellipsoid will be drawn.
    :type modes: :class:`ModeSet`, :class:`PCA`, :class:`ANM` or :class:`NMA`
    
    :arg onto: 3 modes onto which ellipsoid will be projected.
    :type modes: :class:`ModeSet`, :class:`PCA`, :class:`ANM` or :class:`NMA`
       
    :arg n_std: Number of standard deviations to scale the ellipsoid.
    :type n_std: float
    
    :arg scale: Used for scaling the volume of ellipsoid. This can be
        obtaine from :func:`sampleModes`.
    :type scale: float


    .. plot::
       :context:
       :include-source:
        
       # Show projection of subspace spanned by ANM 1-3 onto subspace of PC 1-3 
       plt.figure(figsize=(5,4))
       showEllipsoid(p38_anm[:3], p38_pca[:3])
       
       # Let's compare this with that of ANM modes 18-20
       showEllipsoid(p38_anm[17:], p38_pca[:3], color='red')
       # This ANM subspace appears as a tiny volume at the center
       # since faster ANM modes does not correspond to top ranking PCA modes
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')

    """
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(modes, (NMABase, ModeSet)):
        raise TypeError('modes must be a NMA or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if not modes.is3d():
        raise ValueError('modes must be from a 3-dimensional model')
    if len(modes) != 3:
        raise ValueError('length of modes is not equal to 3')
    if onto is not None:
        if not isinstance(onto, (NMABase, ModeSet)):
            raise TypeError('onto must be a NMA or ModeSet instance, '
                            'not {0:s}'.format(type(onto)))
        if not onto.is3d():
            raise ValueError('onto must be from a 3-dimensional model')
        if len(onto) != 3:
            raise ValueError('length of onto is not equal to 3')
        if onto.getNumOfAtoms() != modes.getNumOfAtoms():
            raise ValueError('modes and onto does not have same number of atoms')
        
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    
    var = modes.getVariances()
    #randn = np.random.standard_normal((1000, 3))
    #coef = ((randn ** 2 * var).sum(1) ** 0.5).mean()
    #scale = float(n_std) * modes.getNumOfAtoms()**0.5 * float(rmsd) / coef * var ** 0.5
    scale = float(n_std) * scale * var ** 0.5
    #scale = float(n_std) * modes.getNumOfAtoms() ** 0.5 * float(rmsd) / var.sum() ** 0.5 * var ** 0.5   

    x = scale[0] * np.outer(np.cos(u), np.sin(v))
    y = scale[1] * np.outer(np.sin(u), np.sin(v))
    z = scale[2] * np.outer(np.ones(np.size(u)), np.cos(v))
    if onto is not None:
        change_of_basis = np.dot(modes.getArray().T, onto.getArray())

        xyz = np.array([x.flatten(), y.flatten(), z.flatten()])
        xyz = np.dot(xyz.T, change_of_basis)
        x = xyz[:,0].reshape((100,100))
        y = xyz[:,1].reshape((100,100))
        z = xyz[:,2].reshape((100,100))

    cf = plt.gcf()
    show = None
    for child in cf.get_children():
        if isinstance(child, Axes3D):
            show = child
            break 
    if show is None:
        show = Axes3D(cf)
    show.plot_wireframe(x, y, z, rstride=6, cstride=6, *args, **kwargs)
    if onto is not None:
        onto = onto.getModes()
        show.set_xlabel('Mode {0:d} coordinate'.format(onto[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(onto[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(onto[2].getIndex()+1))
    else:
        modes = modes.getModes()
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    return show

def traverseMode(mode, atoms, n_steps=10, rmsd=1.5):
    """Generates a trajectory along a given *mode*, which can be used to
    animate fluctuations in an external program. 
    
    :arg mode: Mode along which a trajectory will be generated.
    :type mode: :class:`Mode`   
    
    :arg atoms: Atoms whose active coordinate set will be used as the initial 
        conformation.
    :type atoms: :class:`~prody.atomic.Atomic` 
    
    :arg n_steps: Number of steps to take along each direction. 
        For example, for ``n_steps=10``, 20 conformations will be 
        generated along the first mode. Default is 10.
    :type n_steps: int 
    
    :arg rmsd: The maximum RMSD that the conformations will have with 
        respect to the initial conformation. Default is 1.5 A.
    :type rmsd: float

    :returns: :class:`~prody.ensemble.Ensemble`
    
    For given normal mode :math:`u_i`, its eigenvalue
    :math:`\lambda_i`, number of steps :math:`n`, and maximum :math:`RMSD`  
    conformations :math:`[R_{-n} R_{-n+1} ... R_{-1} R_0 R_1 ... R_n]` are 
    generated.
    
    :math:`R_0` is the active coordinate set of *atoms*. 
    :math:`R_k = R_0 + sk\lambda_iu_i`, where :math:`s` is found using
    :math:`s = ((N (\\frac{RMSD}{n})^2) / \lambda_i^{-1}) ^{0.5}`, where
    :math:`N` is the number of atoms.
    
    
    .. plot::
       :context:
       :include-source:
        
       trajectory = traverseMode( p38_anm[0], p38_structure.select('calpha'), 
                                  n_steps=8, rmsd=1.4 )
       rmsd = calcRMSD(trajectory)
       plt.figure(figsize=(5,4))
       plt.plot(rmsd, '-o')
       plt.xlabel('Frame index')
       plt.ylabel('RMSD (A)')
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    """
    if not isinstance(mode, VectorBase):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0:s}'.format(type(mode)))
    if not mode.is3d():
        raise ValueError('mode must be from a 3-dimensional model.')
    n_atoms = mode.getNumOfAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, Atomic):
            raise TypeError('{0:s} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.getNumOfAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoordinates()

    name = str(mode)
    
    rmsd = float(rmsd) + 0.000004
    LOGGER.info('Parameter: rmsd = {0:.2f} A'.format(rmsd))
    n_steps = int(n_steps)
    LOGGER.info('Parameter: n_steps = {0:d}'.format(n_steps))
    step = rmsd / n_steps
    LOGGER.info('Step size is {0:.2f} A RMSD'.format(step))
    arr = mode.getArrayNx3()
    var = mode.getVariance()
    scale = ((n_atoms * step**2) / var) **0.5
    LOGGER.info('Mode is scaled by {0:g}.'.format(scale))

    array = arr * var**0.5 * scale 
    confs_add = [initial + array]
    for s in range(1, n_steps):
        confs_add.append( confs_add[-1] + array)
    confs_sub = [initial - array]
    for s in range(1, n_steps):
        confs_sub.append( confs_sub[-1] - array)
    confs_sub.reverse()
    ensemble = Ensemble('Conformations along {0:s}'.format(name))
    ensemble.setCoordinates(initial)    
    ensemble.addCoordset(np.array(confs_sub + [initial] + confs_add))
    return ensemble
  
def deform(atoms, mode, rmsd=None):  
    """|new| Generate a new coordinate set for *atoms* along the *mode*.
    
    .. versionadded:: 0.7
    
    *atoms* must be a :class:`~prody.atomic.AtomGroup` instance.
    New coordinate set will be appended to *atoms*. If *rmsd* is provided,
    *mode* will be scaled to generate a coordinate set with given RMSD distance
    to the active coordinate set.
    
    Below example shows how to deform a structure along a normal mode
    or linear combinations of normal modes:
    
    >>> deform(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5)
    >>> deform(p38_structure, -p38_pca[1] * p38_pca[1].getVariance()**0.5)
    >>> deform(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5 + 
    ...                       p38_pca[1] * p38_pca[1].getVariance()**0.5)
    >>> deform(p38_structure, p38_pca[0], rmsd=1.0)
    >>> calcRMSD(p38_structure)
    array([ 0.   ,  0.41 ,  0.308,  0.513,  1.   ])
    
    """

    if not isinstance(atoms, AtomGroup):
        raise TypeError('atoms must be an AtomGroup, not {0:s}'
                        .format(type(atoms)))
    if not isinstance(mode, VectorBase):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0:s}'.format(type(mode)))
    if not mode.is3d():
        raise ValueError('mode must be from a 3-dimensional model.')
    if atoms.getNumOfAtoms() != mode.getNumOfAtoms():
        raise ValueError('number of atoms do not match')
    
    array = mode.getArrayNx3()
    
    if rmsd is not None:
        rmsd = float(rmsd)
        # rmsd = ( ((scalar * array)**2).sum() / n_atoms )**0.5
        scalar = (atoms.getNumOfAtoms() * rmsd**2 / (array**2).sum())**0.5
        LOGGER.info('Mode is scaled by {0:g}.'.format(scalar))
        atoms.addCoordset( atoms.getCoordinates() + array * scalar)
    else:     
        atoms.addCoordset( atoms.getCoordinates() + array)
    
def calcSqFlucts(modes):
    """Return sum of square-fluctuations for given set of normal *modes*.
    
    .. versionchanged:: 0.7.1
       :class:`Vector` instances are accepted as *modes* argument.
    
    >>> calcSqFlucts(p38_pca) # doctest: +SKIP 
    array([  0.94163178,   0.97486815,   0.81056074,   0.59926465,
             0.80505867,   0.93568339,   1.1634173 ,   1.39873827,
             ...,
             0.17928521,   0.12930389,   0.21500264,   0.16070006,   
             0.21077809])
             
    """
    
    if not isinstance(modes, (VectorBase, NMABase, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Vector):
        if self.is3d():
            return (self.getArrayNx3()**2).sum(axis=1)
        else:
            return (self.getArray() ** 2)
    else:
        square_fluctuations = np.zeros(modes.getNumOfAtoms()) 
        if isinstance(modes, VectorBase):
            modes = [modes]
        for mode in modes:
            square_fluctuations += mode.getSqFlucts()
        return square_fluctuations
 
def calcCrossCorrelations(modes, n_cpu=1):
    """Return cross-correlations matrix.
    
    For a 3-d model, cross-correlations matrix is an NxN matrix, where N is the 
    number of atoms. Each element of this matrix is the trace of the 
    submatrix corresponding to a pair of atoms.
    
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.

    For large systems, calculation of cross-correlations matrix may be time 
    consuming. Optionally, multiple processors may be employed to perform
    calculations by passing ``n_cpu=2`` or more. 

    :arg n_cpu: Number of CPUs to use. Default is 1. 
    :type n_cpu: int 
    
    >>> calcCrossCorrelations(p38_anm) # doctest: +SKIP
    array([[ 1.   ,  0.948,  0.874, ...,  0.604,  0.552,  0.32 ],
           [ 0.948,  1.   ,  0.97 , ...,  0.517,  0.532,  0.317],
           [ 0.874,  0.97 ,  1.   , ...,  0.323,  0.353,  0.145],
           ..., 
           [ 0.604,  0.517,  0.323, ...,  1.   ,  0.956,  0.865],
           [ 0.552,  0.532,  0.353, ...,  0.956,  1.   ,  0.953],
           [ 0.32 ,  0.317,  0.145, ...,  0.865,  0.953,  1.   ]])

    """
    
    if not isinstance(n_cpu, int):
        raise TypeError('n_cpu must be an integer')
    elif n_cpu < 1:
        raise ValueError('n_cpu must be equal to or greater than 1')
        
    if not isinstance(modes, (Mode, NMABase, ModeSet)):
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


def calcCumulativeOverlap(modes1, modes2):
    """Return cumulative overlap of modes in *modes2* with those in *modes1*.
    
    Returns a number of *modes1* contains a single :class:`Mode` or a 
    :class:`Vector` instance. If *modes1* contains multiple modes, returns an
    array. Elements of the array correspond to cumulative overlaps for modes 
    in *modes1* with those in *modes2*.
    
    >>> calcCumulativeOverlap(p38_pca[0], p38_anm) # doctest: +SKIP
    0.89357715148440353
    
    """
    
    overlap = calcOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    return cumov

def calcCumulativeOverlapArray(modes1, modes2):
    """Return array of cumulative overlaps.
   
    Returned array has the shape ``(len(modes1), len(modes2))``. Each row
    corresponds to cumulative overlaps calculated for modes in *modes1* with
    those in *modes2*. Each value in a row corresponds to cumulative overlap
    calculated using upto that many number of modes from *modes2*.
    
    >>> calcCumulativeOverlapArray(p38_pca[0], p38_anm) # doctest: +SKIP
    array([ 0.38956674,  0.39199772,  0.81423637,  0.81953805,  0.82296797,
            0.83320235,  0.83947803,  0.85141434,  0.87426872,  0.88529856,
            0.88647576,  0.88878561,  0.88882618,  0.88894772,  0.89121041,
            0.89162571,  0.89276304,  0.89336947,  0.8934872 ,  0.89357715])
    
    """
    
    overlap = calcOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))
    return cumov


def calcSubspaceOverlap(modes1, modes2):
    """Return subspace overlap between two sets of modes (*modes1* and *modes2*).
    
    Also known as the root mean square inner product (RMSIP) of essential 
    subspaces [AA99]_.
    
    This function returns a single number.
    
    >>> calcSubspaceOverlap(p38_pca[:3], p38_anm[:3]) # doctest: +SKIP
    0.75174192334469792
    """
    overlap = calcOverlap(modes1, modes2)
    if isinstance(modes1, Mode):
        length = 1
    else:
        length = len(modes1)
    rmsip = np.sqrt(np.power(overlap, 2).sum() / length)
    return rmsip


def calcCovarianceOverlap(modelA, modelB):
    """Return overlap between covariances of *modelA* and *modelB*.
    
    .. versionadded:: 0.5.3
    
    Overlap between covariances are calculated using normal modes 
    (eigenvectors), hence modes in both models must have been calculated.
    
    This function implements equation 11 in [BH02]_.
    
    """
    if not modelA.is3d() or not modelB.is3d(): 
        raise TypeError('both models must be 3-dimensional') 
    if len(modelA) == 0 or len(modelB) == 0:  
        raise TypeError('modes must be calculated for both models, '
                        'try calcModes method')
    if modelA.getNumOfAtoms() != modelB.getNumOfAtoms(): 
        raise ValueError('modelA and modelB must have same number of atoms')
    arrayA = modelA.getArray()
    varA = modelA.getVariances()
    arrayB = modelB.getArray()
    varB = modelB.getVariances()
    
    dotAB = np.dot(arrayA.T, arrayB)**2
    outerAB = np.outer(varA**0.5, varB**0.5)
    diff = (np.sum(varA.sum() + varB.sum()) - 2 * np.sum(outerAB * dotAB))
    if diff < ZERO:
        diff = 0
    else:
        diff = diff ** 0.5
    return 1 - diff / np.sqrt(varA.sum() + varB.sum())

def calcCovariance(modes):
    """Return covariance matrix calculated for given *modes*."""
    
    return modes.getCovariance()

def calcTempFactors(modes, atoms):
    """Return temperature (β) factors calculated using *modes* from a 
    :class:`ANM` or :class:`GNM` instance scaled according to the experimental 
    β factors from *atoms*."""
    
    model = modes.getModel()
    if not isinstance(model, GNMBase):
        raise TypeError('modes must come from GNM or ANM')
    if model.getNumOfAtoms() != atoms.getNumOfAtoms():
        raise ValueError('modes and atoms must have same number of modes')
    sqf = calcSqFlucts(modes)
    return sqf / ((sqf**2).sum()**0.5) * (
                                        (atoms.getTempFactors()**2).sum()**0.5)
    
def showFractOfVariances(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.bar`.
    
    Note that mode indices are incremented by 1.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showFractOfVariances(p38_pca) 
       showCumFractOfVariances(p38_pca)
      
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(modes, (ModeSet, NMABase)):
        raise TypeError('modes must be NMABase, or ModeSet, not {0:s}'.format(type(modes)))
    
    fracts = [(mode.getIndex(), mode.getFractOfVariance()) for mode in modes]
    fracts = np.array(fracts)
    show = plt.bar(fracts[:,0]+0.5, fracts[:,1], *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    return show

def showCumFractOfVariances(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.plot`.
    
    Note that mode indices are incremented by 1.
    See :func:`showFractOfVariances` for an example.
    
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(modes, (Mode, NMABase, ModeSet)):
        raise TypeError('modes must be a Mode, NMABase, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Mode):
        indices = modes.getIndices() + 0.5
        modes = [modes]
    elif isinstance(modes, ModeSet):
        indices = modes.getIndices() + 0.5
    else:
        indices = np.arange(len(modes)) + 0.5
    fracts = np.array([mode.getFractOfVariance() for mode in modes]).cumsum()
    show = plt.plot(indices, fracts, *args, **kwargs)
    axis = list(plt.axis())
    axis[0] = 0.5
    axis[2] = 0
    axis[3] = 1
    plt.axis(axis)
    plt.xlabel('Mode index')
    plt.ylabel('Fraction of variance')
    return show

def showProjection(ensemble, modes, *args, **kwargs):
    """Show projection of conformational deviations onto given modes.
    
    :arg ensemble: a :class:`~prody.ensemble.Ensemble` instance
    
    Matplotlib function used for plotting depends on the number of modes:
        
      * 1 mode: :func:`~matplotlib.pyplot.hist`
      * 2 modes: :func:`~matplotlib.pyplot.plot`
      * 3 modes: :meth:`mpl_toolkits.mplot3d.Axes3D.plot`
      
    By default ``marker='o', ls='None'`` is passed to the plotting function 
    to disable lines in projections onto 2 or 3-d spaces.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[0]) 
       plt.title('Projection onto PC1')

       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[:2])
       plt.title('Projection onto PC1-2')
       
       plt.figure(figsize=(5,4))
       showProjection(p38_ensemble, p38_pca[:3]) # onto top 3 PCs
       plt.title('Projection onto PC1-3')

       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'
                        .format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    if isinstance(modes, Mode) or (isinstance(modes, (ModeSet, NMABase)) and 
                                   len(modes)==1):
        if not isinstance(modes, Mode):
            modes = modes[0]
        projection = calcProjection(ensemble, modes)
        show = plt.hist(projection.flatten(), *args, **kwargs)
        plt.xlabel('Mode {0:d} coordinate'.format(modes.getIndex()+1))
        plt.ylabel('Number of conformations')
    elif len(modes) == 2:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        projection = calcProjection(ensemble, modes)
        show = plt.plot(projection[:,0], projection[:,1], *args, **kwargs)
        modes = [m for m in modes]
        plt.xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        plt.ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
    elif len(modes) == 3:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        projection = calcProjection(ensemble, modes)
        modes = [m for m in modes]
        cf = plt.gcf()
        show = None
        for child in cf.get_children():
            if isinstance(child, Axes3D):
                show = child
                break 
        if show is None:
            show = Axes3D(cf)
        show.plot(projection[:, 0], projection[:, 1], projection[:, 2], 
                  *args, **kwargs)
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    else:
        raise ValueError('Projection onto upto 3 modes can be shown. '
                         'You have given {0:d} mode.'.format(len(modes)))
    
    return show

def showCrossProjection(ensemble, mode_x, mode_y, scale=None, scalar=None, 
                        *args, **kwargs):
    """Show projection of conformational deviations using 
    :func:`~matplotlib.pyplot.plot`.
    
    This function is differs from :func:`showProjection` by accepting modes
    from two different models.
    
    :arg ensemble: Ensemble for which deviations will be projected
    :type ensemble: :class:`~prody.ensemble.Ensemble`
    :arg mode_x: Projection onto this mode will be shown along x-axis. 
    :type mode_x: :class:`Mode`
    :arg mode_y: Projection onto this mode will be shown along y-axis.
    :type mode_y: :class:`Mode`
    :arg scale: Scale width of the projection onto one of modes. 
                ``x`` and ``y`` are accepted.
    :type scale: str
    
    By default ``marker='o', ls='None'`` is passed to the plotting function 
    to disable lines.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showCrossProjection(p38_ensemble, p38_pca[0], p38_anm[2])
    
    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    |example| See :ref:`pca-xray-plotting` for a more elaborate example.
       
    """
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'
                        .format(type(ensemble)))
    if not isinstance(mode_x, VectorBase):
        raise TypeError('mode_x must be a Mode instance, not {0:s}'
                        .format(type(mode_x)))
    if not mode_x.is3d():
        raise ValueError('mode_x must be 3-dimensional')
    if not isinstance(mode_y, VectorBase):
        raise TypeError('mode_y must be a Mode instance, not {0:s}'
                        .format(type(mode_y)))
    if not mode_y.is3d():
        raise ValueError('mode_y must be 3-dimensional')
    
    xcoords = calcProjection(ensemble, mode_x) 
    ycoords = calcProjection(ensemble, mode_y)
    if isinstance(scale, str) and scale.lower() in ('x', 'y'):
        if scalar is not None:
            scalar = float(scalar)
        else:
            xmin = xcoords.min()
            xmax = xcoords.max()
            ymin = ycoords.min()
            ymax = ycoords.max()
            scalar = ((ymax - ymin) / (xmax - xmin)) 
            scalar = scalar * np.sign(calcOverlap(mode_x, mode_y))
            if scale == 'x':
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_x, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_y, scalar))
        if scale.lower() == 'x':
            xcoords = xcoords * scalar  
        else:
            ycoords = ycoords * scalar
    elif scale is not None:
        LOGGER.warning('{0:s} is not a valid value for scale argument. '
                       'Only "x" or "y" are accepted.'.format(str(scale)))
    if 'ls' not in kwargs:
        kwargs['ls'] = 'None'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'
    show = plt.plot(xcoords, ycoords, *args, **kwargs)
    plt.xlabel('{0:s} coordinate'.format(mode_x))
    plt.ylabel('{0:s} coordinate'.format(mode_y))
    return show

def showOverlapTable(rows, cols, *args, **kwargs):
    """Show overlap table using :func:`~matplotlib.pyplot.pcolor`.
    
    *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the displayed matrix.
    
    Note that mode indices are increased by 1. List of modes should contain
    a set of contiguous modes from the same model. 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showOverlapTable( p38_pca[:6], p38_anm[:6] )
       plt.title('p38 PCA vs ANM')

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all') 
       
    """
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(rows, (NMABase, ModeSet)):
        raise TypeError('rows must be an NMA model or a ModeSet, not {0:s}'
                        .format(type(rows)))
    if not isinstance(rows, (NMABase, ModeSet)):
        raise TypeError('cols must be an NMA model or a ModeSet, not {0:s}'
                        .format(type(cols)))
    overlap = abs(calcOverlap(rows, cols))
    if isinstance(rows, NMABase):
        rows = rows[:]
    if isinstance(cols, NMABase):
        cols = cols[:]
    show = plt.pcolor(overlap, cmap=plt.cm.jet, *args, **kwargs), plt.colorbar()
    x_range = np.arange(1, len(cols)+1)
    plt.xticks(x_range-0.5, x_range)
    plt.xlabel(str(cols))
    y_range = np.arange(1, len(rows)+1)
    plt.yticks(y_range-0.5, y_range)
    plt.ylabel(str(rows))
    plt.axis([0, len(cols), 0, len(rows)])
    return show

def showCrossCorrelations(modes, *args, **kwargs):
    """Show cross-correlations for given modes using 
    :func:`~matplotlib.pyplot.imshow`.
    
    See also :func:`getCrossCorrelations`. 
    
    By default, *origin=lower* and *interpolation=bilinear* keyword
    arguments are passed to imshow function. User can overwrite these
    parameters.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,5))
       # Show cross-correlations for ANM modes 1-3
       showCrossCorrelations( p38_anm[:3] )
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    arange = np.arange(modes.getNumOfAtoms())
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:, 
                       arange[0]+1:] = calcCrossCorrelations(modes)
    if not kwargs.has_key('interpolation'):
        kwargs['interpolation'] = 'bilinear'
    if not kwargs.has_key('origin'):
        kwargs['origin'] = 'lower'
    show = plt.imshow(cross_correlations, *args, **kwargs), plt.colorbar()
    plt.axis([arange[0]+0.5, arange[-1]+1.5, arange[0]+0.5, arange[-1]+1.5])
    plt.title('Cross-correlations for {0:s}'.format(str(modes))) 
    plt.xlabel('Indices')
    plt.ylabel('Indices')
    return show

def showMode(mode, *args, **kwargs):
    """Show mode array using :func:`~matplotlib.pyplot.plot`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,4))
       showMode( p38_anm[0] )
       plt.grid()
       plt.legend(loc='lower right')
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, '
                        'not {0:s}'.format(type(mode)))
    if mode.is3d():
        a3d = mode.getArrayNx3()
        show = plt.plot(a3d[:, 0], *args, label='x-component', **kwargs)
        plt.plot(a3d[:, 1], *args, label='y-component', **kwargs)
        plt.plot(a3d[:, 2], *args, label='z-component', **kwargs)
    else:
        show = plt.plot(mode.getArray(), *args, **kwargs)
    plt.title(str(mode))
    plt.xlabel('Indices')
    return show

def showSqFlucts(modes, *args, **kwargs):
    """Show square fluctuations using :func:`~matplotlib.pyplot.plot`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,4))
       showSqFlucts( p38_anm[0] )
       showSqFlucts( p38_anm[1] )
       plt.legend()
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    sqf = calcSqFlucts(modes)
    if not 'label' in kwargs:
        kwargs['label'] = str(modes) 
    show = plt.plot(sqf, *args, **kwargs)
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    plt.title(str(modes))
    return show

def showScaledSqFlucts(modes, *args, **kwargs):
    """Show scaled square fluctuations using :func:`~matplotlib.pyplot.plot`.
    
    Modes or mode sets given as additional arguments will be scaled to have
    the same mean squared fluctuations as *modes*. 
    
    .. plot::
       :context:
       :include-source:
       
       plt.figure(figsize=(5,4))
       showScaledSqFlucts(p38_pca[0], p38_anm[2])
       plt.legend()

    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    sqf = calcSqFlucts(modes)
    mean = sqf.mean()
    args = list(args)
    modesarg = []
    i = 0
    while i < len(args):
        if isinstance(args[i], (VectorBase, ModeSet, NMABase)):
            modesarg.append(args.pop(i))
        else:
            i += 1
    show = [plt.plot(sqf, *args, label=str(modes), **kwargs)]
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        scalar = mean / sqf.mean()
        show.append(plt.plot(sqf * scalar, *args, 
                           label='{0:s} (x{1:.2f})'.format(str(modes), scalar), 
                           **kwargs))
    return show

def showNormedSqFlucts(modes, *args, **kwargs):
    """Show normalized square fluctuations using 
    :func:`~matplotlib.pyplot.plot`.
    
    .. plot::
       :context:
       :include-source:
       
       plt.figure(figsize=(5,4))
       showNormedSqFlucts(p38_pca[0], p38_anm[2])
       plt.legend()

    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    sqf = calcSqFlucts(modes)
    args = list(args)
    modesarg = []
    i = 0
    while i < len(args):
        if isinstance(args[i], (VectorBase, ModeSet, NMABase)):
            modesarg.append(args.pop(i))
        else:
            i += 1
    show = [plt.plot(sqf/(sqf**2).sum()**0.5, *args, 
                        label='{0:s}'.format(str(modes)), **kwargs)]    
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        show.append(plt.plot(sqf/(sqf**2).sum()**0.5, *args, 
                    label='{0:s}'.format(str(modes)), **kwargs))
    return show

def showContactMap(enm, *args, **kwargs):
    """Show Kirchhoff matrix using :func:`~matplotlib.pyplot.spy`.
    
    .. plot::
       :context:
       :include-source:
        
       p38_gnm = GNM('p38')
       p38_gnm.buildKirchhoff( p38_structure )
       plt.figure(figsize=(4,4))
       showContactMap( p38_gnm )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all') 
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(enm, GNMBase):
        raise TypeError('model argument must be an ENM instance')
    kirchhoff = enm.getKirchhoff()
    if kirchhoff is None:
        LOGGER.warning('kirchhoff matrix is not set')
        return None
    show = plt.spy(kirchhoff, *args, **kwargs)
    plt.title('{0:s} contact map'.format(enm.getName())) 
    plt.xlabel('Residue index')
    plt.ylabel('Residue index')
    return show

def showOverlap(mode, modes, *args, **kwargs):
    """Show overlap :func:`~matplotlib.pyplot.bar`.
    
    :arg mode: a single mode/vector
    :type mode: :class:`Mode`, :class:`Vector` 
    :arg modes: multiple modes
    :type modes: :class:`ModeSet`, :class:`ANM`, :class:`GNM`, or :class:`PCA` 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(4,4))
       showOverlap( p38_pca[0], p38_anm[:6] )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all') 
       
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be Mode or Vector, not {0:s}'
                        .format(type(mode)))
    if not isinstance(modes, (NMABase, ModeSet)):
        raise TypeError('modes must be NMA or ModeSet, not {0:s}'
                        .format(type(modes)))
    overlap = abs(calcOverlap(mode, modes))
    if isinstance(modes, NMABase):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    show = plt.bar(arange, overlap, *args, **kwargs)
    plt.title('Overlap with {0:s}'.format(str(mode)))
    plt.xlabel('{0:s} mode index'.format(modes))
    plt.ylabel('Overlap')
    return show
    
def showCumulativeOverlap(mode, modes, *args, **kwargs):
    """Show cumulative overlap using :func:`~matplotlib.pyplot.plot`.
    
    :type mode: :class:`Mode`, :class:`Vector` 
    :arg modes: multiple modes
    :type modes: :class:`ModeSet`, :class:`ANM`, :class:`GNM`, or :class:`PCA` 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showCumulativeOverlap( p38_pca[0], p38_anm )
       # Let's also show the overlap
       showOverlap( p38_pca[0], p38_anm )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')  
          
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(mode, (Mode, Vector)):
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0:s}'
                        .format(type(mode)))
    if not isinstance(modes, (NMABase, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    cumov = (calcOverlap(mode, modes) ** 2).cumsum() ** 0.5
    if isinstance(modes, NMABase):
        arange = np.arange(0.5, len(modes)+0.5)
    else:
        arange = modes.getIndices() + 0.5
    show = plt.plot(arange, cumov, *args, **kwargs)
    plt.title('Cumulative overlap with {0:s}'.format(str(mode)))
    plt.xlabel('{0:s} mode index'.format(modes))
    plt.ylabel('Cumulative overlap')
    plt.axis((arange[0]-0.5, arange[-1]+0.5, 0, 1))
    return show
    
def resetTicks(x, y=None):
    """Reset X (and Y) axis ticks using values in given *array*.
    
    Ticks in the current figure should not be fractional values for this 
    function to work as expected. 
    
    """
    
    if x is not None:
        try:    
            xticks = plt.xticks()[0]
            xlist = list(xticks.astype(np.int32))
            if xlist[-1] > len(x):
                xlist.pop()
            if xlist:
                xlist = list(x[xlist]) 
                plt.xticks(xticks, xlist + [''] * (len(xticks) - len(xlist)))
        except:
            LOGGER.warning('xticks could not be reset.')
    if y is not None:
        try:    
            yticks = plt.yticks()[0]
            ylist = list(yticks.astype(np.int32))
            if ylist[-1] > len(y):
                ylist.pop()
            if ylist:
                ylist = list(y[ylist]) 
                plt.yticks(yticks, ylist + [''] * (len(yticks) - len(ylist)))
        except:
            LOGGER.warning('xticks could not be reset.')
    
