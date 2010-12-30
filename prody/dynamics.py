# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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
  * :class:`Mode`
  * :class:`ModeSet`
  * :class:`Vector`
  
Base Classes
------------

  * :class:`NMABase`
  * :class:`GNMBase`
  * :class:`VectorBase`

Inheritance Diagram
-------------------

.. inheritance-diagram:: prody.dynamics
   :parts: 1

Functions
---------

Many of the functions documented in this page accepts a *modes* argument (may 
also appear in different names). One of the following may be accepted as this 
agument:

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

  * :func:`calcCovariance`
  * :func:`calcCrossCorrelations`
  * :func:`calcSqFlucts`  
  * :func:`calcProjection`

**Write data**:

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
  * :func:`printOverlapTable`
  
**Sampling**:

  * :func:`sampleModes`
  * :func:`traverseMode`

**Model reduction**:
    
  * :func:`reduceModel`
  * :func:`sliceVector`
  * :func:`sliceMode`

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
    

Examples
--------

This is a preliminary part for examples shown in this section.


>>> from prody import *
>>> import matplotlib.pyplot as plt
>>> import numpy as np
>>> import os

>>> # First let's parse a PDB structure of Cyclin-dependent kinase 2
>>> cdk2 = parsePDB('1hcl')
>>> # Step 1: Instantiate an ANM instance
>>> cdk2_anm = ANM('cdk2_1hcl')
>>> # Step 2: Build Hessian matrix for selected atoms 
>>> cdk2_anm.buildHessian( cdk2.select('calpha') ) # 294 alpha carbons are passed
>>> # Step 3: Calculate modes, using default parameters
>>> cdk2_anm.calcModes() # by default 20 modes will be calculated
>>> cdk2_anm 
<ANM: cdk2_1hcl (20 modes, 294 nodes)>

Following are from :ref:`p38-xray-calculations`, and will be used to
demonstrate comparative analysis functions. 

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
   import os
   
   plt.close('all')
   
   cdk2 = parsePDB('1hcl')
   cdk2_anm = ANM('cdk2_1hcl')
   cdk2_anm.buildHessian( cdk2.select('calpha') )
   cdk2_anm.calcModes()


   p38_pca = loadModel('p38_xray.pca.npz')
   p38_anm = loadModel('1p38.anm.npz') 
   p38_ensemble = loadEnsemble('p38_X-ray.ens.npz')
   p38_structure = parsePDB('p38_ref_chain.pdb')
   
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import os.path
import time
import os
import gzip

import numpy as np
linalg = None
scipyla = None
plt = None
KDTree = None

import prody
from .atomic import *
from .ensemble import *
from prody import ProDyLogger as LOGGER
from prody import ProDyAtomSelect as SELECT

__all__ = ['ANM', 'GNM', 'PCA', 'EDA', 'Mode', 'ModeSet', 'Vector', 
           
           'NMABase', 'GNMBase', 'VectorBase',
           
           'calcANM', 'calcGNM', 
           
           'calcCovariance', 'calcCrossCorrelations', 'calcSqFlucts',
           
           'calcProjection',  
           
           'writeArray', 'writeModes', 'writeNMD', 'writeOverlapTable',
           
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD', 
           
           'calcOverlap', 'calcCumulativeOverlap', 'calcCumulativeOverlapArray', 
           'calcSubspaceOverlap', 'printOverlapTable',
           
           'sampleModes', 'traverseMode',
            
           'reduceModel', 'sliceVector', 'sliceMode',
            
           'showContactMap', 'showCrossCorrelations', 'showCumulativeOverlap', 
           'showCumFractOfVariances', 'showFractOfVariances', 'showMode', 
           'showOverlap', 'showOverlapTable', 'showProjection', 'showCrossProjection',
           'showEllipsoid', 'showSqFlucts', 'showScaledSqFlucts', 'showNormedSqFlucts',
           ]

VMDPATH = '/usr/local/bin/vmd'
ZERO = 1e-8

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
        
        >>> mode = cdk2_anm[0]
        >>> mode
        <Mode: 1 from ANM cdk2_1hcl>
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
        294
        """
        return self._model._n_atoms
    
    def getNumOfDegOfFreedom(self):
        """Return number of degrees of freedom.
        
        >>> mode.getNumOfDegOfFreedom()
        882
        """
        return self._model._dof
    
    def getName(self):
        """A descriptive name for the mode instance.
        
        >>> mode.getName()
        'Mode 1 from ANM cdk2_1hcl'
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
        <ANM: cdk2_1hcl (20 modes, 294 nodes)>
        """
        return self._model
    
    def getArray(self):
        """Return the normal mode array (eigenvector).
        
        >>> mode.getArray() # doctest: +SKIP
        array([ -6.07335804e-02,  -8.41067011e-03,  -9.12056373e-03,
                -1.03352547e-01,   9.21216689e-04,  -2.41896444e-02,
                -1.02868623e-01,   2.57252944e-02,  -3.77778634e-02,
                -1.02059320e-01,  -4.99044459e-04,  -2.25037521e-02,
                 ....
                 3.26339508e-02,  -7.98732544e-03,  -3.70856447e-03,
                 2.40309445e-02,  -7.74481996e-03,  -6.85740998e-03,
                 2.17636226e-02,  -8.72281427e-03,  -1.29909447e-02,
                 6.52858846e-03,  -9.49345509e-03,  -7.46635551e-03])
        """
        return self._model._array[:, self._index].copy()
    getEigenvector = getArray
    
    
    def getEigenvalue(self):
        """Return normal mode eigenvalue.
        
        >>> mode.getEigenvalue() # doctest: +SKIP
        0.50192959268070148
        """
        return self._model._eigvals[self._index]
    
    def getVariance(self):
        """Variance along the mode. 
        
        If the model is not a PCA, inverse of the eigenvalue is returned.
        
        >>> mode.getVariance() # doctest: +SKIP
        1.9923113013903169
        >>> mode.getEigenvalue()**-1 # doctest: +SKIP
        1.9923113013903169
        """
        return self._model._vars[self._index]

    def getFractOfVariance(self):
        """Return fraction of variance explained by the mode.
        
        Fraction of variance is the ratio of the variance along this mode to 
        the trace of the covariance matrix.
        
        See :meth:`getVariance`
        
        >>> mode.getFractOfVariance() # doctest: +SKIP
        0.15330741156510538
        """
        return self.getVariance() / self._model._trace
    
    def getCollectivity(self, masses=None):
        """Return collectivity of the mode.
        
        :arg masses: atomic masses
        :type masses: :class:`numpy.ndarray`
        
        This function implements collectivity as defined in equation 5 of [BR95]_.
        
        If *masses* are provided, they will be incorporated in the calculation.
        Otherwise, atoms are assumed to have uniform masses.
        
        >>> mode.getCollectivity() # doctest: +SKIP
        0.43143628650850246
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
        """Return covariance matrix calculated for this mode instance.
        
        >>> mode.getCovariance() # doctest: +SKIP
        array([[  7.34877529e-03,   1.01769275e-03,   1.10359002e-03, ...,
                 -7.89960500e-04,   1.14870995e-03,   9.03430500e-04],
               [  1.01769275e-03,   1.40934850e-04,   1.52830305e-04, ...,
                 -1.09397423e-04,   1.59078724e-04,   1.25111279e-04],
               [  1.10359002e-03,   1.52830305e-04,   1.65729784e-04, ...,
                 -1.18630995e-04,   1.72505593e-04,   1.35671162e-04],
               ..., 
               [ -7.89960500e-04,  -1.09397423e-04,  -1.18630995e-04, ...,
                  8.49172232e-05,  -1.23481186e-04,  -9.71147411e-05],
               [  1.14870995e-03,   1.59078724e-04,   1.72505593e-04, ...,
                 -1.23481186e-04,   1.79558430e-04,   1.41218035e-04],
               [  9.03430500e-04,   1.25111279e-04,   1.35671162e-04, ...,
                 -9.71147411e-05,   1.41218035e-04,   1.11064312e-04]])
        """
        array = self.getArray()
        return np.outer(array, array) * self.getVariance()
    
    def getSqFlucts(self):
        """Return square fluctuations.
        
        Square fluctuations are obtained by multiplying the squared the mode 
        array with the variance (:meth:`getVariance`) along the mode.
        
        >>> mode.getSqFlucts() # doctest: +SKIP
        array([  7.65543992e-03,   2.24488389e-02,   2.52444001e-02,
                 2.17615637e-02,   3.60212802e-02,   2.70340666e-02,
                 ...
                 2.32777245e-03,   1.86683818e-03,   2.27626659e-03,
                 1.36372235e-03,   1.43149042e-03,   3.75539965e-04])
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
            raise TypeError('array must be a 1-dimentional numpy.ndarray')
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
        
        For a 3-dimentional vector, returns length of the vector divided by 3.
        """
        if self._is3d: 
            return len(self._array)/3
        else:
            return len(self._array)

class NMABase(object):
    """Base class for Normal Mode Analysis Calculaations.
    
    Derived classes are:
        
        * :class:`GNMBase`
        * :class:`NMA`
        * :class:`PCA`
    """
    def __init__(self, name):
        """Initialize a Normal Mode analysis with a given name."""
        self._name = str(name)
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
        
class NMA(NMABase):
    """A class for analysis of externally calculated Hessian matrices and 
    normal modes.
    
    """
    
    def __init__(self):
        pass
    
    def setHessian(self, hessian):
        pass
    
    def calcModes(self, n_modes):
        pass
    
    def setEigens(self, vectors, values):
        pass


class ModeSet(object):
    """A class for providing access to data for a subset of modes.
    
    Instances are obtained by slicing an NMA model (:class:`ANM`, :class:`GNM`, 
    or :class:`PCA`).
    
    ModeSet's contain a reference to the model and a list of mode indices.
    Methods common to NMA models are also defined for mode sets.
    
    >>> modes = cdk2_anm[:3]
    >>> modes
    <ModeSet: 1 to 3 from cdk2_1hcl (3 modes)>
    """
    
    __slots__ = ['_model', '_indices', '_slice']
    
    def __init__(self, model, indices):
        if not isinstance(model, NMABase):
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
        """Return True if mode instance is from a 3-dimensional model.
        
        >>> modes.is3d()
        True
        """
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Return number of nodes.
        
        >>> modes.getNumOfAtoms()
        294
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
        882
        """
        return self._model._dof
        
    def getModes(self):
        """Return all modes in the subset in a list.
        
        >>> modes.getModes()
        [<Mode: 1 from ANM cdk2_1hcl>, <Mode: 2 from ANM cdk2_1hcl>, <Mode: 3 from ANM cdk2_1hcl>]
        """
        getMode = self._model.getMode
        return [getMode(i) for i in self._indices]
    
    def getName(self):
        """Return name.
        
        >>> modes.getName()
        'Modes 1 to 3 from ANM cdk2_1hcl'
        """
        return str(self)
    
    def getModel(self):
        """Return the model that the modes belongs to.
        
        >>> modes.getModel()
        <ANM: cdk2_1hcl (20 modes, 294 nodes)>
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
        
        >>> modes.getEigenvalues()
        array([ 0.50192959,  0.63339708,  0.74119293])
        """
        return self._model._eigvals[self._indices].copy()

    def getEigenvectors(self):
        """Return eigenvectors.
        
        >>> modes.getEigenvectors()
        array([[-0.06073358,  0.08408354,  0.09105849],
               [-0.00841067,  0.06913276, -0.00746136],
               [-0.00912056, -0.06069708,  0.0470834 ],
               ..., 
               [ 0.00652859, -0.03727421, -0.11962175],
               [-0.00949346, -0.00743206, -0.02354374],
               [-0.00746636,  0.025803  ,  0.01264865]])
        """
        return self._model._array[:, self._indices].copy()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues).
        
        >>> modes.getVariances()
        array([ 1.9923113 ,  1.57878846,  1.34917639])
        >>> modes.getEigenvalues()**-1
        array([ 1.9923113 ,  1.57878846,  1.34917639])
        """
        return self._model._vars[self._indices].copy()

    def getArray(self):
        """Return eigenvectors.
        
        >>> modes.getArray()
        array([[-0.06073358,  0.08408354,  0.09105849],
               [-0.00841067,  0.06913276, -0.00746136],
               [-0.00912056, -0.06069708,  0.0470834 ],
               ..., 
               [ 0.00652859, -0.03727421, -0.11962175],
               [-0.00949346, -0.00743206, -0.02354374],
               [-0.00746636,  0.025803  ,  0.01264865]])
        """
        return self._model._array[:, self._indices].copy()
        
    def getCovariance(self):
        """Return covariance matrix calculated for modes in the set.
        
        >>> modes.getCovariance()
        array([[ 0.02969777,  0.00927842, -0.00116957, ..., -0.02043412,
                -0.00273034,  0.00588272],
               [ 0.00927842,  0.00776161, -0.00694599, ..., -0.00297353,
                -0.00041509,  0.00281408],
               [-0.00116957, -0.00694599,  0.00897312, ..., -0.00414555,
                -0.00061088, -0.00153348],
               ..., 
               [-0.02043412, -0.00297353, -0.00414555, ...,  0.02158429,
                 0.00411363, -0.00365695],
               [-0.00273034, -0.00041509, -0.00061088, ...,  0.00411363,
                 0.00101462, -0.00056333],
               [ 0.00588272,  0.00281408, -0.00153348, ..., -0.00365695,
                -0.00056333,  0.00137807]])
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
        
        When available, this method uses Bio.KDTree.
        
        *masses* is not used yet.        """
        if KDTree is None: 
            prody.importBioKDTree()
            if not KDTree:
                LOGGER.debug('Using a slower method for building Kirchhoff matrix.')
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
        start = time.time()
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        if KDTree:
            kdtree = KDTree(3)
            kdtree.set_coords(coords) 
            kdtree.all_search(cutoff)
            for i, j in kdtree.all_get_indices():
                kirchhoff[i, j] = -gamma
                kirchhoff[j, i] = -gamma
                kirchhoff[i, i] += gamma
                kirchhoff[j, j] += gamma
        else:
            cutoff2 = cutoff * cutoff
            for i in range(n_atoms):
                res_i3 = i*3
                res_i33 = res_i3+3
                xyz_i = coords[i, :]
                for j in range(i+1, n_atoms):
                    i2j = coords[j, :] - xyz_i
                    i2j2 = np.dot(i2j, i2j)
                    if i2j2 > cutoff2:
                        continue             
                    kirchhoff[i, j] = -gamma
                    kirchhoff[j, i] = -gamma
                    kirchhoff[i, i] += gamma
                    kirchhoff[j, j] += gamma
            
        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'.format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._cutoff = cutoff
        self._gamma = gamma
        self._n_atoms = n_atoms
        self._dof = n_atoms
        
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize
        Kirchhoff matrix.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if self._kirchhoff is None:
            raise RuntimeError('Kirchhoff matrix is not set')
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
                eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            values, vectors = linalg.eigh(self._kirchhoff, turbo=turbo, eigvals=eigvals)
        else:
            values, vectors = linalg.eigh(self._kirchhoff)
        n_zeros = sum(values < ZERO)
        if n_zeros < 1: 
            LOGGER.warning('Less than 6 zero eigenvalues are calculated.')
            shift = n_zeros - 1
        elif n_zeros > 1: 
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
        
        When available, this method uses Bio.KDTree.
        
        *masses* is not used yet.                       
        """
        if KDTree is None: 
            prody.importBioKDTree()
            if not KDTree:
                LOGGER.debug('Using a slower method for building Hessian matrix.')
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
        start = time.time()
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        hessian = np.zeros((dof, dof), 'd')
        if KDTree:
            kdtree = KDTree(3)
            kdtree.set_coords(coords) 
            kdtree.all_search(cutoff)
            for i, k in kdtree.all_get_indices():
                if k > i: 
                    j = i
                    i = k
                else:
                    j = k
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
        else:
            cutoff2 = cutoff * cutoff 
            for i in range(n_atoms):
                res_i3 = i*3
                res_i33 = res_i3+3
                xyz_i = coords[i, :]
                for j in range(i+1, n_atoms):
                    i2j = coords[j, :] - xyz_i
                    i2j2 = np.dot(i2j, i2j)
                    if i2j2 > cutoff2:
                        continue             
                    res_j3 = j*3
                    res_j33 = res_j3+3
                    super_element = -np.outer(i2j, i2j) / i2j2 * gamma
                    hessian[res_i3:res_i33, res_j3:res_j33] = super_element 
                    hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                    hessian[res_i3:res_i33, res_i3:res_i33] -= super_element
                    hessian[res_j3:res_j33, res_j3:res_j33] -= super_element
                    kirchhoff[i, j] = -gamma
                    kirchhoff[j, i] = -gamma
                    kirchhoff[i, i] += gamma
                    kirchhoff[j, j] += gamma
        LOGGER.info('Hessian was built in {0:.2f}s.'.format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._gamma = gamma
        self._cutoff = cutoff
        self._n_atoms = n_atoms
        self._dof = dof

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize
        Hessian matrix.

        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg zeros: If ``True``, modes with zero eigenvalues will be kept.
        :type zeros: bool, default is ``False``
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if self._hessian is None:
            raise RuntimeError('Hessian matrix is not set')
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
                eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            values, vectors = linalg.eigh(self._hessian, turbo=turbo, eigvals=eigvals)
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
    (also known as Essential Dynamics Analysis (EDA) in [AA93]_)."""

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
        
        *coordsets* must have the following methods: 
            * ``getCoordinates``
            * ``getCoordsets``
            * ``getNumOfCoordsets``
            * ``getNumOfAtoms``
        
        :class:`~prody.ensemble.Ensemble` and :class:`~prody.atomic.Atomic`
        instances are acceptable.
        
        If *weights* is ``None``, but *coordsets* has :meth:`getWeights` method,
        weights from that method will be used. 
        
        """
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
        """Calculate principal (or essential) modes.

        This method uses :func:`scipy.linalg.eigh` function to diagonalize
        covariance matrix.
        
        :arg n_modes: number of non-zero eigenvalues/vectors to calculate. 
                      If ``None`` is given, all modes will be calculated. 
        :type n_modes: int or None, default is 20
        
        :arg turbo: Use a memory intensive, but faster way to calculate modes.
        :type turbo: bool, default is ``True``
        
        """
        if linalg is None:
            prody.importLA()

        start = time.time()
        dof = self._dof
        if scipyla:        
            if n_modes: 
                eigvals = (dof - n_modes, dof - 1)
            else: 
                eigvals = None
                n_modes = dof
            values, vectors = linalg.eigh(self._cov, turbo=turbo, eigvals=eigvals)
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

EDA = PCA

def saveModel(nma, filename=None, matrices=False):
    """Save *nma* model data as :file:`filename.nma.npz`. 
    
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
    return Vector(attr_dict['array'], str(attr_dict['name']), bool(attr_dict['is3d']))

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
    
    NMD file format is described at 
    http://www.csb.pitt.edu/People/abakan/software/NMWiz/nmdformat.html.
    
    .. note:: 
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`Vector` instance is given, it will be normalized before
          it is written. It's length before normalization will be written
          as the scaling factor of the vector.
        
    """
    if not isinstance(modes, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(modes)))
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
        out.write('mode 1 {0:.2f} {1:s}\n'.format(abs(modes), ' '.join(['{0:.3f}'.format(x) for x in modes.getNormed().getArray()])))
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
    
def calcANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Return an ANM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`~prody.atomic.Atomic` instance.  
    
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
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    anm = ANM(name)
    anm.buildHessian(ag.select(selstr), cutoff, gamma)
    anm.calcModes(n_modes)
    return anm

def calcGNM(pdb, selstr='all', cutoff=15., gamma=1., n_modes=20, zeros=False):
    """Return an GNM instance for given PDB identifier or atom data.
    
    By default only alpha carbons are considered, but selection string
    helps selecting a subset of it.
    
    *pdb* can be :class:`~prody.atomic.Atomic` instance.  
    
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
        raise TypeError('pdb must be an atom container, not {0:s}'.format(type(pdb)))
    gnm = GNM(name)
    gnm.buildKirchhoff(ag.select(selstr), cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm

def calcProjection(ensemble, modes):
    """Return projection of conformational deviations onto given modes.

    For K conformations and M modes, a (K,M) matrix is returned.     
    
    >>> calcProjection(p38_ensemble, p38_pca[:3])
    array([[ 11.4104979 ,   2.65096089,   1.11556184],
           [  7.84594842,  -1.85640127,   1.59130996],
           ...
           [  5.51572309,  -5.10756938,  -5.98343901],
           [  5.87457149,  -7.40788451,   2.20797633]])
    """
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
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


def calcOverlap(rows, cols):
    """Return overlap (or correlation) between two sets of modes (*rows* and *cols*).
    
    Returns a matrix whose rows correspond to modes passed as 
    *rows* argument, and columns correspond to those passed as *cols* 
    argument.
    
    >>> calcOverlap(p38_pca[0], p38_anm[2]) # doctest: +SKIP
    -0.71366564906422636
    """
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'.format(type(rows)))
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
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
    """Write table of overlaps (correlations) between two sets of modes to a file.

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

def sliceVector(vector, atoms, selstr):
    """Return a slice of *vector* matching *atoms* specified by *selstr*.
    
    Note that retuned :class:`Vector` instance is not normalized.
    
    :arg vector: vector instance to be sliced
    :type vector: :class:`VectorBase`
    
    :arg atoms: atoms for which *vector* describes a deformation, motion, etc.
    :type atoms: :class:`~prody.atomic.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: :class:`Vector`
    """
    if not isinstance(vector, VectorBase):
        raise TypeError('vector must be a VectorBase instance, not {0:s}'.format(type(vector)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'.format(type(atoms)))
    if atoms.getNumOfAtoms() != vector.getNumOfAtoms(): 
        raise ValueError('number of atoms in *vector* and *atoms* must be equal')
    return Vector(vector.getArrayNx3()[SELECT.getBoolArray(atoms, selstr), :].flatten(),
                  '{0:s} slice "{1:s}"'.format(str(vector), selstr), vector.is3d())

def sliceMode(mode, atoms, selstr):
    """Return a slice of *mode* matching *atoms* specified by *selstr*.
    
    This works sligtly difference from :func:`sliceVector`. Mode array 
    (eigenvector) is multiplied by square-root of the variance along the mode.
    If mode is from an elastic network model, variance is defined as the 
    inverse of the eigenvalue.
    
    Note that retuned :class:`Vector` instance is not normalized.
    
    :arg mode: mode instance to be sliced
    :type mode: :class:`Mode`
    
    :arg atoms: atoms for which *mode* describes a deformation, motion, etc.
    :type atoms: :class:`~prody.atomic.Atomic`
    
    :arg selstr: selection string
    :type selstr: str 
    
    :returns: :class:`Vector`
    """
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, not {0:s}'.format(type(mode)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0:s}'.format(type(atoms)))
    if atoms.getNumOfAtoms() != mode.getNumOfAtoms(): 
        raise ValueError('number of atoms in *mode* and *atoms* must be equal')
    return Vector(mode.getArrayNx3()[SELECT.getBoolArray(atoms, selstr), :].flatten() * mode.getVariance()**0.5,
                  '{0:s} slice "{1:s}"'.format(str(mode), selstr), mode.is3d())
    
def reduceModel(model, atoms, selstr):
    """Returns reduced NMA model.
    
    Reduces a :class:`NMA` model to a subset of *atoms* matching a selection *selstr*.
    
    This function behaves depending on the type of the model.
    
    :arg model: dynamics model
    :type model: :class:`ANM`, :class:`GNM`, or :class:`PCA`
    :arg atoms: atoms that were used to build the model
    :arg selstr: a selection string specifying subset of atoms  

    For ANM and GNM:    
       This function implements the method in [KH00]_. Selected atoms 
       constitute the system and the rest are the environment.
    
    For PCA:
       This function simply takes the sub-covariance matrix for the selected
       atoms.
    """
    if linalg is None:
        prody.importLA()

    #LOGGER.warning('Implementation of this function is not finalized. Use it with caution.')
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
        LOGGER.warning('selection results in same number of atoms, model is not reduced')
        return None

    if model.is3d():
        system = np.tile(system, (3,1)).transpose().flatten()
        other = np.tile(other, (3,1)).transpose().flatten()
    ss = matrix[system,:][:,system]
    if isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(ss)
        return eda
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(linalg.inv(oo), os))
    
    if isinstance(model, GNM):
        gnm = GNM(model.getName()+' reduced')
        gnm.setKirchhoff(matrix)
        return gnm
    elif isinstance(model, ANM):
        anm = ANM(model.getName()+' reduced')
        anm.setHessian(matrix)
        return anm
    elif isinstance(model, PCA):
        eda = PCA(model.getName()+' reduced')
        eda.setCovariance(matrix)
        return eda

def writeModes(filename, modes, format='g', sep=' ', compressed=False):
    """Write *modes* (eigenvectors) into a plain text file with name *filename*.
    
    See also :func:`writeArray`.
        
    >>> writeModes('p38_pca_modes_1-3.txt', p38_pca[:3])
    'p38_pca_modes_1-3.txt'
    >>> os.remove('p38_pca_modes_1-3.txt')
    """
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
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
     
    >>> writeArray('p38_cross-correlations.txt', calcCrossCorrelations(p38_pca))
    'p38_cross-correlations.txt'
    >>> os.remove('p38_cross-correlations.txt')
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

def sampleModes(modes, atoms=None, n_confs=1000, rmsd=1.0):
    """Return an ensemble of randomly sampled conformations along given *modes*.
    
    If *atoms* are provided, sampling will be around its active coordinate set. 
    Otherwise, sampling is around the 0 coordinate set.
    
    :arg modes: Modes along which sampling will be performed.
    :type modes: :class:`Mode`, :class:`ModeSet`, :class:`PCA`, :class:`ANM` or :class:`NMA`   
    
    :arg atoms: Atoms whose active coordinate set will be used as the initial 
        conformation.
    :type atoms: :class:`~prody.atomic.Atomic`  
    
    :arg n_confs: Number of conformations to generate. Default is 1000.
    :type n_steps: int 
    
    :arg rmsd: The maximum RMSD that the conformations will have with 
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
        raise ValueError('modes must be from a 3-dimentional model')
    n_confs = int(n_confs)
    n_atoms = modes.getNumOfAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, (Atomic)):
            raise TypeError('{0:s} is not correct type for atoms'.format(type(modes)))
        if atoms.getNumOfAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = [atoms.getCoordinates()]

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
       # since faster ANM modes does not correspong to top ranking PCA modes
       
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
        raise ValueError('modes must be from a 3-dimentional model')
    if len(modes) != 3:
        raise ValueError('length of modes is not equal to 3')
    if onto is not None:
        if not isinstance(onto, (NMABase, ModeSet)):
            raise TypeError('onto must be a NMA or ModeSet instance, '
                            'not {0:s}'.format(type(onto)))
        if not onto.is3d():
            raise ValueError('onto must be from a 3-dimentional model')
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
    conformations :math:`[R_{-n} R_{-n+1} ... R_{-1} R_0 R_1 ... R_n]` are generated.
    
    :math:`R_0` is the active coordinate set of *atoms*. 
    :math:`R_k = R_0 + sk\lambda_iu_i`, where :math:`s` is found using
    :math:`s = ((N (\\frac{RMSD}{n})^2) / \lambda_i^{-1}) ^{0.5}`, where
    :math:`N` is the number of atoms.
    
    
    .. plot::
       :context:
       :include-source:
        
       trajectory = traverseMode( cdk2_anm[0], cdk2.select('calpha'), n_steps=8, rmsd=1.4 )
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
        raise ValueError('mode must be from a 3-dimentional model.')
    n_atoms = mode.getNumOfAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, Atomic):
            raise TypeError('{0:s} is not correct type for atoms'.format(type(modes)))
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
    
    
def calcSqFlucts(modes):
    """Return sum of square-fluctuations for given set of normal *modes*.
    
    >>> calcSqFlucts(p38_pca)
    array([  0.94163178,   0.97486815,   0.81056074,   0.59926465,
             0.80505867,   0.93568339,   1.1634173 ,   1.39873827,
             ...
             0.13855481,   0.10944998,   0.09378239,   0.14173194,
             0.17928521,   0.12930389,   0.21500264,   0.16070006,   0.21077809])
    """
    if not isinstance(modes, (Mode, NMABase, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    square_fluctuations = np.zeros(modes.getNumOfAtoms()) 
    if isinstance(modes, Mode):
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
    
    >>> calcCrossCorrelations(p38_anm)
    array([[ 1.        ,  0.94782223,  0.8736113 , ...,  0.60437502,
             0.55178945,  0.31997558],
           [ 0.94782223,  1.        ,  0.96980278, ...,  0.51726959,
             0.53243166,  0.31677069],
           [ 0.8736113 ,  0.96980278,  1.        , ...,  0.32284386,
             0.35342037,  0.14512422],
           ..., 
           [ 0.60437502,  0.51726959,  0.32284386, ...,  1.        ,
             0.95600469,  0.86549921],
           [ 0.55178945,  0.53243166,  0.35342037, ...,  0.95600469,
             1.        ,  0.9528413 ],
           [ 0.31997558,  0.31677069,  0.14512422, ...,  0.86549921,
             0.9528413 ,  1.        ]])
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

def calcCovariance(modes):
    """Calculate covariance matrix from given modes and return it."""
    return modes.getCovariance()

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
      
    By default ``marker='o', ls='None'`` is passed to the plotting function to disable 
    lines in projections onto 2 or 3-d spaces.
    
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
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    if isinstance(modes, Mode):
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
        show.plot(projection[:, 0], projection[:, 1], projection[:, 2], *args, **kwargs)
        show.set_xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        show.set_ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
        show.set_zlabel('Mode {0:d} coordinate'.format(modes[2].getIndex()+1))
    else:
        raise ValueError('Projection onto upto 3 modes can be shown. You have given {0:d} mode.'.format(len(modes)))
    
    return show

def showCrossProjection(ensemble, mode_x, mode_y, scale=None, scalar=None, *args, **kwargs):
    """Show projection of conformational deviations using :func:`~matplotlib.pyplot.plot`.
    
    This function is differs from :func:`showProjection` by accepting modes
    from two different models.
    
    :arg ensemble: Ensemble for which deviations will be projected
    :type ensemble: :class:`~prody.ensemble.Ensemble`
    :arg mode_x: Projection onto this mode will be shown along x-axis. 
    :type mode_x: :class:`Mode`
    :arg mode_y: Projection onto this mode will be shown along y-axis.
    :type mode_y: :class:`Mode`
    :arg scale: Scale width of the projection onto one of modes. ``x`` and ``y`` are accepted.
    :type scale: str
    
    By default ``marker='o', ls='None'`` is passed to the plotting function to disable 
    lines.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showCrossProjection(p38_ensemble, p38_pca[0], p38_anm[2])
    
    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    |more| See :ref:`p38-xray-plotting` for a more elaborate example.
       
    """
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(ensemble, prody.Ensemble):
        raise TypeError('ensemble must be an Ensemble, not {0:s}'.format(type(ensemble)))
    if not isinstance(mode_x, Mode):
        raise TypeError('mode_x must be a Mode instance, not {0:s}'.format(type(mode_x)))
    if not mode_x.is3d():
        raise Exception('mode_x must be 3-dimensional')
    if not isinstance(mode_y, Mode):
        raise TypeError('mode_y must be a Mode instance, not {0:s}'.format(type(mode_y)))
    if not mode_y.is3d():
        raise Exception('mode_y must be 3-dimensional')
    
    xcoords = calcProjection(ensemble, mode_x) 
    ycoords = calcProjection(ensemble, mode_y)
    if scale is not None and isinstance(scale, str) and scale.lower() in ('x', 'y'):
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
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'.format(mode_x, scalar))
            else:
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'.format(mode_y, scalar))
        if scale.lower() == 'x':
            xcoords = xcoords * scalar  
        elif scale.lower() == 'y': 
            ycoords = ycoords / scalar
    else:
        LOGGER.warning('{0:s} is not a valid value for scale argument. Only "x" or "y" are accepted.'
                       .format(str(scale)))
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
        raise TypeError('rows must be an NMA model or a ModeSet, not {0:s}'.format(type(rows)))
    if not isinstance(rows, (NMABase, ModeSet)):
        raise TypeError('cols must be an NMA model or a ModeSet, not {0:s}'.format(type(cols)))
    overlap = abs(calcOverlap(rows, cols))
    if isinstance(rows, NMABase):
        rows = rows[:]
    if isinstance(cols, NMABase):
        cols = cols[:]
    Y = rows.getIndices()+0.5#np.arange(modes1.indices[0]+0.5, len(modes1)+1.5)
    X = cols.getIndices()+0.5#np.arange(modes2.indices[0]+0.5, len(modes2)+1.5)
    axis = (X[0], X[-1], Y[0], Y[-1])
    X, Y = np.meshgrid(X, Y)
    show = plt.pcolor(X, Y, overlap, cmap=plt.cm.jet, *args, **kwargs), plt.colorbar()
    plt.xlabel(str(cols))
    plt.ylabel(str(rows))
    plt.axis(axis)
    return show

def showCrossCorrelations(modes, *args, **kwargs):
    """Show cross-correlations for given modes using :func:`~matplotlib.pyplot.imshow`.
    
    See also :func:`getCrossCorrelations`. 
    
    By default, *origin=lower* and *interpolation=bilinear* keyword
    arguments are passed to imshow function. User can overwrite these
    parameters.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,5))
       # Show cross-correlations for ANM modes 1-3
       showCrossCorrelations( cdk2_anm[:3] )
       
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
       showMode( cdk2_anm[0] )
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
        raise TypeError('mode must be a Mode, not {0:s}'.format(type(modes)))
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
       showSqFlucts( cdk2_anm[0] )
       showSqFlucts( cdk2_anm[1] )
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
    plt.title(str(modes.getModel()))
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
    for arg in args:
        if isinstance(arg, (Mode, ModeSet, NMA)):
            modesarg.append(args.pop(0))
    show = [plt.plot(sqf, *args, label=str(modes), **kwargs)]
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        scalar = mean / sqf.mean()
        show.append(plt.plot(sqf * scalar, *args, label='{0:s} (x{1:.2f})'.format(str(modes), scalar), **kwargs))
    return show

def showNormedSqFlucts(modes, *args, **kwargs):
    """Show normalized square fluctuations using :func:`~matplotlib.pyplot.plot`.
    
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
    for arg in args:
        if isinstance(arg, (Mode, ModeSet, NMA)):
            modesarg.append(args.pop(0))
    show = [plt.plot(sqf/(sqf**2).sum()**0.5, *args, label='{0:s}'.format(str(modes)), **kwargs)]    
    plt.xlabel('Indices')
    plt.ylabel('Square fluctuations')
    for modes in modesarg:
        sqf = calcSqFlucts(modes)
        show.append(plt.plot(sqf/(sqf**2).sum()**0.5, *args, label='{0:s}'.format(str(modes)), **kwargs))
    return show

def showContactMap(enm, *args, **kwargs):
    """Show Kirchhoff matrix using :func:`~matplotlib.pyplot.spy`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(4,4))
       showContactMap( cdk2_anm )

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
        raise TypeError('mode must be Mode or Vector, not {0:s}'.format(type(mode)))
    if not isinstance(modes, (NMABase, ModeSet)):
        raise TypeError('modes must be NMA or ModeSet, not {0:s}'.format(type(modes)))
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
        raise TypeError('mode must be NMA, ModeSet, Mode or Vector, not {0:s}'.format(type(mode)))
    if not isinstance(modes, (NMABase, ModeSet)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'.format(type(modes)))
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
    
