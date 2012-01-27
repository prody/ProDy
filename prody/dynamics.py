# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
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
"""This module defines classes and functions for protein dynamics analysis. 

Protein dynamics
===============================================================================

Modeling and analysis 
-------------------------------------------------------------------------------

Following classes are designed for modeling and analysis of protein dynamics:

============  =================================================================
Model         Description
============  =================================================================
:class:`ANM`  Anisotropic network model, for coarse-grained NMA
:class:`GNM`  Gaussian network model, for coarse-grained dynamics analysis
:class:`PCA`  Principal component analysis of conformation ensembles
:class:`EDA`  Essential dynamics analysis of dynamics trajectories
:class:`NMA`  Normal mode analysis, for analyzing data from external programs
============  =================================================================

Following classes are for supporting above dynamics models by allowing
analysis of individual modes:

================  =============================================================
Class             Description
================  =============================================================
:class:`Mode`     analyze individual normal/principal/essential modes
:class:`ModeSet`  analyze subset of modes from one of aforementioned models
:class:`Vector`   analyze modified modes or deformation vectors
================  =============================================================

Customize models 
-------------------------------------------------------------------------------

Following classes allow for using structure or distance based, or other custom 
force constants and cutoff distances in :class:`ANM` and :class:`GNM` 
calculations: 

============================  =================================================
Class                         Description
============================  =================================================
:class:`Gamma`                base class for developing property custom force 
                              constant calculation methods
:class:`GammaStructureBased`  secondary structure based force constants
:class:`GammaVariableCutoff`  atom type based variable cutoff function
============================  =================================================

Function library
===============================================================================

Dynamics of the functions described below accept a *modes* argument (may also 
appear in different names), which may refer to one or more of the following:

  * a dynamics model, :class:`ANM`, :class:`GNM`, :class:`NMA`, :class:`PCA`, 
    or :class:`EDA` 
  * a :class:`Mode` obtained by indexing an NMA model, e.g. ``anm[0]``
  * a :class:`ModeSet` obtained by slicing an NMA model, e.g. ``anm[0:10]``

Some of these functions may also accept :class:`Vector` instances as *mode* 
argument.  These are noted in function documentations. 


Analyze models
-------------------------------------------------------------------------------

Following functions are for calculating atomic properties using modes
from dynamics models:

===========================  ==================================================
Function                     Calculated data or property
===========================  ==================================================
:func:`calcCollectivity`     degree of collectivity of a mode
:func:`calcCovariance`       covariance matrix for given modes
:func:`calcCrossCorr`        cross-correlations of fluctuations
:func:`calcPerturbResponse`  response to perturbations in positions
:func:`calcProjection`       projection of conformations onto modes
:func:`calcSqFlucts`         square-fluctuations
:func:`calcTempFactors`      temperature factors fitted to exp. data
===========================  ==================================================

Compare models
-------------------------------------------------------------------------------

Following functions are for comparing normal modes or dynamics models:

===========================  ==================================================
Function                      Output
===========================  ==================================================
:func:`calcOverlap`          overlap (correlation) between modes
:func:`calcCumOverlap`       cumulative overlap between modes
:func:`calcCumOverlapArray`  incremental cumulative overlap
:func:`calcSubspaceOverlap`  overlap between normal mode subspaces 
:func:`calcCovOverlap`       covariance overlap between models
:func:`printOverlapTable`    formatted overlap table printed on screen
===========================  ==================================================


Generate conformers
-------------------------------------------------------------------------------

Following functions can be used to generate conformers along normal modes:

====================  =========================================================
Function              Method
====================  =========================================================
:func:`deformAtoms`   deform atoms along a mode
:func:`sampleModes`   deform along random combination of a set of modes 
:func:`traverseMode`  traverse a mode along both directions
====================  =========================================================

Edit model
-------------------------------------------------------------------------------

Following functions can be used to reduce, slice, or extrapolate models:

========================  =====================================================
Function                  Description
========================  =====================================================
:func:`sliceMode`         take a slice of the normal mode      
:func:`sliceModel`        take a slice of a model
:func:`sliceVector`       take a slice of a vector
:func:`reduceModel`       reduce a model to a subset of atoms
:func:`extrapolateModel`  extrapolate a coarse-grained model to all-atoms  
========================  =====================================================

Parse/write data
-------------------------------------------------------------------------------

Following functions are for calculating atomic properties using modes
from dynamics models:

=========================  ====================================================
Function                   Input/output
=========================  ====================================================
:func:`parseArray`         numeric arrays, e.g. coordinates, eigenvectors
:func:`parseModes`         normal modes
:func:`parseNMD`           normal mode, coordinate, and atomic data for NMWiz
:func:`parseSparseMatrix`  matrix data in sparse coordinate list format
:func:`writeArray`         numeric arrays, e.g. coordinates, eigenvectors
:func:`writeModes`         normal modes
:func:`writeNMD`           normal mode, coordinate, and atomic data
:func:`writeOverlapTable`  overlap between modes in a formatted table
=========================  ====================================================

Save/load models 
-------------------------------------------------------------------------------

The above models and objects can be efficiently saved and loaded in later
Python sessions using the following functions:

| :func:`loadModel`, :func:`saveModel` - load/save dynamics models
| :func:`loadVector`, :func:`saveVector` - load/save normal modes and vectors
  

Short-hand functions
-------------------------------------------------------------------------------

Following allow for performing some dynamics calculations in one function call:

| :func:`calcANM` - perform ANM calculations
| :func:`calcGNM` - perform GNM calculations

Plotting functions
-------------------------------------------------------------------------------

Plotting functions are called by the name of the plotted data/property 
and are prefixed with "show".  Function documentations refers to the 
:mod:`matplotlib.pyplot` function utilized for actual plotting. 
Arguments and keyword arguments are passed to the Matplotlib functions.  


===========================  ==================================================
Function                     Plotted data
===========================  ==================================================
:func:`showMode`             mode shape
:func:`showOverlap`          overlap between modes
:func:`showSqFlucts`         square-fluctuations
:func:`showEllipsoid`        depict projection of a normal mode space on 
                             another 
:func:`showContactMap`       contact map based on a Kirchhoff matrix
:func:`showProjection`       projection of conformations onto normal modes
:func:`showOverlapTable`     overlaps between two models
:func:`showScaledSqFlucts`   square-fluctuations fitted to experimental data
:func:`showNormedSqFlucts`   normalized square-fluctuations
:func:`showCrossProjection`  project conformations onto modes from different 
                             models
:func:`showFractOfVar`       fraction of variances
:func:`showCrossCorr`        cross-correlations between fluctuations
                             in atomic positions
:func:`showCumOverlap`       cumulative overlap of a mode with multiple modes 
                             from another model
:func:`showCumFractOfVar`    cumulative fraction of variances 
:func:`resetTicks`           change ticks in a plot
===========================  ==================================================

Visualize modes
-------------------------------------------------------------------------------

Finally, normal modes can be visualized and animated using VMD plugin 
:ref:`nmwiz`. Following functions allow for running NMWiz from within Python: 
 
| :func:`viewNMDinVMD` - run VMD and load normal mode data
| :func:`getVMDpath`, :func:`setVMDpath` - get/set path to VMD executable


Examples
-------------------------------------------------------------------------------

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

:mod:`prody.dynamics`
===============================================================================

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
   
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os.path
import time
import os
from types import FunctionType
import sys

import numpy as np
linalg = None
scipyla = None
plt = None
scipy_sparse = None
scipy_sparse_la = None

import prody
from atomic import *
from ensemble import *
from tools import *
import measure
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS
from prody import ProDyAtomSelect as SELECT
from prody import ProDyException


__all__ = ['ANM', 'GNM', 'NMA', 'PCA', 'EDA', 'Mode', 'ModeSet', 'Vector', 
           
           'NMABase', 'GNMBase', 'VectorBase', 
           
           'Gamma', 'GammaStructureBased', 'GammaVariableCutoff',
           
           'calcANM', 'calcGNM', 
           
           'calcCollectivity', 'calcCovariance', 'calcCrossCorr',
           'calcCrossCorrelations', 
           
           'calcSqFlucts', 'calcTempFactors',
           
           'calcProjection', 'calcPerturbResponse', 'scanPerturbationResponse', 
           
           'parseArray', 'parseModes', 'parseNMD',
           
           'writeArray', 'writeModes', 'parseSparseMatrix',
           'writeNMD', 'writeOverlapTable',
           
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD', 
           
           'calcOverlap', 'calcCumOverlap', 'calcCumulativeOverlap', 
           'calcCumOverlapArray', 'calcCumulativeOverlapArray', 
           'calcSubspaceOverlap', 
           'calcCovOverlap', 'calcCovarianceOverlap', 'printOverlapTable',
           
           'deformAtoms', 'deform', 'sampleModes', 'traverseMode',
            
           'extrapolateModel', 'reduceModel', 'sliceVector', 
           'sliceMode', 'sliceModel',
            
           'showContactMap', 'showCrossCorr', 'showCrossCorrelations', 
           'showCumOverlap', 'showCumulativeOverlap', 
           'showFractOfVar', 'showFractOfVariances', 
           'showCumFractOfVar', 'showCumFractOfVariances', 'showMode', 
           'showOverlap', 'showOverlapTable', 'showProjection', 
           'showCrossProjection', 'showEllipsoid', 'showSqFlucts', 
           'showScaledSqFlucts', 'showNormedSqFlucts', 'resetTicks'
           ]

ZERO = 1e-6

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
   
    __slots__ = []
    
    def __abs__(self):
        return np.sqrt((self._getArray()**2).sum())
    
    def __neg__(self):
        return Vector(-self._getArray(), '-({0:s})'.format(str(self)), 
                      self.is3d())
    
    def __div__(self, other):
        if isinstance(other, (int, float, long)):
            return Vector(self._getArray() / other, 
                          '({1:s})/{0}'.format(other, str(self)), self.is3d())
        else:
            raise TypeError('{0} is not a scalar'.format(other))
    
    def __idiv__(self, other):
        return self.__div__(other)
    
    def __mul__(self, other):
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector(other * self._getArray(), 
                          '{0}*({1:s})'.format(other, str(self)), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self._getArray(), other._getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
    
    def __rmul__(self, other):   
        """Return scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector(other * self._getArray(), 
                          '{0}*({1:s})'.format(other, str(self)), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self._getArray(), other._getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
            
    def __imul__(self, other):
        return self.__mul__(other)    

    def __add__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self._getArray() + other._getArray(), 
                          '{0:s} + {1:s}'.format(str(self), str(other)), 
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __radd__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self._getArray() + other._getArray(), 
                          '{0:s} + {1:s}'.format(str(other), str(self)), 
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))
        
    def __iadd__(self, other):
        return self.__add__(other)   

    def __sub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self._getArray() - other._getArray(), 
                          '{0:s} - {1:s}'.format(str(self), str(other)), 
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __rsub__(self, other):
        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return  Vector(other._getArray() - self._getArray(), 
                           '{0:s} - {1:s}'.format(str(other), str(self)), 
                           self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __isub__(self, other):
        return self.__sub__(other)   

    def __pow__(self, other):
        if isinstance(other, (int, float, long)): 
            return Vector(self._getArray() ** other, 
                          '({0:s})**{1}'.format(str(self), other), self.is3d())
        else:
            raise TypeError('{0} is not a scalar'.format(other))

    def getArray(self):
        """Return a copy of array."""
        pass
        
    def getArrayNx3(self):
        """Return a copy of array with shape (N, 3)."""
        
        if self.is3d():
            return self.getArray().reshape((self.numAtoms(), 3))
        else:
            return self.getArray()
    
    def _getArrayNx3(self):
        """Return a copy of array with shape (N, 3)."""
        
        if self.is3d():
            return self._getArray().reshape((self.numAtoms(), 3))
        else:
            return self._getArray()
        

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
        self._index = int(index)
        
    def __len__(self):
        return self._model._dof
    
    def __repr__(self):
        return '<Mode: {0:d} from {1:s}>'.format(self._index+1, str(self._model))

    def __str__(self):
        return 'Mode {0:d} from {1:s}'.format(self._index+1, str(self._model))

    def is3d(self):
        """Return ``True`` if mode instance is from a 3-dimensional model."""
        
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""
        
        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._model._n_atoms
    
    def getNumOfDegOfFreedom(self):
        """Deprecated, use :meth:`numDOF`."""
        
        prody.deprecate('getNumOfDegOfFreedom', 'numDOF')
        return self.numDOF()
        
    def numDOF(self):
        """Return number of degrees of freedom (three times the number of 
        atoms)."""
        
        return self._model._dof
    
    def getName(self):
        """Deprecated, use :meth:`getTitle`."""
        
        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
        
    def getTitle(self):
        """A descriptive title for the mode instance."""
        
        return str(self)
    
    def getIndex(self):
        """Return the index of the mode.  Note that mode indices are 
        zero-based."""
        
        return self._index
    
    def getModel(self):
        """Return the model that the mode instance belongs to."""
        
        return self._model
    
    def getArray(self):
        """Return a copy of the normal mode array (eigenvector)."""
        
        return self._model._array[:, self._index].copy()
    
    getEigenvector = getArray

    def _getArray(self):
        """Return a copy of the normal mode array (eigenvector)."""
    
        return self._model._array[:, self._index]
    
    def getEigenvalue(self):
        """Return normal mode eigenvalue."""
        
        return self._model._eigvals[self._index]
    
    def getVariance(self):
        """Variance along the mode.  If the model is not a PCA, inverse of the
        eigenvalue is returned."""
        
        return self._model._vars[self._index]

    def getFractOfVariance(self):
        """Return fraction of variance explained by the mode.  Fraction of 
        variance is the ratio of the variance along this mode to the trace 
        of the covariance matrix.  See :meth:`getVariance`."""
        
        return self.getVariance() / self._model._trace
    
    def getCollectivity(self, masses=None):
        """Return the degree of collectivity of the mode.  This function 
        implements collectivity as defined in equation 5 of [BR95]_. See 
        also the :func:`calcCollectivity`."""

        return calcCollectivity(self)

    def getCovariance(self):
        """Return covariance matrix calculated for this mode instance."""
        
        array = self._getArray()
        return np.outer(array, array) * self.getVariance()
    
    def getSqFlucts(self):
        """Return square fluctuations.  Square fluctuations are obtained by 
        multiplying the squared the mode array with the variance (:meth:`
        getVariance`) along the mode."""
        
        if self.is3d():
            return (self._getArrayNx3()**2).sum(axis=1) * self.getVariance()
        else:
            return (self._getArray() ** 2)  * self.getVariance()


class Vector(VectorBase):
    
    """A class to provide operations on a modified mode array.  This class 
    holds only mode array (i.e. eigenvector) data, and has no associations 
    with an NMA instance.  Scalar multiplication of :class:`Mode` instance 
    or addition of two :class:`Mode` instances results in a :class:`Vector` 
    instance. 
    
    """
    
    __slots__ = ['_title', '_array', '_is3d']
    
    def __init__(self, array, title='Unknown', is_3d=True):
        """Instantiate with a name, an array, and a 3d flag."""
        
        if not isinstance(array, np.ndarray) or array.ndim != 1:
            raise TypeError('array must be a 1-dimensional numpy.ndarray')
        if not isinstance(is_3d, bool):
            raise TypeError('is_3d must be a boolean')
        self._title = str(title)
        self._array = array
        self._is3d = is_3d
        
    def __len__(self):
        return len(self._array)
    
    def __repr__(self):
        return '<Vector: {0:s}>'.format(self._title)
    
    def __str__(self):
        return self._title 
    
    def is3d(self):
        return self._is3d
    
    def getName(self):
        """Deprecated, use :meth:`getTitle`."""

        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
        
    def getTitle(self):
        """Get the descriptive title for the vector instance."""
        
        return self._title
    
    def setName(self, name):
        """Deprecated, use :meth:`setTitle`."""

        prody.deprecate('setName', 'setTitle')
        return self.setTitle(name)
        
    def setTitle(self, title):
        """Set the descriptive title for the vector instance."""
        
        self._title = str(title) 
    
    def getArray(self):
        """Return a copy of array."""
        
        return self._array.copy()
    
    def _getArray(self):
        """Return array."""
        
        return self._array

    def getNormed(self):
        """Return mode after normalizing it."""
        
        return Vector(self._array/(self._array**2).sum()**0.5, 
                      '({0:s})/||{0:s}||'.format(self._title), self._is3d)

    def getNumOfDegOfFreedom(self):
        """Deprecated, use :meth:`numDOF`."""
        
        prody.deprecate('getNumOfDegOfFreedom', 'numDOF')
        return self.numDOF()
        
    def numDOF(self):

        """Return number of degrees of freedom."""
        return len(self._array)

    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""
        
        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms.  For a 3-dimensional vector, returns length 
        of the vector divided by 3."""
        
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
    
    def __init__(self, title='Unknown'):
        """Initialize a Normal Mode analysis with a given name.
        
        .. versionchanged:: 0.7
           When an empty string is passed as *title* argument, NMA instance 
           is called "Unknown".
        
        """
        
        title = str(title)
        if title == '':
            title = 'Unknown'
        self._title = title
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
           A list or tuple of integers can be used for indexing."""
        
        if self._n_modes == 0:
            raise ProDyException('{0:s} modes are not calculated, try '
                                 'calcModes() method'.format(str(self)))
        if isinstance(index, int):
            return self.getMode(index)
        elif isinstance(index, slice):
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
                self._title, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'NMA {0:s}'.format(self._title)

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
        """Deprecated, use :meth:`numAtoms`."""
        
        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
    
    def getNumOfModes(self):
        """Deprecated, use :meth:`numModes`."""
        
        prody.deprecate('getNumOfModes', 'numModes')
        return self.numModes()
        
    def numModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes)."""
        
        return self._n_modes
    
    def getNumOfDegOfFreedom(self):
        """Deprecated, use :meth:`numDOF`."""
        
        prody.deprecate('getNumOfDegOfFreedom', 'numDOF')
        return self.numDOF()
        
    def numDOF(self):
        """Return number of degrees of freedom."""
        
        return self._dof
        
    def getName(self):
        """Deprecated, use :meth:`getTitle`."""

        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
        
    def getTitle(self):
        """Return title of the model."""
        
        return self._title
    
    def setName(self, name):
        """Deprecated, use :meth:`setTitle`."""

        prody.deprecate('setName', 'setTitle')
        return self.setTitle(name)
        
    def setTitle(self, title):
        """Set title of the model."""
        
        self._title = str(title)
    
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
        return Mode(self, index)

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
        """Return a copy of eigenvectors array."""
        
        if self._array is None: return None
        return self._array.copy()
        

    def _getArray(self):
        """Return eigenvectors array."""
        
        if self._array is None: return None
        return self._array

    def getCovariance(self):
        """Return covariance matrix.  If covariance matrix is not set or yet 
        calculated, it will be calculated using available modes."""
        
        if self._cov is None:
            if self._array is None:
                return None
            self._cov = np.dot(self._array, np.dot(np.diag(self._vars), 
                                                   self._array.T))
        return self._cov
        
    def calcModes(self):
        pass
    
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to the :class:`NMA` 
        instance.  If *eigenvalue* is not given, it will be set to 1. 
        Variances are set as the inverse eigenvalues.
        
        .. versionadded:: 0.5.3"""
        
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
        else:
            if vector.shape[0] != self._array.shape[0]: 
                raise ValueError('shape of vector do not match shape of ' 
                                 'existing vectors')
            self._array = np.concatenate((self._array, vector), 1)
            self._eigvals = np.concatenate((self._eigvals, value))
            self._n_modes += vector.shape[1]            
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
        self._vars = 1 / values


class NMA(NMABase):
    
    """A class for analysis of externally calculated Hessian matrices and 
    normal modes.
    
    """
    
    def __init__(self, name='Unknown'):
        NMABase.__init__(self, name)

        
class ModeSet(object):

    """A class for providing access to data for a subset of modes.
    
    Instances are obtained by slicing an NMA model (:class:`ANM`, :class:`GNM`, 
    or :class:`PCA`).
    
    ModeSet's contain a reference to the model and a list of mode indices.
    Methods common to NMA models are also defined for mode sets."""
    
    __slots__ = ['_model', '_indices']
    
    def __init__(self, model, indices):
        if not isinstance(model, NMABase):
            raise TypeError('model must be an NMA, not {0:s}'
                            .format(type(model)))
        self._model = model
        self._indices = np.array(indices, int)
        
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
        """Return True if mode instance is from a 3-dimensional model."""
        
        return self._model._is3d
    
    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""
        
        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._model._n_atoms
    
    def getNumOfModes(self):
        """Deprecated, use :meth:`numModes`."""
        
        prody.deprecate('getNumOfModes', 'numModes')
        return self.numModes()
        
    def numModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes)."""
        
        return len(self._indices)
    
    def getNumOfDegOfFreedom(self):
        """Deprecated, use :meth:`numDOF`."""
        
        prody.deprecate('getNumOfDegOfFreedom', 'numDOF')
        return self.numDOF()
        
    def numDOF(self):
        """Return number of degrees of freedom."""
        
        return self._model._dof
        
    def getModes(self):
        """Return a list that contains the modes in the mode set."""
        
        getMode = self._model.getMode
        return [getMode(i) for i in self._indices]
    
    def getName(self):
        """Deprecated, use :meth:`getTitle`."""
        
        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
    
    def getTitle(self):
        """Return title of the mode set."""
        
        return str(self)
    
    
    def getModel(self):
        """Return the model that the modes belongs to."""
        
        return self._model
    
    def getIndices(self):
        """Return indices of modes in the mode set."""
        
        return self._indices
    
    def getEigenvalues(self):
        """Return eigenvalues."""
        
        return self._model._eigvals[self._indices]

    def getEigenvectors(self):
        """Return a copy of eigenvectors."""
        
        return self._model._array[:, self._indices]
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues)."""
        
        return self._model._vars[self._indices]

    def getArray(self):
        """Return a copy of eigenvectors array."""
        
        return self._model._array[:, self._indices]

    def _getArray(self):
        """Return a copy of eigenvectors array."""

        return self._model._array[:, self._indices]
        
    def getCovariance(self):
        """Return covariance matrix calculated for modes in the set."""
        
        array = self._getArray()
        return np.dot(array, np.dot(np.diag(self.getVariances()), array.T))
   

class GNMBase(NMABase):

    """Class for Gaussian Network Model analysis of proteins."""

    def __init__(self, name='Unknown'):
        NMABase.__init__(self, name)
        self._is3d = False
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        
    def __repr__(self):
        return '<GNM: {0:s} ({1:d} modes, {2:d} nodes)>'.format(
                                    self._title, self.__len__(), self._n_atoms)
    
    def __str__(self):
        return 'GNM {0:s}'.format(self._title)
    
    def _reset(self):
        NMABase._reset(self)
        self._cutoff = None
        self._gamma = None
        self._kirchhoff = None
        self._is3d = False
    
    def getCutoff(self):
        """Return cutoff distance."""
        
        return self._cutoff
    
    def getGamma(self):
        """Return spring constant (or the gamma function or :class:`Gamma`
        instance)."""
        
        return self._gamma

    def getKirchhoff(self):
        """Return a copy of the Kirchhoff matrix."""
        
        if self._kirchhoff is None: return None
        return self._kirchhoff.copy()

    def _getKirchhoff(self):
        """Return the Kirchhoff matrix."""
        
        return self._kirchhoff

def checkENMParameters(cutoff, gamma):
    """Check type and values of *cutoff* and *gamma*."""

    if not isinstance(cutoff, (float, int)):
        raise TypeError('cutoff must be a float or an integer')
    elif cutoff < 4:
        raise ValueError('cutoff must be greater or equal to 4')
    if isinstance(gamma, Gamma):
        gamma_func = gamma.gamma
    elif isinstance(gamma, FunctionType):
        gamma_func = gamma
    else:
        if not isinstance(gamma, (float, int)):
            raise TypeError('gamma must be a float, an integer, derived '
                             'from Gamma, or a function')
        elif gamma <= 0:
            raise ValueError('gamma must be greater than 0')
        gamma = float(gamma)
        gamma_func = lambda dist2, i, j: gamma 
    return cutoff, gamma, gamma_func 

class GNM(GNMBase):
    
    """A class for Gaussian Network Model (GNM) analysis of proteins 
    ([IB97]_, [TH97]_).
    
    |example| See example :ref:`gnm`.
    
    """
    
    def setKirchhoff(self, kirchhoff):
        """Set Kirchhoff matrix."""
        
        if not isinstance(kirchhoff, np.ndarray):
            raise TypeError('kirchhoff must be a Numpy array')
        elif not kirchhoff.ndim == 2 or \
                 kirchhoff.shape[0] != kirchhoff.shape[1]:
            raise ValueError('kirchhoff must be a square matrix')
        elif kirchhoff.dtype != float:
            try:
                kirchhoff = kirchhoff.astype(float)
            except:
                raise ValueError('kirchhoff.dtype must be float')
                
        self._reset()
        self._kirchhoff = kirchhoff
        self._n_atoms = kirchhoff.shape[0]
        self._dof = kirchhoff.shape[0]
    
    def buildKirchhoff(self, coords, cutoff=10., gamma=1., **kwargs):
        """Build Kirchhoff matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        :type coords: :class:`~numpy.ndarray` or :class:`~prody.atomic.Atomic`
        
        :arg cutoff: cutoff distance (Å) for pairwise interactions
            default is 10.0 Å, , minimum is 4.0 Å
        :type cutoff: float
        
        :arg gamma: spring constant, default is 1.0
        :type gamma: float
        
        :arg sparse: Elect to use sparse matrices. Default is ``False``. If 
            Scipy is not found, :class:`ImportError` is raised.
        :type sparse: bool
        
        .. versionchanged:: 0.6
            Instances of :class:`Gamma` classes and custom functions are
            accepted as *gamma* argument.        

        .. versionchanged:: 0.7.3
           When Scipy is available, user can select to use sparse matrices for
           efficient usage of memory at the cost of computation speed."""
        
        slow = kwargs.get('slow', False)
        try:
            from KDTree import KDTree
        except ImportError:
            KDTree = False
        if not slow and not KDTree: 
            LOGGER.info('Using a slower method for building the Kirchhoff '
                         'matrix.')
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoords()
            except AttributeError:
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoordinates attribute')
        coords = checkCoords(coords, 'coords')
        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        self._reset()
        self._cutoff = cutoff
        self._gamma = g
                    
        n_atoms = coords.shape[0]
        start = time.time()
        if kwargs.get('sparse', False):
            prody.importScipySparse()
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        
        if not slow and KDTree:
            kdtree = measure.getKDTree(coords) 
            kdtree.all_search(cutoff)
            radii = kdtree.all_get_radii()
            r = 0
            for i, j in kdtree.all_get_indices():
                g = gamma(radii[r]**2, i, j)
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] = kirchhoff[i, i] + g
                kirchhoff[j, j] = kirchhoff[j, j] + g
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
                    kirchhoff[i, i] = kirchhoff[i, i] + g
                    kirchhoff[j, j] = kirchhoff[j, j] + g
            
        LOGGER.debug('Kirchhoff was built in {0:.2f}s.'
                     .format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._n_atoms = n_atoms
        self._dof = n_atoms
        
    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh` 
        function to diagonalize the Kirchhoff matrix. When Scipy is not found, 
        :func:`numpy.linalg.eigh` is used.

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
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean'
        assert isinstance(turbo, bool), 'turbo must be a boolean'
        if linalg is None:
            prody.importLA()
        start = time.time()
        shift = 0
        if scipyla:
            if n_modes is None:
                eigvals = None
                n_modes = self._dof 
            else:
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
                try:
                    values, vectors = scipy_sparse_la.eigsh(
                                self._kirchhoff, k=n_modes + 1, which='SA')
                except:
                    values, vectors = scipy_sparse_la.eigen_symmetric(
                                self._kirchhoff, k=n_modes + 1, which='SA')                
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
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                          ''.format(self._n_modes, time.time()-start))


class ANM(GNMBase):

    """Class for Anisotropic Network Model (ANM) analysis of proteins 
    ([PD00]_, [ARA01]_)
    
    |example| See example :ref:`anm`.
    
    """

    def __init__(self, name='Unknown'):
        GNMBase.__init__(self, name)
        self._is3d = True
        self._cutoff = None
        self._kirchhoff = None
        self._gamma = None
        self._hessian = None


    def __repr__(self):
        return '<ANM: {0:s} ({1:d} modes, {2:d} nodes)>'.format(
                                    self._title, self.__len__(), self._n_atoms)
                                    
    def __str__(self):
        return 'ANM {0:s}'.format(self._title)

    def _reset(self):
        GNMBase._reset(self)
        self._hessian = None
        self._is3d = True
        
    def getHessian(self):
        """Return a copy of the Hessian matrix."""
        
        if self._hessian is None: return None
        return self._hessian.copy()
    
    def _getHessian(self):
        """Return the Hessian matrix."""
        
        return self._hessian
    
    def setHessian(self, hessian):
        """Set Hessian matrix.  A symmetric matrix is expected, i.e. not a 
        lower- or upper-triangular matrix."""
        
        if not isinstance(hessian, np.ndarray):
            raise TypeError('hessian must be a Numpy array')
        elif hessian.ndim != 2 or hessian.shape[0] != hessian.shape[1]:
            raise ValueError('hessian must be square matrix')
        elif hessian.shape[0] % 3:
            raise ValueError('hessian.shape must be (3*n_atoms,3*n_atoms)')
        elif hessian.dtype != float:
            try:
                hessian = hessian.astype(float)
            except:
                raise ValueError('hessian.dtype must be float')
        self._reset()
        self._hessian = hessian
        self._dof = hessian.shape[0]
        self._n_atoms = self._dof / 3 

    def buildHessian(self, coords, cutoff=15., gamma=1., **kwargs):
        """Build Hessian matrix for given coordinate set.
        
        :arg coords: a coordinate set or anything with getCoordinates method
        :type coords: :class:`~numpy.ndarray` or :class:`~prody.atomic.Atomic`
        
        :arg cutoff: cutoff distance (Å) for pairwise interactions,
            default is 15.0 Å, minimum is 4.0 Å 
        :type cutoff: float
        
        :arg gamma: spring constant, default is 1.0 
        :type gamma: float, :class:`Gamma`
        
        :arg sparse: slect to use sparse matrices. Default is ``False``. If 
            Scipy is not found, :class:`ImportError` is raised.
        :type sparse: bool
        
        .. versionchanged:: 0.6
            Instances of :class:`Gamma` classes and custom functions are
            accepted as *gamma* argument.        
    
        .. versionchanged:: 0.7.3
           When Scipy is available, user can select to use sparse matrices for
           efficient usage of memory at the cost of computation speed.
                       
        """
        
        slow = kwargs.get('slow', False)
        try:
            from KDTree import KDTree
        except ImportError:
            KDTree = False
        if not slow and not KDTree: 
            LOGGER.info('Using a slower method for building the Hessian '
                         'matrix.')
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoords()
            except AttributeError:
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoordinates attribute')
        coords = checkCoords(coords, 'coords')
        cutoff, g, gamma = checkENMParameters(cutoff, gamma)
        self._reset()
        self._cutoff = cutoff
        self._gamma = g
        n_atoms = coords.shape[0]
        dof = n_atoms * 3
        start = time.time()
        
        if kwargs.get('sparse', False):
            prody.importScipySparse()
            kirchhoff = scipy_sparse.lil_matrix((n_atoms, n_atoms))
            hessian = scipy_sparse.lil_matrix((dof, dof))
        else:
            kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
            hessian = np.zeros((dof, dof), float)
        if not slow and KDTree:
            kdtree = measure.getKDTree(coords) 
            kdtree.all_search(cutoff)
            for i, j in kdtree.all_get_indices():
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
                hessian[res_i3:res_i33, res_i3:res_i33] = \
                    hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                hessian[res_j3:res_j33, res_j3:res_j33] = \
                    hessian[res_j3:res_j33, res_j3:res_j33] - super_element
                kirchhoff[i, j] = -g
                kirchhoff[j, i] = -g
                kirchhoff[i, i] = kirchhoff[i, i] - g
                kirchhoff[j, j] = kirchhoff[j, j] - g
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
                    hessian[res_i3:res_i33, res_i3:res_i33] = \
                        hessian[res_i3:res_i33, res_i3:res_i33] - super_element
                    hessian[res_j3:res_j33, res_j3:res_j33] = \
                        hessian[res_j3:res_j33, res_j3:res_j33] - super_element
                    kirchhoff[i, j] = -g
                    kirchhoff[j, i] = -g
                    kirchhoff[i, i] = kirchhoff[i, i] - g
                    kirchhoff[j, j] = kirchhoff[j, j] - g
        LOGGER.info('Hessian was built in {0:.2f}s.'.format(time.time()-start))
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._n_atoms = n_atoms
        self._dof = dof

    def calcModes(self, n_modes=20, zeros=False, turbo=True):
        """Calculate normal modes.  This method uses :func:`scipy.linalg.eigh` 
        function to diagonalize the Hessian matrix. When Scipy is not found, 
        :func:`numpy.linalg.eigh` is used.

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
        assert n_modes is None or isinstance(n_modes, int) and n_modes > 0, \
            'n_modes must be a positive integer'
        assert isinstance(zeros, bool), 'zeros must be a boolean' 
        assert isinstance(turbo, bool), 'turbo must be a boolean'
        if linalg is None:
            prody.importLA()
        
        start = time.time()
        shift = 5
        if scipyla:
            if n_modes is None:
                eigvals = None
                n_modes = self._dof 
            else: 
                if n_modes >= self._dof:
                    eigvals = None
                    n_modes = self._dof 
                else:
                    eigvals = (0, n_modes + shift)
            if eigvals: 
                turbo = False
            if isinstance(self._hessian, np.ndarray):            
                values, vectors = linalg.eigh(self._hessian, turbo=turbo, 
                                              eigvals=eigvals)
            else:
                prody.importScipySparseLA()
                try:
                    values, vectors = scipy_sparse_la.eigsh(
                            self._hessian, k=n_modes+6, which='SA')
                except:                
                    values, vectors = scipy_sparse_la.eigen_symmetric(
                            self._hessian, k=n_modes+6, which='SA')
        
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
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                          ''.format(self._n_modes, time.time()-start))


class PCA(NMABase):
    
    """A class for Principal Component Analysis (PCA) of conformational 
    ensembles.
    
    |example| See examples in :ref:`pca`.
    """

    def __init__(self, name='Unknown'):
        NMABase.__init__(self, name)
    
    def __repr__(self):
        return '<PCA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._title, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'PCA {0:s}'.format(self._title)
    
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

    def buildCovariance(self, coordsets):
        """Build a weighted covariance matrix for coodsets.  *coordsets* 
        argument may be an instance of one of the following:
          
        * :class:`~prody.atomic.Atomic`
        * :class:`~prody.ensemble.Ensemble`
        * :class:`~prody.ensemble.PDBEnsemble`
        * :class:`~prody.ensemble.Trajectory`
        * :class:`numpy.ndarray`

        .. versionchanged:: 0.8
           :class:`~numpy.ndarray` and :class:`~prody.ensemble.TrajectoryBase` 
           instances are accepted as *coordsets* argument.

        A NumPy array passed as *coordsets* argument must have the shape 
        (n_coordsets, n_atoms, 3).
        
        When a trajectory instance is passed as *coordsets* argument, 
        covariance will be built by aligning frames to the reference 
        coordinates.  Reference coordinates will be considered as the 
        average of coordinate sets in the trajectory. 
         
        .. note::        
           If *coordsets* is a :class:`~prody.ensemble.PDBEnsemble` instance,
           coordinates are treated specially.  Let's say **C**\_ij is the 
           super element of the covariance matrix that corresponds to atoms 
           *i* and *j*.  This super element is divided by number of coordinate
           sets (PDB structures) in which both of these atoms are observed 
           together."""
        
        if not isinstance(coordsets, (prody.Ensemble, prody.Atomic, 
                                      prody.TrajectoryBase, np.ndarray)):
            raise TypeError('coordsets must be an Ensemble, Atomic, Numpy '
                            'array instance')
        LOGGER.timeit()
        weights = None
        if isinstance(coordsets, np.ndarray): 
            if coordsets.ndim != 3 or coordsets.shape[2] != 3 or \
                coordsets.dtype not in (np.float32, float):
                raise ValueError('coordsets is not a valid coordinate array')
        elif isinstance(coordsets, prody.Atomic):
            coordsets = coordsets._getCoordsets()
        elif isinstance(coordsets, prody.Ensemble):
            if isinstance(coordsets, PDBEnsemble):
                weights = coordsets.getWeights() > 0
            coordsets = coordsets._getCoordsets()
        
        if isinstance(coordsets, prody.TrajectoryBase):
            nfi = coordsets.getNextIndex()
            coordsets.reset()
            n_atoms = coordsets.numSelected()
            dof = n_atoms * 3
            cov = np.zeros((dof, dof))
            mean = coordsets._getCoords().flatten()
            n_confs = 0
            n_frames = len(coordsets)
            LOGGER.info('Covariance will be calculated using {0:d} frames.'
                            .format(n_frames))
            coordsum = np.zeros(dof)
            LOGGER.progress('Calculating covariance', n_frames)
            for frame in coordsets:
                frame.superpose()
                coords = frame._getCoords().flatten()
                coordsum += coords
                cov += np.outer(coords, coords)
                n_confs += 1
                LOGGER.update(n_confs)
            LOGGER.clear()
            cov /= n_confs
            coordsum /= n_confs
            cov -= np.outer(coordsum, coordsum)
            coordsets.goto(nfi)
            self._cov = cov
        else:
            n_confs = coordsets.shape[0]
            if n_confs < 3:
                raise ValueError('coordsets must have more than 3 coordinate '
                                 'sets')
            n_atoms = coordsets.shape[1]
            if n_atoms < 3:
                raise ValueError('coordsets must have more than 3 atoms')
            dof = n_atoms * 3
            LOGGER.info('Covariance is calculated using {0:d} coordinate sets.'
                            .format(len(coordsets)))
            if weights is None:
                if coordsets.dtype == float:
                    self._cov = np.cov(coordsets.reshape((n_confs, dof)).T, 
                                       bias=1)
                else:
                    cov = np.zeros((dof, dof))
                    coordsets = coordsets.reshape((n_confs, dof))
                    mean = coordsets.mean(0)
                    LOGGER.progress('Building covariance', n_confs)
                    for i, coords in enumerate(
                                            coordsets.reshape((n_confs, dof))):
                        deviations = coords - mean
                        cov += np.outer(deviations, deviations)
                        LOGGER.update(n_confs)
                    LOGGER.clear()
                    cov /= n_confs 
                    self._cov = cov
            else:
                # PDB ensemble case
                mean = np.zeros((n_atoms, 3))
                for i, coords in enumerate(coordsets):
                    mean += coords * weights[i]
                mean /= weights.sum(0)
                d_xyz = ((coordsets - mean) * weights).reshape((n_confs, dof))
                divide_by = weights.astype(float).repeat(3, 
                                                axis=2).reshape((n_confs, dof))
                self._cov = np.dot(d_xyz.T, d_xyz) / np.dot(divide_by.T, 
                                                            divide_by)
        self._trace = self._cov.trace()
        self._dof = dof
        self._n_atoms = n_atoms
        LOGGER.timing('Covariance matrix was calculated in %2fs.')
        
    def calcModes(self, n_modes=20, turbo=True):
        """Calculate principal (or essential) modes.  This method uses 
        :func:`scipy.linalg.eigh` function to diagonalize covariance matrix. 
        When Scipy is not found, :func:`numpy.linalg.eigh` is used.
        
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
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                     .format(self._n_modes, time.time()-start))

    def performSVD(self, coordsets):
        """Calculate principal modes using singular value decomposition (SVD).
        *coordsets* argument may be a :class:`~prody.atomic.Atomic`, 
        :class:`~prody.ensemble.Ensemble`, or :class:`numpy.ndarray` instance.
        If *coordsets* is a numpy array it must have the shape 
        ``(n_coordsets, n_atoms, 3)``.
        
        .. versionadded:: 0.6.2
        
        .. versionchanged:: 0.8
           :class:`numpy.ndarray` instances are accepted as *coordsets* 
           argument.
        
        This is a considerably faster way of performing PCA calculations 
        compared to eigenvalue decomposition of covariance matrix, but is
        an approximate method when heterogeneous datasets are analyzed. 
        Covariance method should be preferred over this one for analysis of 
        ensembles with missing atomic data.  See :ref:`pca-xray-calculations`
        example for comparison of results from SVD and covariance methods.
        """

        if linalg is None:
            prody.importLA()

        start = time.time()
        if not isinstance(coordsets, (prody.Ensemble, prody.Atomic, 
                                      np.ndarray)):
            raise TypeError('coordsets must be an Ensemble, Atomic, Numpy '
                            'array instance')
        if isinstance(coordsets, np.ndarray):
            if coordsets.ndim != 3 or coordsets.shape[2] != 3 or \
                coordsets.dtype not in (np.float32, float):
                raise ValueError('coordsets is not a valid coordinate array')
            deviations = coordsets - coordsets.mean(0)
        else:
            if isinstance(coordsets, prody.Ensemble):
                deviations = coordsets.getDeviations()
            elif isinstance(coordsets, prody.Atomic):
                deviations = coordsets._getCoordsets() - \
                             coordsets._getCoords()

        n_confs = deviations.shape[0]
        if n_confs < 3:
            raise ValueError('coordsets must have more than 3 coordinate sets')
        n_atoms = deviations.shape[1]
        if n_atoms < 3:
            raise ValueError('coordsets must have more than 3 atoms')

        dof = n_atoms * 3        
        deviations = deviations.reshape((n_confs, dof)).T

        vectors, values, self._temp = linalg.svd(deviations, 
                                                 full_matrices=False)
        values = (values ** 2) / n_confs
        self._dof = dof
        self._n_atoms = n_atoms
        which = values > 1e-18
        self._eigvals = values[which]
        self._array = vectors[:, which]
        self._vars = self._eigvals
        self._trace = self._vars.sum()
        self._n_modes = len(self._eigvals)
        LOGGER.debug('{0:d} modes were calculated in {1:.2f}s.'
                         .format(self._n_modes, time.time()-start))
        
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to :class:`NMA` instance.
        If *eigenvalue* is not given, it will be set to 1.  Eigenvalue is also 
        set as the variance.
        
        .. versionadded:: 0.7"""

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
        and values must have shape ``(M,)``.  Eigenvalues are also set as the 
        variances.
        """
        
        NMABase.setEigens(self, vectors, values)
        self._vars = self._eigvals.copy()

class EDA(PCA):
    
    """A class for Essential Dynamics Analysis (EDA) [AA93]_.
    
    |example| See examples in :ref:`eda`.
    
    """

    def __repr__(self):
        return '<EDA: {0:s} ({1:d} modes, {2:d} atoms)>'.format(
                self._title, self._n_modes, self._n_atoms)

    def __str__(self):
        return 'EDA {0:s}'.format(self._title)

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
    >>> assignSecstr(header, ubi)
    <AtomGroup: 1aar_A_ca (76 atoms; 1 coordinate sets, active set index: 0)>

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
        n_atoms = atoms.numAtoms()
        sstr = atoms.getSecstrs()
        assert sstr is not None, 'secondary structure assignments must be set'
        chid = atoms.getChids()
        assert chid is not None, 'chain identifiers must be set'
        rnum = atoms.getResnums()
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
    >>> print( ca_p.getNames() ) # doctest: +ELLIPSIS
    ['P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P' 'P'
     'P' 'P' 'P' 'P' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA'
     ...
     'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA' 'CA'
     'CA']
         
    We set the radii of atoms: 
     
    >>> variableCutoff = GammaVariableCutoff(ca_p.getNames(), gamma=1,  
    ... default_radius=7.5, debug=True, P=10)
    >>> print( variableCutoff.getRadii() ) # doctest: +ELLIPSIS
    [ 10.   10.   10.   10.   10.   10.   10.   10.   10.   10.   10.   10.
      10.   10.   10.   10.   10.   10.   10.   10.   10.   10.    7.5   7.5
      ...
       7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5   7.5
       7.5]
    
    The above shows that for phosphate atoms radii is set to 10 Å, because
    we passed the ``P=10`` argument.  As for Cα atoms, the default 7.5 Å
    is set as the radius (``default_radius=7.5``).  Note we also passed
    ``debug=True`` argument for demonstration purposes. This argument 
    allows printing debugging information on the screen.
    
    We build :class:`ANM` Hessian matrix as follows:  
        
    >>> anm = ANM('HhaI-DNA')
    >>> anm.buildHessian(ca_p, gamma=variableCutoff, cutoff=20) # doctest: +ELLIPSIS
    P_5 -- CA_111 effective cutoff: 17.5 distance: 18.0730838818 gamma: 0
    CA_257 -- CA_111 effective cutoff: 15.0 distance: 18.948157826 gamma: 0
    CA_277 -- CA_111 effective cutoff: 15.0 distance: 18.2928823863 gamma: 0
    CA_258 -- CA_111 effective cutoff: 15.0 distance: 17.3706104095 gamma: 0
    CA_277 -- CA_104 effective cutoff: 15.0 distance: 17.6521736055 gamma: 0
    CA_258 -- CA_104 effective cutoff: 15.0 distance: 18.2241431897 gamma: 0
    P_7 -- P_5 effective cutoff: 20.0 distance: 12.4763126364 gamma: 1.0
    ...
    CA_306 -- CA_24 effective cutoff: 15.0 distance: 12.2571321279 gamma: 1.0
    CA_24 -- CA_310 effective cutoff: 15.0 distance: 15.4210115751 gamma: 0
    CA_307 -- CA_308 effective cutoff: 15.0 distance: 3.81520169847 gamma: 1.0
    CA_307 -- CA_309 effective cutoff: 15.0 distance: 5.28885564182 gamma: 1.0
    CA_309 -- CA_308 effective cutoff: 15.0 distance: 3.81530260923 gamma: 1.0
    
    Note that we set passed ``cutoff=20.0`` to the :meth:`ANM.buildHessian` 
    method.  This is equal to the largest possible cutoff distance (between 
    two phosphate atoms) for this system, and ensures that all of the 
    potential interactions are evaluated. 
    
    For pairs of atoms for which the actual distance is larger than the 
    effective cutoff, the :meth:`GammaVariableCutoff.gamma` method returns 
    ``0``.  This annuls the interaction between those atom pairs.
    
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


def saveModel(nma, filename=None, matrices=False, **kwargs):
    """Save *nma* model data as :file:`filename.nma.npz`. 
    
    .. versionadded:: 0.5
    
    By default, eigenvalues, eigenvectors, variances, trace of covariance 
    matrix, and name of the model will be saved.  If *matrices* is ``True``,
    covariance, Hessian or Kirchhoff matrices are saved too, whichever are 
    available.  If *filename* is ``None``, name of the NMA instance will be 
    used as the filename, after ``" "`` (white spaces) in the name are 
    replaced with ``"_"`` (underscores).  Extension may differ based on 
    the type of the NMA model.  For ANM models, it is :file:`.anm.npz`.
    Upon successful completion of saving, filename is returned. This 
    function makes use of :func:`numpy.savez` function.
    """
    
    if not isinstance(nma, NMABase):
        raise TypeError('invalid type for nma, {0:s}'.format(type(nma)))
    if len(nma) == 0:
        raise ValueError('nma instance does not contain data')
    
    dict_ = nma.__dict__
    attr_list = ['_title', '_trace', '_array', '_eigvals', '_vars', '_n_atoms',
                 '_dof', '_n_modes']
    if filename is None:
        filename = nma.getTitle().replace(' ', '_')
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
    elif isinstance(nma, EDA):
        type_ = 'EDA'
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
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename



def loadModel(filename):
    """Return NMA instance after loading it from file (*filename*).  
    This function makes use of :func:`numpy.load` function.  See 
    also :func:`saveModel`.
    
    .. versionadded:: 0.5"""
    
    attr_dict = np.load(filename)
    try:
        type_ = attr_dict['type']
    except KeyError:
        raise IOError('{0:s} is not a valid NMA model file'.format(filename))
    try:
        title = str(attr_dict['_title'])
    except KeyError: 
        title = str(attr_dict['_name'])
    if type_ == 'ANM':
        nma = ANM(title)
    elif type_ == 'PCA':
        nma = PCA(title)
    elif type_ == 'EDA':
        nma = EDA(title)
    elif type_ == 'GNM':
        nma = GNM(title)
    elif type_ == 'NMA':
        nma = NMA(title)
    else:
        raise IOError('NMA model type is not recognized'.format(type_))
    dict_ = nma.__dict__ 
    for attr in attr_dict.files:
        if attr in ('type', '_name', '_title'): 
            continue
        elif attr in ('_trace', '_cutoff', '_gamma'):
            dict_[attr] = float(attr_dict[attr])
        elif attr in ('_dof', '_n_atoms', '_n_modes'):
            dict_[attr] = int(attr_dict[attr])
        else:
            dict_[attr] = attr_dict[attr]
    return nma

def saveVector(vector, filename, **kwargs):
    """Save *vector* data as :file:`filename.vec.npz`.  Upon successful 
    completion of saving, filename is returned.  This function makes use 
    of :func:`numpy.savez` function."""
    
    if not isinstance(vector, Vector):
        raise TypeError('invalid type for vector, {0:s}'.format(type(vector)))
    attr_dict = {}
    attr_dict['title'] = vector.getTitle()
    attr_dict['array'] = vector._getArray()
    attr_dict['is3d'] = vector.is3d()
    filename += '.vec.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename

def loadVector(filename):
    """Return :class:`Vector` instance after loading it from file (*filename*).
    This function makes use of :func:`numpy.load` function.  See also
    :func:`saveVector`."""
    
    attr_dict = np.load(filename)
    try:
        title = str(attr_dict['title'])
    except KeyError:
        title = str(attr_dict['name'])
    return Vector(attr_dict['array'], title, bool(attr_dict['is3d']))

def getVMDpath():
    """Return VMD path set by user or one identified automatically."""
    
    path = SETTINGS.get('vmd', None)
    if isExecutable(path):
        return path   
    else:
        LOGGER.warning('VMD path is not set by user, looking for it.')    

        #from types import StringType, UnicodeType
        vmdbin = None
        vmddir = None
        if sys.platform == 'win32': 
            if PY3K:
                import winreg as _winreg
            else:
                import _winreg
            for vmdversion in ('1.8.7', '1.9', '1.9.1'): 
                try:
                    key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                            'Software\\University of Illinois\\VMD\\' + 
                            vmdversion)
                    vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                    vmdbin = os.path.join(vmddir, 'vmd.exe') 
                except:    
                    pass
                try:
                    key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                'Software\\WOW6432node\\University of Illinois\\VMD\\' + 
                vmdversion)
                    vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                    vmdbin = os.path.join(vmddir, 'vmd.exe') 
                except:    
                    pass
        else:
            vmdbin = which('vmd')
            if False:
                pipe = os.popen('which vmd')
                vmdbin = pipe.next().strip()
                vmdfile = open(vmdbin)
                for line in vmdfile:
                    if 'defaultvmddir' in line:
                        exec(line.strip())
                        vmddir = defaultvmddir
                        break
                vmdfile.close()
        if False and \
           isinstance(vmdbin, (StringType, UnicodeType)) and \
           isinstance(vmddir, (StringType, UnicodeType)) and \
           os.path.isfile(vmdbin) and os.path.isdir(vmddir): 
            pass#return vmdbin, vmddir
        if isExecutable(vmdbin):
            setVMDpath(vmdbin)
            return vmdbin
        

def setVMDpath(path):
    """Set path to a VMD executable."""
    
    if isExecutable(path):
        SETTINGS['vmd'] = path
        SETTINGS.save()
        LOGGER.info("VMD path is set to '{0:s}'.".format(path))
    else:
        raise OSError('{0:s} is not executable.'.format(str(path)))

def parseNMD(filename, type=NMA):
    """Returns normal mode and atomic data parsed from an NMD file.
    Normal mode data is returned in an :class:`NMA` instance. Atomic
    data is returned in an :class:`~prody.atomic.AtomGroup` instance. 
    
    .. versionadded:: 0.5.3
    
    .. versionchanged:: 0.7
       User can pass NMA type for the data, eg. :class:`ANM` or :class:`PCA`.
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
        coords = np.fromstring( coords, dtype=float, sep=' ')
        dof = coords.shape[0]
        ag = None
        n_atoms = dof / 3
        coords = coords.reshape((n_atoms, 3))
        ag = AtomGroup(name)
        ag.setCoords(coords)
        data = atomic.pop('atomnames', None)
        if data is not None:
            ag.setNames(data.split())
        data = atomic.pop('resnames', None)
        if data is not None:
            ag.setResnames(data.split())
        data = atomic.pop('chainids', None)
        if data is not None:
            ag.setChids(data.split())
        data = atomic.pop('resnums', None)
        if data is not None:
            ag.setResnums(np.fromstring(data, int, sep=' '))
        data = atomic.pop('resids', None)
        if data is not None:
            ag.setResnums(np.fromstring(data, int, sep=' '))
        data = atomic.pop('bfactors', None)
        if data is not None:
            ag.setBetas(np.fromstring(data, float, sep=' '))
    nma = type(name)
    for mode in modes:
        
        items = mode.split()
        diff = len(items) - dof
        mode = np.array(items[diff:]).astype(float)
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
    *atoms*.  Returns *filename*, if file is successfully written.  NMD file 
    format is described at :ref:`nmd-format`.
    
    .. note:: 
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`Vector` instance is given, it will be normalized before
          it is written. It's length before normalization will be written
          as the scaling factor of the vector.
    """
    
    if not isinstance(modes, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0:s}'.format(type(modes)))
    if modes.numAtoms() != atoms.numAtoms():
        raise Exception('number of atoms do not match')
    out = openFile(filename, 'w')
    
    #out.write('#!{0:s} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0:s}\n'.format(os.path.abspath(filename)))
    name = modes.getTitle()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = os.path.splitext(os.path.split(filename)[1])[0]
    out.write('name {0:s}\n'.format(name))
    try:
        coords = atoms.getCoords()
    except:
        raise ProDyException('coordinates could not be retrived '
                             'from atoms instance')
    if coords is None:
        raise ProDyException('coordinates could not be retrived '
                             'from atoms instance')
    
    try:
        data = atoms.getNames()
        if data is not None:
            out.write('atomnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnames()
        if data is not None:
            out.write('resnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnums()
        if data is not None:
            out.write('resids {0:s}\n'.format(' '.join(data.astype('|S5'))))
    except:
        pass
    try:
        data = atoms.getChids()
        if data is not None:
            out.write('chainids {0:s}\n'.format(' '.join(data)))
    except:
        pass
    
    try:
        data = atoms.getBetas()
        if data is not None:
            out.write('bfactors {0:s}\n'.format(' '.join(
                            ['{0:.3f}'.format(x) for x in data.flatten()])))
    except:
        pass
    
    out.write('coordinates {0:s}\n'.format(
                    ' '.join(['{0:.3f}'.format(x) for x in coords.flatten()])))
    
    count = 0
    if isinstance(modes, Vector):
        out.write('mode 1 {0:.2f} {1:s}\n'.format(abs(modes), ' '.join(
                ['{0:.3f}'.format(x) for x in modes.getNormed()._getArray()])))
        count += 1
    else:
        if isinstance(modes, Mode):
            modes = [modes]
        for mode in modes:
            if mode.getEigenvalue() < ZERO:
                continue
            out.write('mode {0:d} {1:.2f} {2:s}\n'.format(
                       mode.getIndex()+1, mode.getVariance()**0.5, 
                       ' '.join(
                            ['{0:.3f}'.format(x) for x in mode._getArray()])))
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. '
                       'Given modes might have 0 eigenvalues.')
    out.close() 
    return filename  

def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""
    
    vmd = SETTINGS.get('vmd')
    if vmd:
        os.system('{0:s} -e {1:s}'.format(vmd, os.path.abspath(filename)))
    
def calcANM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, 
            zeros=False):
    """Return an :class:`ANM` instance and atoms used for the calculations.
    By default only alpha carbons are considered, but selection string helps 
    selecting a subset of it.  *pdb* can be :class:`~prody.atomic.Atomic` 
    instance.
    
    .. versionchanged:: 0.6
       Returns also the :class:`~prody.atomic.Selection` instance.
    """
    
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        title = ag.getTitle()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if isinstance(pdb, AtomGroup):
            title = ag.getTitle()
        else: 
            title = ag.getAtomGroup().getTitle()
    else:
        raise TypeError('pdb must be an atomic class, not {0:s}'
                        .format(type(pdb)))
    anm = ANM(title)
    sel = ag.select(selstr)
    anm.buildHessian(sel, cutoff, gamma)
    anm.calcModes(n_modes)
    return anm, sel 

def calcGNM(pdb, selstr='calpha', cutoff=15., gamma=1., n_modes=20, 
            zeros=False):
    """Return a :class:`GNM` instance and atoms used for the calculations.
    By default only alpha carbons are considered, but selection string helps 
    selecting a subset of it.  *pdb* can be :class:`~prody.atomic.Atomic` 
    instance.  
    
    .. versionchanged:: 0.6
       Returns also the :class:`~prody.atomic.Selection` instance.
    
    """
    
    if isinstance(pdb, str):
        ag = prody.parsePDB(pdb)
        title = ag.getTitle()
    elif isinstance(pdb, Atomic):
        ag = pdb
        if isinstance(pdb, AtomGroup):
            title = ag.getTitle()
        else: 
            title = ag.getAtomGroup().getTitle()
    else:
        raise TypeError('pdb must be an atom container, not {0:s}'
                        .format(type(pdb)))
    gnm = GNM(title)
    sel = ag.select(selstr)
    gnm.buildKirchhoff(sel, cutoff, gamma)
    gnm.calcModes(n_modes)
    return gnm, sel

def calcCollectivity(mode, masses=None):
    """Return collectivity of the mode.  This function implements collectivity 
    as defined in equation 5 of [BR95]_.  If *masses* are provided, they will 
    be incorporated in the calculation.  Otherwise, atoms are assumed to have 
    uniform masses.
    
    :arg mode: mode or vector
    :type mode: :class:`Mode` or :class:`Vector`
    
    :arg masses: atomic masses
    :type masses: :class:`numpy.ndarray`
    """
    
    is3d = mode.is3d()
    if masses is not None:
        if len(masses) != mode.numAtoms(): 
            raise ValueError('length of massesmust be equal to number of atoms')
        if is3d:
            u2in = (mode.getArrayNx3() ** 2).sum(1) / masses
    else:
        if is3d:
            u2in = (mode.getArrayNx3() ** 2 ).sum(1)
        else:
            u2in = (mode.getArrayNx3() ** 2 )
    u2in = u2in * (1 / u2in.sum() ** 0.5)
    coll = np.exp(-(u2in * np.log(u2in)).sum()) / mode.numAtoms()
    return coll
    

def calcProjection(ensemble, modes, rmsd=True):
    """Return projection of conformational deviations onto given modes.
    For K conformations and M modes, a (K,M) matrix is returned.
    
    .. versionchanged:: 0.8
       By default root-mean-square deviation (RMSD) along the normal mode is 
       calculated. To calculate the projection pass ``rmsd=True``.
       :class:`Vector` instances are accepted as *ensemble* argument to allow
       for projecting a deformation vector onto normal modes.  
    """
    
    if not isinstance(ensemble, (prody.Ensemble, prody.Conformation, 
                                 prody.Vector)):
        raise TypeError('ensemble must be Ensemble, Conformation, or Vector, '
                        'not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    if not modes.is3d(): 
        raise ValueError('modes must be 3-dimensional')
    if isinstance(ensemble, Vector):
        n_atoms = ensemble.numAtoms()
    else:
        n_atoms = ensemble.numSelected()
    if n_atoms != modes.numAtoms():
        raise ValueError('number of atoms are not the same')
    if isinstance(ensemble, Vector):
        if not ensemble.is3d(): 
            raise ValueError('ensemble must be a 3d vector instance')
        deviations = ensemble._getArray()
    else:
        deviations = ensemble.getDeviations()
    if deviations.ndim == 3:
        deviations = deviations.reshape((deviations.shape[0], 
                                         deviations.shape[1] * 3))
    elif deviations.ndim == 2:
        deviations = deviations.reshape((1, deviations.shape[0] * 3))
    else:
        deviations = deviations.reshape((1, deviations.shape[0]))
    projection = np.dot(deviations, modes._getArray())
    if rmsd:
        projection =  (1 / (n_atoms ** 0.5)) * projection
    return projection


def calcOverlap(rows, cols):
    """Return overlap (or correlation) between two sets of modes (*rows* and 
    *cols*).  Returns a matrix whose rows correspond to modes passed as *rows* 
    argument, and columns correspond to those passed as *cols* argument.
    
    .. versionchanged:: 0.7
       Both rows and columns are normalized prior to calculating overlap.       
    """
    
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('rows must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(rows)))
    if not isinstance(rows, (NMABase, ModeSet, Mode, Vector)):
        raise TypeError('cols must be NMA, ModeSet, Mode, or Vector, not {0:s}'
                        .format(type(cols)))
    
    if rows.numDOF() != cols.numDOF(): 
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
    and columns of the printed table.  This function may be used to take 
    a quick look into mode correspondences between two models.
    
    >>> # Compare top 3 PCs and slowest 3 ANM modes
    >>> printOverlapTable(p38_pca[:3], p38_anm[:3]) # doctest: +SKIP   
    Overlap Table
                            ANM 1p38
                        #1     #2     #3
    PCA p38 xray #1   -0.39  +0.04  -0.71
    PCA p38 xray #2   -0.78  -0.20  +0.22
    PCA p38 xray #3   +0.05  -0.57  +0.06
    """
    
    print getOverlapTable(rows, cols)

def writeOverlapTable(filename, rows, cols):
    """Write table of overlaps (correlations) between two sets of modes to a 
    file.  *rows* and *cols* are sets of normal modes, and correspond to rows
    and columns of the overlap table.  See also :func:`printOverlapTable`."""
    
    assert isinstance(filename, str), 'filename must be a string'
    out = openFile(filename, 'w')
    out.write(getOverlapTable(rows, cols))
    out.close()
    return filename
    
def getOverlapTable(rows, cols):
    """Make a formatted string of overlaps between modes in *rows* and *cols*.
    """

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
            if abs(overlap[i, j]).round(2) == 0.00:
                minplus = ' '
            elif overlap[i, j] < 0: 
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
    grained level to all atom level.  For each atom in *nodes* argument *atoms* 
    argument must contain a corresponding residue.  Note that modes in the 
    extrapolated model will not be normalized.  For a usage example see 
    :ref:`extrapolate`."""
    
    if not isinstance(enm, NMABase):
        raise TypeError('enm must be an NMABase instance')
    if not isinstance(nodes, Atomic):
        raise TypeError('nodes must be an Atomic instance')
    if enm.numAtoms() != nodes.numAtoms():
        raise ValueError('enm and nodes must have same number of atoms')
    
    if isinstance(atoms, Atomic):
        is3d = enm.is3d()            
        atom_indices = []
        indices = []
        hierview = atoms.getHierView()
        for i, node in enumerate(nodes):
            res = hierview[node.getChid(), node.getResnum(), node.getIcode()]
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
                          atoms.getACSIndex())
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
    if atoms.numAtoms() != vector.numAtoms(): 
        raise ValueError('number of atoms in *vector* and *atoms* must be '
                         'equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())
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
    if atoms.numAtoms() != mode.numAtoms(): 
        raise ValueError('number of atoms in *mode* and *atoms* must be equal')
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())
    vec = Vector(mode.getArrayNx3()[
                 which,:].flatten() * mode.getVariance()**0.5,
                 '{0:s} slice "{1:s}"'.format(str(mode), selstr), 
                 mode.is3d()) 
    return (vec, sel)

def sliceModel(model, atoms, selstr):
    """Return a slice of *model* matching *atoms* specified by *selstr*.
    
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
    if atoms.numAtoms() != model.numAtoms(): 
        raise ValueError('number of atoms in *model* and *atoms* must be '
                         'equal')
    
    array = model._getArray()
    if isinstance(atoms, AtomGroup):
        sel = atoms.select(selstr)
        which = sel.getIndices()
    else:
        which = SELECT.getIndices(atoms, selstr)
        sel = Selection(atoms.getAtomGroup(), atoms.getIndices()[which],
                        selstr, atoms.getACSIndex())

    nma = type(model)('{0:s} slice "{1:s}"'.format(model.getTitle(), selstr))
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
    *selstr*.  This function behaves differently depending on the type of the 
    *model* argument.  For ANM and GNM or other NMA models, this functions 
    derives the force constant matrix for system of interest (specified by the 
    *selstr*) from the force constant matrix for the *model* by assuming that 
    for any given displacement of the system of interest, the other atoms move 
    along in such a way as to minimize the potential energy.  This is based on 
    the formulation in in [KH00]_.  For PCA models, this function simply takes 
    the sub-covariance matrix for the selected atoms.

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
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not '
                         'built')

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
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(ss)
        return eda, system
    so = matrix[system,:][:,other]
    os = matrix[other,:][:,system]
    oo = matrix[other,:][:,other]
    matrix = ss - np.dot(so, np.dot(linalg.inv(oo), os))
    
    if isinstance(model, GNM):
        gnm = GNM(model.getTitle() + ' reduced')
        gnm.setKirchhoff(matrix)
        return gnm, system
    elif isinstance(model, ANM):
        anm = ANM(model.getTitle() + ' reduced')
        anm.setHessian(matrix)
        return anm, system
    elif isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(matrix)
        return eda, system

def writeModes(filename, modes, format='%.18e', delimiter=' '):
    """Write *modes* (eigenvectors) into a plain text file with name 
    *filename*. See also :func:`writeArray`.
    
    .. versionchanged:: 0.5.3
       A decompressed file is outputted.
    """
    
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    return writeArray(filename, modes._getArray(), format=format, 
                      delimiter=delimiter)

def parseModes(normalmodes, eigenvalues=None, nm_delimiter=None, 
               nm_skiprows=0, nm_usecols=None, ev_delimiter=None, 
               ev_skiprows=0, ev_usecols=None, ev_usevalues=None):
    """Return :class:`NMA` instance with normal modes parsed from *normalmodes*.
    
    .. versionadded:: 0.5.3
    
    In normal mode file *normalmodes*, columns must correspond to modes 
    (eigenvectors).  Optionally, *eigenvalues* can be parsed from a separate 
    file. If eigenvalues are not provided, they will all be set to 1.
    
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
    making some type and value checks.  Default *format* argument is ``"%d"``.
    Default *delimiter* argument is white space, ``" "``.
    
    *filename* will be returned upon successful writing."""
    
    if not isinstance(array, np.ndarray):
        raise TypeError('array must be a Numpy ndarray, not {0:s}'
                        .format(type(array)))
    elif not array.ndim in (1, 2):
        raise ValueError('array must be a 1 or 2-dimensional Numpy ndarray, '
                         'not {0:d}-d'.format(type(array.ndim)))
    np.savetxt(filename, array, format, delimiter)
    return filename

def parseArray(filename, delimiter=None, skiprows=0, usecols=None, 
               dtype=float):
    """Parse array data from a file.
    
    .. versionadded:: 0.5.3
    
    This function is using :func:`numpy.loadtxt` to parse the file.  Each row 
    in the text file must have the same number of values.
    
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
    
    :arg dtype: Data-type of the resulting array, default is :func:`float`. 
    :type dtype: :class:`numpy.dtype`.
    """

    array = np.loadtxt(filename, dtype=dtype, delimiter=delimiter, 
                       skiprows=skiprows, usecols=usecols)
    return array
        
def parseSparseMatrix(filename, symmetric=False, delimiter=None, skiprows=0,
                      irow=0, icol=1, first=1):
    """Parse sparse matrix data from a file.
    
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

    Data-type of the resulting array, default is :func:`float`."""
    
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
    irow = (sparse[:,irow] - first).astype(int)
    icol = (sparse[:,icol] - first).astype(int)
    matrix[irow, icol] = sparse[:,idata]
    if symmetric:
        matrix[icol, irow] = sparse[:,idata]
    return matrix

def sampleModes(modes, atoms=None, n_confs=1000, rmsd=1.0):
    """Return an ensemble of randomly sampled conformations along given 
    *modes*.
    
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
    ensure that the generated ensemble will have user given average *rmsd* 
    value. 
     
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
       showProjection(ensemble, p38_anm[:3], rmsd=True)
       
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
    n_atoms = modes.numAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, (Atomic)):
            raise TypeError('{0:s} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.numAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoords()

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
    array = modes._getArray()
    if array.ndim > 1:
        for i in range(n_confs):
            append( (array * scale * randn[i]).sum(1).reshape((n_atoms, 3)) )
    else:
        for i in range(n_confs):
            append( (array * scale * randn[i]).reshape((n_atoms, 3)) )

    ensemble = Ensemble('Conformations along {0:s}'.format(modes))
    if initial is None:
        ensemble.setCoords(np.zeros((n_atoms, 3)))
        ensemble.addCoordset(np.array(confs))
    else:    
        ensemble.setCoords(initial)
        ensemble.addCoordset(np.array(confs) + initial)
    return ensemble  

def showEllipsoid(modes, onto=None, n_std=2, scale=1., *args, **kwargs):
    """Show an ellipsoid using  :meth:`~mpl_toolkits.mplot3d.Axes3D.
    plot_wireframe`.
    
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
    
    if not plt: prody.importPyPlot()
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
        if onto.numAtoms() != modes.numAtoms():
            raise ValueError('modes and onto does not have same number of atoms')
        
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    
    var = modes.getVariances()
    #randn = np.random.standard_normal((1000, 3))
    #coef = ((randn ** 2 * var).sum(1) ** 0.5).mean()
    #scale = float(n_std) * modes.numAtoms()**0.5 * float(rmsd) / coef * var ** 0.5
    scale = float(n_std) * scale * var ** 0.5
    #scale = float(n_std) * modes.numAtoms() ** 0.5 * float(rmsd) / var.sum() ** 0.5 * var ** 0.5   

    x = scale[0] * np.outer(np.cos(u), np.sin(v))
    y = scale[1] * np.outer(np.sin(u), np.sin(v))
    z = scale[2] * np.outer(np.ones(np.size(u)), np.cos(v))
    if onto is not None:
        change_of_basis = np.dot(modes._getArray().T, onto._getArray())

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
    n_atoms = mode.numAtoms()
    initial = None
    if atoms is not None:
        if not isinstance(atoms, Atomic):
            raise TypeError('{0:s} is not correct type for atoms'
                            .format(type(atoms)))
        if atoms.numAtoms() != n_atoms:
            raise ValueError('number of atoms do not match')
        initial = atoms.getCoords()

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
    ensemble.setCoords(initial)    
    ensemble.addCoordset(np.array(confs_sub + [initial] + confs_add))
    return ensemble
  
def deform(atoms, mode, rmsd=None):
    """Deprecated, use :func:`deformAtoms`."""
    
    prody.deprecate('deform', 'deformAtoms')
    return deformAtoms(atoms, mode, rmsd)

def deformAtoms(atoms, mode, rmsd=None):
    """Generate a new coordinate set for *atoms* along the *mode*.
    
    .. versionadded:: 0.7
    
    *atoms* must be a :class:`~prody.atomic.AtomGroup` instance.
    New coordinate set will be appended to *atoms*. If *rmsd* is provided,
    *mode* will be scaled to generate a coordinate set with given RMSD distance
    to the active coordinate set.
    
    Below example shows how to deform a structure along a normal mode
    or linear combinations of normal modes:
    
    >>> deformAtoms(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5)
    >>> deformAtoms(p38_structure, -p38_pca[1] * p38_pca[1].getVariance()**0.5)
    >>> deformAtoms(p38_structure, p38_pca[0] * p38_pca[0].getVariance()**0.5 + 
    ...                            p38_pca[1] * p38_pca[1].getVariance()**0.5)
    >>> deformAtoms(p38_structure, p38_pca[0], rmsd=1.0)
    >>> print calcRMSD(p38_structure).round(3)
    [ 0.     0.41   0.308  0.513  1.   ]
    
    """

    if not isinstance(atoms, AtomGroup):
        raise TypeError('atoms must be an AtomGroup, not {0:s}'
                        .format(type(atoms)))
    if not isinstance(mode, VectorBase):
        raise TypeError('mode must be a Mode or Vector instance, '
                        'not {0:s}'.format(type(mode)))
    if not mode.is3d():
        raise ValueError('mode must be from a 3-dimensional model.')
    if atoms.numAtoms() != mode.numAtoms():
        raise ValueError('number of atoms do not match')
    
    array = mode.getArrayNx3()
    
    if rmsd is not None:
        rmsd = float(rmsd)
        # rmsd = ( ((scalar * array)**2).sum() / n_atoms )**0.5
        scalar = (atoms.numAtoms() * rmsd**2 / (array**2).sum())**0.5
        LOGGER.info('Mode is scaled by {0:g}.'.format(scalar))
        atoms.addCoordset( atoms.getCoords() + array * scalar)
    else:     
        atoms.addCoordset( atoms.getCoords() + array)

def scanPerturbationResponse(model, atoms=None, repeats=100):
    """See :func:`calcPerturbResponse`."""
    
    prody.deprecate('scanPerturbationResponse', 'calcPerturResponse')
    return calcPerturResponse(model, atoms, repeats) 

def calcPerturbResponse(model, atoms=None, repeats=100):
    """Return a matrix of profiles from scanning of the response of the 
    structure to random perturbations at specific atom (or node) positions. 
    The function implements the perturbation response scanning (PRS) method 
    described in [CA09]_.  Rows of the matrix are the average magnitude of the 
    responses obtained by perturbing the atom/node position at that row index, 
    i.e. ``prs_profile[i,j]`` will give the response of residue/node *j* to 
    perturbations in residue/node *i*.  PRS is performed using the covariance 
    matrix from *model*, e.t. :class:`ANM` instance.  Each residue/node is 
    perturbed *repeats* times with a random unit force vector.  When *atoms* 
    instance is given, PRS profile for residues will be added as an attribute 
    which then can be retrieved as ``atoms.getData('prs_profile')``.  *model* 
    and *atoms* must have the same number of atoms. *atoms* must be an :class:`
    ~prody.atomic.AtomGroup` instance.
    
    .. versionadded:: 0.8.2
    
    The RPS matrix can be save as follows::
        
     prs_matrix = calcPerturbationResponse(p38_anm)
     writeArray('prs_matrix.txt', prs_matrix, format='%8.6f', delimiter='\t')
    """
    
    if not isinstance(model, NMABase): 
        raise TypeError('model must be an NMA instance')
    elif not model.is3d():
        raise TypeError('model must be a 3-dimensional NMA instance')
    elif len(model) == 0:
        raise ValueError('model must have normal modes calculated')
    if atoms is not None:
        if not isinstance(atoms, prody.AtomGroup):
            raise TypeError('atoms must be an AtomGroup instance')
        elif atoms.numAtoms() != model.numAtoms():
            raise ValueError('model and atoms must have the same number atoms')
            
    assert isinstance(repeats, int), 'repeats must be an integer'
    cov = model.getCovariance()
    if cov is None:
        raise ValueError('model did not return a covariance matrix')
    
    n_atoms = model.numAtoms()
    response_matrix = np.zeros((n_atoms, n_atoms))
    LOGGER.progress('Calculating perturbation response', n_atoms)
    i3 = -3
    i3p3 = 0
    for i in range(n_atoms):
        i3 += 3
        i3p3 += 3
        forces = np.random.rand(repeats * 3).reshape((repeats, 3))
        forces /= ((forces**2).sum(1)**0.5).reshape((repeats, 1))
        for force in forces:
            response_matrix[i] += (np.dot(cov[:, i3:i3p3], force) ** 2
                                            ).reshape((n_atoms, 3)).sum(1)
        LOGGER.update(i)

    response_matrix /= repeats
    LOGGER.clear()
    LOGGER.info('Perturbation response scanning completed in %.1fs.')
    if atoms is not None:
        atoms.setData('prs_profile', response_matrix)
    return response_matrix
    
    # save the original PRS matrix
    np.savetxt('orig_PRS_matrix', response_matrix, delimiter='\t', fmt='%8.6f')
    # calculate the normalized PRS matrix
    self_dp = np.diag(response_matrix) # using self displacement (diagonal of
                               # the original matrix) as a
                               # normalization factor     
    self_dp = self_dp.reshape(n_atoms, 1)
    norm_PRS_mat = response_matrix / np.repeat(self_dp, n_atoms, axis=1)
    # suppress the diagonal (self displacement) to facilitate
    # visualizing the response profile
    norm_PRS_mat = norm_PRS_mat - np.diag(np.diag(norm_PRS_mat))
    np.savetxt('norm_PRS_matrix', norm_PRS_mat, delimiter='\t', fmt='%8.6f')
    return response_matrix
    
    
def calcSqFlucts(modes):
    """Return sum of square-fluctuations for given set of normal *modes*.
    
    .. versionchanged:: 0.7.1
       :class:`Vector` instances are accepted as *modes* argument."""
    
    if not isinstance(modes, (VectorBase, NMABase, ModeSet)):
        raise TypeError('modes must be a Mode, NMA, or ModeSet instance, '
                        'not {0:s}'.format(type(modes)))
    if isinstance(modes, Vector):
        if modes.is3d():
            return (modes._getArrayNx3()**2).sum(axis=1)
        else:
            return (modes._getArray() ** 2)
    else:
        square_fluctuations = np.zeros(modes.numAtoms()) 
        if isinstance(modes, VectorBase):
            modes = [modes]
        for mode in modes:
            square_fluctuations += mode.getSqFlucts()
        return square_fluctuations

def calcCrossCorrelations(modes, n_cpu=1):
    """Deprecated, use :func:`calcCrossCorr`."""
    
    prody.deprecate('calcCrossCorrelations', 'calcCrossCorr')
    return calcCrossCorr(modes, n_cpu)
 
def calcCrossCorr(modes, n_cpu=1):
    """Return cross-correlations matrix.  For a 3-d model, cross-correlations 
    matrix is an NxN matrix, where N is the number of atoms.  Each element of 
    this matrix is the trace of the submatrix corresponding to a pair of atoms.
    Covariance matrix may be calculated using all modes or a subset of modes
    of an NMA instance.  For large systems, calculation of cross-correlations 
    matrix may be time consuming.  Optionally, multiple processors may be 
    employed to perform calculations by passing ``n_cpu=2`` or more."""
    
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
            arvar = (array[:, indices]*variances[indices]).T.reshape(
                                                        (n_modes, n_atoms, 3))
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
    """Deprecated, use :func:`calcCumOverlap`."""
    
    prody.deprecate('calcCumulativeOverlap', 'calcCumOverlap')
    return calcCumOverlap(modes1, modes2) 
    
def calcCumOverlap(modes1, modes2):
    """Return cumulative overlap of modes in *modes2* with those in *modes1*.
    Returns a number of *modes1* contains a single :class:`Mode` or a 
    :class:`Vector` instance. If *modes1* contains multiple modes, returns an
    array. Elements of the array correspond to cumulative overlaps for modes 
    in *modes1* with those in *modes2*."""
    
    overlap = calcOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).sum(axis=overlap.ndim-1))
    return cumov

def calcCumulativeOverlapArray(modes1, modes2):
    """Deprecated, use :func:`calcCumOverlapArray`."""

    prody.deprecate('calcCumulativeOverlapArray', 'calcCumOverlapArray')
    return calcCumOverlapArray(modes1, modes2) 


def calcCumOverlapArray(modes1, modes2):
    """Return array of cumulative overlaps. Returned array has the shape 
    ``(len(modes1), len(modes2))``.  Each row corresponds to cumulative 
    overlaps calculated for modes in *modes1* with those in *modes2*.  
    Each value in a row corresponds to cumulative overlap calculated 
    using upto that many number of modes from *modes2*."""
    
    overlap = calcOverlap(modes1, modes2)
    cumov = np.sqrt(np.power(overlap, 2).cumsum(axis=overlap.ndim-1))
    return cumov


def calcSubspaceOverlap(modes1, modes2):
    """Return subspace overlap between two sets of modes (*modes1* and 
    *modes2*).  Also known as the root mean square inner product (RMSIP) 
    of essential subspaces [AA99]_.  This function returns a single number."""
    
    overlap = calcOverlap(modes1, modes2)
    if isinstance(modes1, Mode):
        length = 1
    else:
        length = len(modes1)
    rmsip = np.sqrt(np.power(overlap, 2).sum() / length)
    return rmsip


def calcCovarianceOverlap(modelA, modelB):
    """Deprecated, use :func:`calcCovOverlap`."""

    prody.deprecate('calcCovarianceOverlap', 'calcCovOverlap')    
    return calcCovOverlap(modelA, modelB)
    
def calcCovOverlap(modelA, modelB):
    """Return overlap between covariances of *modelA* and *modelB*.
    
    .. versionadded:: 0.5.3
    
    Overlap between covariances are calculated using normal modes 
    (eigenvectors), hence modes in both models must have been calculated.
    This function implements equation 11 in [BH02]_."""
    
    if not modelA.is3d() or not modelB.is3d(): 
        raise TypeError('both models must be 3-dimensional') 
    if len(modelA) == 0 or len(modelB) == 0:  
        raise TypeError('modes must be calculated for both models, '
                        'try calcModes method')
    if modelA.numAtoms() != modelB.numAtoms(): 
        raise ValueError('modelA and modelB must have same number of atoms')
    arrayA = modelA._getArray()
    varA = modelA.getVariances()
    arrayB = modelB._getArray()
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
    β-factors from *atoms*."""
    
    model = modes.getModel()
    if not isinstance(model, GNMBase):
        raise TypeError('modes must come from GNM or ANM')
    if model.numAtoms() != atoms.numAtoms():
        raise ValueError('modes and atoms must have same number of nodes')
    sqf = calcSqFlucts(modes)
    return sqf / ((sqf**2).sum()**0.5) * (atoms.getBetas()**2).sum()**0.5
    
    
def showFractOfVariances(modes, *args, **kwargs):
    """Deprecated, use :func:`showFractOfVar`."""
    
    prody.deprecate('showFractOfVariances', 'showFractOfVar')
    return showFractOfVar(modes, *args, **kwargs) 
    
def showFractOfVar(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.
    bar`.  Note that mode indices are incremented by 1.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showFractOfVar(p38_pca) 
       showCumFractOfVar(p38_pca)
      
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if not plt: prody.importPyPlot()
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
    """Deprecated, use :func:`showCumFractOfVar`."""
    
    prody.deprecate('showCumFractOfVariances', 'showCumFractOfVar')
    return showCumFractOfVar(modes, *args, **kwargs)

def showCumFractOfVar(modes, *args, **kwargs):
    """Show fraction of variances of *modes* using :func:`~matplotlib.pyplot.
    plot`.
    
    Note that mode indices are incremented by 1.
    See :func:`showFractOfVar` for an example."""
    
    if not plt: prody.importPyPlot()
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
    """Show a projection of conformational deviations onto up to three normal 
    modes from the same model.
    
    :arg ensemble: a :class:`~prody.ensemble.Ensemble` instance
    :arg modes: a :class:`Mode`, :class:`ModeSet`, or :class:`NMABase` instance
    
    .. versionchanged:: 0.8
       The projected values are by default converted to RMSD.  Pass 
       ``rmsd=False`` to use projection itself. :class:`Vector` instances 
       are accepted as *ensemble* argument to allow for projecting a 
       deformation vector onto normal modes.  
    
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
    
    if not plt: prody.importPyPlot()
    if not isinstance(ensemble, (prody.Ensemble, prody.Conformation, 
                                 prody.Vector)):
        raise TypeError('ensemble must be Ensemble, Conformation, or Vector, '
                        'not {0:s}'.format(type(ensemble)))
    if not isinstance(modes, (NMABase, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
                        .format(type(modes)))
    if not modes.is3d(): 
        raise Exception('modes must be 3-dimensional')
    
    if isinstance(modes, Mode) or (isinstance(modes, (ModeSet, NMABase)) and 
                                   len(modes)==1):
        if not isinstance(modes, Mode):
            modes = modes[0]
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True))
        show = plt.hist(projection.flatten(), *args, **kwargs)
        plt.xlabel('Mode {0:d} coordinate'.format(modes.getIndex()+1))
        plt.ylabel('Number of conformations')
    elif len(modes) == 2:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True))
        show = plt.plot(projection[:, 0], projection[:, 1], *args, **kwargs)
        modes = [m for m in modes]
        plt.xlabel('Mode {0:d} coordinate'.format(modes[0].getIndex()+1))
        plt.ylabel('Mode {0:d} coordinate'.format(modes[1].getIndex()+1))
    elif len(modes) == 3:
        if 'ls' not in kwargs:
            kwargs['ls'] = 'None'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        projection = calcProjection(ensemble, modes, kwargs.pop('rmsd', True)) 
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
    """Show a projection of conformational deviations onto modes from
    different models using :func:`~matplotlib.pyplot.plot`.  This function 
    differs from :func:`showProjection` by accepting modes from two different 
    models.
    
    :arg ensemble: Ensemble for which deviations will be projected
    :type ensemble: :class:`~prody.ensemble.Ensemble`
    :arg mode_x: Projection onto this mode will be shown along x-axis. 
    :type mode_x: :class:`Mode`
    :arg mode_y: Projection onto this mode will be shown along y-axis.
    :type mode_y: :class:`Mode`
    :arg scale: Scale width of the projection onto one of modes. 
                ``x`` and ``y`` are accepted.
    :type scale: str
    :arg scalar: Scalar factor for ``x`` or ``y``.  If ``scalar=None`` is 
        passed, best scaling factor will be calculated and printed on the
        console.
    :type scalar: float
    
    .. versionchanged:: 0.8
       The projected values are by default converted to RMSD. 
       Pass ``rmsd=False`` to calculate raw projection values.
       :class:`Vector` instances are accepted as *ensemble* argument to allow
       for projecting a deformation vector onto normal modes.  
    
    By default ``marker='o', ls='None'`` is passed to the plotting function 
    to disable lines.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5.2,4))
       showCrossProjection(p38_ensemble, p38_pca[0], p38_anm[2])
    
    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    |example| See :ref:`pca-xray-plotting` for a more elaborate example.
       
    """
    if not plt: prody.importPyPlot()
    if not isinstance(ensemble, (prody.Ensemble, prody.Conformation, 
                                 prody.Vector)):
        raise TypeError('ensemble must be Ensemble, Conformation, or Vector, '
                        'not {0:s}'.format(type(ensemble)))
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
    if scale is not None:
        assert isinstance(scale, str), 'scale must be a string'
        scale = scale.lower()
        assert scale in ('x', 'y'), 'scale must be x or y'
    if scalar is not None:
        assert isinstance(scalar, float), 'scalar must be a float'
    xcoords = calcProjection(ensemble, mode_x, kwargs.get('rmsd', True))
    ycoords = calcProjection(ensemble, mode_y, kwargs.pop('rmsd', True))
    if scale:
        if scalar is None:
            scalar = ((ycoords.max() - ycoords.min()) / 
                      (xcoords.max() - xcoords.min())) 
            scalar = scalar * np.sign(calcOverlap(mode_x, mode_y))
            if scale == 'x':
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_x, scalar))
            else:
                scalar = 1 / scalar
                LOGGER.info('Projection onto {0:s} is scaled by {1:.2f}'
                            .format(mode_y, scalar))
        if scale == 'x':
            xcoords = xcoords * scalar  
        else:
            ycoords = ycoords * scalar
    if 'ls' not in kwargs:
        kwargs['ls'] = 'None'
    if 'marker' not in kwargs:
        kwargs['marker'] = 'o'
    show = plt.plot(xcoords, ycoords, *args, **kwargs)
    plt.xlabel('{0:s} coordinate'.format(mode_x))
    plt.ylabel('{0:s} coordinate'.format(mode_y))
    return show

def showOverlapTable(rows, cols, *args, **kwargs):
    """Show overlap table using :func:`~matplotlib.pyplot.pcolor`.  *rows* and 
    *cols* are sets of normal modes, and correspond to rows and columns of the 
    displayed matrix.  Note that mode indices are increased by 1.  List of 
    modes should contain a set of contiguous modes from the same model. 
    
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
    
    if not plt: prody.importPyPlot()
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
    """Deprecated, use :func:`showCrossCorr`."""
    
    prody.deprecate('showCrossCorrelations', 'showCrossCorr')
    return showCrossCorr(modes, *args, **kwargs)

def showCrossCorr(modes, *args, **kwargs):
    """Show cross-correlations for given modes using :func:`~matplotlib.pyplot.
    imshow`.  By default, *origin=lower* and *interpolation=bilinear* keyword 
    arguments are passed to imshow function. User can overwrite these 
    parameters.  See also :func:`getCrossCorr`.
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(6,5))
       # Show cross-correlations for ANM modes 1-3
       showCrossCorr( p38_anm[:3] )
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if not plt: prody.importPyPlot()
    arange = np.arange(modes.numAtoms())
    cross_correlations = np.zeros((arange[-1]+2, arange[-1]+2))
    cross_correlations[arange[0]+1:, 
                       arange[0]+1:] = calcCrossCorr(modes)
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
       plt.legend(loc='lower right', prop={'size': 10})
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
       
    """
    
    if not plt: prody.importPyPlot()
    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, '
                        'not {0:s}'.format(type(mode)))
    if mode.is3d():
        a3d = mode.getArrayNx3()
        show = plt.plot(a3d[:, 0], *args, label='x-component', **kwargs)
        plt.plot(a3d[:, 1], *args, label='y-component', **kwargs)
        plt.plot(a3d[:, 2], *args, label='z-component', **kwargs)
    else:
        show = plt.plot(mode._getArray(), *args, **kwargs)
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
       plt.legend(prop={'size': 10})
       
    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')
    
    """
    
    if not plt: prody.importPyPlot()
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
       plt.legend(prop={'size': 10})

    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    """
    
    if not plt: prody.importPyPlot()
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
    """Show normalized square fluctuations using :func:`~matplotlib.pyplot.
    plot`.
    
    .. plot::
       :context:
       :include-source:
       
       plt.figure(figsize=(5,4))
       showNormedSqFlucts(p38_pca[0], p38_anm[2])
       plt.legend(prop={'size': 10})

    .. plot::
       :context:
       :nofigs:

       plt.close('all')
       
    """
    
    if not plt: prody.importPyPlot()
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
    
    if not plt: prody.importPyPlot()
    if not isinstance(enm, GNMBase):
        raise TypeError('model argument must be an ENM instance')
    kirchhoff = enm.getKirchhoff()
    if kirchhoff is None:
        LOGGER.warning('kirchhoff matrix is not set')
        return None
    show = plt.spy(kirchhoff, *args, **kwargs)
    plt.title('{0:s} contact map'.format(enm.getTitle())) 
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
    
    if not plt: prody.importPyPlot()
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
    """Deprecated, use :func:`showCumOverlap`."""
    
    prody.deprecate('showCumulativeOverlap', 'showCumOverlap')
    return showCumOverlap(mode, modes, *args, **kwargs)
    
def showCumOverlap(mode, modes, *args, **kwargs):
    """Show cumulative overlap using :func:`~matplotlib.pyplot.plot`.
    
    :type mode: :class:`Mode`, :class:`Vector` 
    :arg modes: multiple modes
    :type modes: :class:`ModeSet`, :class:`ANM`, :class:`GNM`, or :class:`PCA` 
    
    .. plot::
       :context:
       :include-source:
        
       plt.figure(figsize=(5,4))
       showCumOverlap( p38_pca[0], p38_anm )
       # Let's also show the overlap
       showOverlap( p38_pca[0], p38_anm )

    .. plot::
       :context:
       :nofigs:
        
       plt.close('all')  
          
    """
    
    if not plt: prody.importPyPlot()
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
    """Reset X (and Y) axis ticks using values in given *array*.  Ticks in the 
    current figure should not be fractional values for this function to work as
    expected."""
    
    if x is not None:
        try:    
            xticks = plt.xticks()[0]
            xlist = list(xticks.astype(int))
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
            ylist = list(yticks.astype(int))
            if ylist[-1] > len(y):
                ylist.pop()
            if ylist:
                ylist = list(y[ylist]) 
                plt.yticks(yticks, ylist + [''] * (len(yticks) - len(ylist)))
        except:
            LOGGER.warning('xticks could not be reset.')
    
