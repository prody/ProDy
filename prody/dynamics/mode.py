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

"""This module defines classes for handling mode data.

.. _normalmode-operations:

*******************************************************************************
Normal Mode Algebra 
*******************************************************************************

In this example we will compare modes from two ANMs for the same protein, 
but everything applies to comparison of ANMs and PCAs (as long as they contain 
same number of atoms).

Let's get started by getting ANM models for two related protein structures:

>>> from prody import *
>>> str1 = parsePDB('1p38')
>>> str2 = parsePDB('1r39')

**Find and align matching chains**

>>> matches = matchChains(str1, str2)
>>> match = matches[0]
>>> ch1 = match[0]
>>> ch2 = match[1]

Minimize RMSD by superposing ``ch2`` onto ``ch1``:

>>> ch2, t = superpose(ch2, ch1) 
>>> # t is transformation, which is already applied to ch2
>>> rmsd = calcRMSD(ch1, ch2)
>>> print( '{0:.2f}'.format(rmsd) ) # Print rmsd with some formatting
0.90


**Get ANM models for each chain**

>>> anm1, ch1 = calcANM(ch1)
>>> anm2, ch2 = calcANM(ch2)

>>> print( anm1[0] )
Mode 1 from ANM 1p38

Let's rename these :class:`ANM` instances, so that they print short: 

>>> anm1.setTitle('1p38_anm')
>>> anm2.setTitle('1r39_anm')

This is how they print now:

>>> print( anm1[0] )
Mode 1 from ANM 1p38_anm
>>> print( anm2[0] )
Mode 1 from ANM 1r39_anm

Calculate overlap
===============================================================================

We need Numpy in this part:

>>> import numpy as np
>>> np.set_printoptions(precision=3)

Multiplication of two :class:`Mode` instances returns dot product
of their eigenvectors. This dot product is the overlap or cosine correlation
between modes.

Let's calculate overlap for slowest modes:

>>> overlap = anm1[0] * anm2[0]
>>> print( '{0:.3f}'.format(overlap) ) 
-0.984

This show that the overlap between these two modes is 0.98, which is not 
surprising since ANM modes come from structures of the *same* protein.

To compare multiple modes, convert a list of modes to a :func:`numpy.array`:

>>> print( np.array(list(anm1[:3])) * np.array(list(anm2[:3])) )
[-0.98402119545 -0.98158348545 -0.991357811832]

This shows that slowest three modes are almost identical.

We could also generate a matrix of overlaps using :func:`numpy.outer`:

>>> outer = np.outer( np.array(list(anm1[:3])),  np.array(list(anm2[:3])) )
>>> print( outer.astype(np.float64).round(2) )
[[-0.98 -0.14 -0.  ]
 [ 0.15 -0.98  0.08]
 [ 0.01 -0.08 -0.99]]

This could also be printed in a pretty table format using 
:func:`~.printOverlapTable`:

>>> printOverlapTable(anm1[:3], anm2[:3])
Overlap Table
                      ANM 1r39_anm
                    #1     #2     #3
ANM 1p38_anm #1   -0.98  -0.14   0.00
ANM 1p38_anm #2   +0.15  -0.98  +0.08
ANM 1p38_anm #3   +0.01  -0.08  -0.99
<BLANKLINE>


**Scaling**

:class:`Mode` instances can be scaled, but after this operation they will
become :class:`Vector` instances:

>>> anm1[0] * 10
<Vector: 10*(Mode 1 from ANM 1p38_anm)>

Linear combination
===============================================================================

It is also possible to linearly combine normal modes:

>>> anm1[0] * 3 + anm1[1] + anm1[2] * 2
<Vector: 3*(Mode 1 from ANM 1p38_anm) + Mode 2 from ANM 1p38_anm + 2*(Mode 3 from ANM 1p38_anm)>

Or, we could use eigenvalues for linear combination:

>>> lincomb = anm1[0] * anm1[0].getEigval() + anm1[1] * anm1[1].getEigval()

It is the name of the :class:`Vector` instance that keeps track of operations.

>>> print( lincomb.getTitle() )  
0.148971269751*(Mode 1 from ANM 1p38_anm) + 0.24904210757*(Mode 2 from ANM 1p38_anm)

Approximate a deformation vector
===============================================================================

Let's get the deformation vector between *ch1* and *ch2*:

>>> defvec = calcDeformVector(ch1, ch2)
>>> defvec_magnitude = abs(defvec)
>>> print( '{0:.2f}'.format(defvec_magnitude) )
16.69

Let's see how deformation projects onto ANM modes:

>>> print( np.array(list(anm1[:3])) * defvec )
[-5.60860594784 2.15393365959 -3.13701609199]

We can use these numbers to combine ANM modes:

>>> approximate_defvec = np.sum( (np.array(list(anm1[:3])) * defvec) * np.array(list(anm1[:3])) ) 
>>> print( approximate_defvec )
-5.60860594784*(Mode 1 from ANM 1p38_anm) + 2.15393365959*(Mode 2 from ANM 1p38_anm) + -3.13701609199*(Mode 3 from ANM 1p38_anm)

Let's deform 1r39 chain along this approximate deformation vector and see
how RMSD changes:

>>> ch2.setCoords(ch2.getCoords() - approximate_defvec.getArrayNx3())
>>> rmsd = calcRMSD(ch1, ch2)
>>> print( '{0:.2f}'.format(rmsd) )
0.82

RMSD decreases from 0.89 A to 0.82 A.

""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np
import prody

__all__ = ['Mode', 'Vector']
           
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
        
    def _getArray(self):
        """Return array."""

        pass
    
    def numAtoms(self):
        """Return number of atoms."""
        
        pass
    
    def is3d(self):
        """Return true if vector is 3d."""
        
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
        :type model: :class:`~.ANM`, :class:`~.GNM`, or :class:`~.PCA` 
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
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._model._n_atoms
    
    def numDOF(self):
        """Return number of degrees of freedom (three times the number of 
        atoms)."""
        
        return self._model._dof
    
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
    
    getEigvec = getArray
    
    def getEigenvector(self):
        """Deprecated, use :meth:`getEigvec` instead."""
        
        prody.deprecate('getEigenvector', 'getEigvec')
        return self.getEigvec()

    def _getArray(self):
        """Return a copy of the normal mode array (eigenvector)."""
    
        return self._model._array[:, self._index]
    
    def getEigenvalue(self):
        """Deprecated, use :meth:`getEigval` instead."""
        
        prody.deprecate('getEigenvalue', 'getEigval')
        return self.getEigval()
        
    def getEigval(self):
        """Return normal mode eigenvalue."""
        
        return self._model._eigvals[self._index]
    
    def getVariance(self):
        """Variance along the mode.  If the model is not a PCA, inverse of the
        eigenvalue is returned."""
        
        return self._model._vars[self._index]

    def getCovariance(self):
        """Deprecated, use :func:`~.calcCovariance` instead."""
        
        prody.deprecate('getCovariance', 'calcCovariance')
        return prody.calcCovariance(self)


class Vector(VectorBase):
    
    """A class to provide operations on a modified mode array.  This class 
    holds only mode array (i.e. eigenvector) data, and has no associations 
    with an NMA instance.  Scalar multiplication of :class:`Mode` instance 
    or addition of two :class:`Mode` instances results in a :class:`Vector` 
    instance."""
    
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
    
    def getTitle(self):
        """Get the descriptive title for the vector instance."""
        
        return self._title
    
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

    def numDOF(self):

        """Return number of degrees of freedom."""
        return len(self._array)

    def numAtoms(self):
        """Return number of atoms.  For a 3-dimensional vector, returns length 
        of the vector divided by 3."""
        
        if self._is3d: 
            return len(self._array)/3
        else:
            return len(self._array)
