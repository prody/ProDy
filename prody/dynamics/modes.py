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

""":mod:`modes` module defines class for accessing individual normal modes and 
for defining operations on modes. 

Classes:

  * :class:`VectorBase`
  * :class:`Vector`
  * :class:`Mode`

"""
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np

__all__ = ['Mode', 'Vector']


class VectorBase(object):
    """A base class for :class:`Mode` and :class:`Vector`.
    
    This base class defines some shared methods, such as scalar multiplication
    or addition of mode instances.
    
    Defined functions are:
        
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
        """Returns scaled mode or dot product between modes."""
        if isinstance(other, (int, float, long)): 
            return Vector('{0}*({1:s})'.format(other, str(self)),
                              other * self.getArray(), self.is3d())
        elif isinstance(other, VectorBase):
            return np.dot(self.getArray(), other.getArray())
        else:
            raise TypeError('{0} is not a scalar or a mode'.format(other))
    
    def __rmul__(self, other):   
        """Returns scaled mode or dot product between modes."""
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
        """Returns mode after normalizing it."""
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
