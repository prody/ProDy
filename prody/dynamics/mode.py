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

"""This module defines classes for handling mode data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

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
        * Power (mode**x)"""

    __slots__ = []

    def __abs__(self):

        return np.sqrt((self._getArray()**2).sum())

    def __neg__(self):

        return Vector(-self._getArray(), '-({0})'.format(str(self)),
                      self.is3d())

    def __div__(self, other):

        try:
            result = self._getArray() / other
        except Exception as err:
            raise TypeError('{0} is not a scalar {1}'
                            .format(other, str(err)))
        return Vector(result, '({0})/{1}'.format(str(self), other),
                      self.is3d())

    def __idiv__(self, other):

        return self.__div__(other)

    def __mul__(self, other):
        """Return scaled mode or dot product between modes."""

        try:
            other = other._getArray()
        except AttributeError:
            try:
                result = other * self._getArray()
            except Exception as err:
                raise TypeError('{0} is not a scalar or a mode ({1})'
                                .format(other, str(err)))
            else:
                return Vector(result, '({1})*{0}'.format(other, str(self)),
                              self.is3d())
        else:
            try:
                return np.dot(self._getArray(), other)
            except Exception:
                raise ValueError('{0} and {1} do not have same dimensions '
                                 '({2})'.format(str(self), str(other),
                                                str(err)))

    def __rmul__(self, other):
        """Return scaled mode or dot product between modes."""

        try:
            other = other._getArray()
        except AttributeError:
            try:
                result = other * self._getArray()
            except Exception as err:
                raise TypeError('{0} is not a scalar or a mode ({1})'
                                .format(other, str(err)))
            else:
                return Vector(result, '{0}*({1})'.format(other, str(self)),
                              self.is3d())
        else:
            try:
                return np.dot(self._getArray(), other)
            except Exception:
                raise ValueError('{0} and {1} do not have same dimensions '
                                 '({2})'.format(str(self), str(other),
                                                str(err)))

    def __imul__(self, other):

        return self.__mul__(other)

    def __add__(self, other):

        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self._getArray() + other._getArray(),
                          '({0}) + ({1})'.format(str(self), str(other)),
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __radd__(self, other):

        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(self._getArray() + other._getArray(),
                          '({0}) + ({1})'.format(str(other), str(self)),
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
                          '({0}) - ({1})'.format(str(self), str(other)),
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __rsub__(self, other):

        if isinstance(other, VectorBase):
            if len(self) != len(other):
                raise ValueError('modes do not have the same length')
            return Vector(other._getArray() - self._getArray(),
                          '({0}) - ({1})'.format(str(other), str(self)),
                          self.is3d())
        else:
            raise TypeError('{0} is not a mode instance'.format(other))

    def __isub__(self, other):

        return self.__sub__(other)

    def __pow__(self, other):

        try:
            result = self._getArray() ** other
        except Exception as err:
            raise TypeError('{0} is not a scalar ({0})'
                            .format(other, str(err)))
        else:
            return Vector(result, '({0})**{1}'.format(str(self), other),
                          self.is3d())

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
        """Return **True** if vector is 3d."""

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

    def numModes(self):
        """Return 1."""

        return 1


class Mode(VectorBase):

    """A class to provide access to and operations on mode data.

    """

    __slots__ = ['_model', '_index']

    def __init__(self, model, index):
        """Initialize mode object as part of an NMA model.

        :arg model: a normal mode analysis instance
        :type model: :class:`.NMA`, :class:`.GNM`, or :class:`.PCA`
        :arg index: index of the mode
        :type index: int"""

        self._model = model
        self._index = int(index)

    def __len__(self):

        return self._model.numDOF()

    def __repr__(self):

        return '<Mode: {0} from {1}>'.format(self._index + 1,
                                             str(self._model))

    def __str__(self):

        return 'Mode {0} from {1}'.format(self._index+1, str(self._model))

    def __int__(self):

        return self._index

    def __float__(self):

        return self.getEigval()

    def is3d(self):
        """Return **True** if mode instance is from a 3-dimensional model."""

        return self._model.is3d()

    def numAtoms(self):
        """Return number of atoms."""

        return self._model.numAtoms()

    def numDOF(self):
        """Return number of degrees of freedom (three times the number of
        atoms)."""

        return self._model.numDOF()

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

    def _getArray(self):
        """Return a copy of the normal mode array (eigenvector)."""

        return self._model._array[:, self._index]

    def getEigval(self):
        """Return normal mode eigenvalue.  For :class:`.PCA` and :class:`.EDA`
        models built using coordinate data in Å, unit of eigenvalues is |A2|.
        For :class:`.ANM` and :class:`.GNM`, on the other hand, eigenvalues
        are in arbitrary or relative units but they correlate with stiffness
        of the motion along associated eigenvector."""

        return self._model._eigvals[self._index]

    def getVariance(self):
        """Return variance along the mode.  For :class:`.PCA` and :class:`.EDA`
        models built using coordinate data in Å, unit of variance is |A2|.  For
        :class:`.ANM` and :class:`.GNM`, on the other hand, variance is the
        inverse of the eigenvalue, so it has arbitrary or relative units."""

        return self._model._vars[self._index]


class Vector(VectorBase):

    """A class to provide operations on a modified mode array.  This class
    holds only mode array (i.e. eigenvector) data, and has no associations
    with an NMA instance.  Scalar multiplication of :class:`Mode` instance
    or addition of two :class:`Mode` instances results in a :class:`Vector`
    instance."""

    __slots__ = ['_title', '_array', '_is3d']

    def __init__(self, array, title='Unknown', is3d=True):
        """Instantiate with a name, an array, and a 3d flag."""

        try:
            ndim, shape = array.ndim, array.shape
        except AttributeError:
            array = np.array(array)
            ndim, shape = array.ndim, array.shape
        if ndim != 1:
            raise ValueError('array.ndim must be 1')
        is3d = bool(is3d)
        if is3d and shape[0] % 3 != 0:
            raise ValueError('len(array) must be a multiple of 3')
        self._title = str(title)
        self._array = array
        self._is3d = is3d

    def __len__(self):
        return len(self._array)

    def __repr__(self):
        return '<Vector: {0}>'.format(self._title)

    def __str__(self):
        return self._title

    def is3d(self):
        """Return **True** if vector instance describes a 3-dimensional
        property, such as a deformation for a set of atoms."""

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
                      '({0})/||{0}||'.format(self._title), self._is3d)

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
