# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2013 Ahmet Bakan
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

"""This module defines a class handling normal mode analysis data."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2013 Ahmet Bakan'

import numpy as np

from .mode import Mode
from .modeset import ModeSet

from prody import PY2K

if PY2K:
    range = xrange

__all__ = ['NMA']

class NMA(object):

    """A class for handling Normal Mode Analysis (NMA) data."""

    def __init__(self, title='Unknown'):

        self._title = str(title).strip()
        self._n_modes = 0
        self._cov = None
        self._n_atoms = 0
        self._dof = 0
        self._array = None      # modes/eigenvectors
        self._eigvals = None
        self._vars = None       # evs for PCA, inverse evs for ENM
        self._trace = None
        self._is3d = True       # is set to false for GNM
        self._indices = None

    def __len__(self):

        return self._n_modes

    def __getitem__(self, index):
        """A list or tuple of integers can be used for indexing."""

        if self._n_modes == 0:
            raise ValueError('{0} modes are not calculated, use '
                             'calcModes() method'.format(str(self)))
        if isinstance(index, slice):
            indices = np.arange(*index.indices(len(self)))
            if len(indices) > 0:
                return ModeSet(self, indices)
        elif isinstance(index, (list, tuple)):
            for i in index:
                assert isinstance(i, int), 'all indices must be integers'
            if len(index) == 1:
                return self._getMode(index[0])
            return ModeSet(self, index)
        try:
            index = int(index)
        except Exception:
            raise IndexError('indices must be int, slice, list, or tuple')
        else:
            return self._getMode(index)

    def __iter__(self):

        for i in range(self._n_modes):
            yield self[i]

    def __repr__(self):

        return '<{0}: {1} ({2} modes; {3} atoms)>'.format(
                self.__class__.__name__, self._title, self._n_modes,
                self._n_atoms)

    def __str__(self):

        return self.__class__.__name__ + ' ' + self._title

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

    def _getMode(self, index):

        if self._n_modes == 0:
            raise ValueError('{0} modes are not calculated, use '
                             'calcModes() method'.format(str(self)))
        if index >= self._n_modes or index < -self._n_modes:
            raise IndexError('{0} contains {1} modes, try 0 <= index < '
                             '{1}'.format(str(self), self._n_modes))
        if index < 0:
            index += self._n_modes
        return Mode(self, index)

    def _getTrace(self):
        """Return trace, and emit a warning message if trace is calculated
        using eigenvalues of a subset of variances (eigenvalues or inverse
        eigenvalues)."""

        trace = self._trace
        if trace is None:
            if self._vars is None:
                raise ValueError('variances are not set or calculated')
            trace = self._vars.sum()
            diff = self._dof - self._n_modes
            if self._is3d and diff > 6:
                diff = True
            elif diff > 1:
                diff = True
            else:
                diff = False
            if diff:
                from prody import LOGGER
                LOGGER.warn('Total variance for {0} is calculated using '
                            '{1} available modes out of {2} possible.'
                            .format(str(self), self._n_modes, self._dof))
        return trace

    def getModel(self):
        """Return self."""

        return self

    def is3d(self):
        """Return **True** is model is 3-dimensional."""

        return self._is3d

    def numAtoms(self):
        """Return number of atoms."""

        return self._n_atoms

    def numModes(self):
        """Return number of modes in the instance (not necessarily maximum
        number of possible modes)."""

        return self._n_modes

    def numDOF(self):
        """Return number of degrees of freedom."""

        return self._dof

    def getTitle(self):
        """Return title of the model."""

        return self._title

    def setTitle(self, title):
        """Set title of the model."""

        self._title = str(title)

    def getEigvals(self):
        """Return eigenvalues.  For :class:`.PCA` and :class:`.EDA` models
        built using coordinate data in Å, unit of eigenvalues is |A2|.  For
        :class:`.ANM` and :class:`.GNM`, on the other hand, eigenvalues are
        in arbitrary or relative units but they correlate with stiffness of
        the motion along associated eigenvector."""

        if self._eigvals is None: return None
        return self._eigvals.copy()

    def getVariances(self):
        """Return variances.  For :class:`.PCA` and :class:`.EDA` models
        built using coordinate data in Å, unit of variance is |A2|.  For
        :class:`.ANM` and :class:`.GNM`, on the other hand, variance is the
        inverse of the eigenvalue, so it has arbitrary or relative units."""

        if self._vars is None: return None
        return self._vars.copy()

    def getArray(self):
        """Return a copy of eigenvectors array."""

        if self._array is None: return None
        return self._array.copy()

    getEigvecs = getArray

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
        """"""

        pass

    def addEigenpair(self, vector, value=None):
        """Add eigen *vector* and eigen *value* pair(s) to the instance.
        If eigen *value* is omitted, it will be set to 1.  Inverse
        eigenvalues are set as variances."""

        if self._array is None:
            self.setEigens()

        try:
            ndim, shape = vector.ndim, vector.shape
        except AttributeError:
            raise TypeError('vector must be a numpy array')

        if ndim > 2:
            raise ValueError('vector must be 1- or 2-dimensional')
        elif ndim == 2 and shape[0] < shape[1]:
            raise ValueError('eigenvectors must correspond to vector columns')
        else:
            vector = vector.reshape((shape[0], 1))

        if eigval is None:
            if ndim == 1:
                eigval = np.ones(1)
            else:
                eigval = np.ones(shape[2])
        else:
            try:
                ndim2, shape2 = eigval.ndim, eigval.shape
            except AttributeError:
                try:
                    eigval = np.array([value], float)
                except Exception:
                    raise TypeError('value must be a number or array')
            else:
                if value.ndim > 1:
                    raise ValueError('value must be a 1-dimensional array')
                elif value.shape[0] != value.shape[0]:
                    raise ValueError('number of eigenvectors and eigenvalues '
                                     'must match')

        if vector.shape[0] != self._array.shape[0]:
            raise ValueError('shape of vector do not match shape of '
                             'existing eigenvectors')
        self._array = np.concatenate((self._array, vector), 1)
        self._eigvals = np.concatenate((self._eigvals, value))
        self._n_modes += shape[1]
        self._vars = 1 / self._eigvals

    def setEigens(self, vectors, values=None):
        """Set eigen *vectors* and eigen *values*.  If eigen *values* are
        omitted, they will be set to 1.  Inverse eigenvalues are set as
        variances."""

        try:
            ndim, shape = vectors.ndim, vectors.shape
        except AttributeError:
            raise TypeError('vectors must be a numpy array')

        if ndim != 2:
            raise ValueError('vectors must be a 2-dimensional array')
        else:
            dof = shape[0]
            if self._is3d:
                n_atoms = dof / 3
            else:
                n_atoms = dof
            if self._n_atoms > 0 and n_atoms != self._n_atoms:
                    raise ValueError('vectors do not have the right shape, '
                                 'which is (M,{0})'.format(n_atoms*3))
            n_modes = vectors.shape[1]
        if values is not None:
            if not isinstance(vectors, np.ndarray):
                raise TypeError('values must be a numpy.ndarray, not {0}'
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


