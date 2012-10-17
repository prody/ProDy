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

"""This module defines a class handling normal mode analysis data.""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from .mode import Mode

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
            raise ValueError('{0:s} modes are not calculated, try '
                             'calcModes() method'.format(str(self)))
        indices = self._indices
        if indices is None:
            self._indices = indices = np.arange(self._n_modes)
        try:
            indices = indices[index]
        except IndexError as err:
            raise IndexError(str(err))
        try:
            length = len(indices)
        except TypeError:            
            return Mode(self, indices)
        else:
            return ModeSet(self, indices)
        
    def __iter__(self):
        
        for i in xrange(self._n_modes):
            yield self[i]
    
    def __repr__(self):
        
        return '<{0:s}: {1:s} ({2:d} modes; {3:d} atoms)>'.format(
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
    
    def addEigenpair(self, eigenvector, eigenvalue=None):
        """Add *eigenvector* and *eigenvalue* pair to the :class:`NMA` 
        instance.  If *eigenvalue* is not given, it will be set to 1. 
        Variances are set as the inverse eigenvalues."""
        
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
        
        :arg vectors: eigenvectors
        :type vectors: numpy.ndarray
        
        :arg values: Eigenvalues. When ``None`` is passed (default value), 
            all eigenvalues will be set to ``1``.
        :type values: numpy.ndarray
        
        For M modes and N atoms, *vectors* must have shape ``(3*N, M)``
        and values must have shape ``(M,)``.
        
        Variances are set as the inverse eigenvalues."""
        
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


