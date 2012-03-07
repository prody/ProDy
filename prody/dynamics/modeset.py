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

"""This module defines a pointer class for handling subsets of normal modes."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np
import prody

__all__ = ['ModeSet']

class ModeSet(object):

    """A class for providing access to subset of mode data.  Instances 
    are obtained by slicing an NMA model (:class:`ANM`, :class:`GNM`, or 
    :class:`PCA`).  ModeSet's contain a reference to the model and a list 
    of mode indices.  Methods common to NMA models are also defined for 
    mode sets."""
    
    __slots__ = ['_model', '_indices']
    
    def __init__(self, model, indices):
        if not isinstance(model, NMA):
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
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._model._n_atoms
    
    def numModes(self):
        """Return number of modes in the instance (not necessarily maximum 
        number of possible modes)."""
        
        return len(self._indices)
    
    def numDOF(self):
        """Return number of degrees of freedom."""
        
        return self._model._dof
        
    def getModes(self):
        """Return a list that contains the modes in the mode set."""
        
        getMode = self._model.getMode
        return [getMode(i) for i in self._indices]
    
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
        """Deprecated, use :meth:`getEigvals` instead."""
        
        prody.deprecate('getEigenvalue', 'getEigvals')
        return self.getEigvals()
        
    def getEigvals(self):
        """Return eigenvalues."""
        
        return self._model._eigvals[self._indices]

    def getEigenvectors(self):
        """Deprecated, use :meth:`getEigvecs` instead."""
        
        prody.deprecate('getEigenvector', 'getEigvecs')
        return self.getEigvecs()
    
    def getVariances(self):
        """Return variances (~inverse eigenvalues)."""
        
        return self._model._vars[self._indices]

    def getArray(self):
        """Return a copy of eigenvectors array."""
        
        return self._model._array[:, self._indices]

    getEigvecs = getArray

    def _getArray(self):
        """Return a copy of eigenvectors array."""

        return self._model._array[:, self._indices]
        
    def getCovariance(self):
        """Deprecated, use :func:`~.calcCovariance` instead."""
        
        prody.deprecate('getCovariance', 'calcCovariance')
        return prody.calcCovariance(self)
