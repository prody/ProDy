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

from numpy import array

__all__ = ['ModeSet']

class ModeSet(object):

    """A class for providing access to subset of mode data.  Instances 
    are obtained by slicing an NMA model (:class:`.ANM`, :class:`.GNM`, or 
    :class:`.PCA`).  ModeSet's contain a reference to the model and a list 
    of mode indices.  Methods common to NMA models are also defined for 
    mode sets."""
    
    __slots__ = ['_model', '_indices']
    
    def __init__(self, model, indices):
        self._model = model
        self._indices = array(indices, int)
        
    def __len__(self):
        return len(self._indices)
        
    def __iter__(self):
        for i in self._indices:
            yield self._model[i]
    
    def __repr__(self):
        return '<ModeSet: {0:d} modes from {1:s}>'.format(len(self),
                                                       str(self._model))

    def __str__(self):
        return '{0:d} modes from {1:s}'.format(len(self._indices), 
                                               str(self._model))
    
    def is3d(self):
        """Return **True** is model is 3-dimensional."""
                
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

    def getTitle(self):
        """Return title of the mode set."""
        
        return str(self)
    
    def getModel(self):
        """Return the model that the modes belongs to."""
        
        return self._model
    
    def getIndices(self):
        """Return indices of modes in the mode set."""
        
        return self._indices
        
    def getEigvals(self):
        """Return eigenvalues.  For :class:`.PCA` and :class:`.EDA` models 
        built using coordinate data in Å, unit of eigenvalues is |A2|.  For
        :class:`.ANM` and :class:`.GNM`, on the other hand, eigenvalues are 
        in arbitrary or relative units but they correlate with stiffness of 
        the motion along associated eigenvector."""
        
        return self._model._eigvals[self._indices]

    def getVariances(self):
        """Return variances.  For :class:`.PCA` and :class:`.EDA` models 
        built using coordinate data in Å, unit of variance is |A2|.  For
        :class:`.ANM` and :class:`.GNM`, on the other hand, variance is the 
        inverse of the eigenvalue, so it has arbitrary or relative units."""
        
        return self._model._vars[self._indices]

    def getArray(self):
        """Return a copy of eigenvectors array."""
        
        return self._model._array[:, self._indices]

    getEigvecs = getArray

    def _getArray(self):
        """Return eigenvectors array."""
        
        return self._model._array[:, self._indices]
