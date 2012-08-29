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

from .nma import NMA

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
        self._indices = array(indices, int)
        
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
        
        return self._model._is3d
    
    is3d.__doc__ = NMA.is3d.__doc__
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._model._n_atoms
    
    def numModes(self):
        
        return len(self._indices)

    numModes.__doc__ = NMA.numModes.__doc__
    
    def numDOF(self):
        
        return self._model._dof

    numDOF.__doc__ = NMA.numDOF.__doc__

        
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
        
    def getEigvals(self):
        
        return self._model._eigvals[self._indices]

    getEigvals.__doc__ = NMA.getEigvals.__doc__

    def getVariances(self):
        
        return self._model._vars[self._indices]

    getVariances.__doc__ = NMA.getVariances.__doc__

    def getArray(self):
        
        return self._model._array[:, self._indices]

    getArray.__doc__ = NMA.getArray.__doc__

    getEigvecs = getArray

    def _getArray(self):

        return self._model._array[:, self._indices]

    _getArray.__doc__ = NMA._getArray.__doc__
