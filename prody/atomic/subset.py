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

import numpy as np

from atomic import ATOMIC_DATA_FIELDS
from atomic import wrapGetMethod, wrapSetMethod
from pointer import AtomPointer, MCAtomPointer


class AtomSubsetMeta(type):

    def __init__(cls, name, bases, dict):

        for field in ATOMIC_DATA_FIELDS.values():
            meth = field.meth_pl
            getMeth = 'get' + meth
            setMeth = 'set' + meth
            # Define public method for retrieving a copy of data array
            if field.call:
                def getData(self, var=field.var, call=field.call):
                    for meth in call:
                        getattr(self._ag, meth)()
                    array = self._ag._data[var]
                    return array[self._indices]
            else:
                def getData(self, var=field.var):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    return array[self._indices] 
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            getData.__doc__ = field.getDocstr('get')
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)
            
            if field.readonly:
                continue
            
            # Define public method for setting values in data array
            def setData(self, value, var=field.var, none=field.none):
                array = self._ag._data[var]
                if array is None:
                    raise AttributeError(var + ' data is not set')
                array[self._indices] = value
                if none:
                    self._ag.__setattr__('_'+none,  None)
            setData = wrapSetMethod(setData)
            setData.__name__ = setMeth 
            setData.__doc__ = field.getDocstr('set')  
            setattr(cls, setMeth, setData)

                        
class AtomSubset(AtomPointer):
    
    """A class for manipulating subset of atomic data in an :class:`AtomGroup`.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, and active coordinate set index for the atom group.
    
    """
    
    __metaclass__ = AtomSubsetMeta    
    
    def __init__(self, ag, indices, **kwargs):
        
        AtomPointer.__init__(self, ag)
                
        if not isinstance(indices, np.ndarray):
            indices = np.array(indices, int)
        elif not indices.dtype == int:
            indices = indices.astype(int)
            
        if kwargs.get('unique'):
            self._indices = indices
        else:
            self._indices = np.unique(indices)
            
        self._selstr = kwargs.get('selstr')
    
    def __len__(self):
        return len(self._indices)

    def __invert__(self):
        
        ones = np.ones(self._ag.numAtoms(), bool)
        ones[self._indices] = False
        sel = Selection(self._ag, ones.nonzero()[0], 
                        "not ({0:s}) ".format(self.getSelstr()))
        return sel
    
    def __or__(self, other):
        
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
            
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
            
        if self is other:
            return self

        if isinstance(other, Atom):
            other_indices = np.array([other._index])
        else:
            other_indices = other._indices
            
        indices = np.unique(np.concatenate((self._indices, other_indices)))
        return Selection(self._ag, indices, '({0:s}) or ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()))

    def __and__(self, other):
        
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
            
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
            
        if self is other:
            return self

        indices = set(self._indices)
        if isinstance(other, Atom):
            other_indices = set([other._index])
        else:
            other_indices = set(other._indices)
            
        indices = indices.intersection(other_indices)
        if indices:
            indices = np.unique(indices)
            return Selection(self._ag, indices, '({0:s}) and ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()))
               
    def getCoords(self):
        """Return a copy of coordinates."""
        
        if self._ag._coords is not None:
            # Since this is not slicing, a view is not returned
            return self._ag._coords[self._indices]
    
    _getCoords = getCoords
    
    def setCoords(self, coords):
        """Set coordinates."""
        
        self._ag._coords[self._indices] = coords
        self._ag._setTimeStamp()
    
    def getIndices(self):
        """Return a copy of the indices of atoms."""
        
        return self._indices.copy()
    
    def _getIndices(self):
        """Return indices of atoms."""
        
        return self._indices
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return len(self._indices)

    def getData(self, label):
        """Return a copy of the data associated with *label*, if it exists.
        
        .. versionadded:: 0.7.1"""
        
        if self._ag.isData(label):
            return self._ag._data[label][self._indices]
    
    _getData = getData
    
    def setData(self, label, data):
        """Update *data* with label *label* for the atom subset.
        
        .. versionadded:: 0.7.1
        
        :raise AttributeError: when data associated with *label* is not present
        """
        
        if self._ag.isData(label):
            if label in READONLY:
                raise AttributeError("{0:s} is read-only".format(label))
            self._ag._data[label][self._indices] = data 
        else:
            raise AttributeError("AtomGroup '{0:s}' has no data with label "
                            "'{1:s}'".format(self._ag.getTitle(), label))


class MCAtomSubset(MCAtomPointer, AtomSubset):

    def __init__(self, ag, indices, acsi, **kwargs):
        
        MCAtomPointer.__init__(self, ag, acsi)

        if not isinstance(indices, np.ndarray):
            indices = np.array(indices, int)
        elif not indices.dtype == int:
            indices = indices.astype(int)
        
        if kwargs.get('unique'):
            self._indices = indices
        else:
            self._indices = np.unique(indices)
        
        self._selstr = kwargs.get('selstr')

    def __invert__(self):
        
        ones = np.ones(self._ag.numAtoms(), bool)
        ones[self._indices] = False
        sel = Selection(self._ag, ones.nonzero()[0], 
                        "not ({0:s}) ".format(self.getSelstr()), self._acsi)        
        return sel
    
    def __or__(self, other):
        
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
            
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
            
        if self is other:
            return self
    
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
            
        if isinstance(other, Atom):
            other_indices = np.array([other._index])
        else:
            other_indices = other._indices
            
        indices = np.unique(np.concatenate((self._indices, other_indices)))
        return Selection(self._ag, indices, '({0:s}) or ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()), acsi)

    def __and__(self, other):
        
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
            
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
        
        if self is other:
            return self
    
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
    
        indices = set(self._indices)
        if isinstance(other, Atom):
            other_indices = set([other._index])
        else:
            other_indices = set(other._indices)
    
        indices = indices.intersection(other_indices)
        if indices:
            indices = np.unique(indices)
            return Selection(self._ag, indices, '({0:s}) and ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()), acsi)
               

    def getCoords(self):
        """Return a copy of coordinates from the active coordinate set."""
        
        if self._ag._coords is not None:
            # Since this is not slicing, a view is not returned
            return self._ag._coords[self._acsi, self._indices]
    
    _getCoords = getCoords
    
    def setCoords(self, coords):
        """Set coordinates in the active coordinate set."""
        
        self._ag._coords[self._acsi, self._indices] = coords
        self._ag._setTimeStamp(self._acsi)

    
    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer 
        or a list/array of integers."""
        
        if self._ag._coords is None:
            return None
        if indices is None:
            return self._ag._coords[:, self._indices]
        if isinstance(indices, (int, slice)):
            return self._ag._coords[indices, self._indices]
        if isinstance(indices, (list, np.ndarray)):
            return self._ag._coords[indices, self._indices]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')
                         
    _getCoordsets = getCoordsets

    def iterCoordsets(self):
        """Yield copies of coordinate sets."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coords[i, self._indices]

    _iterCoordsets = iterCoordsets


