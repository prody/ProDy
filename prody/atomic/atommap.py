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

"""This module defines atom map classes that allow for pointing to atoms
in arbitrary order.

.. currentmodule:: prody.atomic"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from atomic import ATOMIC_DATA_FIELDS
from atomic import wrapGetMethod, wrapSetMethod
from atom import Atom, MCAtom
from pointer import AtomPointer, MCAtomPointer

__all__ = ['AtomMap', 'MCAtomMap']

class AtomMapMeta(type):
    
    def __init__(cls, name, bases, dict):
        
        for field in ATOMIC_DATA_FIELDS.values():
            meth = field.meth_pl
            getMeth = 'get' + meth
    
            if field.call:
                def getData(self, var=field.var, call=field.call, 
                            dtype=field.dtype):
                    for meth in call:
                        getattr(self._ag, meth)()
                    data = self._ag._data[var][self._indices]
                    result = np.zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result 
    
            else:
                def getData(self, var=field.var, dtype=field.dtype):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    data = self._ag._data[var][self._indices]
                    result = np.zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result
    
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            if field.dtype in (int, float):
                zero = '0'
            elif field.dtype == bool:
                zero = 'True'
            else:
                zero = '""'
            getData.__doc__ = field.getDocstr('get', selex=False) + \
                   'Entries for unmapped atoms will be ``{0:s}``.'.format(zero) 
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)


class AtomMap(AtomPointer):
    
    """A class for mapping atomic data."""
    
    __metaclass__ = AtomMapMeta
    
    def __init__(self, ag, indices, mapping, unmapped, title='Unnamed'):
        """Instantiate an atom map.        
        
        :arg ag: AtomGroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms
        :arg unmapped: list of indices for unmapped atoms
        :arg title: title of the instance, default is 'Unnamed'
        
        Length of *mapping* must be equal to length of *indices*.  Number of 
        atoms (including unmapped dummy atoms) are determined from the total 
        lengths of *mapping* and *unmapped* arrays."""
        
        AtomPointer.__init__(self, ag)
        
        if not isinstance(indices, np.ndarray):
            self._indices = np.array(indices, int)
        elif not indices.dtype == int:
            self._indices = indices.astype(int)
        else:
            self._indices = indices

        if not isinstance(mapping, np.ndarray):
            self._mapping = np.array(mapping, int)
        elif not mapping.dtype == int:
            self._mapping = mapping.astype(int)
        else:
            self._mapping = mapping

        if not isinstance(unmapped, np.ndarray):
            self._unmapped = np.array(unmapped, int)
        elif not unmapped.dtype == int:
            self._unmapped = unmapped.astype(int)
        else:
            self._unmapped = unmapped
        
        self._title = str(title)
        self._len = len(self._unmapped) + len(self._mapping)
        
    def __repr__(self):

        return ('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; {3:d} mapped; '
                '{4:d} unmapped)>').format( self._title, self._ag.getTitle(), 
                self._len, self.numMapped(), self.numUnmapped())

            
    def __str__(self):
    
        return 'AtomMap {0:s}'.format(self._title)
    
    def __len__(self):
    
        return self._len
    
    def getTitle(self):
        """Return title of the instance."""
        
        return self._title
    
    def setTitle(self, title):
        """Set title of the instance."""
        
        self._title = str(title)

    def numAtoms(self):
        """Return number of atoms."""
        
        return self._len

    def iterAtoms(self):
        """Yield atoms. Note that ``None`` will be yielded for dummy atoms."""
    
        indices = np.zeros(self._len, int)
        indices[self._unmapped] = -1
        indices[self._mapping] = self._indices
        ag = self._ag
        for index in indices:
            if index > -1:
                yield Atom(ag, index)
            else:
                yield None
    
    __iter__ = iterAtoms 
    
    def getCoords(self):
        """Return a copy of coordinates."""
        
        if self._ag._coords is None:
            return None
        coords = np.zeros((self._len, 3), float)
        coords[self._mapping] = self._ag._coords[self._indices] 
        return coords
    
    _getCoords = getCoords
    
    def setCoords(self, coords):
        """Set coordinates.  Length of the *coords* array must match the 
        number of mapped atoms."""
        
        self._ag._coords[self._indices] = coords
    
    def getData(self, label):
        """Return a copy of data associated with *label*, if it exists."""
        
        if self._ag.isData(label):
            data = self._ag._data[label][self._indices]
            result = np.zeros((self._len,) + data.shape[1:], data.dtype)
            result[self._mapping] = data
            return result

    _getData = getData

    def getIndices(self):
        """Return a copy of indices of mapped atoms."""
        
        return self._indices.copy()

    def getMapping(self):
        """Return a copy of mapping of indices."""
        
        return self._mapping.copy()
    
    def _getMapping(self):
        """Return mapping of indices."""
        
        return self._mapping    

    def getUnmappedFlags(self):
        """Return an array with 1s for unmapped atoms."""
        
        flags = np.zeros(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 1
        return flags
    
    def getMappedFlags(self):
        """Return an array with 1s for mapped atoms."""
        
        flags = np.ones(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 0
        return flags

    def numMapped(self):
        """Return number of mapped atoms."""
        
        return len(self._mapping)

    def numUnmapped(self):
        """Return number of unmapped atoms."""
        
        return len(self._unmapped)


class MCAtomMap(MCAtomPointer, AtomMap):
    
    """A class for mapping atomic data with multiple coordinate sets."""
    
    def __init__(self, ag, indices, mapping, unmapped, title='Unnamed',
                 acsi=None):
        """Instantiate an atom map.        
        
        :arg ag: MCAtomGroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms
        :arg unmapped: list of indices for unmapped atoms
        :arg title: title of the instance, default is 'Unnamed'
        :arg acsi: active coordinate set index, defaults is that of *ag*
        
        Length of *mapping* must be equal to length of *indices*.  Number of 
        atoms (including unmapped dummy atoms) are determined from the total 
        lengths of *mapping* and *unmapped* arrays."""
        
        MCAtomPointer.__init__(self, ag, acsi)
        
        if not isinstance(indices, np.ndarray):
            self._indices = np.array(indices, int)
        elif not indices.dtype == int:
            self._indices = indices.astype(int)
        else:
            self._indices = indices

        if not isinstance(mapping, np.ndarray):
            self._mapping = np.array(mapping, int)
        elif not mapping.dtype == int:
            self._mapping = mapping.astype(int)
        else:
            self._mapping = mapping

        if not isinstance(unmapped, np.ndarray):
            self._unmapped = np.array(unmapped, int)
        elif not unmapped.dtype == int:
            self._unmapped = unmapped.astype(int)
        else:
            self._unmapped = unmapped
        
        self._title = str(title)
        self._len = len(self._unmapped) + len(self._mapping)
    
    def __repr__(self):
        
        n_csets = self._ag.numCoordsets()
        if n_csets:
            return ('<MCAtomMap: {0:s} (from {1:s}; {2:d} atoms; {3:d} mapped;'
                    ' {4:d} unmapped; {5:d} coordsets, active {6:d})>').format(
                    self._title, self._ag.getTitle(), self._len, 
                    self.numMapped(), self.numUnmapped(), n_csets, self._acsi)
        else:
            return ('<MCAtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                    '{3:d} mapped; {4:d} unmapped; 0 coordsets)>').format(
                    self._title, self._ag.getTitle(), self._len, 
                    self.numMapped(), self.numUnmapped())
            
    def __str__(self):
        return 'MCAtomMap {0:s}'.format(self._title)

    def iterAtoms(self):
        """Yield atoms. Note that ``None`` will be yielded for dummy atoms."""
    
        indices = np.zeros(self._len, int)
        indices[self._unmapped] = -1
        indices[self._mapping] = self._indices
        ag = self._ag
        acsi = self._acsi
        for index in indices:
            if index > -1:
                yield MCAtom(ag, index, acsi)
            else:
                yield None
    
    __iter__ = iterAtoms 
    
    def getCoords(self):
        """Return a copy of coordinates from the active coordinate set."""
        
        if self._ag._coords is None:
            return None
        coords = np.zeros((self._len, 3), float)
        coords[self._mapping] = self._ag._coords[self._acsi, self._indices] 
        return coords
    
    _getCoords = getCoords
    
    def setCoords(self, coords):
        """Set coordinates in the active coordinate set.  Length of the 
        *coords* array must match the number of mapped atoms."""
        
        self._ag._coords[self._acsi, self._indices] = coords
    

    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer 
        or a list/array of integers."""
        
        if self._ag._coords is None:
            return None
        
        if indices is None:
            indices = np.arange(self._ag._n_csets)
        elif isinstance(indices, (int, long)):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(indices.indices(self._ag._n_csets))
        
        try:
            coordsets = np.zeros((len(indices), self._len, 3))
            coordsets[:, self._mapping] = self._ag._coords[indices][:, 
                                                                self._indices]  
            return coordsets
        
        except IndexError:
            raise IndexError('indices may be an integer or a list/array '
                             'of integers')

    _getCoordsets = getCoordsets
    
    def iterCoordsets(self):
        """Yield copies of coordinate sets."""
        
        for i in range(self._ag._n_csets):
            coords = np.zeros((self._len, 3), float)
            coords[self._mapping] = self._ag._coords[i, self._indices] 
            yield coords
    
    _iterCoordsets = iterCoordsets
