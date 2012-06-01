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

"""This module defines :class:`AtomMap` class that allows for pointing atoms in 
arbitrary order.

.. _atommaps:

How AtomMap's work
===============================================================================
    
:class:`AtomMap` class adds great flexibility to manipulating atomic data.

First let's see how an instance of :class:`.Selection` (:class:`.Chain`, or 
:class:`.Residue`) works.  Below table shows indices for a selection of atoms 
in an :class:`~.AtomGroup` and values returned when 
:meth:`~.Selection.getNames`, :meth:`~.Selection.getResnames` and
:meth:`~.Selection.getResnums` methods are called.

.. csv-table:: **Atom Subset** 
   :header: "Indices", "Names", "Resnames", "Resnums"

   0, N, PHE, 1
   1, CA, PHE, 1
   2, C, PHE, 1
   3, O, PHE, 1
   4, CB, PHE, 1
   5, CG, PHE, 1
   6, CD1, PHE, 1
   7, CD2, PHE, 1
   8, CE1, PHE, 1
   9, CE2, PHE, 1
   10, CZ, PHE, 1

:class:`~.Selection` instances keep indices ordered and do not allow duplicate 
values, hence their use is limited. In an :class:`AtomMap`, indices do not need
to be sorted, duplicate indices may exist, even "DUMMY" atoms are allowed.

Let's say we instantiate the following AtomMap::
    
    amap = AtomMap(atomgroup, indices=[0, 1, 3, 8, 8, 9, 10], 
                   mapping=[5, 6, 7, 0, 1, 2, 3])


The size of the :class:`AtomMap` based on this mapping is 8, since the larger 
mapping is 7.

Calling the same functions for this AtomMap instance would result in the 
following:

.. csv-table:: **Atom Map**
   :header: "Mapping", "Indices", "Names", "Resnames", "Resnums", \
            "MappedFlags", "DummyFlags"

   0, 8, CE1, PHE, 1, 1, 0
   1, 8, CE1, PHE, 1, 1, 0
   2, 9, CE2, PHE, 1, 1, 0
   3, 10, CZ, PHE, 1, 1, 0
   4, , , , 0, 0, 1
   5, 0, N, PHE, 1, 1, 0
   6, 1, CA, PHE, 1, 1, 0
   7, 3, O, PHE, 1, 1, 0
   
For unmapped atoms, numeric attributes are set to 0, others to empty string,
i.e. ``""``.

.. seealso::
   :class:`AtomMap` are used by :mod:`.proteins` module functions that 
   match or map protein chains.  :ref:`pca-xray` and :ref:`pca-dimer` 
   examples that make use of these functions and :class:`AtomMap` class.


"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import arange, array, ndarray, ones, zeros 

from prody.utilities import rangeString

from atom import Atom
from fields import ATOMIC_FIELDS
from fields import wrapGetMethod, wrapSetMethod
from pointer import AtomPointer

__all__ = ['AtomMap']

class AtomMapMeta(type):
    
    def __init__(cls, name, bases, dict):
        
        for field in ATOMIC_FIELDS.values():
            
            if field.private:
                continue
            
            meth = field.meth_pl
            getMeth = 'get' + meth
    
            if field.call:
                def getData(self, var=field.var, call=field.call, 
                            dtype=field.dtype):
                    for meth in call:
                        getattr(self._ag, meth)()
                    data = self._ag._data[var][self._indices]
                    result = zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result 
    
            else:
                def getData(self, var=field.var, dtype=field.dtype):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    data = self._ag._data[var][self._indices]
                    result = zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result
    
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            if field.dtype in (int, float):
                zero = 0
            elif field.dtype == bool:
                zero = True
            else:
                zero = ''
            getData.__doc__ = (field.getDocstr('get', selex=False) +
                               ' Entries for dummy atoms will be ``{0:s}``.'
                               .format(repr(zero))) 
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)


class AtomMap(AtomPointer):
    
    """A class for mapping atomic data."""
    
    __metaclass__ = AtomMapMeta
    
    __slots__ = ['_ag', '_indices', '_acsi', '_mapping', '_dummies', '_title',
                 '_len', '_idarray']
    
    def __init__(self, ag, indices, mapping, dummies, title='Unnamed',
                 acsi=None, **kwargs):
        """Instantiate an atom map.        
        
        :arg ag: AtomGroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms
        :arg dummies: list of indices for dummy atoms
        :arg title: title of the instance, default is 'Unnamed'
        :arg acsi: active coordinate set index, defaults is that of *ag*
        
        Length of *mapping* must be equal to length of *indices*.  Number of 
        atoms (including dummy atoms) are determined from the total lengths 
        of *mapping* and *dummies* arrays.
        
        Following built-in functions are customized for this class:
        
        * :func:`len` returns the number of atoms in the instance.
        * :func:`iter` yields :class:`~.Atom` instances.
        * Indexing is not available."""
        
        AtomPointer.__init__(self, ag, acsi)
        
        if kwargs.get('intarrays'):
            self._indices = indices
            self._mapping = mapping
            if dummies is None:
                self._dummies = array([], int)
            else:
                self._dummies = dummies
        else:            
            self._indices = array(indices, int)
            self._mapping = array(mapping, int)
            if dummies is None:
                self._dummies = array([], int)
            else:
                self._dummies = array(dummies, int)
        
        self._idarray = None
        self._title = str(title)
        self._len = len(self._dummies) + len(self._mapping)
        
    
    def __repr__(self):
        
        rep = '<AtomMap: {0:s} from {1:s} ({2:d} atoms'.format(
                self._title, self._ag.getTitle(), self._len) 
        if self.numDummies():
            rep += ', {0:d} mapped, {1:d} dummy'.format(self.numMapped(), 
                                                         self.numDummies())
        
        n_csets = self._ag.numCoordsets()
        if n_csets > 1:
            rep += '; active {0:d} of {1:d} coordsets)>'.format(
                    self.getACSIndex(), n_csets)
        elif n_csets == 0:
            rep += '; no coordinates'
        return rep + ')>'
        
    def __str__(self):
    
        return 'AtomMap {0:s}'.format(self._title)
    
    def __len__(self):
    
        return self._len
    
    def __getitem__(self, index):

        if self._idarray is None:        
            idarray = zeros(self._len, int)
            idarray[self._mapping] = self._indices
            if len(self._dummies):
                idarray[self._dummies] = -1
            self._idarray = idarray
            indices = idarray[index]
        else:
            indices = self._idarray[index]
        
        try:
            n_sel = len(indices)
        except TypeError: 
            if indices > -1:
                return self._ag[indices]
        else:
            mapping = (indices > -1).nonzero()[0]
            if len(mapping) < n_sel:
                return AtomMap(self._ag, indices, mapping, 
                       (indices == -1).nonzero()[0], 
                       '({0:s})[{1:s}]'.format(self._title, repr(index)), 
                       self._acsi, intarrays=True)
            else:
                return AtomMap(self._ag, indices, mapping, None,
                       '({0:s})[{1:s}]'.format(self._title, repr(index)), 
                       self._acsi, intarrays=True)
        
    
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
    
        acsi = self.getACSIndex()
        indices = zeros(self._len, int)
        indices[self._dummies] = -1
        indices[self._mapping] = self._indices
        ag = self._ag
        for index in indices:
            if index > -1:
                yield Atom(ag, index, acsi)
            else:
                yield None
    
    __iter__ = iterAtoms 
    
    def getCoords(self):
        """Return a copy of coordinates from the active coordinate set."""
        
        coords = self._ag._getCoordsets()
        if coords is not None:
            xyz = zeros((self._len, 3), float)
            xyz[self._mapping] = coords[self.getACSIndex(), self._indices] 
            return xyz
    
    _getCoords = getCoords
    
    def setCoords(self, coords):
        """Set coordinates in the active coordinate set.  Length of the 
        *coords* array must match the number of mapped atoms."""
        
        coordsets = self._ag._getCoordsets()
        if coordsets is not None: 
            coordsets[self.getACSIndex(), self._indices] = coords
    

    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer 
        or a list/array of integers."""
        
        coords = self._ag._getCoordsets()
        if coords is not None:
            n_csets = self._ag.numCoordsets()
            if indices is None:
                indices = arange(n_csets)
            elif isinstance(indices, (int, long)):
                indices = array([indices])
            elif isinstance(indices, slice):
                indices = arange(indices.indices(n_csets))
            
            try:
                coordsets = zeros((len(indices), self._len, 3))
                coordsets[:, self._mapping] = coords[indices][:, self._indices]  
                return coordsets
            
            except IndexError:
                raise IndexError('indices may be an integer or a list/array '
                                 'of integers')

    _getCoordsets = getCoordsets
    
    def iterCoordsets(self):
        """Yield copies of coordinate sets."""
        
        coords = self._ag._getCoordsets()
        if coords is not None:
            mapping = self._mapping
            n_atoms = self._len
            indices = self._indices
            for i in range(self._ag.numCoordsets()):
                xyz = zeros((n_atoms, 3), float)
                xyz[mapping] = coords[i, indices] 
                yield xyz
    
    _iterCoordsets = iterCoordsets

    def getData(self, label):
        """Return a copy of data associated with *label*, if it is present."""
        
        try:
            data = self._ag._data[label]
        except KeyError:
            pass
        else:
            result = zeros((self._len,) + data.shape[1:], data.dtype)
            result[self._mapping] = data[self._indices]
            return result

    _getData = getData

    def getIndices(self):
        """Return a copy of indices of mapped atoms."""
        
        return self._indices.copy()

    def _getIndices(self):
        """Return indices of mapped atoms."""
        
        return self._indices

    def getMapping(self):
        """Return a copy of mapping of indices."""
        
        return self._mapping.copy()
    
    def _getMapping(self):
        """Return mapping of indices."""
        
        return self._mapping    

    def getDummyFlags(self):
        """Return an array with 1s for dummy atoms."""
        
        flags = zeros(self._len, bool)
        if len(self._dummies):
            flags[self._dummies] = 1
        return flags
    
    def getMappedFlags(self):
        """Return an array with 1s for mapped atoms."""
        
        flags = ones(self._len, bool)
        if len(self._dummies):
            flags[self._dummies] = 0
        return flags

    def numMapped(self):
        """Return number of mapped atoms."""
        
        return len(self._mapping)

    def numDummies(self):
        """Return number of dummy atoms."""
        
        return len(self._dummies)

    def getSelstr(self):
        """Return selection string that selects mapped atoms."""
        
        return 'index ' + rangeString(self._indices)
