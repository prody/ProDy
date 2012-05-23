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

"""This module defines base class :class:`Atomic` that all other 
:mod:`~prody.atomic` classes are derived from."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import zeros

from prody import LOGGER

from bond import trimBonds
from fields import ATOMIC_ATTRIBUTES, READONLY

__all__ = ['Atomic']

isMacro = lambda none: None
isKeyword = lambda none: None
SELECT = None

class Atomic(object):
    
    """Base class for all atomic classes.  This class can be used for type
    checking:
    
    >>> from prody import *
    >>> ag = parsePDB('1aar')
    >>> isinstance(ag, Atomic)
    True
    >>> prot = ag.select('protein')
    >>> isinstance(prot, Atomic)
    True"""
    
    __slots__ = []
    
    def __getattribute__(self, name):
        
        try:
            return object.__getattribute__(self, name)

        except AttributeError:
            selstr = name
            items = name.split('_')
            word = items[0]
            if (isKeyword(word) or word == 'not' or isMacro(word) or
                self.isData(word)):
                selstr = ' '.join(items)
                return SELECT.select(self, selstr)

        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))

    def copy(self):
        """Return a copy of atoms (and atomic data) in a new :class:`AtomGroup`
        instance."""
        
        atommap = False
        indices = None
        try:
            ag = self.getAtomGroup()
        except AttributeError:
            ag = self
            new = AtomGroup(ag.getTitle())
        else:
            indices = self._getIndices()
            if hasattr(self, 'getDummyFlags'):
                new = AtomGroup(ag.getTitle() + ' mapping ' + self.getTitle())
                new._n_atoms = len(self)
                new.setData('dummy', self.getDummyFlags())
                new.setData('mapped', self.getMappedFlags())
                atommap = True
            else:
                new = AtomGroup(ag.getTitle() + ' selection ' + 
                                repr(self.getSelstr()))

        coords = ag._getCoordsets()
        if coords is not None:
            if indices is None or atommap:
                new.setCoords(self.getCoordsets())
            else:
                # for an Atom, array from getCoordsets will have wrong shape  
                new.setCoords(coords[:, indices])
            new._cslabels = ag.getCSLabels()
        
        for key, array in ag._data.iteritems():
            if key in READONLY or array is None:
                continue
            if key in ATOMIC_ATTRIBUTES:
                if indices is None:
                    new._data[key] = array.copy()
                elif atommap:
                    new._data[key] = getattr(self, 'get' + 
                                             ATOMIC_ATTRIBUTES[key].meth_pl)()
                else:
                    new._data[key] = array[indices]
            else:
                new._data[key] = self.getData(key)
                
        bonds = ag._bonds
        bmap = ag._bmap
        if bonds is not None and bmap is not None:
            if indices is None:
                new._bonds = bonds.copy()
                new._bmap = bmap.copy()
                new._data['numbonds'] = ag._data['numbonds'].copy()
                if ag._data['fragindices'] is not None:
                    new._data['fragindices'] = ag._data['fragindices'].copy()
            elif atommap:
                if len(set(indices)) == len(indices): 
                    bonds = trimBonds(bonds, indices)
                    if bonds is not None:
                        revmap = zeros(len(ag), int)
                        revmap[indices] = self._getMapping()
                        new.setBonds(revmap[bonds])
                else:
                    LOGGER.warn('Bond information is not copied to {0:s}.'
                                .format(new.getTitle()))
            else:
                bonds = trimBonds(bonds, indices)
                if bonds is not None:
                    new.setBonds(bonds)
        return new

    __copy__ = copy

    def select(self, selstr, **kwargs):
        """Return atoms matching *selstr* criteria.  See :mod:`~.select` module
        documentation for details and usage examples."""
        
        return SELECT.select(self, selstr, **kwargs)
