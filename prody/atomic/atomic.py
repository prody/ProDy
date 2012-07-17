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

from numpy import all, zeros

from prody import LOGGER

from bond import trimBonds
from fields import ATOMIC_FIELDS, READONLY

__all__ = ['Atomic']

SELECT = None
isKeyword = None
isSelectionMacro = None

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
            if name.startswith('is') and self.isFlagLabel(name[2:]):
                return all(self._getFlags(name[2:]))
            else:
                if self.isFlagLabel(name):
                    try:
                        ag = self.getAtomGroup()
                    except AttributeError:
                        ag = self
                        selstr = name
                    else:
                        selstr = '({0:s}) and ({1:s})'.format(name, 
                                                              self.getSelstr())
                    try:
                        dummies = self.numDummies()
                    except AttributeError:
                        return Selection(ag, self._getSubset(name), selstr, 
                                         self.getACSIndex(), unique=True)
                    else:
                        return AtomMap(ag, self._getSubset(name), 
                                       self.getACSIndex(), intarrays=True,
                                       dummies=dummies)
                else:
                    selstr = name
                    items = name.split('_')
                    word = items[0]
                    if (isKeyword(word) or word == 'not' or 
                        isSelectionMacro(word) or self.isDataLabel(word)):
                        selstr = ' '.join(items)
                        return SELECT.select(self, selstr)

        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))

    def copy(self):
        """Return a copy of atoms (and atomic data) in a new :class:`AtomGroup`
        instance."""
        
        dummies = None
        indices = None
        try:
            ag = self.getAtomGroup()
        except AttributeError:
            ag = self
            new = AtomGroup(ag.getTitle())
        else:
            indices = self.getIndices()
            new = AtomGroup(ag.getTitle() + ' ' + str(self))
            try:
                dummies = self.numDummies()
            except AttributeError:
                pass
            else:
                if dummies:
                    dummy = self.getDummyFlags()
                    mapped = self.getMappedFlags()

        try:
            self.getIndex()
        except AttributeError:
            this = self
        else:
            this = self.all

        if self.numCoordsets():
            new.setCoords(this.getCoordsets(), label=ag.getCSLabels())
        
        for key, array in ag._data.iteritems():
            if key in READONLY:
                continue
            try:
                field = ATOMIC_FIELDS[key]
            except KeyError:
                new._data[key] = this.getData(key)
            else:
                meth = field.meth_pl
                getattr(new, 'set' + meth)(getattr(this, 'get' + meth)())
                
        for label in ag.getFlagLabels():
            new.setFlags(label, this.getFlags(label))
        if dummies:
            new.setFlags('dummy', dummy)
            new.setFlags('mapped', mapped)
            
        bonds = ag._bonds
        bmap = ag._bmap
        if bonds is not None and bmap is not None:
            if indices is None:
                new._bonds = bonds.copy()
                new._bmap = bmap.copy()
                new._data['numbonds'] = ag._data['numbonds'].copy()
                if ag._data['fragindex'] is not None:
                    new._data['fragindex'] = ag._data['fragindex'].copy()
            elif dummies is not None:
                if dummies:
                    indices = indices[self._getMapping()]
                if len(set(indices)) == len(indices): 
                    new.setBonds(trimBonds(bonds, indices))
                else:
                    LOGGER.warn('Duplicate atoms in mapping, bonds are '
                                'not copied.')
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
