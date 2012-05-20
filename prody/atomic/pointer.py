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

"""This module defines atom pointer base class."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from atomic import Atomic
from bond import Bond

from prody import LOGGER

__all__ = ['AtomPointer']


class AtomPointer(Atomic):
    
    """A base for classes pointing to atoms in :class:`~.AtomGroup` instances.
    Derived classes are:
        
      * :class:`~.Atom`
      * :class:`~.AtomSubset`
      * :class:`~.AtomMap`"""
    
    __slots__ = ['_ag', '_acsi']
    
    def __init__(self, ag, acsi):
        
        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance, not {0:s}'
                            .format(type(ag)))
        self._ag = ag
        if acsi is None:
            self._acsi = ag.getACSIndex()
        else: 
            self._acsi = int(acsi)

    def __contains__(self, item):
        
        return (isinstance(item, AtomPointer) and
                self._ag == item.getAtomGroup() and len(item) <= len(self) and 
                set(item._getIndices()).issubset(set(self._getIndices())))
    
    def __eq__(self, other):
        
        return (isinstance(other, AtomPointer) and 
                self._ag == other.getAtomGroup() and 
                len(other) == len(self) and 
                np.all(self._getIndices() == other._getIndices()))
    
    def __ne__(self, other):
        
        return not self.__eq__(other)
    
    def __invert__(self):
        
        ones = np.ones(self._ag.numAtoms(), bool)
        ones[self._getIndices()] = False
        return Selection(self._ag, ones.nonzero()[0], 
                         "not ({0:s}) ".format(self.getSelstr()), 
                         self.getACSIndex(), unique=True)

    def __or__(self, other):
        
        if self is other:
            return self
    
        if not isinstance(other, AtomPointer):
            raise TypeError('other must be an AtomPointer')
            
        if self._ag != other.getAtomGroup():
            raise ValueError('both selections must be from the same AtomGroup')
            
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warn('Active coordinate set indices do not match, it will '
                        'be set to zero.')
            acsi = 0
            
        indices = np.unique(np.concatenate((self._getIndices(), 
                                            other._getIndices())))
        return Selection(self._ag, indices, '({0:s}) or ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()), 
                                    acsi, unique=True)

    def __and__(self, other):
        
        if self is other:
            return self
    
        if not isinstance(other, AtomPointer):
            raise TypeError('other must be an AtomPointer')
            
        if self._ag != other.getAtomGroup():
            raise ValueError('both selections must be from the same AtomGroup')
    
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0

        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warn('Active coordinate set indices do not match, it will '
                        'be set to zero.')
            acsi = 0
            
        indices = set(self._getIndices())
    
        indices = indices.intersection(other.getIndices())
        if indices:
            indices = np.unique(indices)
            return Selection(self._ag, indices, '({0:s}) and ({1:s})'.format(
                                    self.getSelstr(), other.getSelstr()), acsi)

    def __add__(self, other):
        """Returns an :class:`~.AtomMap` instance. Order of pointed atoms are
        preserved."""
        
        if not isinstance(other, AtomPointer):
            raise TypeError('unsupported operand type(s) for +: {0:s} and '
                            '{1:s}'.format(repr(type(self).__name__), 
                                           repr(type(other).__name__)))
                    
        ag = self._ag
        if ag != other._ag:
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('Active coordset indices of atoms are not the same.'
                           ' Result will have ACSI {0:d}.'.format(acsi))
        
        title = '({0:s}) + ({1:s})'.format(str(self), str(other))
        indices = np.concatenate([self._getIndices(), other._getIndices()])
        length = len(self)
        
        if isinstance(self, AtomMap):
            mapping = [self._getMapping()]
            unmapped = [self._dummies]
        else:
            mapping = [np.arange(length)]
            unmapped = [np.array([])]
        
        if isinstance(other, AtomMap):
            mapping.append(other._getMapping() + length)
            unmapped.append(other._dummies + length) 
        else:
            mapping.append(np.arange(length, length + len(other)))
            unmapped.append(np.array([]))
            
        return AtomMap(ag, indices, np.concatenate(mapping), 
                           np.concatenate(unmapped), title, acsi)
                       
    def _getTimeStamp(self, index=None):
        
        if index is None:
            index = self.getACSIndex()
        return self._ag._getTimeStamp(index)
    
    def _getKDTree(self):
        """Return KDTree for the active coordinate set from the atom group."""
        
        return self._ag._getKDTree(self.getACSIndex())

    def getAtomGroup(self):
        """Return associated atom group."""
        
        return self._ag
    
    def numCoordsets(self):
        """Return number of coordinate sets."""
        
        return self._ag._n_csets
    
    def getACSIndex(self):
        """Return index of the coordinate set."""
        
        acsi = self._acsi
        if acsi >= self._ag._n_csets:
            raise ValueError('{0:s} has fewer coordsets than assumed by {1:s}'
                             .format(str(self._ag), str(self)))
        return acsi

    def setACSIndex(self, index):
        """Set coordinates at *index* active."""
        
        if self._ag._coords is None:
            raise AttributeError('coordinates are not set')
            
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
            
        n_csets = self._ag._n_csets
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
            
        if index < 0:
            index += n_csets
            
        self._acsi = index
        
    def getACSLabel(self):
        """Return active coordinate set label."""
        
        if self._ag._n_csets:
            return self._ag._cslabels[self.getACSIndex()]
    
    def isData(self, label):
        """Return ``True`` if data *label* is present."""
        
        return self._ag.isData(label)

    def getDataType(self, label):
        """Return type of data, or ``None`` if data *label* is not present."""
        
        return self._ag.getDataType(label)

    def _iterBonds(self):
        """Yield pairs of indices for bonded atoms that are within the pointer. 
        Use :meth:`setBonds` for setting bonds."""
        
        if self._ag._bonds is None:
            raise ValueError('bonds are not set, use `AtomGroup.setBonds`')
            
        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) * 0.5 >= len(self): 
            for a, b in self._ag._iterBonds():
                if a in iset and b in iset:
                    yield a, b
        else:
            for a, bmap in zip(indices, self._ag._bmap[indices]):
                for b in bmap:
                    if b > -1 and b in iset:
                        yield a, b
                iset.remove(a)
            
