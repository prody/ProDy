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

"""This module defines atom pointer base class.

.. currentmodule:: prody.atomic"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from atomic import Atomic, MultiCoordset

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class AtomPointer(Atomic):
    
    """A base for classes pointing to atoms in :class:`~atomgroup.AtomGroup` 
    instances."""
    
    
    def __init__(self, ag):
        
        if not isinstance(ag, pkg.AtomGroup):
            raise TypeError('ag must be an AtomGroup instance, not {0:s}'
                            .format(type(ag)))
        self._ag = ag

    def __contains__(self, item):
        
        return (isinstance(item, AtomPointer) and
                self._ag == item.getAtomGroup() and len(item) <= len(self) and 
                set(item._getIndices()).issubset(set(self._getIndices())))
    
    def __eq__(self, other):
        
        return (isinstance(other, AtomPointer) and 
                self._ag == other.getAtomGroup() and 
                len(other) == len(self) and 
                np.all(self._indices == other._getIndices()))
    
    def __ne__(self, other):
        
        return not self.__eq__(other)

    def __add__(self, other):
        """Returns an :class:`AtomMap` instance. Order of pointed atoms are
        preserved."""
        
        if not isinstance(other, AtomPointer):
            raise TypeError('an AtomPointer instance cannot be added to a '
                            '{0:s} instance'.format(type(other)))
        ag = self._ag
        if ag != other.getAtomGroup():
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        
        title = '({0:s}) + ({1:s})'.format(str(self), str(other))
        indices = np.concatenate([self._getIndices(), other._getIndices()])
        length = len(self)
        
        if isinstance(self, pkg.AtomMap):
            mapping = [self.getMapping()]
            unmapped = [self._unmapped]
        else:
            mapping = [np.arange(length)]
            unmapped = [np.array([])]
        
        if isinstance(other, pkg.AtomMap):
            mapping.append(other.getMapping() + length)
            unmapped.append(other._unmapped + length) 
        else:
            mapping.append(np.arange(length, length+len(other)))
            unmapped.append(np.array([]))
            
        return pkg.AtomMap(ag, indices, np.concatenate(mapping), 
                           np.concatenate(unmapped), title)
                       
    def _getTimeStamp(self):
        
        return self._ag._getTimeStamp()
    
    def _getKDTree(self):
        """Return KDTree for from the atom group."""
        
        return self._ag._getKDTree()

    def getAtomGroup(self):
        """Return associated atom group."""
        
        return self._ag
    
    def copy(self):
        """Make a copy of atoms."""
        
        return self._ag.copy(self)
    
    __copy__ = copy

    def isData(self, label):
        """Return ``True`` if *label* is a user data."""
        
        return self._ag.isData(label)

    def getDataType(self, label):
        """Return type of the user data, ``None`` if data label is not present.
        """
        
        return self._ag.getDataType(label)


class MCAtomPointer(MultiCoordset, AtomPointer):     
                
    """A base for classes pointing to atoms in :class:`~atomgroup.MCAtomGroup` 
    instances."""


    def __init__(self, ag, acsi=None):
        
        if not isinstance(ag, pkg.MCAtomGroup):
            raise TypeError('ag must be an MCAtomGroup instance, not {0:s}'
                            .format(type(ag)))
        self._ag = ag
        if acsi is None:
            self._acsi = ag.getACSIndex()
        else: 
            self._acsi = int(acsi)
    
    def __add__(self, other):
        """Returns an :class:`AtomMap` instance. Order of pointed atoms are
        preserved."""
        
        if not isinstance(other, MCAtomPointer):
            raise TypeError('an MCAtomPointer instance cannot be added to a '
                            '{0:s} instance'.format(type(other)))
        ag = self._ag
        if ag != other._ag:
            raise ValueError('MCAtomPointer instances must point to the same '
                             'MCAtomGroup instance')
        acsi = self._acsi
        if self._acsi != other._acsi:
            LOGGER.warning('Active coordinate set indices of operands are not '
                           'the same.  Result will have {0:d}'.format(acsi))
        
        title = '({0:s}) + ({1:s})'.format(str(self), str(other))
        indices = np.concatenate([self._getIndices(), other._getIndices()])
        length = len(self)
        
        if isinstance(self, pkg.MCAtomMap):
            mapping = [self._getMapping()]
            unmapped = [self._unmapped]
        else:
            mapping = [np.arange(length)]
            unmapped = [np.array([])]
        
        if isinstance(other, pkg.MCAtomMap):
            mapping.append(other._getMapping() + length)
            unmapped.append(other._unmapped + length) 
        else:
            mapping.append(np.arange(length, length + len(other)))
            unmapped.append(np.array([]))
            
        return pkg.MCAtomMap(ag, indices, np.concatenate(mapping), 
                         np.concatenate(unmapped), title, acsi)
     
    def _getTimeStamp(self, index=None):
        
        if index is None:
            index = self._acsi
        return self._ag._getTimeStamp(index)
    
    def _getKDTree(self):
        """Return KDTree for the active coordinate set from the atom group."""
        
        return self._ag._getKDTree(self._acsi)
     
    def numCoordsets(self):
        """Return number of coordinate sets."""
        
        return self._ag._n_csets

    def setACSIndex(self, index):
        """Set coordinates at *index* active."""
        
        if self._ag._coords is None:
            raise AttributeError('coordinates are not set')
            
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
            
        if self._ag._n_csets <= index or  self._ag._n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
            
        if index < 0:
            index += self._ag._n_csets
            
        self._acsi = index
        
    def getACSLabel(self):
        """Return active coordinate set label."""
        
        if self._ag._n_csets:
            return self._ag._cslabels[self._acsi]
