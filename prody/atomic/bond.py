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
"""This module defined classes for dealing with bond information provided
by the user.  Bonds can be set using :meth:`~atomgroup.AtomGroup.setBonds` 
method."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

__all__ = ['Bond']

class Bond(object):
    
    """A pointer class for bonded atoms."""
    
    __slots__ = ['_ag', '_acsi', '_indices']
    
    def __init__(self, ag, indices, acsi=None):
        
        self._ag = ag
        self._indices = indices
        if acsi is None:
            self._acsi = ag.getACSIndex()
        else:
            self._acsi = acsi

    def __repr__(self):
        
        one, two = self._indices
        names = self._ag._getNames()
        return '<Bond: {0:s}({1:d})--{2:s}({3:d}) from {4:s}>'.format(
                            names[one], one, names[two], two, str(self._ag))
    
    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0:s}({1:d})--{2:s}({3:d})'.format(
                                            names[one], one, names[two], two)

    def getAtomGroup(self):
        """Return atom group."""
        
        return self._ag

    def getAtoms(self):
        """Return bonded atoms."""
        
        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Return indices of bonded atoms."""
        
        return self._indices.copy()
    
    def getLength(self):
        
        vector = self.getVector()
        return np.multiply(vector, vector, vector).sum() ** 0.5

    def getVector(self):
        """Return bond vector that originates from the first atom."""
        
        one, two = self._indices
        acsi = self.getACSIndex()
        return self._ag._coords[acsi, two] - self._ag._coords[acsi, one]
    
    def getACSIndex(self):
        """Return index of the coordinate set."""
        
        acsi = self._acsi
        if acsi >= self._ag._n_csets:
            raise ValueError('{0:s} has fewer coordsets than assumed by {1:s}'
                             .format(str(self._ag), str(self)))
        return acsi

    def setACSIndex(self, index):
        """Set the coordinate set at *index* active."""

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

def evalBonds(bonds, n_atoms):
    """Return an array mapping atoms to their bonded neighbors and an array
    that stores number of bonds made by each atom."""
    
    numbonds = np.bincount(bonds.reshape((bonds.shape[0] * 2)))
    bmap = np.zeros((n_atoms, numbonds.max()), int)
    bmap.fill(-1)
    index = np.zeros(n_atoms, int)
    for bond in bonds:
        a, b = bond
        bmap[a, index[a]] = b
        bmap[b, index[b]] = a
        index[bond] += 1
    return bmap, numbonds

def trimBonds(bonds, indices):
    """Return bonds made by atoms at given indices."""
    
    newindices = np.zeros(indices.max()+1, int)
    newindices[indices] = np.arange(len(indices))
    iset = set(indices)
    newbonds = []
    for i, bond in enumerate(bonds):
        if bond[0] in iset and bond[1] in iset:
            newbonds.append(newindices[bond])

    newbonds = np.array(newbonds)
    if len(newbonds) > 0:
        return newbonds
