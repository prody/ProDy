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

""" This module defines a class and function for identifying contacts."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import array, ndarray

from prody.atomic import Atomic, AtomGroup, AtomSubset, Selection
from prody.kdtree import KDTree
from prody.utilities import rangeString

__all__ = ['Contacts', 'iterNeighbors', 'findNeighbors']

class Contacts(object):
    
    """A class for contact identification.    
    Contacts are identified using the coordinates of atoms at the time
    of instantiation."""
    
    def __init__(self, atoms, unitcell=None):
        """*atoms* must be an :class:`.Atomic` instance.  When an orthorhombic
        *unitcell* array is given"""

        try:
            self._acsi = atoms.getACSIndex()
        except AttributeError:
            try:
                self._ag = atoms.getAtoms()
                unitcell = unitcell or atoms.getUnitcell()[:3]
                self._indices = atoms.getSelection()                
            except AttributeError:
                try:
                    ndim, shape = atoms.ndim, atoms.shape
                except AttributeError:
                    raise TypeError('atoms must be an Atomic or Frame instance'
                                    ', not a {0:s}'.format(type(atoms)))
                else:
                    if not (ndim == 2 and shape[1] == 3):
                        raise ValueError('atoms.shape must be (n_atoms, 3) or '
                                         '(3,).')
                    self._ag = None
                    self._indices = None
                    self._kdtree = KDTree(atoms, unitcell=unitcell)
            else:
                if self._ag is not None:
                    self._acsi = self._ag.getACSIndex()
                    if self._indices is not None:
                        self._indices = self._indices.getIndices()
                else:
                    self._acsi = None
                self._kdtree = KDTree(self._atoms._getCoords(), 
                                      unitcell=unitcell)
        else:
            try:        
                self._ag = atoms.getAtomGroup()
            except AttributeError:
                self._ag = atoms 
                self._indices = None
                self._kdtree = KDTree(atoms._getCoords(), unitcell=unitcell)
            else:
                self._indices = atoms._getIndices()
                self._kdtree = KDTree(atoms._getCoords(), unitcell=unitcell)
        self._unitcell = unitcell
        self._atoms = atoms

    def __repr__(self):
        
        return '<Contacts: {0:s} (active coordset index: {1:d})>'.format(
                                                str(self._atoms), self._acsi)

    def __str__(self):
        
        return 'Contacts ' + str(self._atoms)

    def __call__(self, radius, center):
        """Select atoms radius *radius* (Å) of *center*, which can be point(s)
        in 3-d space (:class:`numpy.ndarray` with shape ``(n_atoms, 3)``) or a
        set of atoms, e.g. :class:`.Selection`."""
        
        try:
            center = center._getCoords()
        except AttributeError:
            try:
                ndim, shape = center.ndim, center.shape
            except AttributeError:
                raise TypeError('center must be an Atomic instance or a'
                                'coordinate array')
            else:        
                if shape == (3,):
                    center = [center]
                elif not ndim == 2 and shape[1] == 3:
                    raise ValueError('center.shape must be (n_atoms, 3) or'
                                     '(3,)')
        else:
            if center is None:
                raise ValueError('center does not have coordinate data')

        search = self._kdtree.search
        get_indices = self._kdtree.getIndices
        get_count = self._kdtree.getCount
        indices = set()
        update = indices.update
        radius = float(radius)
        for xyz in center:
            search(radius, xyz)
            if get_count():
                update(get_indices())
        indices = list(indices)
        if indices:
            indices.sort()
            if self._ag is None:
                return array(indices)
            else:
                if self._indices is not None:        
                    indices = self._indices[indices]
                return Selection(self._ag, array(indices), 'index ' + 
                                 rangeString(indices), acsi=self._acsi, 
                                 unique=True)

    select = __call__

    def getAtoms(self):
        """Return atoms, or coordinate array, provided at instantiation.."""
        
        return self._atoms
    
    def getUnitcell(self):
        """Return unitcell array, or **None** if one was not provided."""
        
        return self._unitcell.copy()


def iterNeighbors(atoms, radius, atoms2=None, unitcell=None):
    """Yield pairs of *atoms* that are within *radius* of each other and the 
    distance between them.  If *atoms2* is also provided, one atom from *atoms*
    and another from *atoms2* will be yielded."""
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')
    elif not isinstance(radius, (float, int)):
        raise TypeError('radius must be a number')
    elif radius <= 0:
        raise ValueError('radius must be a positive number')
        
    if atoms2 is None:
        if len(atoms) <= 1:
            raise ValueError('length of atoms must be more than 1')
        ag = atoms
        if not isinstance(ag, AtomGroup):
            ag = ag.getAtomGroup()
            indices = atoms._getIndices()
            index = lambda i: indices[i]
        else:
            index = lambda i: i
        kdtree = KDTree(atoms._getCoords(), none=list)
        
        _dict = {}
        for (i, j), r in zip(*kdtree(radius)): 
            a1 = _dict.get(i)
            if a1 is None:      
                a1 = ag[index(i)]
                _dict[i] = a1
            a2 = _dict.get(j)
            if a2 is None:      
                a2 = ag[index(j)]
                _dict[j] = a2
            yield (a1, a2, r)   
    else:
        if len(atoms) >= len(atoms2): 
            ag = atoms
            if not isinstance(ag, AtomGroup):
                ag = ag.getAtomGroup()
                indices = atoms._getIndices()
                index = lambda i: indices[i]
            else:
                index = lambda i: i
            kdtree = KDTree(atoms._getCoords(), none=list)
            
            _dict = {}
            for a2 in atoms2.iterAtoms():
                for i, r in zip(*kdtree(radius, a2._getCoords())): 
                    a1 = _dict.get(i)
                    if a1 is None:      
                        a1 = ag[index(i)]
                        _dict[i] = a1
                    yield (a1, a2, r)   
        else:    
            ag = atoms2
            if not isinstance(ag, AtomGroup):
                ag = ag.getAtomGroup()
                indices = atoms2._getIndices()
                index = lambda i: indices[i]
            else:
                index = lambda i: i
            kdtree = KDTree(atoms2._getCoords(), none=list)
            
            _dict = {}
            for a1 in atoms.iterAtoms():
                for i, r in zip(*kdtree(radius, a1._getCoords())): 
                    a2 = _dict.get(i)
                    if a2 is None:      
                        a2 = ag[index(i)]
                        _dict[i] = a2
                    yield (a1, a2, r)   


def findNeighbors(atoms, radius, atoms2=None):
    """Return list of pairs of *atoms* that are within *radius* of each other 
    and the distance between them.  If *atoms2* is also provided, one atom 
    from *atoms* and another from *atoms2* will be returned."""

    return list(iterNeighbors(atoms, radius, atoms2))
