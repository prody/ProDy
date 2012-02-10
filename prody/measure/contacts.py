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

""" This module defines a class for identifying contacts.

.. currentmodule:: prody
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

import prody
from measure import getKDTree

__all__ = ['Contacts']


class Contacts(object):
    
    """A class for identification of contacts in or between atom groups.  
    Contacts are identified using the coordinates of atoms at the time
    of instantiation."""
    
    def __init__(self, atoms):
        """*atoms* for which contacts will be identified. *atoms* can be
        instances of one of :class:`~atomic.atomgroup.AtomGroup` or 
        :class:`~atomic.subset.AtomSubset`."""

        if not isinstance(atoms, (prody.AtomGroup, prody.AtomSubset)):                
            raise TypeError('{0:s} is not a valid type for atoms'
                            .format(type(atoms)))
        self._atoms = atoms
        if isinstance(atoms, prody.AtomGroup):
            self._ag = atoms 
            self._indices = None
            self._kdtree = atoms._getKDTree()
        else:
            self._ag = atoms.getAtomGroup()
            self._indices = atoms.getIndices()
            self._kdtree = getKDTree(self._atoms._getCoords())
        self._acsi = atoms.getACSIndex()

    def __repr__(self):
        
        return '<Contacts: {0:s} (active coordset index: {1:d})>'.format(
                                                str(self._atoms), self._acsi)

    def __str__(self):
        
        return 'Contacts ' + str(self._atoms)


    def select(self, within, what):
        """Select atoms *within* of *what*.  *within* is distance in Å and 
        *what* can be point(s) in 3-d space (:class:`~numpy.ndarray` with 
        shape N,3) or a set of atoms, i.e. :class:`~atomic.bases.Atomic` 
        instances."""
        
        if isinstance(what, np.ndarray):
            if what.ndim == 1 and len(what) == 3:
                what = [what]
            elif not (what.ndim == 2 and what.shape[1] == 3):
                raise ValueError('*what* must be a coordinate array, '
                                 'shape (N, 3) or (3,).')
        else:
            try:
                what = what._getCoords()
            except:
                raise TypeError('*what* must have a getCoords() method.')
            if not isinstance(what, np.ndarray):
                raise ValueError('what.getCoords() method must '
                                 'return a numpy.ndarray instance.')

        search = self._kdtree.search
        get_indices = self._kdtree.get_indices
        indices = []
        append = indices.append
        for xyz in what:
            search(xyz, float(within))
            append(get_indices())
        indices = np.unique(np.concatenate(indices))
        if len(indices) != 0:
            if self._indices is not None:        
                indices = self._indices[indices]
            return prody.Selection(self._ag, np.array(indices), 
                    'index {0:s}'.format(' '.join(np.array(indices, '|S'))), 
                                         acsi=self._acsi, unique=True)
