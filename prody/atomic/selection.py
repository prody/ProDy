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

"""This module defines classes for handling arbitrary atom selections.

.. currentmodule:: prody.atomic"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from subset import AtomSubset, MCAtomSubset

__all__ = ['Selection', 'MCSelection']

ellipsis = lambda s: s[:15] + '...' + s[-15:] if len(s) > 33 else s

class Selection(AtomSubset):
    
    """A class for accessing and manipulating attributes of select of atoms 
    in an :class:`AtomGroup` instance.
    
    """
    
    def __init__(self, ag, indices, selstr, **kwargs):
        
        kwargs['selstr'] = selstr
        AtomSubset.__init__(self, ag, indices, **kwargs)
    
    def __repr__(self):

        return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms)>').format(
                ellipsis(self._selstr), self._ag.getTitle(), len(self))

    def __str__(self):

        return 'Selection "{0:s}" from {1:s}'.format(
                ellipsis(self._selstr), self._ag.getTitle())
    
    def iterAtoms(self):
        """Yield atoms."""

        ag = self._ag
        for index in self._indices:
            yield Atom(ag=ag, index=index)

    __iter__ = iterAtoms
    
    def getSelstr(self):
        """Return selection string that selects this atom subset."""
        
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom selection."""
        
        return HierView(self)

    def update(self):    
        """Update selection.
        
        .. versionadded:: 0.9.3"""
        
        self._indices = prody.ProDyAtomSelect.getIndices(self._ag, 
                                                         self._selstr)

class MCSelection(MCAtomSubset, Selection):
    
    def __init__(self, ag, indices, selstr, acsi, **kwargs):
        
        kwargs['selstr'] = selstr
        MCAtomSubset.__init__(self, ag, indices, acsi, **kwargs)


    def __repr__(self):
        
        n_csets = self._ag.numCoordsets()
        selstr = ellipsis(self._selstr)
        if n_csets:
            return ('<MCSelection: "{0:s}" from {1:s} ({2:d} atoms; '
                    '{3:d} coordsets, active {4:d})>').format(selstr, 
                    self._ag.getTitle(), len(self), n_csets, self._acsi)
        else:
            return ('<MCSelection: "{0:s}" from {1:s} ({2:d} atoms; 0 '
                    'coordsets)>').format(selstr, self._ag.getTitle(), 
                    len(self), n_csets)

    def __str__(self):
        
        return 'MCSelection "{0:s}" from {1:s}'.format(
                ellipsis(self._selstr), self._ag.getTitle())
    
    def iterAtoms(self):
        """Yield atoms."""

        ag = self._ag
        acsi = self._acsi
        for index in self._indices:
            yield Atom(ag=ag, index=index, acsi=acsi)

    __iter__ = iterAtoms
