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

"""This module defines classes for handling arbitrary atom selections."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from atom import Atom
from subset import AtomSubset

__all__ = ['Selection']

SELECT = None

ellipsis = lambda s: s[:15] + '...' + s[-15:] if len(s) > 33 else s

class Selection(AtomSubset):
    
    """A class for accessing and manipulating attributes of selection of atoms 
    in an :class:`AtomGroup` instance."""
    
    __slots__ = ['_ag', '_indices', '_acsi', '_selstr']
    
    def __init__(self, ag, indices, selstr, acsi=None, **kwargs):
        
        kwargs['selstr'] = selstr
        AtomSubset.__init__(self, ag, indices, acsi, **kwargs)

    def __repr__(self):
        
        n_csets = self._ag.numCoordsets()
        selstr = ellipsis(self._selstr)
        if n_csets:
            if n_csets == 1:
                return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms)>'
                    ).format(selstr, self._ag.getTitle(), len(self), n_csets)
            else:
                return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                        'active of {3:d} coordsets)>').format(selstr, 
                        self._ag.getTitle(), len(self), self.getACSIndex(), 
                        n_csets)
        else:
            return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; no '
                    'coordinates)>').format(selstr, self._ag.getTitle(), 
                    len(self))

    def __str__(self):

        return 'Selection "{0:s}" from {1:s}'.format(
                ellipsis(self._selstr), self._ag.getTitle())
    
    def getSelstr(self):
        """Return selection string that selects this atom subset."""
        
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom selection."""
        
        return HierView(self)

    def update(self):    
        """Update selection.
        
        .. versionadded:: 0.9.3"""
        
        self._indices = SELECT.getIndices(self._ag, self._selstr)
