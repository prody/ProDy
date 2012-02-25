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

"""This module defines a class for handling trajectory frames."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from prody import measure
from prody.measure import superpose

__all__ = ['Frame']

class Frame(object):
    
    """A class for storing trajectory frame coordinates and provide methods 
    acting on them."""
    
    __slots__ = ['_traj', '_index', '_coords', '_indices', '_ag', '_sel', 
                 '_unitcell', '_velocs']

    def __init__(self, traj, index, coords, unitcell=None, velocs=None):

        self._coords = coords
        self._unitcell = unitcell
        self._velocs = velocs
        
        self._traj = traj
        self._index = index
        self._ag = traj.getAtoms()
        self._sel = traj.getSelection()
        self._indices = traj._getSelIndices()
        
    def __repr__(self):

        sel = ''
        if self._sel:
            sel = 'selected {0:d} of '.format(self.numSelected())

        return ('<Frame: {0:d} from {1:s} ({3:d} atoms)>'
               ).format(self._index, self._traj.getTitle(), 
                        sel, self._traj.numAtoms())
            
    def __str__(self):
        
        return 'Frame {0:d} from {1:s}'.format(self._index, 
                                              self._traj.getTitle())

    def numAtoms(self):
        """Return number of atoms."""
        
        return self._traj.numAtoms()
    
    def numSelected(self):
        """Return number of selected atoms."""
        
        if self._sel is None:
            return self._traj.numAtoms()
        else:
            return len(self._indices)
    
    def getIndex(self):
        """Return index."""
        
        return self._index
    
    def getWeights(self):
        """Return coordinate weights for selected atoms."""
        
        if self._sel is None:
            return self._traj.getWeights()
        else:
            return self._traj.getWeights()[self._indices]

    def _getWeights(self):
        
        if self._sel is None:
            return self._traj._getWeights()
        else:
            return self._traj._getWeights()[self._indices]

    def getTrajectory(self):
        """Return the trajectory."""
        
        return self._traj
    
    def getCoords(self):
        """Return a copy of coordinates of (selected) atoms."""
        
        if self._ag is None:
            coords = self._coords
        else:
            coords = self._ag._getCoords()
        
        if self._indices is None:
            return coords.copy()
        else:
            return coords[self._indices]
    
    def _getCoords(self):
        """Return coordinates of (selected) atoms."""
        
        if self._ag is None:
            coords = self._coords
        else:
            coords = self._ag._getCoords()
        
        if self._indices is None:
            return coords
        else:
            return coords[self._indices]
    
    def getUnitcell(self):
        """Return a copy of unitcell array."""
        
        if self._unitcell is not None:
             return self._unitcell.copy()
    
    def _getUnitcell(self):
    
        return self._unitcell
    
    def getDeviations(self):
        """Return deviations from the trajectory reference coordinates."""

        indices = self._indices
        coords = self._getCoords()
        if indices is None:
            return coords - self._traj._coords
        else:
            return coords - self._traj._coords[indices]

    def getRMSD(self):
        """Return RMSD from the trajectory reference coordinates.  If weights 
        for the trajectory are set, weighted RMSD will be returned."""
        
        indices = self._indices 
        traj = self._traj
        coords = self._getCoords()
        if indices is None:
            return measure._calcRMSD(coords, traj._coords, traj._weights)
        else:
            if traj._weights is None:
                return measure._calcRMSD(coords, traj._coords[indices])
            else:
                return measure._calcRMSD(coords, traj._coords[indices], 
                                                 traj._weights[indices])

    def superpose(self):
        """Superpose frame onto the trajectory reference coordinates.  Note 
        that transformation matrix is calculated based on selected atoms and 
        applied to all atoms. If weights for the trajectory are set, they will 
        be used to calculate to transformation."""
        
        indices = self._indices 
        traj = self._traj
        coords = self._getCoords()
        if indices is None:
            self._coords, t = superpose(coords, traj._coords, traj._weights)
        else:
            if traj._weights is None:
                measure._superpose(coords, traj._coords[indices])
            else:
                measure._superpose(coords, traj._coords[indices], 
                                           traj._weights[indices])
