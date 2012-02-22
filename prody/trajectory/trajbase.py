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

"""This module defines base class for trajectory handling."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody.atomic import AtomGroup
from prody.ensemble import Ensemble, checkWeights
from prody.tools import checkCoords

__all__ = ['TrajBase']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class TrajBase(object):
    
    """Base class for :class:`~.Trajectory` and :class:`~.TrajFile`.
    Derived classes must implement functions described in this class."""

    def __init__(self, title='Unknown'):
    
        self._title = str(title).strip()
        self._coords = None         # reference
        self._n_atoms = 0
        self._n_csets = 0 # number of conformations/frames/coordinate sets
        self._weights = None
        self._ag = None
        self._sel = None
        self._indices = None # indices of selected atoms
        self._nfi = 0
        self._closed = False

    def __iter__(self):
        if self._closed:
            raise ValueError('I/O operation on closed file')
        while self._nfi < self._n_csets: 
            yield self.next()
    
    def __getitem__(self, index):
    
        if self._closed: 
            raise ValueError('I/O operation on closed file')
            
        if isinstance(index, int):
            return self.getFrame(index)
            
        elif isinstance(index, (slice, list, np.ndarray)):
            if isinstance(index, slice):
                ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                    self._title, index.indices(len(self))))
            else:
                ens = Ensemble('{0:s} slice'.format(self._title))
            ens.setCoords(self.getCoords())
            if self._weights is not None:
                ens.setWeights(self._weights.copy())
            ens.addCoordset(self.getCoordsets(index))
            return ens
            
        else:
            raise IndexError('invalid index')

    def __len__(self):
    
        return self._n_csets
    
    def __enter__(self):
        
        return self
    
    def __exit__(self, type, value, tb):
        
        self.close()

    def getTitle(self):
        """Return title of the ensemble."""
        
        return self._title

    def setTitle(self, title):
        """Set title of the ensemble."""
        
        self._title = str(title)
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
   
    def numCoordsets(self):
        """Return number of coordinate sets, i.e conformations or frames."""
        
        return self._n_csets
    
    def numSelected(self):  
        """Return number of selected atoms."""
        
        if self._sel is None:
            return self._n_atoms
        return len(self._indices) 
    
    def getAtoms(self):
        """Return associated atom group."""
        
        return self._ag
    
    def setAtoms(self, ag, setref=True):
        """Associate the instance with an :class:`~.AtomGroup`.  Note that at 
        association, active coordinate set of the :class:`~.AtomGroup`, if it 
        has one, will be set as the reference coordinates for the ensemble or 
        trajectory. Changes in :class:`~.AtomGroup` active coordinate set will 
        not be reflected to the reference coordinates. If you want to preserve 
        the present reference coordinates, pass ``setref=False``."""
        
        if ag is None:
            self._ag = None
        else:
            if not isinstance(ag, AtomGroup):
                raise TypeError('ag must be an AtomGroup instance')
            if self._n_atoms != 0 and ag.numAtoms() != self._n_atoms:
                raise ValueError('AtomGroup must have same number of atoms')
            self._ag = ag
            if setref:
                coords = ag.getCoords()
                if coords is not None:
                    self._coords = coords 
                    LOGGER.info('Coordinates of {0:s} is set as the reference '
                               'for {1:s}.'.format(ag.getTitle(), self._title))
        self._sel = None
        self._indices = None
        
    def getSelection(self):
        """Return the current selection. If ``None`` is returned, it means
        that all atoms are selected."""
        
        return self._sel
    
    def _getSelIndices(self):
        return self._indices 
    
    def select(self, selstr):
        """Select a subset atoms. When a subset of atoms are selected, their 
        coordinates will be evaluated in, for example, RMSD calculations. 
        If *selstr* results in selecting no atoms, all atoms will be 
        considered. For more information on atom selections see 
        :ref:`selections`."""
        
        if selstr is None:
            self._sel = None
            self._indices = None
        else:
            if self._ag is None:
                raise AttributeError(self.__class__.__name__ + ' must be be '
                                     'associated with an AtomGroup, use '
                                     '`setAtomGroup` method.')
            sel = self._ag.select(selstr)
            if sel is not None:
                self._indices = sel.getIndices()
            self._sel = sel
            return sel
    
    def getCoords(self):
        """Return a copy of reference coordinates for selected atoms."""
        
        if self._coords is None:
            return None
        if self._sel is None:
            return self._coords.copy()
        return self._coords[self._indices]
    
    def _getCoords(self):
        """Return a view of reference coordinates for selected atoms."""

        if self._coords is None:
            return None
        if self._sel is None:
            return self._coords
        return self._coords[self._indices]

    def setCoords(self, coords):
        """Set reference coordinates."""

        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoords()
            except AttributeError:
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoordinates attribute')
        self._coords = checkCoords(coords, arg='coords', 
                                   n_atoms=self._n_atoms, cset=False)
        
    def getWeights(self):
        """Return a copy of weights of selected atoms."""
        
        if self._weights is None:
            return None
        if self._sel is None:
            return self._weights.copy()
        if self._weights.ndim == 2:
            return self._weights[self._indices]
        else:
            return self._weights[:, self._indices]
    
    def _getWeights(self):

        if self._weights is None:
            return None
        if self._sel is None:
            return self._weights
        if self._weights.ndim == 2:
            return self._weights[self._indices]
        else:
            return self._weights[:, self._indices]

    
    def setWeights(self, weights):
        """Set atomic weights."""
        
        if self._n_atoms == 0:
            raise AttributeError('first set reference coordinates')
        self._weights = checkWeights(weights, self._n_atoms, None)
        
    def numFrames(self):
        """Return number of frames."""
        
        return self._n_csets
    
    def getNextIndex(self):
        """Return the index of the next frame."""
        
        return self._nfi
    
    def iterCoordsets(self):
        """Iterate over coordinate sets. Coordinates for selected atoms are 
        returned. Reference coordinates are not included. Iteration starts
        from the current position in the trajectory."""

        if self._closed:
            raise ValueError('I/O operation on closed file')        
        while self._nfi < self._n_csets: 
            yield self.nextCoordset()
    
    def getCoordsets(self, indices=None):
        """Returns coordinate sets at given *indices*. *indices* may be an 
        integer, a list of ordered integers or ``None``. ``None`` returns all 
        coordinate sets. If a list of indices is given, unique numbers will
        be selected and sorted. That is, this method will always return unique 
        coordinate sets in the order they appear in the trajectory file.
        Shape of the coordinate set array is (n_sets, n_atoms, 3)."""
        
        pass
    
    def getFrame(self, index):
        """Return frame at given *index*."""
        
        pass

    def nextCoordset(self):
        """Return next coordinate set."""
        
        pass
    
    def next(self):
        """Return next frame."""
        
        pass

    def goto(self, n):
        """Go to the frame at index *n*. ``n=0`` will rewind the trajectory
        to the beginning. ``n=-1`` will go to the last frame."""
        
        pass

    def skip(self, n):
        """Skip *n* frames. *n* must be a positive integer."""
        
        pass

    def reset(self):
        """Go to first frame whose index is 0."""
        
        pass

    def close(self):
        """Close trajectory file."""
        
        pass
    
    def hasUnitcell(self):
        """Return ``True`` if trajectory has unitcell data."""
        
        pass
