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

from prody import LOGGER
from prody.atomic import AtomGroup
from prody.ensemble import Ensemble, checkWeights
from prody.utilities import checkCoords

from .frame import Frame

__all__ = ['TrajBase']


class TrajBase(object):

    """Base class for :class:`.Trajectory` and :class:`.TrajFile`.  Derived
    classes must implement functions described in this class."""

    def __init__(self, title='Unknown'):

        self._title = str(title).strip()
        self._coords = None         # reference
        self._n_atoms = 0
        self._n_csets = 0 # number of conformations/frames/coordinate sets
        self._weights = None
        self._ag = None
        self._atoms = None
        self._indices = None # indices of selected atoms
        self._frame = None # if atoms are set, always return the same frame
        self._nfi = 0
        self._closed = False

    def __iter__(self):

        if self._closed:
            raise ValueError('I/O operation on closed file')
        while self._nfi < self._n_csets:
            yield self.next()

    def __str__(self):

        return '{0:s} {1:s}'.format(self.__class__.__name__, self._title)

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

    def numFrames(self):
        """Return number of frames."""

        return self._n_csets

    numCoordsets = numFrames

    def numSelected(self):
        """Return number of selected atoms.  A subset of atoms can be selected 
        by passing a selection to :meth:`setAtoms`."""

        return self._n_atoms if self._indices is None else len(self._indices)

    def getAtoms(self):
        """Return associated/selected atoms."""

        return self._atoms

    def setAtoms(self, atoms):
        """Set *atoms* or specify a selection of atoms to be considered in 
        calculations and coordinate requests.  When a selection is set, 
        corresponding subset of coordinates will be considered in, for 
        example, alignments and RMSD calculations.  Setting atoms also 
        allows some functions to access atomic data when needed.  For 
        example, :class:`.Trajectory` and :class:`.Frame` instances become 
        suitable arguments for :func:`.writePDB`.  Passing **None** as *atoms* 
        argument will deselect atoms."""

        if atoms is None:
            self._atoms = self._indices = None
            return
        
        try:
            atoms.getACSIndex()
        except AttributeError:
            raise TypeError('atoms must be an Atomic instance')
        
        n_atoms = self._n_atoms
        if n_atoms:
            
            if atoms.numAtoms() > n_atoms:
                raise ValueError('atoms must be same size or smaller than '
                                   'the trajectory')
            
            elif atoms.numAtoms() == n_atoms:
                self._atoms = atoms
                self._indices = None
            
            else:
                try:
                    ag = atoms.getAtomGroup()
                except AttributeError:
                    raise ValueError('atoms must indicate a subset or must '
                                     'match the trajectory size')
                else:
                    if ag.numAtoms() != n_atoms:
                        raise ValueError('atoms must point to an AtomGroup '
                                         'of the same size as the trajectory')
                    self._atoms = atoms
                    self._indices = atoms.getIndices()
        
        else:
            self._n_atoms = atoms.numAtoms()
            self._atoms = atoms

    def link(self, ag):
        """Link :class:`.AtomGroup` instance *ag* to the trajectory.  When a 
        new frame is parsed from the trajectory file, coordinates of *ag* and 
        of all selections and atom subsets pointing to it, will be updated. 
        To break an established link, pass **None** argument.

        .. warning::

           Every time a frame is parsed from the trajectory, all coordinate 
           sets present in the linked :class:`.AtomGroup` will be overwritten.
        """

        if ag is None:
            self._ag = None
        else:
            try:
                ag.getACSIndex()
            except AttributeError:
                raise TypeError('ag must be an AtomGroup instance')
            try:
                ag.getAtomGroup()
            except AttributeError:
                pass
            else:
                raise TypeError('ag must be an AtomGroup instance')
            if self.numAtoms() != ag.numAtoms():
                raise ValueError('atom group and trajectory must '
                                   'have same number of atoms')
            self._ag = ag
            self._frame = Frame(None, None, None)
            self._atoms = self._atoms or ag

    def getLinked(self):
        """Return linked :class:`.AtomGroup` instance, or **None** if a link
        is not established."""
        
        return self._ag
    
    def isLinked(self):
        """Return **True** if trajectory is linked to an :class:`.AtomGroup`
        instance."""
        
        return self._ag is not None

    def getCoords(self):
        """Return a copy of reference coordinates for (selected) atoms."""

        if self._coords is None:
            return None
        if self._indices is None:
            return self._coords.copy()
        return self._coords[self._indices]

    def _getCoords(self):
        """Return a view of reference coordinates for (selected) atoms."""

        if self._coords is None:
            return None
        if self._indices is None:
            return self._coords
        return self._coords[self._indices]

    def setCoords(self, coords):
        """Set *coords* as the trajectory reference coordinate set.  *coords*
        must be an object with :meth:`getCoords` method, or a Numpy array with
        suitable data type, shape, and dimensionality."""

        atoms = coords
        try:
            coords = atoms.getCoords()
        except AttributeError:
            pass
        else:
            if coords is None:
                raise ValueError('coordinates of {0:s} are not set'
                                   .format(str(atoms)))

        try:
            checkCoords(coords, natoms=self._n_atoms)
        except TypeError:
            raise TypeError('coords must be a numpy array or an object '
                            'with `getCoords` method')

        self._coords = coords

    def setWeights(self, weights):
        """Set atomic weights."""

        if self._n_atoms == 0:
            raise AttributeError('coordinates must be set first')
        self._weights = checkWeights(weights, self._n_atoms, None)

    def getWeights(self):
        """Return a copy of weights of (selected) atoms."""

        if self._weights is not None:
            if self._indices is None:
                return self._weights.copy()
            else:
                return self._weights[self._indices]

    def _getWeights(self):

        if self._weights is not None:
            if self._indices is None:
                return self._weights
            else:
                return self._weights[self._indices]

    def nextIndex(self):
        """Return the index of the next frame."""

        return self._nfi

    def iterCoordsets(self):
        """Yield coordinate sets for (selected) atoms. Reference coordinates
        are not included. Iteration starts from the next frame in line."""

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
        """Return next coordinate set in a :class:`.Frame` instance.  Note that
        when atoms are set for the trajectory, this method will return the same
        frame instance after updating its coordinates."""

        pass

    def goto(self, n):
        """Go to the frame at index *n*. ``n=0`` will rewind the trajectory
        to the beginning, same as calling :meth:`reset` method. ``n=-1``
        will go to the last frame.  Frame *n* will not be parsed until one
        of :meth:`next` or :meth:`nextCoordset` methods is called."""

        pass

    def skip(self, n):
        """Skip *n* frames.  *n* must be a positive integer.  Skipping some
        frames will only change the next frame index (:meth:`nextIndex`)
        Next frame will not be parsed until one of :meth:`next` or
        :meth:`nextCoordset` methods is called."""

        pass

    def reset(self):
        """Go to first frame at index 0.  First frame will not be parsed until
        one of :meth:`next` or :meth:`nextCoordset` methods is called."""

        pass

    def close(self):
        """Close trajectory file."""

        pass

    def hasUnitcell(self):
        """Return ``True`` if trajectory has unitcell data."""

        pass
