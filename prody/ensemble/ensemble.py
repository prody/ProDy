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

"""This module defines a class for handling ensembles of conformations."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup
from prody.measure import getRMSD
from prody.utilities import importLA, checkCoords

from conformation import *

__all__ = ['Ensemble']

class Ensemble(object):
    
    """A class for analysis of arbitrary conformational ensembles. 
    
    Indexing (e.g. ``ens[0]``) returns a :class:`.Conformation` instance that 
    points to a coordinate set in the ensemble.  Slicing (e.g. ``ens[0:10]``) 
    returns an :class:`Ensemble` instance that contains a copy of the subset 
    of conformations (coordinate sets). """

    def __init__(self, title='Unknown'):
        """Instantiate with a *title* or a :class:`.Atomic` instance.  All 
        coordinate sets from atomic instances will be added to the ensemble."""
        
        self._title = str(title).strip()

        self._coords = None         # reference
        self._n_atoms = 0
        self._n_csets = 0 # number of conformations/frames/coordinate sets
        self._weights = None
        self._atoms = None
        self._indices = None # indices of selected atoms

        self._confs = None       # coordinate sets
        
        if isinstance(title, (Atomic, Ensemble)):
            self.setCoords(title.getCoords())
            self.addCoordset(title)
    
    def __repr__(self):
    
        if self._indices is None:
            return ('<Ensemble: {0:s} ({1:d} conformations; {2:d} atoms)>'
                    ).format(self._title, len(self), self._n_atoms)
        else:
            return ('<Ensemble: {0:s} ({1:d} conformations; selected {2:d} of '
                    '{3:d} atoms)>').format(self._title, len(self), 
                                           self.numSelected(), self._n_atoms)

    def __str__(self):
    
        return 'Ensemble {0:s}'.format(self._title)

    def __len__(self):
    
        return self._n_csets
    
    def __getitem__(self, index):
        """Return a conformation at given index."""
        
        if self._confs is None:
            return None
        
        # use sth like follows
        # which = np.arange(self._n_csets)[index].nonzero()[0]
        # if len(which) == 1:
        #       return getConf...
        # else: 
        #     return SubEnsemble
        
        if isinstance(index, int):
            return self.getConformation(index) 
            
        elif isinstance(index, slice):
            ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                self._title, index.indices(len(self))))
            ens.setCoords(self.getCoords())
            ens.addCoordset(self.getCoordsets(index))
            if self._weights is not None:
                ens.setWeights(self.getWeights())
            return ens
            
        elif isinstance(index, (list, np.ndarray)):
            ens = Ensemble('Conformations of {0:s}'.format(self._title))
            ens.setCoords(self.getCoords())
            ens.addCoordset(self.getCoordsets(index))
            if self._weights is not None:
                ens.setWeights(self.getWeights())
            return ens
            
        else:
            raise IndexError('invalid index')
    
    def __add__(self, other):
        """Concatenate ensembles. The reference coordinates and weights of 
        *self* is used in the resulting ensemble."""
        
        if not isinstance(other, Ensemble):
            raise TypeError('an Ensemble instance cannot be added to an {0:s} '
                            'instance'.format(type(other)))
        elif self.numAtoms() != other.numAtoms():
            raise ValueError('Ensembles must have same number of atoms.')
    
        ensemble = Ensemble('{0:s} + {1:s}'.format(self.getTitle(), 
                                                   other.getTitle()))
        ensemble.setCoords(self._coords.copy())
        ensemble.addCoordset(self._confs.copy())
        ensemble.addCoordset(other.getCoordsets())
        if self._weights is not None: 
            LOGGER.info('Atom weights from {0:s} are used in {1:s}.'
                        .format(repr(self._title), repr(ensemble.getTitle())))
            ensemble.setWeights(self._weights)
        return ensemble
    
    def __iter__(self):
    
        n_csets = self._n_csets
        for i in range(n_csets):
            if n_csets != self._n_csets:
                raise RuntimeError('number of conformations in the ensemble '
                                   'changed during iteration')
            yield Conformation(self, i)

    def getTitle(self):
        """Return title of the ensemble."""
        
        return self._title

    def setTitle(self, title):
        """Set title of the ensemble."""
        
        self._title = str(title)
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
   
    def numConfs(self):
        """Return number of conformations."""
        
        return self._n_csets
    
    numCoordsets = numConfs
    
    def numSelected(self):  
        """Return number of selected atoms.  Number of all atoms will be 
        returned if a selection is not made.  A subset of atoms can be 
        selected by passing a selection to :meth:`setAtoms`."""
        
        return self._n_atoms if self._indices is None else len(self._indices) 
    
    def getAtoms(self):
        """Return associated atom group."""
        
        return self._atoms
    
    def setAtoms(self, atoms):
        """Set *atoms* or specify a selection of atoms to be considered in 
        calculations and coordinate requests.  When a selection is set, 
        corresponding subset of coordinates will be considered in, for 
        example, alignments and RMSD calculations.  Setting atoms also 
        allows some functions to access atomic data when needed.  For 
        example, :class:`.Ensemble` and :class:`.Conformation` instances 
        become suitable arguments for :func:`.writePDB`.  Passing **None** 
        as *atoms* argument will deselect atoms."""
        
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
                                 'the ensemble')
            
            elif atoms.numAtoms() == n_atoms:
                self._atoms = atoms
                self._indices = None
            
            else:
                try:
                    ag = atoms.getAtomGroup()
                except AttributeError:
                    raise ValueError('atoms must indicate a subset or must '
                                     'match the ensemble size')
                else:
                    if ag.numAtoms() != n_atoms:
                        raise ValueError('atoms must point to an AtomGroup '
                                         'of the same size as the ensemble')
                    self._atoms = atoms
                    self._indices = atoms.getIndices()
        
        else:
            self._n_atoms = atoms.numAtoms()
            self._atoms = atoms
            
    def getCoords(self):
        """Return a copy of reference coordinates for selected atoms."""
        
        if self._coords is None:
            return None
        if self._indices is None:
            return self._coords.copy()
        return self._coords[self._indices]
    
    def _getCoords(self):
        """Return a view of reference coordinates for selected atoms."""

        if self._coords is None:
            return None
        if self._indices is None:
            return self._coords
        return self._coords[self._indices]

    def setCoords(self, coords):
        """Set *coords* as the ensemble reference coordinate set.  *coords* 
        may be an array with suitable data type, shape, and dimensionality, or
        an object with :meth:`getCoords` method."""

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
        self._n_atoms = coords.shape[0]
        
    def getWeights(self):
        """Return a copy of weights of selected atoms."""
        
        if self._weights is None:
            return None
        if self._indices is None:
            return self._weights.copy()
        if self._weights.ndim == 2:
            return self._weights[self._indices]
        else:
            return self._weights[:, self._indices]
    
    def _getWeights(self):

        if self._weights is None:
            return None
        if self._indices is None:
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

    def addCoordset(self, coords):
        """Add coordinate set(s) to the ensemble.  *coords* must be a Numpy
        array with suitable data type, shape and dimensionality, or an object 
        with :meth:`getCoordsets` method."""
        
        n_atoms = self._n_atoms
        try:
            if self._coords is not None and hasattr(coords, '_getCoordsets'): 
                coords = coords._getCoordsets()
            else:
                coords = coords.getCoordsets()
                
        except AttributeError:
            pass
        else:
            if coords is None:
                raise ValueError('coordinates are not set')
            
        try:
            checkCoords(coords, csets=True, natoms=n_atoms)
        except TypeError:
            raise TypeError('coords must be a numpy array or an object '
                            'with `getCoords` method')
        
        if not n_atoms:
            self._n_atoms = n_atoms = coords.shape[-2]
            
        if coords.ndim == 2:
            coords = coords.reshape((1, n_atoms, 3))
            n_confs = 1
        else:
            n_confs = coords.shape[0]
            
        if self._confs is None: 
            self._confs = coords
        else:
            self._confs = np.concatenate((self._confs, coords), axis=0)
        self._n_csets += n_confs

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate set(s) at given *indices*, which may be
        an integer, a list of integers or ``None``. ``None`` returns all 
        coordinate sets.  For reference coordinates, use :meth:`getCoordinates`
        method."""
        
        if self._confs is None:
            return None
        if self._indices is None:
            if indices is None:
                return self._confs.copy()
            if isinstance(indices, (int, long, slice)): 
                return self._confs[indices].copy()
            if isinstance(indices, (list, np.ndarray)):        
                return self._confs[indices]
        else:
            selids = self._indices
            if indices is None:
                return self._confs[:,selids]
            if isinstance(indices, (int, long, slice)): 
                return self._confs[indices, selids]
            if isinstance(indices, (list, np.ndarray)):        
                return self._confs[indices, selids]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')

    def _getCoordsets(self, indices=None):

        if self._confs is None:
            return None
        if self._indices is None:
            if indices is None:
                return self._confs
            if isinstance(indices, (int, long, slice)): 
                return self._confs[indices]
            if isinstance(indices, (list, np.ndarray)):        
                return self._confs[indices]
        else:
            selids = self._indices
            if indices is None:
                return self._confs[:,selids]
            if isinstance(indices, (int, long, slice)): 
                return self._confs[indices, selids]
            if isinstance(indices, (list, np.ndarray)):        
                return self._confs[indices, selids]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')
    
    def delCoordset(self, index):
        """Delete a coordinate set from the ensemble."""
        
        if isinstance(index, int):
            index = [index]
        else:
            index = list(index)
        length = self._n_csets
        which = np.ones(length, np.bool)
        which[index] = False
        if which.sum() == 0:
            self._confs = None
            self._weights = None
        else:
            self._confs = self._confs[which]
            if self._weights is not None:
                self._weights = self._weights[which]
        self._n_csets -= len(index)

    def iterCoordsets(self):
        """Iterate over coordinate sets. A copy of each coordinate set for
        selected atoms is returned. Reference coordinates are not included."""
        
        if self._indices is None:
            for conf in self._confs:
                yield conf.copy()
        else:
            indices = self._indices        
            for conf in self._confs:
                yield conf[indices].copy()
    
    def getConformation(self, index):
        """Return conformation at given index."""
        
        if self._confs is None:
            raise AttributeError('conformations are not set')
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        n_confs = self._n_csets
        if -n_confs <= index < n_confs:
            if index < 0:
                index = n_confs + index
            return Conformation(self, index)
        else:
            raise IndexError('conformation index out of range')
            
    def superpose(self):
        """Superpose the ensemble onto the reference coordinates."""
        
        if self._coords is None:
            raise ValueError('coordinates are not set, use `setCoords`')
        if self._confs is None or len(self._confs) == 0: 
            raise ValueError('conformations are not set, use `addCoordset`')
        LOGGER.timeit()
        self._superpose(trans=True) # trans kwarg is used by PDBEnsemble
        LOGGER.timing('Superposition completed in %.2f seconds.')
        
    def _superpose(self, **kwargs):
        """Superpose conformations and update coordinates."""
        
        indices = self._indices
        weights = self._weights
        mobs = self._confs
        if indices is None:
            idx = False
            tar = self._coords
            movs = None
        else:
            idx = True
            if self._weights is not None:
                weights = weights[indices]
            tar = self._coords[indices]
            movs = self._confs
                           
        linalg = importLA()
        svd = linalg.svd
        det = linalg.det
        dot = np.dot
        add = np.add
        subtract = np.subtract
        array = np.array
        sign = np.sign
        
        if weights is None:
            tar_com = tar.mean(0)
            tar_org = (tar - tar_com)
            mob_org = np.zeros(tar_org.shape, dtype=mobs.dtype)
            tar_org = tar_org.T
        else:
            weights_sum = weights.sum()
            weights_dot = dot(weights.T, weights)
            tar_com = (tar * weights).sum(axis=0) / weights_sum
            tar_org = (tar - tar_com)
            mob_org = np.zeros(tar_org.shape, dtype=mobs.dtype)

        LOGGER.progress('Superposing ', len(mobs))
        for i, mob in enumerate(mobs):        
            if idx:            
                mob = mob[indices]
            if weights is None:
                mob_com = mob.mean(0)        
                matrix = dot(tar_org, subtract(mob, mob_com, mob_org))
            else:
                mob_com = (mob * weights).sum(axis=0) / weights_sum
                subtract(mob, mob_com, mob_org)
                matrix = np.dot((tar_org * weights).T, 
                                (mob_org * weights)) / weights_dot
                
            U, s, Vh = svd(matrix)
            Id = array([ [1, 0, 0], [0, 1, 0], [0, 0, sign(det(matrix))] ])
            rotation = dot(Vh.T, dot(Id, U.T))

            if movs is None:
                mobs[i] = dot(mob_org, rotation) 
                add(mobs[i], tar_com, mobs[i]) 
            else:
                add(dot(movs[i], rotation), 
                    (tar_com - dot(mob_com, rotation)), movs[i])
            LOGGER.update(i)
        LOGGER.clear()
            
    def iterpose(self, rmsd=0.0001):
        """Iteratively superpose the ensemble until convergence.  Initially, 
        all conformations are aligned with the reference coordinates.  Then 
        mean coordinates are calculated, and are set as the new reference 
        coordinates.  This is repeated until reference coordinates do not 
        change.  This is determined by the value of RMSD between the new and 
        old reference coordinates.  Note that at the end of the iterative 
        procedure the reference coordinate set will be average of conformations
        in the ensemble. 
        
        :arg rmsd: change in reference coordinates to determine convergence,
            default is 0.0001 Å RMSD
        :type rmsd: float"""
        
        if self._coords is None:
            raise AttributeError('coordinates are not set, use `setCoords`')
        if self._confs is None or len(self._confs) == 0: 
            raise AttributeError('conformations are not set, use' 
                                    '`addCoordset`')
        LOGGER.info('Starting iterative superposition:')
        LOGGER.timeit()
        rmsdif = 1
        step = 0
        weights = self._weights
        if weights is not None and weights.ndim == 3:
            weightsum = weights.sum(axis=0)
        length = len(self)
        while rmsdif > rmsd:
            self._superpose()
            if weights is None:
                newxyz = self._confs.sum(0) / length
            else:
                newxyz = (self._confs * weights).sum(0) / weightsum
            rmsdif = getRMSD(self._coords, newxyz)
            self._coords = newxyz
            step += 1
            LOGGER.info(('Step #{0:d}: RMSD difference = '
                               '{1:.4e}').format(step, rmsdif))
        LOGGER.timing('Iterative superposition completed in %.2fs.')

    def getMSFs(self):
        """Return mean square fluctuations (MSFs) for selected atoms.  
        Conformations can be aligned using one of :meth:`superpose` or 
        :meth:`iterpose` methods prior to MSF calculation."""
        
        if self._confs is None: 
            return
        indices = self._indices
        if indices is None:
            mean = self._confs.mean(0)
            ssqf = np.zeros(mean.shape)
            for conf in self._confs:
                ssqf += (conf - mean) ** 2
        else:
            mean = self._confs[:, indices].mean(0)
            ssqf = np.zeros(mean.shape)
            for conf in self._confs[:, indices]:
                ssqf += (conf - mean) ** 2
        return ssqf.sum(1) / self._n_csets
    
    def getRMSFs(self):
        """Return root mean square fluctuations (RMSFs) for selected atoms.  
        Conformations can be aligned using one of :meth:`superpose` or 
        :meth:`iterpose` methods prior to RMSF calculation."""

        return self.getMSFs() ** 0.5
            
    def getDeviations(self):
        """Return deviations from reference coordinates for selected atoms.  
        Conformations can be aligned using one of :meth:`superpose` or 
        :meth:`iterpose` methods prior to calculating deviations."""
        
        if not isinstance(self._confs, np.ndarray):
            LOGGER.warning('Conformations are not set.')
            return None
        if not isinstance(self._coords, np.ndarray):
            LOGGER.warning('Coordinates are not set.')
            return None
        
        return self._getCoordsets() - self._coords 
        
    def getRMSDs(self):
        """Return root mean square deviations (RMSDs) for selected atoms.  
        Conformations can be aligned using one of :meth:`superpose` or 
        :meth:`iterpose` methods prior to RMSD calculation."""
        
        if self._confs is None or self._coords is None: 
            return None
        if self._indices is None:
            return getRMSD(self._coords, self._confs, self._weights)
        else:
            indices = self._indices
            if self._weights is None:
                return getRMSD(self._coords[indices], self._confs[:,indices])
            else:
                return getRMSD(self._coords[indices], self._confs[:,indices],
                               self._weights[indices])
