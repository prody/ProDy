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

"""This module defines classes for handling ensembles of conformations."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from time import time

import numpy as np
from prody import Atomic, AtomGroup
from prody.measure import _calcRMSD, _calcTransformation, _applyTransformation
from prody.measure import _superposeTraj
from prody.tools import checkCoords

from conformation import *

__all__ = ['Ensemble', 'PDBEnsemble']

pkg = __import__(__package__) 
LOGGER = pkg.LOGGER

class Ensemble(object):
    
    """A class for analysis of arbitrary conformational ensembles. 
    
    Indexing (e.g. ``ens[0]``) returns a :class:`~.Conformation` instance that 
    points to a coordinate set in the ensemble. Slicing (e.g. ``ens[:10]``) 
    returns an :class:`Ensemble` instance that contains subset of conformations
    (coordinate sets). The ensemble obtained by slicing will have a copy of 
    the reference coordinates."""

    def __init__(self, title='Unknown'):
        """Instantiate with a *title* or a :class:`~.Atomic` instance.  All 
        coordinate sets from atomic instances will be added to the ensemble 
        automatically."""
        
        self._title = str(title).strip()

        self._coords = None         # reference
        self._n_atoms = 0
        self._n_csets = 0 # number of conformations/frames/coordinate sets
        self._weights = None
        self._ag = None
        self._sel = None
        self._indices = None # indices of selected atoms

        self._confs = None       # coordinate sets
        
        if isinstance(title, (Atomic, Ensemble)):
            self.setCoords(title.getCoords())
            self.addCoordset(title)

    
    def __repr__(self):
    
        if self._sel is None:
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
            LOGGER.info('Atom weights from "{0:s}" are used in "{1:s}".'
                        .format(self._title, ensemble.getTitle()))
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
        trajectory.  Changes in :class:`~.AtomGroup` active coordinate set will
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
                                     '`setAtoms` method.')
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

    def numConfs(self):  
        """Return number of conformations."""

        return self._n_csets

    def addCoordset(self, coords, allcoordsets=True):
        """Add coordinate set(s) as conformation(s).  :class:`~.Atomic` 
        instances are accepted as *coords* argument.  If *allcoordsets* is 
        ``True``, all coordinate sets from the :class:`~.Atomic` instance 
        will be appended to the ensemble.  Otherwise, only the active 
        coordinate set will be appended."""
        
        assert isinstance(allcoordsets, bool), 'allcoordsets must be boolean'
        if not isinstance(coords, np.ndarray):
            if isinstance(coords, (Atomic, Ensemble)):
                atoms = coords
                if allcoordsets:
                    coords = atoms.getCoordsets()
                else:
                    coords = atoms.getCoords()
                if coords is None:
                    raise ValueError('{0:s} must contain coordinate data'
                                     .format(atoms))
            else:
                raise TypeError('coords must be a Numpy array or '
                                'ProDy Atomic or Ensemble instance')
        
        coords = checkCoords(coords, arg='coords', cset=True, 
                             n_atoms=self._n_atoms, reshape=True)
        if self._n_atoms == 0:
            self._n_atoms = coords.shape[-2]
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
        
        if self._sel is None:
            for conf in self._confs:
                yield conf.copy()
        else:
            indices = self._getSelIndices()        
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
            raise AttributeError('coordinates are not set, '
                                 'use setCoordinates() method to set it.')
        if self._confs is None or len(self._confs) == 0: 
            raise AttributeError('conformations are not set, '
                                 'use addCoordset() method to set it.')
        LOGGER.info('Superposing structures.')
        start = time()
        self._superpose()
        LOGGER.info(('Superposition is completed in {0:.2f} '
                           'seconds.').format((time() - start)))
        
    def _superpose(self):
        """Superpose conformations and update coordinates."""
        
        indices = self._indices
        if indices is None:
            _superposeTraj(self._confs, self._coords, 
                           self._weights, self._confs)
        else:
            weights = None
            if self._weights is not None:
                weights = self._weights[indices]
            _superposeTraj(self._confs[:,indices], self._coords[indices], 
                           weights, self._confs)
            
    def iterpose(self, rmsd=0.0001):
        """Iteratively superpose the ensemble until convergence.
        
        Initially, all conformations are aligned with the reference 
        coordinates. Then mean coordinates are calculated, and are set
        as the new reference coordinates. This is repeated until 
        reference coordinates do not change. This is determined by
        the value of RMSD between the new and old reference coordinates.        
        
        :arg rmsd: RMSD (A) between old and new reference coordinates 
                     to converge
        :type rmsd: float, default is 0.0001
            
        """
        
        if self._coords is None:
            raise AttributeError('coordinates are not set, '
                                 'use setCoordinates() method to set it.')
        if self._confs is None or len(self._confs) == 0: 
            raise AttributeError('conformations are not set, '
                                 'use addCoordset() method to set it.')
        LOGGER.info('Starting iterative superposition')
        start = time()
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
            rmsdif = _calcRMSD(self._coords, newxyz)
            self._coords = newxyz
            step += 1
            LOGGER.info(('Step #{0:d}: RMSD difference = '
                               '{1:.4e}').format(step, rmsdif))
        LOGGER.info('Iterative superposition completed in {0:.2f}s.'
                    .format((time() - start)))

    def getMSFs(self):
        """Calculate and return mean square fluctuations (MSFs). 
        Note that you might need to align the conformations using 
        :meth:`superpose` or :meth:`iterpose` before calculating MSFs."""
        
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
        """Calculate and return root mean square fluctuations (RMSFs). 
        Note that you might need to align the conformations using 
        :meth:`superpose` or meth:`iterpose` before calculating RMSFs."""

        return self.getMSFs() ** 0.5
            
    def getDeviations(self):
        """Return deviations from reference coordinates. Note that you
        might need to align the conformations using :meth:`superpose` or 
        :meth:`iterpose` before calculating deviations."""
        
        if not isinstance(self._confs, np.ndarray):
            LOGGER.warning('Conformations are not set.')
            return None
        if not isinstance(self._coords, np.ndarray):
            LOGGER.warning('Coordinates are not set.')
            return None
        
        return self._getCoordsets() - self._coords 
        
    def getRMSDs(self):
        """Calculate and return root mean square deviations (RMSDs). Note that 
        you might need to align the conformations using :meth:`superpose` or 
        :meth:`iterpose` before calculating RMSDs."""
        
        if self._confs is None or self._coords is None: 
            return None
        if self._indices is None:
            return _calcRMSD(self._coords, self._confs, self._weights)
        else:
            indices = self._indices
            if self._weights is None:
                return _calcRMSD(self._coords[indices], self._confs[:,indices])
            else:
                return _calcRMSD(self._coords[indices], 
                                 self._confs[:,indices],
                                 self._weights[indices])
                

class PDBEnsemble(Ensemble):
    
    """This class enables handling coordinates for heterogeneous structural 
    datasets and stores identifiers for individual conformations.
    
    |example| See usage usage in :ref:`pca-xray`, :ref:`pca-dimer`, and 
    :ref:`pca-blast`.
    
    .. note:: This class is designed to handle conformations with missing
       coordinates, e.g. atoms that are note resolved in an X-ray structure.
       For unresolved atoms, the coordinates of the reference structure is
       assumed in RMSD calculations and superpositions."""

    def __init__(self, title='Unknown'):
        self._labels = []
        Ensemble.__init__(self, title)
        
    def __repr__(self):
        return '<PDB' + Ensemble.__repr__(self)[1:]
    
    def __str__(self):
        return 'PDB' + Ensemble.__str__(self)
    
    def __add__(self, other):
        """Concatenate two ensembles. The reference coordinates of *self* is 
        used in the result."""
        
        if not isinstance(other, Ensemble):
            raise TypeError('an Ensemble instance cannot be added to an {0:s} '
                            'instance'.format(type(other)))
        elif self.numAtoms() != other.numAtoms():
            raise ValueError('Ensembles must have same number of atoms.')
    
        ensemble = PDBEnsemble('{0:s} + {1:s}'.format(self.getTitle(), 
                                                   other.getTitle()))
        ensemble.setCoords(self._coords.copy())
        ensemble.addCoordset(self._confs.copy(), self._weights.copy())
        other_weights = other.getWeights()
        if other_weights is None:
            ensemble.addCoordset(other.getCoordsets())
        else:
            ensemble.addCoordset(other._confs.copy(), other_weights)
        return ensemble
    
    def __iter__(self):
        n_confs = self._n_csets
        for i in range(n_confs):
            if n_confs != self._n_csets:
                raise RuntimeError('number of conformations in the ensemble '
                                   'changed during iteration')
            yield PDBConformation(self, i)
    
    def __getitem__(self, index):
        """Return a conformation at given index."""
        
        if isinstance(index, int):
            return self.getConformation(index) 
        elif isinstance(index, slice):
            ens = PDBEnsemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                self._title, index.indices(len(self))))
            ens.setCoords(self.getCoords())
            ens.addCoordset(self._confs[index].copy(), 
                            self._weights[index].copy())
            return ens
        elif isinstance(index, (list, np.ndarray)):
            ens = PDBEnsemble('Conformations of {0:s}'.format(self._title))
            ens.setCoords(self.getCoords())
            ens.addCoordset(self._confs[index].copy(), 
                            self._weights[index].copy())
            return ens
        else:
            raise IndexError('invalid index')
            
    def _superpose(self):
        """Superpose conformations and return new coordinates."""

        calcT = _calcTransformation
        applyT = _applyTransformation
        if self._sel is None:
            weights = self._weights
            coords = self._coords
            confs = self._confs
            for i, conf in enumerate(confs):
                confs[i] = applyT(calcT(conf, coords, weights[i]), confs[i])
        else:            
            indices = self._getSelIndices()
            weights = self._weights[:, indices]
            coords = self._coords[indices]
            confs_selected = self._confs[:,indices]
            confs = self._confs
            for i, conf in enumerate(confs_selected):
                confs[i] = applyT(calcT(conf, coords, weights[i]), confs[i]) 

    def addCoordset(self, coords, weights=None, allcoordsets=True):
        """Add coordinate set(s) as conformation(s).

        :class:`~.Atomic` instances are accepted as *coords* argument.  
        If *allcoordsets* is ``True``, all coordinate sets from the 
        :class:`~.Atomic` instance will be appended to the ensemble. 
        Otherwise, only the active coordinate set will be appended.

        *weights* is an optional argument. If provided, its length must
        match number of atoms. Weights of missing (not resolved) atoms 
        must be equal to ``0`` and weights of those that are resolved
        can be anything greater than ``0``. If not provided, weights of 
        atoms in this coordinate set will be set equal to ``1``."""
        
        assert isinstance(allcoordsets, bool), 'allcoordsets must be boolean'
        if weights is not None:
            assert isinstance(weights, np.ndarray), 'weights must be ndarray'
            
        ag = None
        if isinstance(coords, Atomic):
            atoms = coords
            if isinstance(coords, AtomGroup):
                ag = atoms
            else:
                ag = atoms.getAtomGroup()
            if allcoordsets:
                coords = atoms.getCoordsets()
            else: 
                coords = atoms.getCoords()
            title = ag.getTitle() 
        elif isinstance(coords, np.ndarray):
            title = 'Unknown'
        else:
            title = str(coords)
            try:
                if allcoordsets:
                    coords = coords.getCoordsets()
                else: 
                    coords = coords.getCoords()
            except AttributeError:            
                raise TypeError('coords must be a Numpy array or must have '
                                'getCoordinates attribute')

        coords = checkCoords(coords, 'coords', cset=True, 
                             n_atoms=self._n_atoms, reshape=True)
        n_csets, n_atoms, _ = coords.shape
        if self._n_atoms == 0:
            self._n_atoms = n_atoms
        if weights is None:
            weights = np.ones((n_csets, n_atoms, 1), dtype=float)
        else:
            weights = checkWeights(weights, n_atoms, n_csets)

        while '  ' in title:
            title = title.replace('  ', ' ')
        title = title.replace(' ', '_')
        
        if n_csets > 1:
            self._labels += ['{0:s}_{1:d}'
                                  .format(title, i+1) for i in range(n_csets)]
        else:
            if ag is not None and ag.numCoordsets() > 1:
                self._labels.append('{0:s}_{1:d}'.format(title, 
                                         atoms.getACSIndex()))
            else:                
                self._labels.append(title)
        if self._confs is None and self._weights is None:
            self._confs = coords
            self._weights = weights
            self._n_csets = n_csets
        elif self._confs is not None and self._weights is not None:
            self._confs = np.concatenate((self._confs, coords), axis=0)
            self._weights = np.concatenate((self._weights, weights), axis=0)
            self._n_csets += n_csets
        else:
            raise RuntimeError('_confs and _weights must be set or None at '
                               'the same time')

    def getLabels(self):
        """Return identifiers of the conformations in the ensemble."""
        
        return list(self._labels)
    
    def getCoordsets(self, indices=None):
        """Return a copy of coordinate set(s) at given *indices* for selected 
        atoms. *indices* may be an integer, a list of integers or ``None``. 
        ``None`` returns all coordinate sets. 
    
        .. warning:: When there are atoms with weights equal to zero (0),
           their coordinates will be replaced with the coordinates of the
           ensemble reference coordinate set.
        """
        
        if self._confs is None:
            return None
    
        if indices is None:
            indices = slice(None)
        elif isinstance(indices, (int, long)): 
            indices = np.array([indices])
        coords = self._coords
        if self._sel is None:
            confs = self._confs[indices].copy()
            for i, w in enumerate(self._weights[indices]):
                which = w.flatten()==0
                confs[i, which] = coords[which]
        else:
            selids = self._getSelIndices()
            coords = coords[selids]
            confs = self._confs[indices, selids]
            for i, w in enumerate(self._weights[indices]):
                which = w[selids].flatten()==0
                confs[i, which] = coords[which]
        return confs
    
    _getCoordsets = getCoordsets
    
    def iterCoordsets(self):
        """Iterate over coordinate sets. A copy of each coordinate set for
        selected atoms is returned. Reference coordinates are not included."""
        
        conf = PDBConformation(self, 0)
        for i in range(self._n_csets):
            conf._index = i
            yield conf.getCoords()
   
    def delCoordset(self, index):
        """Delete a coordinate set from the ensemble."""
        
        Ensemble.delCoordset(self, index)
        if isinstance(index, int):
            index = [index]
        else:
            index = list(index)
        index.sort(reverse=True)
        for i in index:
            self._labels.pop(i)
            
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
            return PDBConformation(self, index)
        else:
            raise IndexError('conformation index out of range')

    def getMSFs(self):
        """Calculate and return mean square fluctuations (MSFs). 
        Note that you might need to align the conformations using 
        :meth:`superpose` or :meth:`iterpose` before calculating MSFs."""
        
        if self._confs is None: 
            return
        indices = self._indices
        if indices is None:
            coords = self._coords
            confs = self._confs
            weights = self._weights > 0
        else:
            coords = self._coords[indices]
            confs = self._confs[:,indices]
            weights = self._weights[:,indices] > 0
        weightsum = weights.sum(0)
        mean = np.zeros(coords.shape)
        for i, conf in enumerate(confs):
            mean += conf * weights[i]
        mean /= weightsum
        ssqf = np.zeros(mean.shape)
        for i, conf in enumerate(confs):
            ssqf += ((conf - mean) * weights[i]) ** 2
        return ssqf.sum(1) / weightsum.flatten()
    
    def getRMSDs(self):
        """Calculate and return root mean square deviations (RMSDs). Note that 
        you might need to align the conformations using :meth:`superpose` or 
        :meth:`iterpose` before calculating RMSDs."""
        
        if self._confs is None or self._coords is None: 
            return None
    
        if self._sel is None:
            return _calcRMSD(self._coords, self._confs, self._weights)
        else:
            indices = self._indices
            return _calcRMSD(self._coords[indices], 
                             self._confs[:,indices],
                             self._weights[:, indices])

    def setWeights(self, weights):
        """Set atomic weights."""
        
        if self._n_atoms == 0:
            raise AttributeError('coordinates are not set')
        elif not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape[:2] != (self._n_csets, self._n_atoms):
            raise ValueError('shape of weights must (n_confs, n_atoms[, 1])')
        if weights.dtype not in (np.float32, float):
            try:
                weights = weights.astype(float)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(float))
        if np.any(weights < 0):
            raise ValueError('weights must greater or equal to 0')
            
        if weights.ndim == 2:
            weights = weights.reshape((self._n_csets, self._n_atoms, 1))
        self._weights = weights


