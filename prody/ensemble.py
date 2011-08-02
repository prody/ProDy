# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module defines a class for working on conformational ensembles and 
also arbitrary coordinate data sets.

Classes
-------

  * :class:`Ensemble`
  * :class:`Conformation`
  * :class:`PDBEnsemble`
  * :class:`PDBConformation`
  * :class:`Trajectory`
  * :class:`DCDFile`
  * :class:`Frame`
  
Inheritance Diagram
-------------------

.. inheritance-diagram:: prody.ensemble
   :parts: 1

Functions
---------

    * :func:`parseDCD`
    * :func:`saveEnsemble`
    * :func:`loadEnsemble`
    * :func:`calcSumOfWeights`
    * :func:`showSumOfWeights`
    * :func:`trimEnsemble`
    
    
Examples
--------

Results from the example :ref:`pca-xray-calculations` will be used to
illustrate class methods and functions in the module.

>>> from prody import *

>>> ensemble = loadEnsemble('p38_X-ray.ens.npz')

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

from time import time
import os.path
import numpy as np
import prody
from prody import ProDyLogger as LOGGER
from prody import measure

__all__ = ['Ensemble', 'Conformation', 'PDBEnsemble', 'PDBConformation',
           'Trajectory', 'DCDFile', 'Frame', 
           'saveEnsemble', 'loadEnsemble', 
           'calcSumOfWeights', 'showSumOfWeights', 'trimEnsemble',
           'parseDCD']
        
class EnsembleError(Exception):
    pass

plt = None

class EnsembleBase(object):
    
    """Base class for :class:`Ensemble` and :class:`TrajectoryBase`. 
    This class provides methods for associating instances with an 
    :class:`~prody.atomic.AtomGroup` and handling reference coordinates."""

    def __init__(self, name):
        self._name = str(name)
        if self._name == '':
            self._name = 'Unnamed'
        self._coords = None         # reference
        self._n_atoms = 0
        self._n_confs = 0 # number of conformations/frames/coordinate sets
        self._weights = None
        self._ag = None
        self._sel = None
        self._indices = None # indices of selected atoms

    def __len__(self):
        return self._n_confs

    def getName(self):
        """Return name of the instance."""
        
        return self._name

    def setName(self, name):
        """Set name of the instance."""
        
        self._name = str(name)
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
   
    def getNumOfCoordsets(self):
        """Return number of coordinate sets, i.e conformations or frames."""
        
        return self._n_confs
    
    def getNumOfSelected(self):
        """Return number of selected atoms."""
        
        if self._sel is None:
            return self._n_atoms
        return len(self._indices) 
    
    def getAtomGroup(self):
        """Return associated atom group instance."""
        
        return self._ag
    
    def setAtomGroup(self, ag, setref=True):
        """Associate the instance with an :class:`~prody.atomic.AtomGroup`.
        
        Note that at association, active coordinate set of the 
        :class:`~prody.atomic.AtomGroup`, if it has one, will be set as 
        the reference coordinates for the ensemble or trajectory. Changes in 
        :class:`~prody.atomicAtomGroup` active coordinate set will not be
        reflected to the reference coordinates. If you want to preserve the 
        present reference coordinates, pass ``setref=False``. 
        
        """
        
        if ag is None:
            self._ag = None
        else:
            if not isinstance(ag, prody.AtomGroup):
                raise TypeError('ag must be an AtomGroup instance')
            if self._n_atoms != 0 and ag.getNumOfAtoms() != self._n_atoms:
                raise ValueError('AtomGroup must have same number of atoms')
            self._ag = ag
            if setref:
                coords = ag.getCoordinates()
                if coords is not None:
                    self._coords = coords 
                    LOGGER.info('Coordinates of {0:s} is set as the reference '
                                'for {1:s}.'.format(ag.getName(), self._name))
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
                raise AttributeError('instance needs to be associated with an '
                                     'AtomGroup for making selections')
            sel = self._ag.select(selstr)
            if sel is not None:
                self._indices = sel.getIndices()
            self._sel = sel
            return sel
    
    def getCoordinates(self):
        """Return a copy of reference coordinates of selected atoms."""
        
        if self._coords is None:
            return None
        if self._sel is None:
            return self._coords.copy()
        return self._coords[self._indices].copy()

    def setCoordinates(self, coords):
        """Set reference coordinates."""

        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoordinates()
            except AttributeError:
                raise TypeError('coords must be an ndarray instance or '
                                'must contain getCoordinates as an attribute')

        elif coords.ndim != 2:
            raise ValueError('coordinates must be a 2d array')
        elif coords.shape[1] != 3:
            raise ValueError('shape of coordinates must be (n_atoms,3)')
        elif not coords.dtype in (np.float32, np.float64):
            try:
                coords = coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        
        if self._n_atoms == 0:
            self._n_atoms = coords.shape[0]    
        elif coords.shape[0] != self._n_atoms:
            raise ValueError('length of coords does not match number of atoms')
        self._coords = coords
        
    def getWeights(self):
        """Return a copy of weights of selected atoms."""
        
        if self._weights is None:
            return None
        if self._sel is None:
            return self._weights.copy()
        if self.weights.ndim == 2:
            return self._weights[self._indices].copy()
        else:
            return self._weights[:, self._indices].copy()
    
    def setWeights(self, weights):
        """Set atomic weights."""
        
        if self._n_atoms == 0:
            raise AttributeError('first set reference coordinates')
        elif not isinstance(weights, np.ndarray):
            raise TypeError('weights must be an ndarray instance')
        elif not weights.ndim in (1, 2): 
            raise ValueError('weights.ndim must be equal to 1 or 2')        
        elif len(weights) != self._n_atoms:
            raise ValueError('length of weights must be equal to number of '
                             'atoms')
        if weights.dtype in (np.float32, np.float64):
            try:
                weights = weights.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
            
        if weights.ndim == 1:
            weights = weights.reshape((self._n_atoms, 1))
        self._weights = weights


class Ensemble(EnsembleBase):
    
    """A class for analysis of arbitrary conformational ensembles. 
    
    Indexing returns a :class:`Conformation` instance, which corresponds to
    a coordinate set. Slicing returns an :class:`Ensemble` instance that 
    contains subset of conformations (coordinate sets). The ensemble obtained 
    by slicing will have a copy of the reference coordinates.

    """

    def __init__(self, name='Unnamed'):
        """Instantiate with a name.
        
        .. versionchanged:: 0.6
           At instantiation, :class:`~prody.atomic.Atomic` instances are 
           accepted as *name* argument. All coordinate sets from *name* will be
           added to the ensemble automatically. 
           
        .. versionchanged:: 0.7
           When an empty string is passed as *name* argument, Ensemble instance
           is called "Unnamed". 
          
        :arg name: A name (:class:`str`) or an :class:`~prody.atomic.Atomic`
            instance.
        
        """
        EnsembleBase.__init__(self, name)
        self._confs = None       # coordinate sets
        
        if isinstance(name, prody.Atomic):
            self.setCoordinates(name.getCoordinates())
            self.addCoordset(name)
        
    def __getitem__(self, index):
        """Return a conformation at given index."""
        if self._confs is None:
            return None
        if isinstance(index, int):
            return self.getConformation(index) 
        elif isinstance(index, slice):
            ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                self._name, index.indices(len(self))))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self.getCoordsets(index))
            ens.setWeights(self.getWeights())
            return ens
        elif isinstance(index, (list, np.ndarray)):
            ens = Ensemble('Conformations of {0:s}'.format(self._name))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self.getCoordsets(index))
            ens.setWeights(self.getWeights())
            return ens
        else:
            raise IndexError('invalid index')
            
    
    def __add__(self, other):
        """Return ensemble that contains conformations in *self* and *other*.
        
        The reference coordinates of *self* is used in the result.
        
        """
        raise NotImplemented('addition for ensembles is not implemented')
        
        if not isinstance(other, Ensemble):
            raise TypeError('an Ensemble instance cannot be added to an {0:s} '
                            'instance'.format(type(other)))
        elif self.getNumOfAtoms() != other.getNumOfAtoms():
            raise ValueError('Ensembles must have same number of atoms.')
    
        ensemble = Ensemble('{0:s} + {1:s}'.format(self.getName(), 
                                                   other.getName()))
        return None
    
    def __iter__(self):
        n_confs = self._n_confs
        for i in range(n_confs):
            if n_confs != self._n_confs:
                raise RuntimeError('number of conformations in the ensemble '
                                   'changed during iteration')
            yield Conformation(self, i)
    
    def __repr__(self):
        return ('<Ensemble: {0:s} ({1:d} conformations, {2:d} atoms, {3:d} '
               'selected)>').format(self._name, len(self), self._n_atoms, 
                                    self.getNumOfSelected())

    def __str__(self):
        return 'Ensemble {0:s}'.format(self._name)

    def getNumOfConfs(self):
        """Return number of conformations."""

        return self._n_confs

    def addCoordset(self, coords, allcoordsets=True):
        """Add coordinate set(s) as conformation(s).
        
        :class:`~prody.atomic.Atomic` instances are accepted as *coords* 
        argument. If *allcoordsets* is ``True``, all coordinate sets from
        the :class:`~prody.atomic.Atomic` instance will be appended to the 
        ensemble. Otherwise, only the active coordinate set will be appended.
        
        """
        
        if not isinstance(coords, np.ndarray):
            if isinstance(coords, prody.Atomic):
                if allcoordsets:
                    coords = coords.getCoordsets()
                else:
                    coords = coords.getCoordinates()
            else:
                raise TypeError('coords must be a numpy ndarray or '
                                'prody Atomic instance')
            
        if not coords.ndim in (2, 3):
            raise ValueError('coords must be a 2d or a 3d array')
        elif not coords.dtype in (np.float32, np.float64):
            try:
                coords = coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        
        if self._n_atoms == 0:
            self._n_atoms = coords.shape[-2]
        else:
            if coords.shape[-2] != self._n_atoms:
                raise ValueError('shape of conf must be (n_atoms,3)')
        n_atoms = self._n_atoms

        if coords.ndim == 2:
            coords = coords.reshape((1, self._n_atoms, 3))
        n_confs = coords.shape[0]
            
        if self._confs is None: 
            self._confs = coords
        else:
            self._confs = np.concatenate((self._confs, coords), axis=0)
        self._n_confs += n_confs

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate sets at given indices.
        
        *indices* may be an integer, a list of integers or ``None``. ``None``
        returns all coordinate sets. 
    
        For reference coordinates, use getCoordinates method.

        """
        
        if self._confs is None:
            return None
        if indices is None:
            indices = slice(None)
        elif isinstance(indices, (int, long)): 
            indices = np.array([indices])
        if self._sel is None:
            coords = self._confs[indices].copy()
            if self._weights is not None and self._weights.ndim == 3:
                for i, w in enumerate(self._weights[indices]):
                    which = w.flatten()==0
                    coords[i, which] = self._coords[which]
        else:
            selids = self._indices
            coords = self._confs[indices, selids].copy()
            if self._weights is not None and self._weights.ndim == 3:
                ref_coords = self._coords[selids]
                for i, w in enumerate(self._weights[indices, selids]):
                    which = w.flatten()==0
                    coords[i, which] = ref_coords[which]
        return coords
    
    def delCoordset(self, index):
        """Delete a coordinate set from the ensemble."""
        
        if isinstance(index, int):
            index = [index]
        else:
            index = list(index)
        length = self._n_confs
        which = np.ones(length, np.bool)
        which[index] = False
        if which.sum() == 0:
            self._confs = None
            self._weights = None
        else:
            self._confs = self._confs[which]
            if self._weights is not None:
                self._weights = self._weights[which]
        self._n_confs -= len(index)

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
        n_confs = self._n_confs
        if -n_confs <= index < n_confs:
            if index < 0:
                index = n_confs - index
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
        LOGGER.info('Superimposing structures.')
        start = time()
        self._superpose()
        LOGGER.info(('Superposition is completed in {0:.2f} '
                           'seconds.').format((time() - start)))
        
    def _superpose(self):
        """Superpose conformations and return new coordinates."""
        
        if self._sel is None:
            weights = self._weights
            coords = self._coords
            confs = self._confs
            for i, conf in enumerate(confs):
                confs[i] = measure._calcTransformation(conf, coords, weights
                                                            ).apply(confs[i]) 
        else:            
            indices = self._indices
            weights = None
            if self._weights is not None:
                weights = self._weights[indices]
            coords = self._coords[indices]
            confs_selected = self._confs[:,indices]
            confs = self._confs
            for i, conf in enumerate(confs_selected):
                confs[i] = measure._calcTransformation(conf, coords, weights
                                                            ).apply(confs[i]) 
            
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
            rmsdif = measure._calcRMSD(self._coords, newxyz)
            self._coords = newxyz
            step += 1
            LOGGER.info(('Step #{0:d}: RMSD difference = '
                               '{1:.4e}').format(step, rmsdif))
        LOGGER.info('Iterative superposition completed in {0:.2f}s.'
                    .format((time() - start)))
#left here        
    def getMSF(self):
        """Calculate and return Mean-Square-Fluctuations."""
        
        if self._confs is None: 
            return
        
        xyzmean = (self._confs * 
                   self._weights).sum(0) / self._weights.sum(0)
        xyzdiff2 = np.power(self._confs - xyzmean, 2).sum(2)
        weightsum = self._weights.sum(2)
        msf = (xyzdiff2 * weightsum).sum(0) / weightsum.sum(0)  
        return msf
            
    def getDeviations(self):
        """Return deviations from reference coordinates."""
        
        if not isinstance(self._confs, np.ndarray):
            LOGGER.warning('Conformations are not set.')
            return None
        if not isinstance(self._coords, np.ndarray):
            LOGGER.warning('Coordinates are not set.')
            return None
        
        return self.getCoordsets() - self._coords 
        
    def getRMSDs(self):
        """Calculate and return Root Mean Square Deviations."""
        
        if self._confs is None or self._coords is None: 
            return None
        return measure._calcRMSD(self._coords, self._confs, self._weights)

class PDBEnsemble(Ensemble):
    
    """This class enables handling coordinates for heterogeneous structural 
    datasets and stores identifiers for individual conformations.
    
    |example| See usage usage in :ref:`pca-xray`, :ref:`pca-dimer`, and 
    :ref:`pca-blast`.
    
    .. versionadded:: 0.8
    
    .. note:: This class is designed to handle conformations with missing
       coordinates, e.g. atoms that are note resolved in an X-ray structure.
       For unresolved atoms, the coordinates of the reference structure is
       assumed in RMSD calculations and superpositions.
    
    >>> ensemble
    <PDB Ensemble: p38 X-ray (75 conformations, 321 atoms, 321 selected)>

    """

    def __init__(self, name='Unnamed'):
        self._identifiers = []
        Ensemble.__init__(self, name)
        
    def __repr__(self):
        return '<PDB ' + Ensemble.__repr__(self)[1:]
    
    def __str__(self):
        return 'PDB ' + Ensemble.__str__(self)
    
    def __iter__(self):
        n_confs = self._n_confs
        for i in range(n_confs):
            if n_confs != self._n_confs:
                raise RuntimeError('number of conformations in the ensemble '
                                   'changed during iteration')
            yield PDBConformation(self, i)
    
    def __getitem__(self, index):
        """Return a conformation at given index."""
        
        if isinstance(index, int):
            return self.getConformation(index) 
        elif isinstance(index, slice):
            ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                self._name, index.indices(len(self))))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self._confs[index].copy(), 
                            self._weights[index].copy())
            return ens
        elif isinstance(index, (list, np.ndarray)):
            ens = Ensemble('Conformations of {0:s}'.format(self._name))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self._confs[index].copy(), 
                            self._weights[index].copy())
            return ens
        else:
            raise IndexError('invalid index')
            
    def _superpose(self):
        """Superpose conformations and return new coordinates."""

        if self._sel is None:
            weights = self._weights
            coords = self._coords
            confs = self._confs
            for i, conf in enumerate(confs):
                confs[i] = measure._calcTransformation(conf, coords, weights[i]
                                                            ).apply(confs[i])
        else:            
            indices = self._getSelIndices()
            weights = self._weights[:, indices]
            coords = self._coords[indices]
            confs_selected = self._confs[:,indices]
            confs = self._confs
            for i, conf in enumerate(confs_selected):
                confs[i] = measure._calcTransformation(conf, coords, weights[i]
                                                            ).apply(confs[i]) 

    def addCoordset(self, coords, weights=None, allcoordsets=True):
        """Add coordinate set(s) as conformation(s).
        
        :class:`~prody.atomic.Atomic` instances are accepted as *coords* 
        argument. If *allcoordsets* is ``True``, all coordinate sets from
        the :class:`~prody.atomic.Atomic` instance will be appended to the 
        ensemble. Otherwise, only the active coordinate set will be appended.

        
        *weights* is an optional argument. If provided, its length must
        match number of atoms.
        
        """
        ag = None
        if isinstance(coords, prody.AtomGroup):
            name = coords.getName()
            ag = coords
        elif isinstance(coords, np.ndarray):
            name = 'Unnamed'
        else:
            name = str(coords)
        n_confs = self._n_confs
        n_atoms = self._n_atoms
        
        if weights is not None:
            if not isinstance(weights, np.ndarray):
                raise TypeError('weights must be an ndarray')
            elif not weights.ndim in (1, 2, 3):
                raise ValueError('weights must be a 1d, 2d, or 3d array')
            elif weights.ndim in (2, 3) and weights.shape[-1] != 1:
                raise ValueError('shape of weights must be '
                                 '([n_coordsets,] number_of_atoms, 1)')
            elif weights.dtype != np.float64:
                try:
                    weights = weights.astype(np.float64)
                except ValueError:
                    raise ValueError('weights array cannot be assigned type '
                                     '{0:s}'.format(np.float64))
            
            if weights.ndim < 3:
                weights = weights.reshape((1, n_atoms, 1))
                
        Ensemble.addCoordset(self, coords, allcoordsets)
        diff = self._n_confs - n_confs

        if weights.shape[0] != diff:  
            if weights.shape[0] == 1:
                LOGGER.warning('Shape of coords and weights did not match. '
                       'First set of weights are used for all coordinate set.')
                weights = weights[0]          
            weights = np.tile(weights, (diff, 1, 1))
                
        while '  ' in name:
            name = name.replace('  ', ' ')
        name = name.replace(' ', '_')
        if diff > 1:
            self._identifiers += ['{0:s}_{1:d}'
                                  .format(name, i+1) for i in range(diff)]
        else:
            if ag is not None and ag.getNumOfCoordsets() > 0:
                self._identifiers.append('{0:s}_{1:d}'.format(name, 
                                                ag.getActiveCoordsetIndex()))
            else:                
                self._identifiers.append(name)
        if self._weights is None: 
            self._weights = weights
        else:
            if weights is not None:
                self._weights = np.concatenate((self._weights, weights),
                                               axis=0)
            else:
                if self._weights is not None:
                    self._weights = np.concatenate((self._weights, 
                                    np.ones((diff, n_atoms, 1))), axis=0)

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate sets at given *indices* for selected 
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
        if self._sel is None:
            coords = self._confs[indices].copy()
            for i, w in enumerate(self._weights[indices]):
                which = w.flatten()==0
                coords[i, which] = self._coords[which]
            return coords 
        else:
            selids = self._getSelIndices()
            coords = self._confs[indices, selids].copy()
            for i, w in enumerate(self._weights[indices]):
                which = w[selids].flatten()==0
                coords[i, which] = self._coords[selids][which]
            return coords 
    
    def iterCoordsets(self):
        """Iterate over coordinate sets. A copy of each coordinate set for
        selected atoms is returned. Reference coordinates are not included."""
        
        conf = PDBConformation(self, 0)
        for i in range(self._n_confs):
            conf._index = i
            yield conf.getCoordinates()
   
    def delCoordset(self, index):
        """Delete a coordinate set from the ensemble."""
        
        Ensemble.delCoordset(self, index)
        if isinstance(index, int):
            index = [index]
        else:
            index = list(index)
        index.sort(reverse=True)
        for i in index:
            self._identifiers.pop(i)
            
    def getConformation(self, index):
        """Return conformation at given index."""

        if self._confs is None:
            raise AttributeError('conformations are not set')
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        n_confs = self._n_confs
        if -n_confs <= index < n_confs:
            if index < 0:
                index = n_confs - index
            return PDBConformation(self, index)
        else:
            raise IndexError('conformation index out of range')
            
    def getRMSDs(self):
        """Calculate and return Root Mean Square Deviations.
        
        >>> rmsd = ensemble.getRMSDs().round(2)
        >>> print rmsd[0]
        0.74
        >>> print( rmsd )
        [ 0.74  0.53  0.58  0.6   0.61  0.72  0.62  0.74  0.69  0.65  0.48  0.54
          0.48  0.75  0.56  0.76  0.84  0.49  0.69  0.74  0.69  0.49  0.69  0.61
          0.73  0.64  0.52  0.51  0.65  0.84  0.77  0.69  0.82  1.19  0.6   1.12
          0.71  0.61  0.73  0.57  0.99  0.94  0.85  0.98  0.58  0.61  0.64  0.89
          0.82  0.95  0.88  0.86  1.09  0.7   0.72  0.86  0.76  0.82  0.88  0.95
          0.63  0.92  1.08  0.44  0.43  0.49  0.64  0.88  0.72  0.9   0.96  1.23
          0.58  0.66  0.83]
          """
        
        if self._confs is None: 
            return None
        weights = self._weights
        wsum_axis_2 = weights.sum(2)
        wsum_axis_1 = wsum_axis_2.sum(1)
        rmsd = np.sqrt((np.power(self.getDeviations(), 2).sum(2) * 
                        wsum_axis_2).sum(1) / wsum_axis_1)
        return rmsd

    def setWeights(self, weights):
        """Set atomic weights."""
        
        if self._n_atoms == 0:
            raise AttributeError('coordinates are not set')
        elif not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape[:2] != (self._n_confs, self._n_atoms):
            raise ValueError('shape of weights must (n_confs, n_atoms[, 1])')
        if weights.dtype != np.float64:
            try:
                weights = weights.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
            
        if weights.ndim == 2:
            weights = weights.reshape((self._n_confs, self._n_atoms, 1))
        self._weights = weights


def saveEnsemble(ensemble, filename=None):
    """Save *ensemble* model data as :file:`filename.ens.npz`. 
    
    If *filename* is ``None``, name of the *ensemble* will be used as 
    the filename, after " " (blank spaces) in the name are replaced with "_" 
    (underscores).  
    
    Extension is :file:`.ens.npz`.
    
    Upon successful completion of saving, filename is returned.
    
    This function makes use of :func:`numpy.savez` function.
    
    """
    
    if not isinstance(ensemble, Ensemble):
        raise TypeError('invalid type for ensemble, {0:s}'
                        .format(type(ensemble)))
    if len(ensemble) == 0:
        raise ValueError('ensemble instance does not contain data')
    
    dict_ = ensemble.__dict__
    attr_list = ['_name', '_confs', '_weights', '_coords', '_n_atoms', 
                 '_n_confs']
    if isinstance(ensemble, PDBEnsemble):
        attr_list.append('_identifiers')
    if filename is None:
        filename = ensemble.getName().replace(' ', '_')
    attr_dict = {}
    for attr in attr_list:
        value = dict_[attr]
        if value is not None:
            attr_dict[attr] = value
    filename += '.ens.npz'
    np.savez(filename, **attr_dict)
    return filename

def loadEnsemble(filename):
    """Return ensemble instance after loading it from file (*filename*).
    
    .. seealso: :func:`saveEnsemble`
    
    This function makes use of :func:`numpy.load` function.
    
    """
    
    attr_dict = np.load(filename)
    if '_identifiers' in attr_dict.files:
        ensemble = PDBEnsemble(str(attr_dict['_name']))
    else:
        ensemble = Ensemble(str(attr_dict['_name']))
    dict_ = ensemble.__dict__
    for attr in attr_dict.files:
        if attr == '_name': 
            continue
        elif attr in ('_n_confs', '_n_atoms'):
            dict_[attr] = int(attr_dict[attr])
        else:
            dict_[attr] = attr_dict[attr]
    if isinstance(ensemble, PDBEnsemble):
        ensemble._ensemble = [None] * len(ensemble)
    return ensemble

class ConformationBase(object):

    """Base class for :class:`Conformation` and :class:`Frame`."""

    __slots__ = ['_ensemble', '_index', '_sel', '_indices']

    def __init__(self, ensemble, index):
        self._ensemble = ensemble
        self._index = index
        self._sel = ensemble.getSelection()
        self._indices = ensemble._getSelIndices()
        
    def __repr__(self):
        return ('<{0:s}: {1:d} from {2:s} (selected {3:d} of {4:d} atoms)>'
               ).format(self.__class__.__name__, self._index, 
                        self._ensemble.getName(), self.getNumOfSelected(), 
                        self._ensemble.getNumOfAtoms())

    def __str__(self):
        return '{0:s} {1:d} from {2:s}'.format(
                    self._index, self._ensemble.getName())

    def getNumOfAtoms(self):
        """Return number of atoms."""
        
        return self._ensemble.getNumOfAtoms()
    
    def getNumOfSelected(self):
        """Return number of selected atoms."""
        
        if self._sel is None:
            return self._ensemble.getNumOfAtoms()
        else:
            return len(self._indices)
    
    def getIndex(self):
        """Return index."""
        
        return self._index
    
    def getWeights(self):
        """Return coordinate weights for selected atoms."""
        
        if self._sel is None:
            return self._ensemble.getWeights()
        else:
            return self._ensemble.getWeights()[self._indices]


class Conformation(ConformationBase):
    
    """A class to provide methods on a conformation in an ensemble.  
    Instances of this class do not keep coordinate and weights data.
    
    """
    
    def __init__(self, ensemble, index):
        """Instantiate with an ensemble instance, and an index."""
        ConformationBase.__init__(self, ensemble, index)
        
    def getEnsemble(self):
        """Return the ensemble that this conformation belongs to."""
        
        return self._ensemble
    
    def getCoordinates(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned.
        
        .. warning:: When there are atoms with weights equal to zero (0),
           their coordinates will be replaced with the coordinates of the
           ensemble reference coordinate set.

        """
        
        if self._ensemble._confs is None:
            return None
        indices = self._indices
        if indices is None:
            return self._ensemble._confs[self._index].copy()
        return self._ensemble._confs[self._index, indices].copy()
    
    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates. 
        Deviations are calculated for selected atoms."""
        
        ensemble = self._ensemble 
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        if indices is None:
            return ensemble._coords - ensemble._confs[self._index].copy()
        else:
            return (ensemble._confs[self._index, indices].copy() - 
                    ensemble._coords[indices])

    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates. RMSD is
        calculated for selected atoms."""
        
        ensemble = self._ensemble 
        indices = self._indices
        if indices is None:
            return measure._calcRMSD(ensemble._coords,
                                     ensemble._confs[self._index], 
                                     self.getWeights())[0]
        else:
            return measure._calcRMSD(ensemble._coords[indices],
                                     ensemble._confs[self._index, indices], 
                                     self.getWeights())[0]
    
class PDBConformation(Conformation):
    
    """This class is the same as :class:`Conformation`, except that the 
    conformation has a name (or identifier), e.g. PDB identifier.
    
    .. versionadded:: 0.8

    >>> conf = ensemble[0] 
    >>> conf
    <PDB Conformation: AtomMap_Chain_A_from_1a9u_->_Chain_A_from_p38_reference from p38 X-ray (index: 0; selected 321 of 321 atoms)>

    """
    
    def __repr__(self):
        return ('<PDB Conformation: {0:s} from {1:s} (index: {2:d}; '
                'selected {3:d} of {4:d} atoms)>').format(
                    self._ensemble._identifiers[self._index], 
                    self._ensemble.getName(), self._index, 
                    self.getNumOfSelected(), self.getNumOfAtoms())
    
    def __str__(self):
        return 'PDB Conformation {0:s} from {1:s}'.format(
                    self._ensemble._identifiers[self._index], 
                    self._ensemble.getName())
    
    def getIdentifier(self):
        """Return the identifier of the conformation instance.
        
        >>> print( conf.getIdentifier() )
        AtomMap_Chain_A_from_1a9u_->_Chain_A_from_p38_reference
        
        """
        
        return self._ensemble._identifiers[self._index]
    
    def setIdentifier(self, identifier):
        """Set the identifier of the conformation instance.
        
        >>> conf.setIdentifier('1a9u')
        >>> print conf.getIdentifier()
        1a9u
        
        """
        
        self._ensemble._identifiers[self._index] = str(identifier)
        
    def getCoordinates(self):
        """Return a copy of the coordinates of the conformation. If a subset
        of atoms are selected in the ensemble, coordinates for selected
        atoms will be returned.
        
        .. warning:: When there are atoms with weights equal to zero (0),
           their coordinates will be replaced with the coordinates of the
           ensemble reference coordinate set.

        """
        
        ensemble = self._ensemble
        if ensemble._confs is None:
            return None
        indices = ensemble._indices
        index = self._index
        if indices is None:
            coords = ensemble._confs[index].copy()
            weights = ensemble._weights[index].flatten()
            which = weights==0
            coords[which] = ensemble._coords[which]
        else: 
            coords = ensemble._confs[index, indices].copy()
            weights = ensemble._weights[index, indices].flatten()
            which = weights==0
            coords[which] = ensemble._coords[indices][which]
        return coords
    
    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates. 
        Deviations are calculated for selected atoms."""
        
        indices = self._indices
        if indices is None:
            return self.getCoordinates() - self._ensemble._coords
        else:
            return self.getCoordinates() - self._ensemble._coords[indices]
    
    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates. RMSD is
        calculated for selected atoms.
        
        >>> print( conf.getRMSD().round(2) )
        0.74

        """
        
        ensemble = self._ensemble
        index = self._index
        indices = ensemble._indices
        if indices is None:
            return measure._calcRMSD(ensemble._coords,
                                     ensemble._confs[index], 
                                     ensemble._weights[index])
        else:
            return measure._calcRMSD(ensemble._coords[indices],
                                     ensemble._confs[index, indices], 
                                     ensemble._weights[index, indices])
    
    def getWeights(self):
        """Return coordinate weights for selected atoms."""
        
        ensemble = self._ensemble
        indices = ensemble._indices
        if indices is None:
            return ensemble._weights[self._index]
        else:
            return ensemble._weights[self._index, indices]
    
class Frame(ConformationBase):
    
    """A class to provide methods on a frame in a trajectory.
    
    """
    
    __slots__ = ['_ensemble', '_index', '_coords', '_indices', '_sel']

    def __init__(self, trajectory, index, coords):
        """Instantiate with an trajectory instance, index, and a name."""
        ConformationBase.__init__(self, trajectory, index)
        self._coords = coords
        
    def getTrajectory(self):
        """Return the trajectory that this frame belongs to."""
        
        return self._ensemble
    
    def getCoordinates(self):
        """Return coordinates for selected atoms."""
        
        return self._coords.copy()
    
    def getDeviations(self):
        """Return deviations from the trajectory reference coordinates."""

        indices = self._indices 
        if indices is None:
            return self._coords - ensemble._coords
        else:
            return self._coords - ensemble._coords[indices]

    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates."""
        
        indices = self._indices 
        ensemble = self._ensemble
        ensemble = self._ensemble
        if indices is None:
            return measure._calcRMSD(self._coords, ensemble._coords, 
                                                    ensemble._weights)
        else:
            if ensemble._weights is None:
                return measure._calcRMSD(self._coords[indices], 
                                         ensemble._coords[indices])
            else:
                return measure._calcRMSD(self._coords[indices], 
                                         ensemble._coords[indices], 
                                         ensemble._weights[indices])
        

    def superpose(self):
        """Superpose frame coordinates onto the trajectory reference 
        coordinates."""
        
        indices = self._indices 
        ensemble = self._ensemble
        if indices is None:
            self._coords, t = measure.superpose(self._coords, 
                                                ensemble._coords, 
                                                ensemble._weights)
        else:
            if ensemble._weights is None:
                self._coords = measure.superpose(self._coords, 
                                                 ensemble._coords[indices])
            else:
                self._coords = measure.superpose(self._coords, 
                                                 ensemble._coords[indices], 
                                                 ensemble._weights[indices])
    
def trimEnsemble(ensemble, **kwargs):
    """Return a PDB ensemble obtained by trimming given *ensemble*.
    
    .. versionadded:: 0.5.3
    
    This function helps selecting atoms in an ensemble based on one of the 
    following criteria, and returns them in a new :class:`PDBEnsemble` instance.
        
    **Occupancy**
    
    If weights in the ensemble correspond to atomic occupancies, then this 
    option may be used to select atoms with average occupancy values higher
    than the argument *occupancy*. Function will calculate average occupancy 
    by dividing the sum of weights for each atom to the number of conformations
    in the ensemble.
    
    :arg occupancy: average occupancy for selection atoms in an ensemble, 
                    must be 0 < *occupancy* <= 1
    :type occupancy: float
    
    """
    
    if not isinstance(ensemble, PDBEnsemble):
        raise TypeError('trimEnsemble() argument must be an Ensemble instance')
    if ensemble.getNumOfConfs() == 0 and ensemble.getNumOfAtoms() == 0:
        raise ValueError('coordinates or conformations must be set for '
                         'ensemble')
    
    if 'occupancy' in kwargs:
        occupancy = float(kwargs['occupancy'])
        assert 0 < occupancy <=1, ('occupancy is not > 0 and <= 1: '
                                   '{0:s}'.format(repr(occupancy)))
        n_confs = ensemble.getNumOfConfs()
        assert n_confs > 0, 'ensemble does not contain anyconformations'
        weights = calcSumOfWeights(ensemble)
        assert weights is not None, 'weights must be set for ensemble'
        weights = weights.flatten()
        mean_weights = weights / n_confs
        
        torf = mean_weights >= occupancy
        
    else:
        return None
    
    trimmed = PDBEnsemble(ensemble.getName())
    coords = ensemble.getCoordinates()
    if coords is not None:
        trimmed.setCoordinates( coords[torf] )
    confs = ensemble.getCoordsets()
    if confs is not None:
        weights = ensemble.getWeights()
        trimmed.addCoordset( confs[:, torf], weights[:, torf] )
    return trimmed

def calcSumOfWeights(ensemble):
    """Return sum of weights from a PDB ensemble.
    
    Weights are summed for each atom over conformations in the ensemble.
    Size of the plotted array will be equal to the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    
    >>> calcSumOfWeights(ensemble) # doctest: +ELLIPSIS
    array([ 74.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  73.,  73.,
            74.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,
            ...
            75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,  75.,
            75.,  75.])
            
    Each number in the above example corresponds to a residue (or atoms) and 
    shows the number of structures in which the corresponding residue is 
    resolved. 
    
    """
    
    if not isinstance(ensemble, PDBEnsemble):
        raise TypeError('ensemble must be an Ensemble instance')
    
    weights = ensemble.getWeights()
    
    if weights is None:
        return None
    
    return weights.sum(0).flatten()
    
    
def showSumOfWeights(ensemble, *args, **kwargs):
    """Show sum of weights for a PDB ensemble using 
    :func:`~matplotlib.pyplot.plot`.
    
    Weights are summed for each atom over conformations in the ensemble.
    Size of the plotted array will be equal to the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    
    """
    """
    *indices*, if given, will be used as X values. Otherwise, X axis will
    start from 0 and increase by 1 for each atom. 
    
    """
    
    if plt is None: prody.importPyPlot()
    if not plt: return None
    if not isinstance(ensemble, PDBEnsemble):
        raise TypeError('ensemble must be an Ensemble instance')
    weights = calcSumOfWeights(ensemble)
    if weights is None:
        return None
    show = plt.plot(weights, *args, **kwargs)
    axis = list(plt.axis())
    axis[2] = 0
    axis[3] += 1
    plt.axis(axis)
    plt.xlabel('Atom index')
    plt.ylabel('Sum of weights')
    return show

# Translation of dcdplugin.c
from struct import calcsize, unpack

RECSCALE32BIT = 1
RECSCALE64BIT = 2

class TrajectoryBase(EnsembleBase):
    
    def __init__(self, name):
        EnsembleBase.__init__(self, name)
        self._nfi = 0
        self._closed = False
    
    def __iter__(self):
        if not self._closed:
            while self._nfi < self._n_confs: 
                yield self.next()
    
    def getNumOfFrames(self):
        """Return number of frames."""
        
        return self._n_confs
    
    def getNextFrameIndex(self):
        """Return the index of the next frame."""
        
        return self._nfi
    
    def iterCoordsets(self):
        """Iterate over coordinate sets. Coordinates for selected atoms are 
        returned. Reference coordinates are not included. Iteration starts
        from the current position in the trajectory."""

        if not self._closed:
            while self._nfi < self._n_confs: 
                yield self.nextCoordset()
    
    def getFrame(self):
        pass


    def next(self):
        pass

    def goto(self, n):
        pass

    def skip(self, n=1):
        pass

    def reset(self):
        pass

    def close(self):
        pass


class TrajectoryFile(TrajectoryBase):
    
    """A base class for providing a consistent interface over trajectories in 
    different formats. 
    
    .. versionadded:: 0.8
    
    """
    
    
    def __init__(self, filename):
        """Instantiate with a trajectory *filename*. The file will be 
        automatically opened for reading at instantiation. 
        Use :meth:`addFile` for appending files to the trajectory."""
        
        if not os.path.isfile(filename):
            raise IOError("[Errno 2] No such file or directory: '{0:s}'"
                          .format(filename))
        self._filename = filename
        self._file = open(filename, 'rb')
        name = os.path.splitext(os.path.split(filename)[1])[0]
        TrajectoryBase.__init__(self, name)
        self._bytes_per_frame = None
        self._first_byte = None
        self._dtype = np.float32

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.getFrame(index)
        elif isinstance(index, (slice, list, np.ndarray)):
            if isinstance(index, slice):
                ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                    self._name, index.indices(len(self))))
            else:
                ens = Ensemble('{0:s} slice'.format(self._name))
            ens.setCoordinates(self.getCoordinates())
            if self._weights is not None:
                ens.setWeights(self._weights.copy())
            ens.addCoordset(self.getCoordsets(index))
            return ens
        else:
            raise IndexError('invalid index')
    
    def __repr__(self):
        if self._closed:
            return ('<{0:s}: {1:s} (closed)>').format(
                        self.__class__.__name__, self._name)
        else:
            return ('<{0:s}: {1:s} (at {2:d} of {3:d} frames, '
                    'selected {4:d} of {5:d} atoms)>').format(
                    self.__class__.__name__, self._name, 
                    self._nfi, self._n_confs, self.getNumOfSelected(),
                    self._n_atoms)
    
    def __str__(self):
        return '{0:s} {1:s}'.format(self.__class__.__name__, self._name)
    
       
    def getFilename(self, absolute=False):
        """Return relative path to the current file. For absolute path,
        pass ``absolute=True`` argument."""
        
        if absolute:
            return os.path.abspath(self._filename)    
        return os.path.relpath(self._filename)
    
    def getFrame(self, index):
        """Return frame at given *index*."""
        
        if closed:
            return None
        if not isinstance(index, (int, long)):
            raise IndexError('index must be an integer')
        if not 0 <= index < self._n_confs:
            raise IndexError('index must be greater or equal to 0 and less '
                             'than number of frames')
        nfi = self._nfi
        if index > nfi:
            self.skip(index - nfi)
        elif index < nfi:
            self.reset()
            if index > 0:
                self.skip(index)
        return self.next()
    
            
    def getCoordsets(self, indices=None):
        """Returns coordinate sets at given *indices*. *indices* may be an 
        integer, a list of ordered integers or ``None``. ``None`` returns all 
        coordinate sets."""
        
        if indices is None:
            indices = np.arange(self._n_confs)
        elif isinstance(indices, (int, long)):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(*indices.indices(self._n_confs))
        elif isinstance(indices, (list, np.ndarray)):
            if isinstance(indices, list):
                indices = np.array(indices)
            indices.sort()
        else:
            raise TypeError('indices must be an integer or a list of integers')

        nfi = self._nfi
        self.reset()

        n_atoms = self.getNumOfSelected() 
        coords = np.zeros((len(indices), n_atoms, 3), self._dtype)

        prev = 0
        next = self.nextCoordset
        for i, index in enumerate(indices):
            diff = index - prev
            if diff > 0:
                self.skip(diff)
            coords[i] = next()

        self.goto(nfi)
        return coords
    
    def skip(self, n=1):
        """Skip *n* frames. *n* must be a positive integer, default is 1."""
        
        if not isinstance(n, (int, long)):
            raise ValueError('n must be an integer')
        if not self._closed and n > 0:
            left = self._n_confs - self._nfi
            if n > left:
                n = left
            self._file.seek(n * self._bytes_per_frame, 1)                
            self._nfi += n
    
    def goto(self, n):
        """Go to the frame at index *n*. ``n=0`` will rewind the trajectory
        to the beginning. ``n=-1`` will go to the last frame."""
        
        if self._closed:
            return None
        if not isinstance(n, (int, long)):
            raise ValueError('n must be an integer')
        n_confs = self._n_confs
        if n > 0:
            if n > n_confs: 
                n = n_confs
            self._file.seek(self._first_byte + n * self._bytes_per_frame)
            self._nfi = n
        elif n < 0:
            n = n_confs + n
            if n < 0:
                n = 0
            self._file.seek(self._first_byte + n * self._bytes_per_frame)
            self._nfi = n
        else:
            self.reset()
        
    
    def reset(self):
        """Go to first frame whose index is 0."""

        if not self._closed:
            self._file.seek(self._first_byte)
            self._nfi = 0

    
    def close(self):
        """Close trajectory file."""
        
        self._file.close()
        self._nfi = 0
        self._closed = True
    
    
class DCDFile(TrajectoryFile):
    
    """A class for reading DCD files."""
    
    def __init__(self, filename):
        """Instantiate with a DCD filename. DCD header and first frame is 
        parsed at instantiation. Coordinates from the first frame is set 
        as the reference coordinates."""
        
        TrajectoryFile.__init__(self, filename)
        self._parseHeader()
        
    def _parseHeader(self):
        """Read the header information from a dcd file.
        Input: fd - a file struct opened for binary reading.
        Output: 0 on success, negative error code on failure.
        Side effects: *natoms set to number of atoms per frame
                      *nsets set to number of frames in dcd file
                      *istart set to starting timestep of dcd file
                      *nsavc set to timesteps between dcd saves
                      *delta set to value of trajectory timestep
                      *nfixed set to number of fixed atoms 
                      *freeind may be set to heap-allocated space
                      *reverse set to one if reverse-endian, zero if not.
                      *charmm set to internal code for handling charmm data.
        """
        
        dcd = self._file
        endian = '' #'=' # native endian
        rec_scale = RECSCALE32BIT
        charmm = None
        dcdcordmagic = unpack(endian+'i', 'CORD')[0]
        # Check magic number in file header and determine byte order
        bits = dcd.read(calcsize('ii'))
        temp = unpack(endian+'ii', bits)
        if temp[0] + temp[1] == 84:
            LOGGER.info('Detected CHARMM -i8 64-bit DCD file of native '
                        'endianness.')
            rec_scale = RECSCALE64BIT
        elif temp[0] == 84 and temp[1] == dcdcordmagic:
            LOGGER.info('Detected standard 32-bit DCD file of native '
                        'endianness.')
        else:
            if unpack('>ii', bits) == temp:
                endian = '>'
            else:
                endian = '<'
            temp = unpack(endian+'ii', bits)
            if temp[0] + temp[1] == 84:
                rec_scale = RECSCALE64BIT
                LOGGER.info('Detected CHARMM -i8 64-bit DCD file of opposite '
                            'endianness.')
            else:
                endian = ''
                temp = unpack(endian+'ii', bits)
                if temp[0] == 84 and temp[1] == dcdcordmagic:
                    LOGGER.info('Detected standard 32-bit DCD file of '
                                'opposite endianness.')
                else:
                    raise IOError('Unrecognized DCD header or unsupported '
                                  'DCD format.')
                    
        
        # check for magic string, in case of long record markers
        if rec_scale == RECSCALE64BIT:
            raise IOError('CHARMM 64-bit DCD files are not yet supported.');
            temp = unpack('I', dcd.read(calcsize('I')))
            if temp[0] != dcdcordmagic: 
                raise IOError('Failed to find CORD magic in CHARMM -i8 64-bit '
                              'DCD file.');
        
        # Buffer the entire header for random access
        bits = dcd.read(80)
        # CHARMm-genereate DCD files set the last integer in the
        # header, which is unused by X-PLOR, to its version number.
        # Checking if this is nonzero tells us this is a CHARMm file
        # and to look for other CHARMm flags.

        temp = unpack(endian + 'i'*20 , bits)
        if temp[-1] != 0:
            charmm = True

        if charmm:
            LOGGER.info('CHARMM format DCD file (also NAMD 2.1 and later).')
            temp = unpack(endian + 'i'*9 + 'f' + 'i'*10 , bits)
        else:
            LOGGER.info('X-PLOR format DCD file (also NAMD 2.0 and earlier) '
                        'is not supported.')
            return None
        
        # Store the number of sets of coordinates (NSET)
        self._n_confs = temp[0]
        # Store ISTART, the starting timestep
        self._first_ts = temp[1]
        # Store NSAVC, the number of timesteps between dcd saves
        self._save_freq = temp[2]
        # Store NAMNF, the number of fixed atoms
        self._n_fixed = temp[8]
        
        if self._n_fixed > 0:
            raise IOError('DCD files with fixed atoms is not yet supported.')
        
        # Read in the timestep, DELTA
        # Note: DELTA is stored as double with X-PLOR but as float with CHARMm
        self._timestep = temp[9]
        self._unitcell = temp[10]
        
        # Get the end size of the first block
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 84:
            raise IOError('Unrecognized DCD format.')
        
        # Read in the size of the next block
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 164:
            raise IOError('Unrecognized DCD format.')

        # Read NTITLE, the number of 80 character title strings there are
        temp = unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))
        self._title = dcd.read(80).strip()
        self._remarks = dcd.read(80).strip()
        # Get the ending size for this block
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 164:
            raise IOError('Unrecognized DCD format.')


        # Read in an integer '4'
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
            raise IOError('Unrecognized DCD format.')

        # Read in the number of atoms
        self._n_atoms = unpack(endian+'i',dcd.read(rec_scale*calcsize('i')))[0]
        # Read in an integer '4'
        if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
            raise IOError('Bad DCD format.')

        self._is64bit = rec_scale == RECSCALE64BIT
        self._endian = endian
        self._n_floats = (self._n_atoms + 2) * 3
        
        if self._is64bit:
            self._bytes_per_frame = 56 + self._n_floats * 8
            self._dtype = np.float64
        else: 
            self._bytes_per_frame = 56 + self._n_floats * 4
            self._dtype = np.float32
        
        self._first_byte = self._file.tell()
        
        self._coords = self.nextCoordset()
        self._file.seek(self._first_byte)
        self._nfi = 0
        
    def next(self):
        """Return next frame."""
        
        if not self._closed and self._nfi < self._n_confs:
            frame = Frame(self, self._nfi, self.nextCoordset())
            return frame
        
    def nextCoordset(self):
        """Return next coordinate set."""
        
        if not self._closed and self._nfi < self._n_confs:
            #Skip extended system coordinates (unit cell data)
            self._file.seek(56, 1)
            n_floats = self._n_floats
            n_atoms = self._n_atoms
            xyz = np.fromfile(self._file, dtype=self._dtype, count=n_floats)
            if len(xyz) != n_floats:
                return None
            xyz = xyz.reshape((3, n_atoms+2)).T[1:-1,:]
            xyz = xyz.reshape((n_atoms, 3))
            self._nfi += 1
            if self._sel is None:
                return xyz
            else:
                return xyz[self._indices]

    def getCoordsets(self, indices=None):
        """Returns coordinate sets at given *indices*. *indices* may be an 
        integer, a list of integers or ``None``. ``None`` returns all 
        coordinate sets."""
        
        if indices is None or indices == slice(None):
            nfi = self._nfi
            self.reset()
            n_floats = self._n_floats
            n_atoms = self._n_atoms
            n_confs = self._n_confs
            data = np.fromfile(self._file, self._dtype, 
                               self._bytes_per_frame * n_confs)
            data = data.reshape((n_confs, 14 + n_floats))
            data = data[:, 14:]
            data = data.reshape((n_confs, 3, n_atoms+2))
            data = data[:, :, 1:-1]
            data = data.transpose(0, 2, 1)
            self.goto(nfi)
            return data
        else:            
            return TrajectoryFile.getCoordsets(self, indices)


class Trajectory(TrajectoryBase):
    
    """A class for handling trajectories in multiple files."""
        
    def __init__(self, name):
        """Trajectory can be instantiated with a *name* or a filename. When
        name is a valid path to a trajectory file it will be opened for 
        reading."""
        
        TrajectoryBase.__init__(self, name)
        self._trajectory = None
        self._trajectories = []
        self._filenames = set()
        self._n_files = 0
        self._cfi = 0 # current file index
        if os.path.isfile(name):
            self.addFile(name)
        
        
    def __repr__(self):
        if self._closed:
            return ('<Trajectory: {0:s} (closed)>').format(self._name)
        else:
            return ('<Trajectory: {0:s} ({1:d} files, at {2:d} of {3:d} '
                    'frames, selected {4:d} of {5:d} atoms)>').format(
                    self._name, self._n_files, self._nfi, 
                    self._n_confs, self.getNumOfSelected(), self._n_atoms)
    
    def __str__(self):
        return '{0:s} {1:s}'.format(self.__class__.__name__, self._name)

    def _nextFile(self):
        self._cfi += 1
        if self._cfi < self._n_files: 
            self._trajectory = self._trajectories[self._cfi]
        
    def addFile(self, filename):
        """Add a file to the trajectory instance. Currently only DCD files
        are supported."""
        
        if os.path.abspath(filename) in self._filenames:        
            raise IOError('{0:s} is already added to the trajectory'
                          .format(filename))
        ext = os.path.splitext(filename)[1]
        if ext.lower() == '.dcd':
            traj = DCDFile(filename)
        else: 
            raise ValueError('Trajectory file type "{0:s}" is not recognized.'
                             .format(ext))
        
        n_atoms = self._n_atoms
        if n_atoms != 0 and n_atoms != traj.getNumOfAtoms():
            raise IOError('{0:s} must have same number of atoms as previously '
                          'loaded files'.format(traj.getName()))
         
        if self._n_files == 0:
            self._trajectory = traj                
            self._name = traj.getName()
        if n_atoms == 0:
            self._n_atoms = traj.getNumOfAtoms()
            self._coords = traj._coords
        self._trajectories.append(traj)
        self._n_confs += traj.getNumOfFrames()
        self._n_files += 1
   
    def getNumOfFiles(self):
        """Return number of open trajectory files."""
        
        return self._n_files

    def getFilenames(self, absolute=False):
        """Return list of filenames opened for reading."""
        
        return [traj.getFilename(absolute) for traj in self._trajectories]
        

    def getCoordsets(self, indices=None):
        pass
    
    def next(self):
        """Return next frame."""

        if not self._closed and self._nfi < self._n_confs:
            return Frame(self, self._nfi, self.nextCoordset())
    
    def nextCoordset(self):
        """Return next coordinate set."""
        
        if not self._closed and self._nfi < self._n_confs:
            traj = self._trajectory
            while traj._nfi == traj._n_confs:
                self._nextFile()
                traj = self._trajectory
            self._nfi += 1
            return traj.nextCoordset()
    
    def goto(self, n):
        pass
    
    def skip(self, n=1):
        """Skip *n* frames. *n* must be a positive integer, default is 1."""
        
        if not isinstance(n, (int, long)):
            raise ValueError('n must be an integer')
        left = self._n_confs - self._nfi
        if n > left:
            n = left
        while not self._closed and self._nfi < self._n_frames and n > 0:
            traj = self._trajectory
            skip = min(n, traj.getNumOfFrames() - traj.getNextFrameIndex())
            traj.skip(skip)
            if n > skip:
                self._nextFile()
            self._nfi += skip
            n -= skip
    
    def reset(self):
        """Go to the first frame of the first file."""

        if not self._closed and self._trajectories:
            for traj in self._trajectories:
                traj.reset()
            self._trajectory = self._trajectories[0]
            self._cfi = 0
            self._nfi = 0

    def close(self):
        """Close all open trajectory files."""
        
        for traj in self._trajectories:
            traj.close()
        self._closed = True
    
def parseDCD(filename, first=None, last=None, step=None):
    """|new| Parse CHARMM format DCD files (also NAMD 2.1 and later).
    
    .. versionadded:: 0.7.2
    
    .. versionchanged:: 0.8
       Returns an :class:`Ensemble` instance. Conformations in the ensemble
       will be ordered as they appear in the trajectory file. 
       *indices* argument is removed. Use :class:`DCDFile` class for parsing 
       coordinates of a subset of atoms.
    
    :arg filename: DCD filename
    :type filename: str
    
    :arg first: index of first frame to read
    :type first: int
        
    :arg last: index of last frame to read
    :type last: int
        
    :arg step: steps between reading frames, default is 1 meaning every frame
    :type step: int
    
    """
    
    dcd = DCDFile(filename)
    start = time()
    n_frames = dcd.getNumOfFrames()
    LOGGER.info('DCD file contains {0:d} coordinate sets for {1:d} atoms.'
                .format(n_frames, dcd.getNumOfAtoms()))
    ensemble = dcd[slice(first,last,step)]    
    dcd.close()
    time_ = time() - start
    dcd_size = 1.0 * dcd.getNumOfFrames() * dcd._bytes_per_frame / (1024*1024)
    LOGGER.info('DCD file was parsed in {0:.2f} seconds.'.format(time_))
    LOGGER.info('{0:.2f} MB parsed at input rate {1:.2f} MB/s.'
                .format(dcd_size, dcd_size/time_))
    LOGGER.info('{0:d} coordinate sets parsed at input rate {1:d} frame/s.'
                .format(n_frames, int(n_frames/time_)))
    return ensemble
  
if __name__ == '__main__':
    ag = prody.parsePDB('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi.pdb')
    dcd = DCDFile('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/sim.dcd')
    dcd.setAtomGroup( ag )
    traj = Trajectory('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/eq1.dcd')
    traj.addFile('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/sim.dcd')
    traj.setAtomGroup( ag )
    ens = parseDCD('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/sim.dcd')
    #dcd = parseDCD('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/eq1.dcd')
    #dcd = parseDCD('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/sim.dcd', indices=np.arange(1000), stride=10)
    #dcd = parseDCD('/home/abakan/research/mkps/dynamics/mkp3/MKP3.dcd', indices=np.arange(1000), stride=10)
    
    
    
