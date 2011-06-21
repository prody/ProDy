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
  
Functions
---------

    * :func:`parseDCD`
    * :func:`saveEnsemble`
    * :func:`loadEnsemble`
    * :func:`calcSumOfWeights`
    * :func:`showSumOfWeights`
    * :func:`trimEnsemble`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import time
import os.path
import numpy as np
import prody
from prody import ProDyLogger as LOGGER
from prody import measure

__all__ = ['Ensemble', 'Conformation', 'saveEnsemble', 'loadEnsemble', 
           'calcSumOfWeights', 'showSumOfWeights', 'trimEnsemble',
           'parseDCD']
        
class EnsembleError(Exception):
    pass

plt = None

class Ensemble(object):
    
    """A class for ensemble analysis. This class enables handling coordinates
    for heterogeneous structural datasets. 
    
    |example| See usage usage in :ref:`pca-xray`, :ref:`pca-dimer`, and 
    :ref:`pca-blast`.
    
    Indexing returns a :class:`Conformation` instance, whereas slicing returns 
    an :class:`Ensemble` instance that contains subset of conformations.
    The ensemble obtained by slicing has the same reference coordinates.

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
        
        self._name = str(name)
        if self._name == '':
            self._name = 'Unnamed'
        
        self._ensemble = []
        self._confs = None       # coordinate data
        self._weights = None
        self._coords = None         # reference
        self._n_atoms = None
        self._n_confs = 0
        self._transformations = []    # from last superimposition
        
        if isinstance(name, prody.Atomic):
            self.setCoordinates(name.getCoordinates())
            self.addCoordset(name.getCoordsets())
        
    def __getitem__(self, index):
        """Return a conformation at given index."""
        
        if isinstance(index, int):
            return self._getConformation(index) 
        elif isinstance(index, slice):
            ens = Ensemble('{0:s} ({1[0]:d}:{1[1]:d}:{1[2]:d})'.format(
                                self._name, index.indices(len(self))))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self._confs[index].copy(), self._weights[index].copy())
            return ens
        elif isinstance(index, (list, np.ndarray)):
            ens = Ensemble('{0:s} slice'.format(self._name))
            ens.setCoordinates(self.getCoordinates())
            ens.addCoordset(self._confs[index].copy(), self._weights[index].copy())
            return ens
        else:
            raise IndexError('invalid index')
            
    
    def __add__(self, other):
        """Return ensemble that contains conformations in *self* and *other*.
        
        The reference coordinates of *self* is used in the result.
        
        """
        
        if not isinstance(other, Ensemble):
            raise TypeError('an Ensemble instance cannot be added to an {0:s} '
                            'instance'.format(type(other)))
        elif self.getNumOfAtoms() != other.getNumOfAtoms():
            raise ValueError('Ensembles must have same number of atoms.')
    
        ensemble = Ensemble('{0:s} + {1:s}'.format(self.getName(), 
                                                   other.getName()))
        
        return None
    
    def __len__(self):
        return self._n_confs
        
    def __iter__(self):
        return self._ensemble.__iter__()
    
    def __repr__(self):
        return '<Ensemble: {0:s} ({1:d} conformations, {2:d} atoms)>'.format(
                                        self._name, len(self), self._n_atoms)

    def __str__(self):
        return 'Ensemble {0:s}'.format(self._name)

    def _getConformation(self, index):
        conf = self._ensemble[index]
        if conf is None:
            conf = Conformation(self, index, str(index))
            self._ensemble[index] = conf
        return conf
    
    def getName(self):
        """Return name of the conformation ensemble."""
        return self._name

    def setName(self, name):
        """Set name of the ensemble instance."""
        self._name = name
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        return self._n_atoms
    
    def getNumOfConfs(self):
        """Return number of conformations."""
        return len(self._ensemble)

    def getCoordinates(self):
        """Return reference coordinates of the ensemble."""
        
        if self._coords is None:
            return None
        return self._coords.copy()

    def setCoordinates(self, coords):
        """Set reference coordinates of the ensemble."""

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
        elif coords.dtype != np.float64:
            try:
                coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        
        if self._n_atoms is not None:
            if coords.shape[0] != self._n_atoms:
                raise ValueError('shape of coords must be (n_atoms,3)')
        else:
            self._n_atoms = coords.shape[0]    
        self._coords = coords
        
   
    def addCoordset(self, coords, weights=None, allcoordsets=True):
        """Add a NumPy array or coordinates from atoms as a conformation.
        
        :class:`~prody.atomic.Atomic` instances are acceptable as *coords* 
        argument. If *allcoordsets* is ``True``, all coordinate sets will be 
        appended to the ensemble.
        
        *weights* is an optional argument. If provided, its length must
        match number of atoms.
        
        """
        
        ag = None
        if not isinstance(coords, np.ndarray):
            try:
                ag = coords
                if allcoordsets:
                    coords = ag.getCoordsets()
                else:
                    coords = ag.getCoordinates()
            except AttributeError:
                raise TypeError('coords must be an ndarray instance or '
                                'must contain getCoordinates as an attribute')
            
        if not coords.ndim in (2, 3):
            raise ValueError('coords must be a 2d or a 3d array')
        elif coords.dtype != np.float64:
            try:
                coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        
        if self._n_atoms is None:
            self._n_atoms = coords.shape[-2]
        else:
            if coords.shape[-2] != self._n_atoms:
                raise ValueError('shape of conf must be (n_atoms,3)')
        n_atoms = self._n_atoms

        if coords.ndim == 2:
            coords = coords.reshape((1, self._n_atoms, 3))
        n_confs = coords.shape[0]
            
        
        if weights is not None:
            if not isinstance(weights, np.ndarray):
                raise TypeError('weights must be an ndarray')
            elif not weights.ndim in (1, 2, 3):
                raise ValueError('weights must be a 1d or a 2d array')
            elif weights.ndim in (2, 3) and weights.shape[-1] != 1:
                raise ValueError('shape of weights must be ([n_coordsets,] number_of_atoms, 1)')
            elif weights.dtype != np.float64:
                try:
                    weights.astype(np.float64)
                except ValueError:
                    raise ValueError('weights array cannot be assigned type '
                                     '{0:s}'.format(np.float64))
            if weights.ndim == 1:
                weights = weights.reshape((1, n_atoms, 1))
                if n_confs > 1:
                    weights = np.tile(weights, (n_confs, 1, 1))
        if self._confs is None: 
            self._confs = coords
            self._weights = weights
        else:
            self._confs = np.concatenate((self._confs, coords), axis=0)
            if weights is not None:
                self._weights = np.concatenate((self._weights, weights), axis=0)
            else:
                if self._weights is not None:
                    self._weights = np.concatenate((self._weights, 
                                        np.ones(n_confs, n_atoms, 1)), axis=0)
        if ag is None:
            self._ensemble += [None] * n_confs
            
        else:
            name = ag.getName()
            if n_confs == 1:
                if ag.getNumOfCoordsets() > 0:
                    name +=  '_' + str(ag.getActiveCoordsetIndex())
                self._ensemble.append(Conformation(self, len(self._ensemble), name))
            else:                
                for i in range(0, n_confs):
                    self._ensemble.append(Conformation(self, len(self._ensemble), name + '_' + str(i)))
        self._n_confs += n_confs
        self._transformations += [None] * n_confs

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

        coords = self._confs[indices].copy()
        if self._weights is not None and self._weights.ndim == 3:
            for i, w in enumerate(self._weights[indices]):
                which = w.flatten()==0
                coords[i, which] = self._coords[which]
        return coords 

    
    def delCoordset(self, index):
        """Delete a coordinate set from the ensemble."""
        
        length = len(self._ensemble)
        which = np.ones(length, np.bool)
        which[index] = False
        if which.sum() == 0:
            self._confs = None
            self._weights = None
        else:
            self._confs = self._confs[which]
            if self._weights is not None:
                self._weights = self._weights[which]
            
        if isinstance(index, int):
            index = [index]
        else:
            index = list(index)
        self._n_confs -= len(index)
        index.sort(reverse=True)
        for i in index:
            conf = self._ensemble.pop(i)
            if conf is not None:
                conf._index = None
                conf._ensemble = None
            self._transformations.pop(i)

    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        return len(self._ensemble)
    
    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate set.
        
        Reference coordinates are not included.
        
        """
        
        for i in xrange(self._confs.shape[0]):
            yield self._confs[i].copy()
    
    def getWeights(self):
        """Return weights."""
        if self._weights is None:
            return None
        else:
            return self._weights.copy()
    
    def getConformation(self, index):
        """Return conformation at given index."""
        return self._getConformation(index)
        
        
    def superpose(self):
        """Superpose the ensemble onto the reference coordinates."""
        
        if self._coords is None:
            raise AttributeError('coordinates are not set, '
                                 'use setCoordinates() method to set it.')
        if self._confs is None or len(self._confs) == 0: 
            raise AttributeError('conformations are not set, '
                                 'use addCoordset() method to set it.')
        LOGGER.info('Superimposing structures.')
        start = time.time()
        self._superpose()
        LOGGER.info(('Superimposition is completed in {0:.2f} '
                           'seconds.').format((time.time() - start)))
        
    def _superpose(self):
        """Superpose conformations and return new coordinates.
        
        This functions saves transformations in self._tranformations.
        
        """

        weights = self._weights
        coords = self._coords
        confs = self._confs
        if not isinstance(self._transformations, list) or \
            len(self._transformations) != self._n_confs:
                self._transformations = [None] * self._n_confs
        transformations = self._transformations
        if weights is None:
            for i in xrange(len(self)):
                conf, t = measure.superpose(confs[i], coords)
                confs[i] = conf
                transformations[i] = t
        else:         
            for i in xrange(len(self)):
                conf, t = measure.superpose(confs[i], coords, weights[i])
                confs[i] = conf
                transformations[i] = t
        
    def transform(self):
        """Apply transformations from previous superimposition step.
        
        A potential use of this method is that superimposition may be performed
        by a core set of atoms. Then set coordinates may be used to select all
        atoms to apply the transformations calculated for the core atoms.
        
        """
        
        if self._transformations.count(None) > 0:
            raise EnsembleError('A transformation for each '
                                    'conformation is not calculated.')
        for i in range(self.conformations.shape[0]):
            self.conformations[i] = GF.transform(self.conformations[i], 
                                                 self._transformations[i][0],
                                                 self._transformations[i][1], 
                                                 self.weights[i])
    
    def iterpose(self, rmsd=0.0001):
        """Iteratively superpose the ensemble until convergence.
        
        Initially, all conformations are aligned with the reference 
        coordinates. Then mean coordinates are calculated, and are set
        as the new reference coordinates. This is repeated until 
        reference coordinates do not change. This is determined by
        the value of RMSD between the new and old reference coordinates.        
        
        :arg cutoff: RMSD (A) between old and new reference coordinates 
                     to converge
        :type cutoff: float, default is 0.0001
            
        """
        
        if self._coords is None:
            raise AttributeError('coordinates are not set, '
                                 'use setCoordinates() method to set it.')
        if self._confs is None or len(self._confs) == 0: 
            raise AttributeError('conformations are not set, '
                                 'use addCoordset() method to set it.')
        LOGGER.info('Starting iterative superimposition')
        start = time.time()
        rmsdif = 1
        step = 0
        weights = self._weights
        if weights is not None:
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
        LOGGER.info('Iterative superimposition completed in {0:.2f}s.'
                    .format((time.time() - start)))
        
    def getMSF(self):
        """Calculate and return Mean-Square-Fluctuations."""
        
        if self.conformations is None: 
            return
        
        xyzmean = (self.conformations * 
                   self.weights).sum(0) / self.weights.sum(0)
        xyzdiff2 = np.power(self.conformations - xyzmean, 2).sum(2)
        weightsum = self.weights.sum(2)
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
        
        if self._confs is None: 
            return None
        weights = self._weights
        if weights is None:
            wsum_axis_2 = 1
            wsum_axis_1 = 1
        else:
            wsum_axis_2 = weights.sum(2)
            wsum_axis_1 = wsum_axis_2.sum(1)
        rmsd = np.sqrt((np.power(self.getDeviations(), 2).sum(2) * 
                        wsum_axis_2).sum(1) / wsum_axis_1)
        return rmsd


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
        raise TypeError('invalid type for ensemble, {0:s}'.format(type(ensemble)))
    if len(ensemble) == 0:
        raise ValueError('ensemble instance does not contain data')
    
    dict_ = ensemble.__dict__
    attr_list = ['_name', '_confs', '_weights', '_coords', '_n_atoms', '_n_confs']
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
    ensemble = Ensemble(str(attr_dict['_name']))
    dict_ = ensemble.__dict__ 
    for attr in attr_dict.files:
        if attr == '_name': 
            continue
        elif attr in ('_n_atoms', '_n_confs'):
            dict_[attr] = int(attr_dict[attr])
        else:
            dict_[attr] = attr_dict[attr]
    ensemble._ensemble = [None] * len(ensemble)
    return ensemble

class Conformation(object):
    
    """A class to provide methods on a conformation in a an ensemble.
    
    Instances of this class do not keep coordinate and weights data.
    
    """
    
    __slots__ = ['_ensemble', '_index', '_name']

    def __init__(self, ensemble, index, name):
        """Instantiate with an ensemble instance, index, and a name."""
        self._name = name
        self._ensemble = ensemble
        self._index = index
        
    def __repr__(self):
        return '<Conformation: {0:s} from {1:s} (index {2:d})>'.format(
                    self._name, self._ensemble._name, self._index)

    def getEnsemble(self):
        """Return the ensemble that this conformation belongs to."""
        return self._ensemble
    
    def getIndex(self):
        """Return the index of the conformation."""
        return self._index
    
    def getName(self):
        """Return name of the conformation instance."""
        return self._name
    
    def setName(self, name):
        """Set name of the conformation instance."""
        self._name = name

    def getCoordinates(self):
        """Return coordinate set for this conformation."""
        if self._ensemble._confs is None:
            return None
        if self._ensemble._weights is None:
            return self._ensemble._confs[self._index].copy()
        else:
            coords = self._ensemble._confs[self._index].copy()
            weights = self._ensemble._weights[self._index].flatten()
            which = weights==0
            coords[which] = self._ensemble._coords[which]
            return coords 
    
    def getWeights(self):
        """Return coordinate weights, eg. occupancy or mass."""
        if self._ensemble._weights is None:
            return None
        else:
            return self._ensemble._weights[self._index].copy()

    def getDeviations(self):
        """Return deviations from the ensemble reference coordinates."""
        return self.getCoordinates() - self._ensemble._coords

    def getRMSD(self):
        """Return RMSD from the ensemble reference coordinates."""
        return measure._calcRMSD(self._ensemble._coords,
                                 self._ensemble._confs[self._index], 
                                 self.getWeights())
    
    def getTransformation(self):
        """Return the transformation from the last superimposition."""
        return self._ensemble._transformations[self._index].copy()
    
def trimEnsemble(ensemble, **kwargs):
    """Return an ensemble obtained by trimming given *ensemble*.
    
    .. versionadded:: 0.5.3
    
    This function helps selecting atoms in an ensemble based on one of the 
    following criteria, and returns them in a new :class:`Ensemble` instance.
        
    **Occupancy**
    
    If weights in the ensemble correspond to atomic occupancies, then this 
    option may be used to select atoms with average occupancy values higher
    than the argument *occupancy*. Function will calculate average occupancy 
    by dividing the sum of weights for each atom to the number of conformations
    in the ensemble.
    
    :arg occupancy: average occupancy for selection atoms in an ensemble, 
                    must be 0 < *occupancy* <= 1
    :type occupancy: float
    
    **Atom selection**
    
    Will be implemented.
    
    """
    
    if not isinstance(ensemble, Ensemble):
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
    
    trimmed = Ensemble(ensemble.getName())
    coords = ensemble.getCoordinates()
    if coords is not None:
        trimmed.setCoordinates( coords[torf] )
    confs = ensemble.getCoordsets()
    if confs is not None:
        weights = ensemble.getWeights()
        if weights is not None:
            if weights.ndim == 1:
                trimmed.addCoordset( confs[:, torf], weights[torf] )
            else:
                trimmed.addCoordset( confs[:, torf], weights[:, torf] )
        else:
            trimmed.addCoordset( confs[:, torf] )
    return trimmed

def calcSumOfWeights(ensemble):
    """Return sum of weights from an ensemble.
    
    Weights are summed for each atom over conformations in the ensemble.
    Size of the plotted array will be equal to the number of atoms.
    
    When analyzing an ensemble of X-ray structures, this function can be used 
    to see how many times a residue is resolved.
    
    """
    
    if not isinstance(ensemble, Ensemble):
        raise TypeError('ensemble must be an Ensemble instance')
    
    weights = ensemble.getWeights()
    
    if weights is None:
        return None
    
    return weights.sum(0)
    
    
def showSumOfWeights(ensemble, *args, **kwargs):
    """Show sum of weights for an ensemble using :func:`~matplotlib.pyplot.plot`.
    
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
    if not isinstance(ensemble, Ensemble):
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

def parseDCD(filename, indices=None, first=None, last=None, stride=None):
    """|new| Parse CHARMM format DCD files (also NAMD 2.1 and later).
    
    .. versionadded:: 0.7.2
    
    The type of coordinate data array returned by this function is 
    :class:`~numpy.float32`. 
    
    :arg filename: DCD filename
    :type filename: str
    
    :arg indices: indices of atoms to read coordinates for, default is all
    :type indices: list
    
    :arg first: index of first frame to read
    :type first: int
        
    :arg last: index of last frame to read
    :type last: int
        
    :arg stride: steps between reading frames, default is 1
    :type stride: int
    
    """
    if not os.path.isfile(filename):
        raise IOError('file not found, {0:s} is not a valid path'
                      .format(filename))
    
    if indices is not None:
        try:
            indices = np.array(indices, dtype=np.int32)
        except: 
            raise TypeError('indices argument is not a list of Numpy array')
        if indices.ndim != 1:
            raise ValueError('indices argument must be one dimensional')

    start = time.time() 
    dcd = open(filename, 'rb')
    header = _parseDCDHeaader(dcd)
    if header is None:
        return None
    endian = header['endian']
    unitcell = header['unitcell']
    bit64 = header['is64bit'] == 2
    n_atoms = header['n_atoms']
    if np.max(indices) >= n_atoms:
        raise ValueError('maximum index exceeds number of atoms in DCD file')
    n_frames = header['n_frames']
    LOGGER.info('DCD file contains {0:d} coordinate sets for {1:d} atoms.'
                .format(n_frames, n_atoms))
    if first is None:
        first = 0
    else:
        first = int(first)
    if stride is None:
        stride = 1
    else:
        stride = int(stride)
    if last is None:
        last = n_frames
    else:
        last = int(last)
    n_floats = (n_atoms + 2) * 3 
    dtype = np.dtype(endian+'f')
    n_frames = 0
    coords = []
    for i in range(first, last, stride):
        if stride > 1:
            for j in range(1, stride):
                dcd.seek(56, 1)
                dcd.seek(n_floats*4, 1)
        if unitcell:
            _parseDCDXSC(dcd)
        xyz = _parseDCDFrame(dcd, n_atoms, n_floats, dtype) 
        if xyz is None:
            LOGGER.error('DCD file ended unexpectedly, returning parsed data.')            
            break
        if indices is not None:
            xyz = xyz[:,indices,:]
        coords.append(xyz)#.astype('d'))
        n_frames += 1
        
    dcd.close()
    coords = np.concatenate(coords)
    time_ = time.time() - start
    dcd_size = 1.0 * n_frames * ((n_atoms + 2) * 4 * 3) / (1024 * 1024) 
    LOGGER.info('DCD file was parsed in {0:.2f} seconds.'.format(time_))
    LOGGER.info('{0:.2f} MB parsed at input rate {1:.2f} MB/s.'
                .format(dcd_size, dcd_size/time_))
    LOGGER.info('{0:d} coordinate sets parsed at input rate {1:d} frame/s.'
                .format(n_frames, int(n_frames/time_)))
    
    return coords

def _parseDCDXSC(dcd):
    """For now, skip extended system coordinates (unit cell data)."""
    dcd.seek(56, 1)

def _parseDCDFrame(dcd, n_atoms, n_floats, dtype):
    xyz = np.fromfile(dcd, dtype=dtype, count=n_floats)
    if len(xyz) != n_floats:
        return None
    xyz = xyz.reshape((3, n_atoms+2)).T[1:-1,:]
    xyz = xyz.reshape((1, n_atoms, 3))
    return xyz
    
def _parseDCDHeaader(dcd):
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
    
    header = {}
    endian = '' #'=' # native endian
    rec_scale = RECSCALE32BIT
    charmm = None
    dcdcordmagic = unpack(endian+'i', 'CORD')[0]
    # Check magic number in file header and determine byte order
    bits = dcd.read(calcsize('ii'))
    temp = unpack(endian+'ii', bits)
    if temp[0] + temp[1] == 84:
        LOGGER.info('Detected CHARMM -i8 64-bit DCD file of native endianness.')
        rec_scale = RECSCALE64BIT
    elif temp[0] == 84 and temp[1] == dcdcordmagic:
        LOGGER.info('Detected standard 32-bit DCD file of native endianness.')
    else:
        if unpack('>ii', bits) == temp:
            endian = '>'
        else:
            endian = '<'
        temp = unpack(endian+'ii', bits)
        if temp[0] + temp[1] == 84:
            rec_scale = RECSCALE64BIT
            LOGGER.info("Detected CHARMM -i8 64-bit DCD file of opposite endianness.")
        else:
            endian = ''
            temp = unpack(endian+'ii', bits)
            if temp[0] == 84 and temp[1] == dcdcordmagic:
                LOGGER.info('Detected standard 32-bit DCD file of opposite endianness.')
            else:
                LOGGER.error('Unrecognized DCD header or unsupported DCD format.')
                return None
                
    
    # check for magic string, in case of long record markers
    if rec_scale == RECSCALE64BIT:
        LOGGER.error("CHARMM 64-bit DCD files are not yet supported.");
        return None
        temp = unpack('I', dcd.read(calcsize('I')))
        if temp[0] != dcdcordmagic: 
            LOGGER.error("Failed to find CORD magic in CHARMM -i8 64-bit DCD file.");
            return None
    
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
        LOGGER.info("CHARMM format DCD file (also NAMD 2.1 and later).")
        temp = unpack(endian + 'i'*9 + 'f' + 'i'*10 , bits)
    else:
        LOGGER.info("X-PLOR format DCD file (also NAMD 2.0 and earlier) is not supported.")
        return None
    
    # Store the number of sets of coordinates (NSET)
    header['n_frames'] = temp[0]
    # Store ISTART, the starting timestep
    header['first_timestep'] = temp[1]
    # Store NSAVC, the number of timesteps between dcd saves
    header['save_frequency'] = temp[2]
    # Store NAMNF, the number of fixed atoms
    header['n_fixed_atoms'] = temp[8]
    
    if header['n_fixed_atoms'] > 0:
        LOGGER.error('DCD files with fixed atoms is not yet supported.')
        return None
    
    # Read in the timestep, DELTA
    # Note: DELTA is stored as a double with X-PLOR but as a float with CHARMm
    header['timestep'] = temp[9]
    header['unitcell'] = temp[10]
    
    # Get the end size of the first block
    if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 84:
        LOGGER.error('Unrecognized DCD format.')
        return None
    
    # Read in the size of the next block
    if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 164:
        LOGGER.error('Unrecognized DCD format.')
        return None
    # Read NTITLE, the number of 80 character title strings there are
    temp = unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))
    header['title'] = dcd.read(80).strip()
    header['remarks'] = dcd.read(80).strip()
    # Get the ending size for this block
    if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 164:
        LOGGER.error('Unrecognized DCD format.')
        return None

    # Read in an integer '4'
    if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
        LOGGER.error('Unrecognized DCD format.')
        return None
    # Read in the number of atoms
    header['n_atoms'] = unpack(endian+'i', dcd.read(rec_scale*calcsize('i')))[0]
    # Read in an integer '4'
    if unpack(endian+'i', dcd.read(rec_scale * calcsize('i')))[0] != 4:
        LOGGER.error('Bad DCD format.')
        return None

    header['is64bit'] = rec_scale == RECSCALE64BIT
    header['endian'] = endian
    
    return header


    

if __name__ == '__main__':
    dcd = parseDCD('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/eq1.dcd')
    dcd = parseDCD('/home/abakan/research/bcianalogs/mdsim/nMbciR/mkp3bcirwi_sim/sim.dcd', indices=np.arange(1000), stride=10)
    #dcd = parseDCD('/home/abakan/research/mkps/dynamics/mkp3/MKP3.dcd', indices=np.arange(1000), stride=10)
    
    
    
