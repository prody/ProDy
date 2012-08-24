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

""" This module defines a class for identifying contacts."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from prody import LOGGER
from prody.atomic import Atomic, AtomGroup, AtomSubset, AtomMap, AtomPointer
from prody.utilities import importLA

from .measure import calcCenter

__all__ = ['Transformation', 'applyTransformation', 'alignCoordsets',
           'calcRMSD', 'calcTransformation', 'superpose', 
           'moveAtoms', 'wrapAtoms',
           'printRMSD']
           

class Transformation(object):
    
    """A class for storing a transformation matrix."""
    
    __slots__ = ['_matrix']
    
    def __init__(self, *args):
        """Either 4x4 transformation *matrix*, or *rotation* matrix and 
        *translation* vector must be provided at instantiation."""

        nargs = len(args)
        if nargs == 1:
            self.setMatrix(args[0])
        elif nargs == 2:
            self._matrix = np.eye(4)
            self.setRotation(args[0])
            self.setTranslation(args[1])
        else:
            raise ValueError('one or two array must be provided as arguments')            
            
    def getRotation(self): 
        """Return rotation matrix."""

        if self._matrix is not None:
            return self._matrix[:3, :3]

    def setRotation(self, rotation):
        """Set rotation matrix."""
        
        if not isinstance(rotation, np.ndarray):
            raise TypeError('rotation must be an ndarray')
        elif rotation.shape != (3,3):
            raise ValueError('rotation must be a 3x3 array')
        self._matrix[:3, :3] = rotation

    def getTranslation(self): 
        """Return translation vector."""
        
        if self._matrix is not None:
            return self._matrix[:3, 3]
    
    def setTranslation(self, translation): 
        """Set translation vector."""
        
        if not isinstance(translation, np.ndarray):
            raise TypeError('translation must be an ndarray')
        elif translation.shape != (3,):
            raise ValueError('translation must be an ndarray of length 3')
        self._matrix[:3, 3] = translation
    
    def getMatrix(self):
        """Returns a copy of the 4x4 transformation matrix whose top left is 
        rotation matrix and last column is translation vector."""

        if self._matrix is not None:        
            return self._matrix.copy()
    
    def setMatrix(self, matrix):
        
        if not isinstance(matrix, np.ndarray):
            raise TypeError('matrix must be an ndarray')
        elif matrix.shape != (4, 4):
            raise ValueError('matrix must have shape (4,4)')
        self._matrix = matrix
    
    def apply(self, atoms):
        """Apply transformation to *atoms*, see :func:`applyTransformation`
        for details."""
        
        return applyTransformation(self, atoms)


def calcTransformation(mobile, target, weights=None):
    """Returns a :class:`Transformation` instance which, when applied to the 
    atoms in *mobile*, minimizes the weighted RMSD between *mobile* and 
    *target*.  *mobile* and *target* may be NumPy coordinate arrays, or 
    :class:`.Atomic` instances, e.g. :class:`.AtomGroup`, :class:`.Chain`, 
    or :class:`.Selection`."""
    
    name = ''
    if not isinstance(mobile, np.ndarray): 
        try:
            mob = mobile._getCoords()
        except AttributeError:
            raise TypeError('mobile must be a numpy array or an object '
                            'with getCoords method')
    else:
        mob = mobile
    if not isinstance(target, np.ndarray): 
        try:
            tar = target._getCoords()
        except AttributeError:
            raise TypeError('target must be a numpy array or an object '
                            'with getCoords method')
    else:
        tar = target
    
    if mob.shape != tar.shape:
        raise ValueError('reference and target coordinate arrays '
                         'must have same number of atoms')
    
    if mob.shape[1] != 3:
        raise ValueError('reference and target must be coordinate arrays')
    
    if weights is not None:
        if not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape != (mob.shape[0], 1):
            raise ValueError('weights must have shape (n_atoms, 1)')

    return Transformation(*getTransformation(mob, tar, weights))


def getTransformation(mob, tar, weights=None):
    
    linalg = importLA()
    
    if weights is None:
        mob_com = mob.mean(0)
        tar_com = tar.mean(0)
        mob = mob - mob_com
        tar = tar - tar_com
        matrix = np.dot(tar.T, mob)
    else:
        weights_sum = weights.sum()
        weights_dot = np.dot(weights.T, weights)
        mob_com = (mob * weights).sum(axis=0) / weights_sum
        tar_com = (tar * weights).sum(axis=0) / weights_sum
        mob = mob - mob_com
        tar = tar - tar_com
        matrix = np.dot((tar * weights).T, (mob * weights)) / weights_dot

    U, s, Vh = linalg.svd(matrix)
    Id = np.array([ [1, 0, 0], 
                    [0, 1, 0], 
                    [0, 0, np.sign(linalg.det(matrix))] ])
    rotation = np.dot(Vh.T, np.dot(Id, U.T))

    return rotation, tar_com - np.dot(mob_com, rotation)


def applyTransformation(transformation, atoms):
    """Return *atoms* after applying *transformation*.  If *atoms* 
    is a :class:`.Atomic` instance, it will be returned after 
    *transformation* is applied to its active coordinate set.  If 
    *atoms* is an :class:`.AtomPointer` instance, *transformation* 
    will be applied to the corresponding coordinate set in the 
    associated :class:`.AtomGroup`."""
    
    coords = None
    ag = None
    if isinstance(atoms, np.ndarray): 
        if atoms.shape[1] != 3:
            raise ValueError('atoms must be a 3-d coordinate array')
        coords = atoms
        atoms = None
    else:
        if isinstance(atoms, AtomPointer):
            ag = atoms.getAtomGroup()
            acsi = ag.getACSIndex()
            ag.setACSIndex(atoms.getACSIndex())
            coords = ag._getCoords()
        else:
            try:
                coords = atoms._getCoords()
            except AttributeError:
                raise TypeError('atoms must be a Atomic instance')
                
    if atoms is None:
        return _applyTransformation(transformation, coords)
    else:
        if ag is None:
            atoms.setCoords(_applyTransformation(transformation, coords))
        else: 
            ag.setCoords(_applyTransformation(transformation, coords))
            ag.setACSIndex(acsi)
        return atoms


def _applyTransformation(t, coords):
    return t.getTranslation() + np.dot(coords, t.getRotation())

    
def superpose(mobile, target, weights=None):
    """Return *mobile*, after its RMSD minimizing superposition onto *target*, 
    and the transformation that minimizes the RMSD."""
    
    t = calcTransformation(mobile, target, weights)
    result = applyTransformation(t, mobile)
    return (result, t)


def moveAtoms(atoms, array):
    """Move or transform *atoms*. *array* must be :class:`numpy.ndarray`.  
    If shape of *array* is one of ``(natoms, 3)``, ``(1, 3)``, or ``(3,)``,
    *atoms* will be translated. If *array* is a ``(4,4)`` matrix, coordinates
    will be transformed."""
    
    try:
        coords = atoms._getCoords()
    except AttributeError: 
        raise TypeError("atoms doesn't have a valid type: " + str(type(atoms)))
    if not isinstance(array, np.ndarray):
        raise TypeError('offset must be a NumPy array')
    if array.shape[-1] == 3 and array.ndim in (1,2):
        coords += array
    elif array.shape == (4,4):
        coords = np.dot(coords, array[:3,:3])
        coords += array[3,:3]
    else:
        raise ValueError('array does not have right shape')
    atoms.setCoords(coords)
    
    
def moveAtoms(atoms, **kwargs):
    """Move *atoms* *to* a new location or *by* an offset.  This method will
    change the active coordinate set of the *atoms*.  Note that only one of 
    *to* or *by* keyword arguments is expected.
    
    ::
       
       # move protein so that its centroid is at the origin, [0., 0., 0.]
       moveAtoms(protein, to=np.zeros(3)) 
       
       # move protein so that its mass center is at the origin
       moveAtoms(protein, to=np.zeros(3), weights=protein.getMasses())
    
       # move protein so that centroid of Cα atoms is at the origin
       moveAtoms(protein.ca, to=np.zeros(3), ag=True) 

       # move protein by 10 A along each direction
       moveAtoms(protein, by=np.ones(3) * 10)


    :arg by: an offset array with shape ``([1,] 3)`` or ``(n_atoms, 3)`` or
        a transformation matrix with shape ``(4, 4)``   
    :type by: :class:`numpy.ndarray`
    
    :arg to: a point array with shape ``([1,] 3)``   
    :type to: :class:`numpy.ndarray`
    
    :arg ag: when *atoms* is a :class:`.AtomSubset`, apply translation vector
        (*to*) or transformation matrix to the :class:`.AtomGroup`, 
        default is **False** 
    :type ag: bool
    
    :arg weights: array of atomic weights with shape ``(n_atoms[, 1])``
    :type weights: :class:`numpy.ndarray`
    
    When *to* argument is passed, :func:`.calcCenter` function is used to 
    calculate centroid or mass center."""
    
    point = kwargs.pop('to', None)
    if point is None:
        offset = kwargs.pop('by', None)
        if offset is None:
            raise TypeError('moveAtoms() expects one of {0:s} or {1:s} '
                            'arguments'.format(repr('to'), repr('by')))

        try:
            shape = offset.shape            
        except AttributeError: 
            raise TypeError('by must be a numpy array')


        if shape == (4, 4):
            if kwargs.pop('ag', False):
                try:
                    atoms = atoms.getAtomGroup()
                except AttributeError:
                    # itself must be an AtomGroup
                    pass
            try:
                coords = atoms._getCoords()
            except AttributeError:
                try:
                    coords = atoms.getCoords()
                except AttributeError:
                    raise TypeError('atoms must be an Atomic instance')

            coords = np.dot(coords, offset[:3,:3])
            coords += offset[3,:3]

            atoms.setCoords(coords)
        else:

            try:
                coords = atoms._getCoords()
            except AttributeError:
                try:
                    coords = atoms.getCoords()
                except AttributeError:
                    raise TypeError('atoms must be an Atomic instance')


            try:
                natoms = atoms.numAtoms()
            except AttributeError: 
                raise TypeError('atoms must be an Atomic instance')
            if shape == (3,) or shape == (1, 3) or shape == (natoms, 3):
                atoms.setCoords(coords + offset)
            else:
                raise ValueError('by.shape is not valid')
        
    else:
        try:
            shape = point.shape            
        except AttributeError: 
            raise TypeError('to must be a numpy array')
        if shape != (3,) and shape != (1, 3): 
            raise ValueError('to.shape must be ([1,] 3)')
    
        center = calcCenter(atoms, weights=kwargs.pop('weights', None))
        offset = point - center
        if kwargs.pop('ag', False):
            try:
                atoms = atoms.getAtomGroup()
            except AttributeError:
                # itself must be an AtomGroup
                pass

        try:
            coords = atoms._getCoords()
        except AttributeError:
            try:
                coords = atoms.getCoords()
            except AttributeError:
                raise TypeError('atoms must be an Atomic instance')
        atoms.setCoords(coords + offset)
    
def calcRMSD(reference, target=None, weights=None):
    """Return root-mean-square deviation(s) (RMSD) between reference and target 
    coordinates.
    
    >>> ens = loadEnsemble('p38_X-ray.ens.npz')
    >>> print ens.getRMSDs().round(2) # doctest: +ELLIPSIS
    [ 0.74  0.53  0.58  0.6   0.61  0.72  0.62  0.74  0.69  0.65  0.48  0.54
      ...
      0.58  0.66  0.83]
    >>> print calcRMSD(ens).round(2) # doctest: +ELLIPSIS
    [ 0.74  0.53  0.58  0.6   0.61  0.72  0.62  0.74  0.69  0.65  0.48  0.54
      ...
      0.58  0.66  0.83]
    >>> print calcRMSD(ens.getCoords(), ens.getCoordsets(), \
ens.getWeights()).round(2) # doctest: +ELLIPSIS
    [ 0.74  0.53  0.58  0.6   0.61  0.72  0.62  0.74  0.69  0.65  0.48  0.54
      ...
      0.58  0.66  0.83]"""
    
    if isinstance(reference, np.ndarray): 
        ref = reference
    else:
        try:
            ref = reference._getCoords()
        except AttributeError:
            raise TypeError('reference must be a numpy array or an object '
                            'with getCoords method')
        if target is None:
            try:
                target = reference._getCoordsets()
            except AttributeError:
                pass
        if weights is None:
            try:
                weights = reference._getWeights()
            except AttributeError:
                pass
    if ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError('reference must have shape (n_atoms, 3)')
    
    if isinstance(target, np.ndarray): 
        tar = target
    else:
        try:
            tar = target._getCoords()
        except AttributeError:
            raise TypeError('target must be a numpy array or an object '
                            'with getCoords method')
    if tar.ndim not in (2, 3) or tar.shape[-1] != 3:
        raise ValueError('target must have shape ([n_confs,] n_atoms, 3)')

    if ref.shape != tar.shape[-2:]:
        raise ValueError('reference and target arrays must have the same '
                         'number of atoms')
    
    if weights is not None:
        if not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif not ((weights.ndim == 2 and len(weights) == len(ref)) or
            (weights.ndim == 3 and weights.shape[:2] == target.shape[:2])) or \
             weights.shape[-1] != 1:
            raise ValueError('weights must have shape ([n_confs,] n_atoms, 1)')
    return getRMSD(ref, tar, weights)
    
    
def getRMSD(ref, tar, weights=None):
    if weights is None:
        divByN = 1.0 / ref.shape[0]
        if tar.ndim == 2:
            return np.sqrt(((ref-tar) ** 2).sum() * divByN)
        else:
            rmsd = np.zeros(len(tar))
            for i, t in enumerate(tar):
                rmsd[i] = ((ref-t) ** 2).sum() 
            return np.sqrt(rmsd * divByN)
    else:
        if tar.ndim == 2:
            return np.sqrt((((ref-tar) ** 2) * weights).sum() * 
                                                    (1 / weights.sum()))
        else:
            rmsd = np.zeros(len(tar))
            if weights.ndim == 2:
                for i, t in enumerate(tar):
                    rmsd[i] = (((ref-t) ** 2) * weights).sum() 
                return np.sqrt(rmsd * (1 / weights.sum()))
            else:
                for i, t in enumerate(tar):
                    rmsd[i] = (((ref-t) ** 2) * weights[i]).sum()
                return np.sqrt(rmsd / weights.sum(1).flatten())
    
    
def printRMSD(reference, target=None, weights=None, log=True, msg=None):
    """Print RMSD to the screen.  If *target* has multiple coordinate sets, 
    minimum, maximum and mean RMSD values are printed.  If *log* is **True** 
    (default), RMSD is written to the standard error using package logger, 
    otherwise standard output is used.  When *msg* string is given, it is
    printed before the RMSD value.  See also :func:`calcRMSD` function. """
    
    if log:
        write = LOGGER.info
    else:
        import sys
        write = lambda line: sys.stdout.write(line + '\n')
    msg = msg or ''
    if msg and msg[-1] != ' ':
        msg += ' '
    rmsd = calcRMSD(reference, target, weights)
    if isinstance(rmsd, np.ndarray) and len(rmsd) > 1:
        write(msg + 'RMSD: min={0:.2f}, max={1:.2f}, mean={2:.2f}'
                    .format(rmsd.min(), rmsd.max(), rmsd.mean()))
    else:
        if isinstance(rmsd, np.ndarray):
            rmsd = rmsd[0]
        write(msg + 'RMSD: {0:.2f}'.format(rmsd))
        
    
def alignCoordsets(atoms, weights=None):
    """Return *atoms* after superposing coordinate sets onto its active 
    coordinate set.  Transformations will be calculated for *atoms* and 
    applied to its :class:`.AtomGroup`, when applicable.  Optionally, 
    atomic *weights* can be passed for weighted superposition."""
    
    try:
        acsi, n_csets = atoms.getACSIndex(), atoms.numCoordsets()
    except AttributeError:
        raise TypeError('atoms must have type Atomic, not {0:s}'
                        .format(type(atoms)))
        if n_csets < 2:
            LOGGER.warning('{0:s} contains fewer than two coordinate sets, '
                           'alignment was not performed.'.format(str(atoms)))
            return
    
    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
    agacsi = ag.getACSIndex()

    tar = atoms._getCoords()
    for i in range(n_csets):
        if i == acsi:
            continue
        atoms.setACSIndex(i)
        ag.setACSIndex(i)
        calcTransformation(atoms, tar, weights).apply(ag)
    atoms.setACSIndex(acsi)
    ag.setACSIndex(agacsi)
    return atoms


def wrapAtoms(frame, unitcell=None, center=np.array([0., 0., 0.])):
    """Wrap atoms into an image of the system simulated under periodic boundary
    conditions. When *frame* is a :class:`.Frame`, unitcell information will be
    retrieved automatically.  
    
    .. note::
       This function will wrap all atoms into the specified periodic image, so 
       covalent bonds will be broken.
    
    :arg frame: a frame instance or a coordinate set
    :type frame: :class:`.Frame`, :class:`.AtomGroup`, :class:`numpy.ndarray`
    
    :arg unitcell: orthorhombic unitcell array with shape (3,)
    :type unitcell: :class:`numpy.ndarray`
    
    :arg center: coordinates of the center of the wrapping cell, default is 
        the origin of the Cartesian coordinate system
    :type center: :class:`numpy.ndarray`"""
    
    try:    
        coords = frame._getCoords()
    except AttributeError:
        coords = frame
    else:
        try:
            frame.getAtomGroup()
        except AttributeError:
            pass
        else:
            raise TypeError('frame must be a Frame, AtomGroup, or numpy array,'
                            ' not a ' + str(type(frame)))            

    if unitcell is None:
        try:
            unitcell = frame.getUnitcell()[:3]
        except AttributeError:
            raise TypeError('unitcell information must be provided')
    
    half = unitcell / 2
    ucmin = center - half
    ucmax = center + half
    for axis in range(3):
        xyz = coords[:, axis]
        which = (xyz < ucmin[axis]).nonzero()[0]
        while len(which):
            coords[which, axis] += unitcell[axis]
            which = which[coords[which, axis] < ucmin[axis]]
        which = (xyz > ucmax[axis]).nonzero()[0]
        while len(which):
            coords[which, axis] -= unitcell[axis]
            which = which[coords[which, axis] > ucmax[axis]]
    return frame
