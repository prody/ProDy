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

from prody.atomic import Atomic, AtomGroup, AtomSubset, AtomMap, AtomPointer
from prody.tools import importLA

__all__ = ['Transformation', 'applyTransformation', 'alignCoordsets',
           'calcRMSD', 'calcTransformation', 'superpose', 'moveAtoms']
           
pkg = __import__(__package__)
LOGGER = pkg.LOGGER

class Transformation(object):
    
    __slots__ = ['_rotation', '_translation']
    
    def __init__(self, rotation, translation):
        self._rotation = None
        self.setRotation(rotation)
        self._translation = None
        self.setTranslation(translation)
    
    def getRotation(self): 
        """Returns rotation matrix."""
        
        return self._rotation.copy()

    def setRotation(self, rotation):
        """Set rotation matrix."""
        
        if not isinstance(rotation, np.ndarray):
            raise TypeError('rotation must be an ndarray')
        elif rotation.shape != (3,3):
            raise ValueError('rotation must be a 3x3 array')
        self._rotation = rotation

    def getTranslation(self): 
        """Returns translation vector."""
        
        return self._translation.copy()
    
    def setTranslation(self, translation): 
        """Set translation vector."""
        
        if not isinstance(translation, np.ndarray):
            raise TypeError('translation must be an ndarray')
        elif translation.shape != (3,):
            raise ValueError('translation must be an ndarray of length 3')
        self._translation = translation
    
    def get4x4Matrix(self):
        """Returns 4x4 transformation matrix whose top left is rotation matrix
        and last column is translation vector."""
        
        fourby4 = np.eye(4)
        fourby4[:3, :3] = self._rotation
        fourby4[:3, 3] = self._translation
        return fourby4
    
    def apply(self, atoms):
        """Applies transformation to given atoms or coordinate set.  ProDy 
        class instances from :mod:`~prody.atomic` are accepted.  Instance is 
        returned after its active coordinate set is transformed.  If a 
        :class:`~.AtomPointer` is passsed, the :class:`~.AtomGroup` that it 
        points to is transformed. 
        
        If an :class:`~numpy.ndarray` instance is given, transformed array 
        is returned."""
        
        return applyTransformation(self, atoms)
    

def calcTransformation(mobile, target, weights=None):
    """Returns a :class:`Transformation` instance which, when applied to the 
    atoms in *mobile*, minimizes the weighted RMSD between *mobile* and 
    *target*.
    
    *mobile* and *target* may be NumPy coordinate arrays, or istances of 
    Molecule, AtomGroup, Chain, or Residue."""
    
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
                         'must have same shape')
    

    if mob.shape[1] != 3:
        raise ValueError('reference and target must be 3-d coordinate arrays')
    
    if weights is not None:
        if not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape != (mob.shape[0], 1):
            raise ValueError('weights must have shape (n_atoms, 1)')

    return _calcTransformation(mob, tar, weights)

def _calcTransformation(mob, tar, weights=None):
    
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

    return Transformation(rotation, tar_com - np.dot(mob_com, rotation))

def _superposeTraj(mobs, tar, weights=None, movs=None):
    # mobs.ndim == 3 and movs.ndim == 3
    # mobs.shape[0] == movs.shape[0]
    linalg = importLA()
    svd = linalg.svd
    det = linalg.det
    dot = np.dot
    add = np.add
    subtract = np.subtract
    array = np.array
    sign = np.sign
    
    tar_com = tar.mean(0)
    tar_org_T = (tar - tar_com).T
    mob_org = np.zeros(mobs.shape[-2:])

    LOGGER.progress('Superposing ', len(mobs))
    for i, mob in enumerate(mobs):      
        mob_com = mob.mean(0)        
        matrix = dot(tar_org_T, subtract(mob, mob_com, mob_org))
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

def _superpose(mob, tar, weights=None, mov=None):
    tar_com = tar.mean(0)
    tar_org = tar - tar_com

    linalg = importLA()
    mob_com = mob.mean(0)
    mob_org = mob - mob_com
    matrix = np.dot(tar_org.T, mob_org)

    U, s, Vh = linalg.svd(matrix)
    Id = np.array([ [1, 0, 0], 
                    [0, 1, 0], 
                    [0, 0, np.sign(linalg.det(matrix))] ])
    rotation = np.dot(Vh.T, np.dot(Id, U.T))

    if mov is None:
        np.add(np.dot(mob_org, rotation), tar_com, mob) 
    else:
        np.add(np.dot(mov, rotation), 
               (tar_com - np.dot(mob_com, rotation)), mov)


def applyTransformation(transformation, atoms):
    """Return *atoms* after applying *transformation*. If *atoms* is a 
    class instance from :mod:`~.atomic`, it will be returned after 
    *transformation* is applied to its active coordinate set. If *atoms*
    is an :class:`~.AtomPointer` instance, *transformation* will
    be applied to the corresponding coordinate set in the associated
    :class:`~.AtomGroup` instance."""
    
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
    return t._translation + np.dot(coords, t._rotation)
    

def calcRMSD(reference, target=None, weights=None):
    """Returns Root-Mean-Square-Deviations between reference and target 
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
    >>> print calcRMSD(ens.getCoords(), ens.getCoordsets(), ens.getWeights()).round(2) # doctest: +ELLIPSIS
    [ 0.74  0.53  0.58  0.6   0.61  0.72  0.62  0.74  0.69  0.65  0.48  0.54
      ...
      0.58  0.66  0.83]
    
    """
    
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
    return _calcRMSD(ref, tar, weights)
    
def _calcRMSD(ref, tar, weights=None):
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
            
    
def superpose(mobile, target, weights=None):
    """Superpose *mobile* onto *target* to minimize the RMSD distance.
    Return *target*, after superposition, and the transformation."""
    
    t = calcTransformation(mobile, target, weights)
    result = applyTransformation(t, mobile)
    return (result, t)

def alignCoordsets(atoms, selstr='calpha', weights=None):
    """Superpose coordinate sets onto the active coordinate set.
    
    Atoms matching *selstr* will be used for calculation of transformation 
    matrix. Transformation matrix will be applied to all atoms in *atoms*,
    or its :class:`~.AtomGroup` if *atoms* is an :class:`~.AtomPointer`.
    
    By default, alpha carbon atoms are used to calculate the transformations.
    
    Optionally, atomic *weights* can be passed for weighted superposition."""
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must have type Atomic, not {0:s}'
                        .format(type(atoms)))
    if not isinstance(selstr, str):
        raise TypeError('selstr must have type str, not {0:s}'
                        .format(type(selstr)))
    n_csets = atoms.numCoordsets()
    if n_csets < 2:
        LOGGER.warning('{0:s} contains only one coordinate set, '
                       'superposition not performed.'.format(str(atoms)))
        return None
    
    acsi = atoms.getACSIndex()
    if isinstance(atoms, AtomGroup):
        ag = atoms
    else: 
        ag = atoms.getAtomGroup()
    agacsi = ag.getACSIndex()
    tar = atoms.select(selstr)
    if tar is None:
        raise ValueError("selstr '{0:s}' did not match any atoms"
                         .format(selstr))
    mob = AtomSubset(ag, tar.getIndices(), 0)
    assert tar.getACSIndex() == acsi
    for i in range(n_csets):
        if i == acsi:
            continue
        mob.setACSIndex(i)
        ag.setACSIndex(i)
        calcTransformation(mob, tar, weights).apply(ag)
    ag.setACSIndex(agacsi)

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
    
