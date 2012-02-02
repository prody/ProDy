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

"""This module defines a class and methods and for comparing coordinate data 
and measuring quantities.

Classes
-------
  
  * :class:`Transformation`
    
Functions
---------
    
  * :func:`alignCoordsets`
  * :func:`applyTransformation`
  * :func:`buildADPMatrix`
  * :func:`buildKDTree`
  * :func:`calcADPAxes`
  * :func:`calcADPs`
  * :func:`calcAngle`
  * :func:`calcCenter`
  * :func:`calcDeformVector`
  * :func:`calcDihedral`
  * :func:`calcDistance`
  * :func:`calcGyradius`
  * :func:`calcOmega`
  * :func:`calcPhi`
  * :func:`calcPsi`
  * :func:`calcRMSD`
  * :func:`calcTransformation`
  * :func:`iterNeighbors`
  * :func:`moveAtoms`
  * :func:`superpose`

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np

from tools import *
import prody
LOGGER = prody.LOGGER

RAD2DEG = 180 / np.pi

__all__ = ['Transformation', 'applyTransformation', 'alignCoordsets',
           'buildADPMatrix', 'buildKDTree', 'iterNeighbors', 
           'calcADPAxes', 'calcADPs',  
           'calcDeformVector', 'calcDistance', 'calcCenter',
           'calcAngle', 'calcDihedral', 'calcOmega', 'calcPhi', 'calcPsi',
           'calcGyradius', 'calcRadiusOfGyration', 
           'calcRMSD', 'calcTransformation', 
           'moveAtoms', 'superpose']
           
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
        """Applies transformation to given atoms or coordinate set.
        
        ProDy class instances from :mod:`~prody.atomic` are accepted. Instance
        is returned after its active coordinate set is transformed.
        If a :class:`~prody.atomic.AtomPointer` is passsed, the 
        :class:`~prody.atomic.AtomGroup` that it points to is transformed. 
        
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
    class instance from :mod:`~prody.atomic`, it will be returned after 
    *transformation* is applied to its active coordinate set. If *atoms*
    is an :class:`~prody.atomic.AtomPointer` instance, *transformation* will
    be applied to the corresponding coordinate set in the associated
    :class:`~prody.atomic.AtomGroup` instance."""
    
    coords = None
    ag = None
    if isinstance(atoms, np.ndarray): 
        if atoms.shape[1] != 3:
            raise ValueError('atoms must be a 3-d coordinate array')
        coords = atoms
        atoms = None
    else:
        if isinstance(atoms, prody.AtomPointer):
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
    

def calcDeformVector(from_atoms, to_atoms):
    """Returns deformation from *from_atoms* to *atoms_to* as a 
    :class:`~prody.dynamics.Vector` instance."""
    
    name = '"{0:s}" => "{1:s}"'.format(str(from_atoms), str(to_atoms))
    if len(name) > 30: 
        name = 'Deformation'
    array = (to_atoms.getCoords() - from_atoms.getCoords()).flatten()
    return prody.Vector(array, name)

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

def calcDistance(one, two):
    """Return the Euclidean distance between *one* and *two*.
    
    Arguments may be :class:`~prody.atomic.Atomic` instances or NumPy arrays. 
    Shape of numpy arrays must be ([M,]N,3), where M is number of coordinate 
    sets and N is the number of atoms."""
    
    if not isinstance(one, np.ndarray):
        try:
            one = one.getCoords()
        except AttributeError:
            raise ValueError('one must be Atom instance or a coordinate array')
    if not isinstance(two, np.ndarray):
        try:
            two = two.getCoords()
        except AttributeError:
            raise ValueError('one must be Atom instance or a coordinate array')
    if one.shape[-1] != 3 or two.shape[-1] != 3:
        raise ValueError('one and two must have shape ([M,]N,3)')
    
    return np.sqrt(np.power(one - two, 2).sum(axis=-1))
    
def alignCoordsets(atoms, selstr='calpha', weights=None):
    """Superpose coordinate sets onto the active coordinate set.
    
    .. versionadded:: 0.5
    
    Atoms matching *selstr* will be used for calculation of transformation 
    matrix. Transformation matrix will be applied to all atoms in *atoms*,
    or its :class:`~prody.atomics.AtomGroup` if *atoms* is an 
    :class:`~prody.atomics.AtomPointer`.
    
    By default, alpha carbon atoms are used to calculate the transformations.
    
    Optionally, atomic *weights* can be passed for weighted superposition."""
    
    if not isinstance(atoms, prody.Atomic):
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
    if isinstance(atoms, prody.AtomGroup):
        ag = atoms
    else: 
        ag = atoms.getAtomGroup()
    agacsi = ag.getACSIndex()
    tar = atoms.select(selstr)
    if tar is None:
        raise ValueError("selstr '{0:s}' did not match any atoms"
                         .format(selstr))
    mob = prody.AtomSubset(ag, tar.getIndices(), 0)
    assert tar.getACSIndex() == acsi
    for i in range(n_csets):
        if i == acsi:
            continue
        mob.setACSIndex(i)
        ag.setACSIndex(i)
        calcTransformation(mob, tar, weights).apply(ag)
    ag.setACSIndex(agacsi)

def calcAngle(atoms1, atoms2, atoms3, radian=False):
    """Return the angle between atoms in degrees."""
    
    if not isinstance(atoms1, prody.Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, prody.Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, prody.Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not atoms1.numAtoms() == atoms2.numAtoms() == atoms3.numAtoms():
        raise ValueError('all arguments must have same number of atoms')
    
    coords2 = atoms2._getCoords()
    v1 = atoms1._getCoords() - coords2
    v2 = atoms3._getCoords() - coords2
    
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if radian:    
        return rad
    else:
        return rad * RAD2DEG

def calcDihedral(atoms1, atoms2, atoms3, atoms4, radian=False):
    """Return the dihedral angle between atoms in degrees."""
    
    if not isinstance(atoms1, prody.Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, prody.Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, prody.Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not isinstance(atoms4, prody.Atomic):
        raise TypeError('atoms4 must be an Atomic instance')
    if not atoms1.numAtoms() == atoms2.numAtoms() == \
           atoms3.numAtoms() == atoms4.numAtoms():
        raise ValueError('all arguments must have same number of atoms')
    
    coords2 = atoms2._getCoords()
    coords3 = atoms3._getCoords()
    a1 = coords2 - atoms1._getCoords()
    a2 = coords3 - coords2
    a3 = atoms4._getCoords() - coords3
    
    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5  
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5  
    sign = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if radian:    
        return sign * rad
    else:
        return sign * rad * RAD2DEG

def calcOmega(residue, radian=False, dist=4.1):
    """Return ω (omega) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues.  Set *dist* to none, to 
    avoid this check."""

    if not isinstance(residue, prody.Residue):
        raise TypeError('{0:s} must be a Residue instance')
    next = residue.getNext()
    if not isinstance(next, prody.Residue):
        raise ValueError('{0:s} is a terminal residue'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0:s} does not have C atom'.format(str(residue)))
    _N = next['N']
    if _N is None:
        raise ValueError('{0:s} does not have N atom'.format(str(residue)))
    _CA = next['CA']
    if _CA is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(next)))
    
    if dist and dist < calcDistance(CA, _CA):
        raise ValueError('{0:s} and {1:s} does not seem to be connected'
                         .format(str(residue), str(next)))
    
    return calcDihedral(CA, C, _N, _CA, radian)

def calcPhi(residue, radian=False, dist=4.1):
    """Return φ (phi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues.  Set *dist* to none, to 
    avoid this check."""

    if not isinstance(residue, prody.Residue):
        raise TypeError('{0:s} must be a Residue instance')
    prev = residue.getPrev()
    if not isinstance(prev, prody.Residue):
        raise ValueError('{0:s} is a terminal residue'.format(str(residue)))

    C_ = prev['C']
    if C_ is None:
        raise ValueError('{0:s} does not have C atom'.format(str(prev)))
    N = residue['N']
    if N is None:
        raise ValueError('{0:s} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0:s} does not have C atom'.format(str(residue)))
    CA_ = prev['CA']
    if C_ is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(prev)))
    if dist and dist < calcDistance(CA, CA_):
        raise ValueError('{0:s} and {1:s} does not seem to be connected'
                         .format(str(residue), str(prev)))
    
    return calcDihedral(C_, N, CA, C, radian)

def calcPsi(residue, radian=False, dist=4.1):
    """Return ψ (psi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues.  Set *dist* to none, to 
    avoid this check."""

    if not isinstance(residue, prody.Residue):
        raise TypeError('{0:s} must be a Residue instance')
    next = residue.getNext()
    if not isinstance(next, prody.Residue):
        raise ValueError('{0:s} is a terminal residue'.format(str(residue)))
    N = residue['N']
    if N is None:
        raise ValueError('{0:s} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0:s} does not have C atom'.format(str(residue)))
    _N = next['N']
    if _N is None:
        raise ValueError('{0:s} does not have N atom'.format(str(next)))
    _CA = next['CA']
    if _CA is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(next)))
    if dist and dist < calcDistance(CA, _CA):
        raise ValueError('{0:s} and {1:s} does not seem to be connected'
                         .format(str(residue), str(next)))
    
    return calcDihedral(N, CA, C, _N, radian)

def calcCenter(atoms, weights=None):
    """Return geometric center of *atoms*.  If *weights* is given it must 
    be a flat array with length equal to number of atoms.  Mass center
    of atoms can be calculated by setting weights equal to mass, i.e.
    ``weights=atoms.getMasses()``."""
    
    try: 
        coords = atoms._getCoords()
    except AttributeError:
        coords = checkCoords(atoms, 'atoms')
    except Exception as err:
        raise type(err)(err)
    
    if weights is None:
        return coords.mean(0) 
    else:
        if not isinstance(weights, np.ndarray):
            raise TypeError('weights must be a numpy array')
        elif weights.ndim != 1:
            raise ValueError('weights must be a 1 dimensional array')
        elif weights.shape[0] != coords.shape[0]:
            raise ValueError('weights length must be equal to number of atoms')
        return (coords * weights).mean(0) / weights.sum()

def moveAtoms(atoms, array):
    """Move or transform *atoms*. *array* must be :class:`numpy.ndarray`.  
    If shape of *array* is one of ``(natoms, 3)``, ``(1, 3)``, or ``(3,)``,
    *atoms* will be translated. Ff *array* is a ``(4,4)`` matrix, coordinates
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
    

        
def buildKDTree(atoms):
    """Return a KDTree built using coordinates of *atoms*.  *atoms* must be
    a ProDy object or a :class:`numpy.ndarray` with shape ``(n_atoms,3)``.  
    This function uses Biopython KDTree module."""
    
    if isinstance(atoms, np.ndarray):
        coords = checkCoords(atoms, 'atoms')
        return getKDTree(coords)
    else:
        try:
            coords = atoms._getCoords()
        except AttributeError:
            raise TypeError('invalid type for atoms')
        finally:
            return getKDTree(coords)

def getKDTree(coords):
    """Internal function to get KDTree for coordinates without any checks."""

    from KDTree import KDTree
    return KDTree(coords)
    
def iterNeighbors(atoms, radius, atoms2=None):
    """Yield pairs of *atoms* that are those within *radius* of each other,
    with the distance between them.  If *atoms2* is also provided, one atom 
    from *atoms* and another from *atoms2* will be yielded."""
    
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be an Atomic instance')
    elif not isinstance(radius, (float, int)):
        raise TypeError('radius must be an Atomic instance')
    elif radius <= 0:
        raise ValueError('radius must have a positive value')
        
    if atoms2 is None:
        if len(atoms) <= 1:
            raise ValueError('length of atoms must be more than 1')
        ag = atoms
        if not isinstance(ag, prody.AtomGroup):
            ag = ag.getAtomGroup()
            indices = atoms._getIndices()
            index = lambda i: indices[i]
        else:
            index = lambda i: i
        kdtree = getKDTree(atoms._getCoords())
        kdtree.all_search(radius)
        
        _dict = {}
        for (i, j), r in zip(kdtree.all_get_indices(), kdtree.all_get_radii()): 
             
            a1 = _dict.get(i)
            if a1 is None:      
                a1 = ag[index(i)]
                _dict[i] = a1
            a2 = _dict.get(j)
            if a2 is None:      
                a2 = ag[index(j)]
                _dict[j] = a2
            yield (a1, a2, r)   
    else:
        if len(atoms) >= len(atoms2): 
            ag = atoms
            if not isinstance(ag, prody.AtomGroup):
                ag = ag.getAtomGroup()
                indices = atoms._getIndices()
                index = lambda i: indices[i]
            else:
                index = lambda i: i
            kdtree = getKDTree(atoms._getCoords())
            
            _dict = {}
            for a2 in atoms2.iterAtoms():
                kdtree.search(a2._getCoords(), radius)
                for i, r in zip(kdtree.get_indices(), kdtree.get_radii()): 
                    a1 = _dict.get(i)
                    if a1 is None:      
                        a1 = ag[index(i)]
                        _dict[i] = a1
                    yield (a1, a2, r)   
        else:    
            ag = atoms2
            if not isinstance(ag, prody.AtomGroup):
                ag = ag.getAtomGroup()
                indices = atoms2._getIndices()
                index = lambda i: indices[i]
            else:
                index = lambda i: i
            kdtree = getKDTree(atoms2._getCoords())
            
            _dict = {}
            for a1 in atoms.iterAtoms():
                kdtree.search(a1._getCoords(), radius)
                for i, r in zip(kdtree.get_indices(), kdtree.get_radii()): 
                    a2 = _dict.get(i)
                    if a2 is None:      
                        a2 = ag[index(i)]
                        _dict[i] = a2
                    yield (a1, a2, r)   



def calcRadiusOfGyration(coords, weights=None):
    """Deprecated, use :meth:`calcGyradius`."""
    
    prody.deprecate('calcRadiusOfGyration', 'calcGyradius')
    return calcGyradius(coords, weights)
    
def calcGyradius(atoms, weights=None):
    """Calculate radius of gyration of *atoms*."""
    
    if not isinstance(atoms, np.ndarray):
        try:
            coords = atoms._getCoords()
        except AttributeError:
            raise TypeError('atoms must have atomic coordinate data')
    else:
        coords = atoms
        if not coords.ndim in (2, 3):
            raise ValueError('coords may be a 2 or 3 dimentional array')
        elif coords.shape[-1] != 3:
            raise ValueError('coords must have shape ([n_coordsets,]n_atoms,3)')
    if weights is not None:
        weights = weights.flatten()
        if len(weights) != coords.shape[-2]:
            raise ValueError('length of weights must match number of atoms')
        wsum = weights.sum()
    else:
        wsum = coords.shape[-2]
        
    if coords.ndim == 2:
        if weights is None:
            com = coords.mean(0)
            d2sum = ((coords - com)**2).sum()
        else:
            
            com = (coords * weights).mean(0) / wsum
            d2sum = (((coords - com)**2).sum(1) * weights).sum()
    else:
        rgyr = []
        for coords in coords:        
            if weights is None:
                com = coords.mean(0)
                d2sum = ((coords - com)**2).sum()
                rgyr.append(d2sum)
            else:
                
                com = (coords * weights).mean(0) / wsum
                d2sum = (((coords - com)**2).sum(1) * weights).sum()
                rgyr.append(d2sum)
        d2sum = np.array(rgyr)
    return (d2sum / wsum) ** 0.5
            
def calcADPAxes(atoms, **kwargs):
    """Return a 3Nx3 array containing principal axes defining anisotropic 
    displacement parameter (ADP, or anisotropic temperature factor) ellipsoids.
    
    .. versionadded:: 0.5.3
    
    .. versionchanged:: 0.7
       *ratio2* optional keyword argument is added.
    
    .. versionchanged:: 0.7.1
       The interpretation of *ratio*, *ratio2*, and *ratio3* keywords
       have changed (see the new definitions below).

    :arg atoms: a ProDy object for handling atomic data
    :type atoms: prody.atomic.Atomic

    :kwarg fract: For an atom, if the fraction of anisotropic displacement 
        explained by its largest axis/eigenvector is less than given value, 
        all axes for that atom will be set to zero. Values
        larger than 0.33 and smaller than 1.0 are accepted. 
    :type fract: float

    :kwarg ratio2: For an atom, if the ratio of the second-largest eigenvalue 
        to the largest eigenvalue axis less than or equal to the given value, 
        all principal axes for that atom will be returned. 
        Values less than 1 and greater than 0 are accepted.  
    :type ratio2: float

    :kwarg ratio3: For an atom, if the ratio of the smallest eigenvalue 
        to the largest eigenvalue is less than or equal to the given value, 
        all principal axes for that atom will be returned. 
        Values less than 1 and greater than 0 are accepted.  
    :type ratio3: float

    :kwarg ratio: Same as *ratio3*.  
    :type ratio: float

    
    Keyword arguments *fract*, *ratio3*, or *ratio3* can be used to set 
    principal axes to 0 for atoms showing relatively lower degree of 
    anisotropy.
    
    3Nx3 axis contains N times 3x3 matrices, one for each given atom. Columns
    of these 3x3 matrices are the principal axes which are weighted by
    square root of their eigenvalues. The first columns correspond to largest
    principal axes.
    
    The direction of the principal axes for an atom is determined based on the 
    correlation of the axes vector with the principal axes vector of the 
    previous atom.  
    
    >>> from prody import *
    >>> protein = parsePDB('1ejg')  
    >>> calphas = protein.select('calpha')
    >>> adp_axes = calcADPAxes( calphas )
    >>> adp_axes.shape
    (138, 3)
    
    These can be written in NMD format as follows:
        
    >>> nma = NMA('ADPs')
    >>> nma.setEigens(adp_axes)
    >>> nma
    <NMA: ADPs (3 modes, 46 atoms)>
    >>> writeNMD( 'adp_axes.nmd', nma, calphas )
    'adp_axes.nmd'
    
    """
    
    linalg = importLA()
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()

    axes = np.zeros((n_atoms*3, 3))    
    variances = np.zeros((n_atoms, 3))
    stddevs = np.zeros((n_atoms, 3))
    anisou = anisous[0]
    element = np.zeros((3,3))
    element[0,0] = anisou[0]
    element[1,1] = anisou[1]
    element[2,2] = anisou[2]
    element[0,1] = element[1,0] = anisou[3]
    element[0,2] = element[2,0] = anisou[4]
    element[1,2] = element[2,1] = anisou[5]
    vals, vecs = linalg.eigh(element)
    # If negative eigenvalues are found (when ADP matrix is not positive 
    # definite) set them to 0   
    vals[ vals < 0 ] = 0
    variances[0] = vals
    vals = vals**0.5
    stddevs[0] = vals
    axes[0:3,:] = vals * vecs

    for i in range(1, n_atoms):
        anisou = anisous[i]
        element[0,0] = anisou[0]
        element[1,1] = anisou[1]
        element[2,2] = anisou[2]
        element[0,1] = element[1,0] = anisou[3]
        element[0,2] = element[2,0] = anisou[4]
        element[1,2] = element[2,1] = anisou[5]
        vals, vecs = linalg.eigh(element)
        # If negative eigenvalues are found (when ADP matrix is not positive 
        # definite) set them to 0   
        vals[ vals < 0 ] = 0
        variances[i] = vals
        vals = vals**0.5
        stddevs[i] = vals
        # Make sure the direction that correlates with the previous atom
        # is selected 
        vals = vals * np.sign((vecs * axes[(i-1)*3:(i)*3,:]).sum(0))
        axes[i*3:(i+1)*3,:] = vals * vecs
    # Resort the columns before returning array
    axes = axes[:, [2,1,0]]
    torf = None
    if 'fract' in kwargs:
        fract = float(kwargs['fract'])
        assert 0.33 < fract < 1.0, 'fract must be > 0.33 and < 1.0'
        variances = variances[:, [2,1,0]]
        torf = variances[:,0] / variances.sum(1) > fract
    elif 'ratio' in kwargs or 'ratio3' in kwargs or 'ratio2' in kwargs: 
        if 'ratio2' in kwargs:
            ratio = float(kwargs['ratio2'])
            assert 0 < ratio < 1.0, 'ratio2 must be > 0 and < 1.0'
            dim = 1
        else:
            ratio = float(kwargs.get('ratio', kwargs.get('ratio3')))
            assert 0 < ratio < 1.0, 'ratio or ratio3 must be > 0 and < 1.0'
            dim = 2
        variances = variances[:, [2,1,0]]
        torf = variances[:,dim] / variances[:,0] <= ratio
    if torf is not None:
        torf =np.tile(torf.reshape((n_atoms,1)), (1,3)).reshape((n_atoms*3, 1))
        axes = axes * torf
    return axes
        
def calcADPs(atom):
    """Calculate anisotropic displacement parameters (ADPs) from 
    anisotropic temperature factors (ATFs).
    
    .. versionadded:: 0.8
    
    *atom* must have ATF values set for ADP calculation. ADPs are returned
    as a tuple, i.e. (eigenvalues, eigenvectors)."""
    
    linalg = importLA()
    if not isinstance(atom, prody.Atom):
        raise TypeError('atom must be of type Atom, not {0:s}'
                        .format(type(atom)))
    anisou = atom.getAnisou()
    if anisou is None:
        raise ValueError('atom does not have anisotropic temperature factors')
    element = np.zeros((3,3))
    element[0,0] = anisou[0]
    element[1,1] = anisou[1]
    element[2,2] = anisou[2]
    element[0,1] = element[1,0] = anisou[3]
    element[0,2] = element[2,0] = anisou[4]
    element[1,2] = element[2,1] = anisou[5]
    vals, vecs = linalg.eigh(element)
    return vals[[2,1,0]], vecs[:, [2,1,0]] 
   
def buildADPMatrix(atoms):
    """Return a 3Nx3N symmetric matrix containing anisotropic displacement
    parameters (ADPs) along the diagonal as 3x3 super elements.
    
    .. versionadded:: 0.5.3
    
    >>> from prody import *
    >>> protein = parsePDB('1ejg')  
    >>> calphas = protein.select('calpha')
    >>> adp_matrix = buildADPMatrix( calphas )
    """
    
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()
    n_dof = n_atoms * 3
    adp = np.zeros((n_dof, n_dof))
    
    for i in range(n_atoms):
        anisou = anisous[i]
        element = np.zeros((3,3))
        element[0,0] = anisou[0]
        element[1,1] = anisou[1]
        element[2,2] = anisou[2]
        element[0,1] = element[1,0] = anisou[3]
        element[0,2] = element[2,0] = anisou[4]
        element[1,2] = element[2,1] = anisou[5]
        adp[i*3:(i+1)*3, i*3:(i+1)*3] = element
    return adp
        
    

