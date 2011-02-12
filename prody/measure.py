# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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

"""This module defines class and methods and for calculating comparing coordinate 
data and measuring quantities.

Classes
-------
  
  * :class:`Transformation`
    
Functions
---------
    
  * :func:`alignCoordsets`
  * :func:`applyTransformation`
  * :func:`buildADPMatrix`
  * :func:`calcADPAxes`
  * :func:`calcDeformVector`
  * :func:`calcDistance`
  * :func:`calcRadiusOfGyration`
  * :func:`calcRMSD`
  * :func:`calcTransformation`
  * :func:`superpose` 

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

import numpy as np
linalg = None

import prody
from prody import ProDyLogger as LOGGER


__all__ = ['Transformation', 'applyTransformation', 'alignCoordsets',
           'buildADPMatrix', 'calcADPAxes', 
           'calcDeformVector', 'calcDistance', 'calcRadiusOfGyration', 
           'calcRMSD', 'calcTransformation', 'superpose']
           
class Transformation(object):
    
    __slots__ = ['_rotation', '_translation']
    
    def __init__(self, rotation, translation):
        self._rotation = rotation
        self._translation = translation
    
    #def __repr__(self):
    #    if self._name is not None: 
    #        return '<Transformation: {0:s}>'.format(self._name)
    #    else:
    #        return object.__repr__(self)

    #def __str__(self):
    #    if self._name is not None: 
    #        return 'Transformation {0:s}'.format(self._name)
    #    else:    
    #        return ''

    #def getName(self): 
    #    """Return name of the translation"""
    #    return self._name
    #def setName(self, name):
    #    """Set name of the translation"""
    #    self._name = name
    
    def getRotation(self): 
        """Returns rotation matrix."""
        return self._rotation.copy()

    def setRotation(self, rotation):
        """Set rotation matrix."""
        if not isinstance(rotation, np.ndarray):
            raise TypeError('rotation must be an ndarray')
        elif rotation.shape != (3,3):
            raise TypeError('rotation must be a 3x3 array')
        self._rotation = rotation

    def getTranslation(self): 
        """Returns translation vector."""
        return self._translation.copy()
    
    def setTranslation(self): 
        """Set translation vector."""
        if not isinstance(translation, np.ndarray):
            raise TypeError('translation must be an ndarray')
        elif translation.shape != (3,):
            raise TypeError('translation must be an ndarray of length 3')
        self._translation = translation
    
    def ge4x4Matrix(self):
        """Returns 4x4 transformation matrix whose top left is rotation matrix
        and last column is translation vector."""
        
        fourby4 = np.eye(4)
        fourby4[:3, :3] = self._rotation
        fourby4[:3, 3] = self._translation
        return fourby4
    
    def apply(self, atoms):
        """Applies transformation to given atoms or coordinate set.
        
        :class:`AtomGroup`, :class:`Chain`, :class:`Residue`, :class:`Atom`, 
        and :class:`Selection` instances are accepted.
        If an instance of one of these is provided, it is returned after
        its active coordinate set is transformed.
        
        If a NumPy array is provided, transformed array is returned.
        
        """
        return applyTransformation(self, atoms)
    

def calcTransformation(mobile, target, weights=None):
    """Returns a :class:`Transformation` instance which, when applied to the 
    atoms in *mobile*, minimizes the weighted RMSD between *mobile* and 
    *target*.
    
    *mobile* and *target* may be NumPy coordinate arrays, or istances of 
    Molecule, AtomGroup, Chain, or Residue.
    
    """
    name = ''
    if not isinstance(mobile, np.ndarray): 
        try:
            mob = mobile.getCoordinates()
        except AttributeError:
            raise TypeError('mobile is not a coordinate array '
                            'and do not contain a coordinate set')
    else:
        mob = mobile
    if not isinstance(target, np.ndarray): 
        try:
            tar = target.getCoordinates()
        except AttributeError:
            raise TypeError('target is not a coordinate array '
                            'and do not contain a coordinate set')
    else:
        tar = target
    
    if mob.shape != tar.shape:
        raise ValueError('reference and target coordinate arrays must have same shape')
    

    if mob.shape[1] != 3:
        raise ValueError('reference and target must be 3-d coordinate arrays')
        
    t = _calcTransformation(mob, tar, weights)
    return t

def _calcTransformation(mob, tar, weights=None):
    if linalg is None: 
        prody.importLA()
    n_atoms = mob.shape[0]
    
    if weights is None:
        weights = 1
        weights_sum = n_atoms
        weights_dot = 1
    else:
        if not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape[0] != n_atoms:
            raise ValueError('lenth of weights array and coordinate arrays must be the same')
        
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
    
    # optalign
    # http://www.pymolwiki.org/index.php/Kabsch
    #E0 = np.sum( np.sum(ref_centered * ref_centered,axis=0),axis=0) + np.sum( np.sum(tar_centered * tar_centered,axis=0),axis=0)
    #reflect = float(str(float(linalg.det(U) * linalg.det(Vh))))
    #if reflect == -1.0:
    #    s[-1] = -s[-1]
    #    U[:,-1] = -U[:,-1]
    #RMSD = E0 - (2.0 * sum(s))
    #RMSD = np.sqrt(abs(RMSD / n_atoms))
    #print RMSD
    #transformation._rotation = np.dot(U, Vh)
    t = Transformation(rotation, tar_com - np.dot(mob_com, rotation))
    return t

def applyTransformation(transformation, coordinates):
    """Applies a transformation to a given coordinate set."""
    if not isinstance(coordinates, np.ndarray): 
        molecule = coordinates
        try:
            coordinates = molecule.getCoordinates()
        except AttributeError:
            raise TypeError('coordinates is not an array of coordinates '
                            'and do not contain a coordinate set')
    else:
        molecule = None
    
    if coordinates.shape[1] != 3:
        raise ValueError('coordinates must be a 3-d coordinate array')
    
    transformed = transformation._translation + np.dot(coordinates, 
                                                       transformation._rotation)
    if molecule is not None:
        molecule.setCoordinates(transformed) 
        return molecule
    else:
        return transformed

def calcDeformVector(from_atoms, to_atoms):
    """Returns deformation :class:`~prody.dynamics.Vector` from *from_atoms* 
    to *atoms_to*."""
    
    name = '"{0:s}" => "{1:s}"'.format(str(from_atoms), str(to_atoms))
    if len(name) > 30: 
        name = 'Deformation'
    array = (to_atoms.getCoordinates() - from_atoms.getCoordinates()).flatten()
    return prody.Vector(array, name)

def calcRMSD(reference, target=None, weights=None):
    """Returns Root-Mean-Square-Deviations between reference and target coordinates."""
    if not isinstance(reference, np.ndarray): 
        try:
            ref = reference.getCoordinates()
            if target is None:
                target = reference.getCoordsets()
        except AttributeError:
            raise TypeError('reference is not an array of coordinates '
                            'and do not contain a coordinate set')
    else:
        ref = reference
        
    if not isinstance(target, np.ndarray): 
        try:
            tar = target.getCoordinates()
        except AttributeError:
            raise TypeError('target is not an array of coordinates '
                            'and do not contain a coordinate set')
    else:
        tar = target
    
    if ref.shape != tar.shape[-2:]:
        raise ValueError('reference and target coordinate arrays must have same shape')
    return _calcRMSD(ref, tar, weights)
    
def _calcRMSD(ref, tar, weights=None):
    n_atoms = ref.shape[0]
    if weights is None:
        weights = 1
        weights_sum = n_atoms
    else:
        if not isinstance(weights, np.ndarray): 
            raise TypeError('weights must be an ndarray instance')
        elif weights.shape[0] != n_atoms:
            raise ValueError('lenth of weights array and coordinate arrays must be the same')
        weights_sum = weights.sum()
    if tar.ndim == 2:
        return np.sqrt(((ref-tar) ** 2).sum() / weights_sum)
    else:
        return np.sqrt(((ref-tar) ** 2).sum(2).sum(1) / weights_sum)
    
def superpose(mobile, target, weights=None):
    """Superpose *mobile* onto *target* to minimize the RMSD distance."""
    t = calcTransformation(mobile, target, weights)
    result = applyTransformation(t, mobile)
    return (result, t) 

def calcDistance(one, two):
    """Return the Euclidean distance between *one* and *two*.
    
    Arguments may be :class:`Atom` instances or NumPy arrays. Shape 
    of numpy arrays must be ([M,]N,3), where M is number of coordinate sets
    and N is the number of atoms.
    
    """
    if not isinstance(one, np.ndarray):
        one = one.getCoordinates()
    if not isinstance(two, np.ndarray):
        two = two.getCoordinates()
    if one.shape != two.shape:
        raise ValueError('shape of coordinates must be the same')
    #if one.shape[-2:] != (1,3):
    #    raise ValueError('shape of coordinates must be ([M,]1,3)')
    
    return np.sqrt(np.power(one - two, 2).sum(axis=-1))
    
def alignCoordsets(atoms, selstr='calpha', weights=None):
    """Superpose coordinate sets onto the active coordinate set.
    
    .. versionadded:: 0.5
    
    Atoms matching *selstr* will be used for calculation of transformation 
    matrix. Transformation matrix will be applied to all atoms in *atoms*,
    or its :class:`~prody.atomics.AtomGroup` if *atoms* is an 
    :class:`~prody.atomics.AtomPointer`.
    
    By default, alpha carbon atoms are used to calculate the transformations.
    
    Optionally, atomic *weights* can be passed for weighted superposition.
        
    """
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must have type Atomic, not {0:s}'.format(type(atoms)))
    if not isinstance(selstr, str):
        raise TypeError('selstr must have type str, not {0:s}'.format(type(selstr)))
    n_coordsets = atoms.getNumOfCoordsets()
    if n_coordsets < 2:
        LOGGER.warning('{0:s} contains only one coordinate set, superposition not performed.'.format(str(atoms)))
        return None
    
    acsi = atoms.getActiveCoordsetIndex()
    if isinstance(atoms, prody.AtomGroup):
        ag = atoms
    else: 
        ag = atoms.getAtomGroup()
    agacsi = ag.getActiveCoordsetIndex()
    tar = atoms.select(selstr)
    mob = prody.AtomSubset(ag, tar.getIndices(), 0)
    assert tar.getActiveCoordsetIndex() == acsi
    for i in range(n_coordsets):
        if i == acsi:
            continue
        mob.setActiveCoordsetIndex(i)
        ag.setActiveCoordsetIndex(i)
        calcTransformation(mob, tar, weights).apply(ag)
    ag.setActiveCoordsetIndex(agacsi)

    
def calcAngle():
    pass

def calcDihedral():
    pass

def calcRadiusOfGyration(coords, weights=None):
    """Calculate radius of gyration for a set of coordinates or atoms."""
    if isinstance(coords, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap)):
        coords = coords.getCoordinates()
    if not isinstance(coords, np.ndarray):
        raise TypeError('coords must be a array or atomic')
    elif not coords.ndim in (2, 3):
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
    
    :arg atoms: a ProDy object for handling atomic data
    :type atoms: prody.atomic.Atomic
    :kwarg fract: For an atom, if the fraction of anisotropic displacement 
        explained by its largest axis is less than given value, 
        all axes for that atom will be set to zero. Values
        larger than 0.33 and smaller than 1.0 are accepted. 
    :type fract: float
    :kwarg ratio: For an atom, if the ratio of the largest principal axis to 
        the smallest principal axis is smaller than given value, 
        all principal axes for that atom will be set to zero. Values
        greater than 1 are accepted.  
    :type ratio: float
    
    Keyword arguments *fract* or *ratio* can be used to set principal axes
    for atoms showing relatively lower degree of anisotropy.
    
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
    
    if linalg is None: 
        prody.importLA()
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'.type(atoms))
    anisous = atoms.getAnisoTempFactors()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.getNumOfAtoms()

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
    elif 'ratio' in kwargs:  
        ratio = float(kwargs['ratio'])
        assert ratio > 1.0, 'ratio must be > 1.0'
        stddevs = stddevs[:, [2,1,0]]
        stddevs[ stddevs[:,2].flatten() == 0, 2 ] = 0.00000001
        torf = stddevs[:,0] / stddevs[:,2] > ratio
    if torf is not None:
        torf =np.tile(torf.reshape((n_atoms,1)), (1,3)).reshape((n_atoms*3, 1))
        axes = axes * torf
    return axes
        
   
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
        raise TypeError('atoms must be of type Atomic, not {0:s}'.type(atoms))
    anisous = atoms.getAnisoTempFactors()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.getNumOfAtoms()
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
        
    

