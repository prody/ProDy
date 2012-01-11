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
  * :func:`calcADPAxes`
  * :func:`calcADPs`
  * :func:`calcDeformVector`
  * :func:`calcDistance`
  * :func:`calcRadiusOfGyration`
  * :func:`calcRMSD`
  * :func:`calcTransformation`
  * :func:`superpose` 

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np
linalg = None

import prody
LOGGER = prody.LOGGER


__all__ = ['Transformation', 'applyTransformation', 'alignCoordsets',
           'buildADPMatrix', 'calcADPAxes', 'calcADPs',  
           'calcDeformVector', 'calcDistance', 'calcRadiusOfGyration', 
           'calcRMSD', 'calcTransformation', 'superpose']
           
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
    if linalg is None:
        prody.importLA()
    
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
    tar_com = tar.mean(0)
    tar_org_T = (tar - tar_com).T
    if linalg is None:
        prody.importLA()
    svd = linalg.svd
    det = linalg.det
    dot = np.dot
    array = np.array
    sign = np.sign
    LOGGER.progress('Superposing ', len(mobs))
    for i, mob in enumerate(mobs):      
        mob_com = mob.mean(0)
        mob_org = mob - mob_com

        matrix = dot(tar_org_T, mob_org)
        U, s, Vh = svd(matrix)
        Id = array([ [1, 0, 0], [0, 1, 0], [0, 0, sign(det(matrix))] ])
        rotation = dot(Vh.T, dot(Id, U.T))

        if movs is None:
            mobs[i] = dot(mob_org, rotation) 
            mobs[i] += tar_com 
        else:
            movs[i] = dot(movs[i], rotation) 
            movs[i] += (tar_com - dot(mob_com, rotation))
    
        LOGGER.update(i)
    LOGGER.clear()

def _superpose(mob, tar, weights=None, mov=None):
    tar_com = tar.mean(0)
    tar_org = tar - tar_com

    if linalg is None:
        prody.importLA()
    mob_com = mob.mean(0)
    mob_org = mob - mob_com
    matrix = np.dot(tar_org.T, mob_org)

    U, s, Vh = linalg.svd(matrix)
    Id = np.array([ [1, 0, 0], 
                    [0, 1, 0], 
                    [0, 0, np.sign(linalg.det(matrix))] ])
    rotation = np.dot(Vh.T, np.dot(Id, U.T))

    if mov is None:
        mob[:] = np.dot(mob_org, rotation) + tar_com 
    else:
        mov[:] = np.dot(mov, rotation) + (tar_com - np.dot(mob_com, rotation))


def applyTransformation(transformation, coords):
    """Returns *coords* after applying *transformation*. If *coords* is a 
    class instance from :mod:`~prody.atomic`, it will be returned after 
    *transformation* is applied to its active coordinate set. If *coords*
    is an :class:`~prody.atomic.AtomPointer` instance, *transformation* will
    be applied to the corresponding coordinate set in the associated
    :class:`~prody.atomic.AtomGroup` instance."""
    
    atoms = None
    ag = None
    if isinstance(coords, np.ndarray): 
        if coords.shape[1] != 3:
            raise ValueError('coordinates must be a 3-d coordinate array')
    else:
        atoms = coords
        if isinstance(atoms, prody.AtomPointer):
            ag = atoms.getAtomGroup()
            acsi = ag.getACSI()
            ag.setACSI(atoms.getACSI())
            coords = ag._getCoords()
        else:
            try:
                coords = atoms._getCoords()
            except AttributeError:
                raise TypeError('coords is not an array of coordinates '
                                'and do not contain a coordinate set')
    if atoms is None:
        return _applyTransformation(transformation, coords)
    else:
        if ag is None:
            atoms.setCoords(_applyTransformation(transformation, coords))
        else: 
            ag.setCoords(_applyTransformation(transformation, coords))
            ag.setACSI(acsi)
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
    
    acsi = atoms.getACSI()
    if isinstance(atoms, prody.AtomGroup):
        ag = atoms
    else: 
        ag = atoms.getAtomGroup()
    agacsi = ag.getACSI()
    tar = atoms.select(selstr)
    if tar is None:
        raise ValueError("selstr '{0:s}' did not match any atoms"
                         .format(selstr))
    mob = prody.AtomSubset(ag, tar.getIndices(), 0)
    assert tar.getACSI() == acsi
    for i in range(n_csets):
        if i == acsi:
            continue
        mob.setACSI(i)
        ag.setACSI(i)
        calcTransformation(mob, tar, weights).apply(ag)
    ag.setACSI(agacsi)

    
def calcAngle():
    pass

def calcDihedral():
    pass

def calcRadiusOfGyration(coords, weights=None):
    """Calculate radius of gyration for a set of coordinates or atoms."""
    
    if isinstance(coords, (prody.AtomGroup, prody.AtomSubset, prody.AtomMap, 
                           prody.ConformationBase)):
        coords = coords._getCoords()
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
    
    if linalg is None: 
        prody.importLA()
    if not isinstance(atoms, prody.Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'.type(atoms))
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
    
    if linalg is None: 
        prody.importLA()
    if not isinstance(atom, prody.Atom):
        raise TypeError('atom must be of type Atom, not {0:s}'.type(atom))
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
        raise TypeError('atoms must be of type Atomic, not {0:s}'.type(atoms))
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
        
    

