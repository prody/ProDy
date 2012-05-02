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
and measuring quantities."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import abs, mod, ndarray, power, sqrt, array, zeros, arccos, cross
from numpy import sign, tile, concatenate, pi

from prody.atomic import Atomic, Residue, Atom, AtomGroup
from prody.tools import importLA, checkCoords
from prody import LOGGER
import prody

__all__ = ['buildDistMatrix', 'calcDistance', 
           'calcCenter', 'calcGyradius', 'calcAngle', 
           'calcDihedral', 'calcOmega', 'calcPhi', 'calcPsi',
           'calcDeformVector', 
           'buildADPMatrix', 'calcADPAxes', 'calcADPs',
           'pickCentral']
           
RAD2DEG = 180 / pi

DISTMAT_FORMATS = set(['sym', 'rcd', 'dist'])

def buildDistMatrix(atoms1, atoms2=None, **kwargs):
    """Return distance matrix.  When *atoms2* is given, a distance matrix 
    with shape ``(len(atoms1), len(atoms2))`` is built.  When *atoms2* is 
    **None**, a symmetric matrix with shape ``(len(atoms1), len(atoms1))``
    is built.  If *unitcell* array is provided, periodic boundary conditions
    will be taken into account.
    
    :arg atoms1: atom or coordinate data
    :type atoms1: :class:`.Atomic`, :class:`.Frame`, :class:`numpy.ndarray`
    
    :arg atoms2: atom or coordinate data
    :type atoms2: :class:`.Atomic`, :class:`.Frame`, :class:`numpy.ndarray`
    
    :arg unitcell: orthorhombic unitcell dimension array with shape ``(3,)``
    :type unitcell: :class:`numpy.ndarray`
    
    :arg format: format of the resulting array, one of ``'sym'`` (symmetric 
        matrix, default), ``'rcd'`` (arrays of row indices, column indices, 
        and distances), or ``'dist'`` (only array of distances) 
    :type format: bool"""
    
    
    if not isinstance(atoms1, ndarray):
        try:
            atoms1 = atoms1._getCoords()
        except AttributeError:
            raise TypeError('atoms1 must be Atomic instance or an array')
    if atoms2 is None:    
        symmetric = True
        atoms2 = atoms1
    else:
        symmetric = False
        if not isinstance(atoms2, ndarray):
            try:
                atoms2 = atoms2._getCoords()
            except AttributeError:
                raise TypeError('atoms2 must be Atomic instance or an array')
    if atoms1.shape[-1] != 3 or atoms2.shape[-1] != 3:
        raise ValueError('one and two must have shape ([M,]N,3)')
        
    unitcell = kwargs.get('unitcell')
    if unitcell is not None:
        if not isinstance(unitcell, ndarray): 
            raise TypeError('unitcell must be an array')
        elif unitcell.shape != (3,):
            raise ValueError('unitcell.shape must be (3,)')

    dist = zeros((len(atoms1), len(atoms2)))
    if symmetric:
        format = kwargs.get('format', 'sym')
        if format not in DISTMAT_FORMATS:
            raise ValueError('format must be one of sym, rcd, or arr')
        if format == 'sym':
            for i, xyz in enumerate(atoms1[:-1]):
                dist[i, i+1:] = dist[i+1:, i] = getDistance(xyz, atoms2[i+1:], 
                                                            unitcell)
        else:
            dist = concatenate([getDistance(xyz, atoms2[i+1:]) 
                                for i, xyz in enumerate(atoms1)])
            if format == 'rcd':        
                n_atoms = len(atoms1)
                rc = array([(i, j) for i in xrange(n_atoms) 
                                   for j in xrange(i + 1, n_atoms)])
                row, col = rc.T
                dist = (row, col, dist) 
                            
    else:
        for i, xyz in enumerate(atoms1):
            dist[i] = getDistance(xyz, atoms2, unitcell)
    return dist    


def calcDistance(atoms1, atoms2, **kwargs):
    """Return the Euclidean distance between *atoms1* and *atoms2*.  Arguments 
    may be :class:`~.Atomic` instances or NumPy arrays.  Shape of numpy arrays 
    must be ``([M,]N,3)``, where *M* is number of coordinate sets and *N* is 
    the number of atoms.  If *unitcell* array is provided, periodic boundary 
    conditions will be taken into account.
    
    :arg unitcell: orthorhombic unitcell dimension array with shape ``(3,)``
    :type unitcell: :class:`numpy.ndarray`"""
    
    if not isinstance(atoms1, ndarray):
        try:
            atoms1 = atoms1._getCoords()
        except AttributeError:
            raise TypeError('atoms1 must be Atomic instance or an array')
    if not isinstance(atoms2, ndarray):
        try:
            atoms2 = atoms2._getCoords()
        except AttributeError:
            raise TypeError('atoms2 must be Atomic instance or an array')
    if atoms1.shape[-1] != 3 or atoms2.shape[-1] != 3:
        raise ValueError('one and two must have shape ([M,]N,3)')
    
    unitcell = kwargs.get('unitcell')
    if unitcell is not None:
        if not isinstance(unitcell, ndarray): 
            raise TypeError('unitcell must be an array')
        elif unitcell.shape != (3,):
            raise ValueError('unitcell.shape must be (3,)')
    
    return getDistance(atoms1, atoms2, unitcell)

    
def getDistance(coords1, coords2, unitcell=None):
    
    if unitcell is None:
        return sqrt(power(coords1 - coords2, 2).sum(axis=-1))
    else:
        diff = coords1 - coords2
        diff = abs(diff, diff)
        mod(diff, unitcell, diff)
        return sqrt(power(diff, 2).sum(axis=-1))
    
    
def calcAngle(atoms1, atoms2, atoms3, radian=False):
    """Return the angle between atoms in degrees."""
    
    if not isinstance(atoms1, Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not atoms1.numAtoms() == atoms2.numAtoms() == atoms3.numAtoms():
        raise ValueError('all arguments must have same number of atoms')
    
    return getAngle(atoms1._getCoords(), atoms2._getCoords(), 
                    atoms3._getCoords(), radian)


def getAngle(coords1, coords2, coords3, radian):
    """Return bond angle in degrees."""

    v1 = coords1 - coords2
    v2 = coords3 - coords2
    
    rad = arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if radian:    
        return rad
    else:
        return rad * RAD2DEG
    

def calcDihedral(atoms1, atoms2, atoms3, atoms4, radian=False):
    """Return the dihedral angle between atoms in degrees."""
    
    if not isinstance(atoms1, Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not isinstance(atoms4, Atomic):
        raise TypeError('atoms4 must be an Atomic instance')
    if not atoms1.numAtoms() == atoms2.numAtoms() == \
           atoms3.numAtoms() == atoms4.numAtoms():
        raise ValueError('all arguments must have same number of atoms')
    
    return getDihedral(atoms1._getCoords(), atoms2._getCoords(), 
                       atoms3._getCoords(), atoms4._getCoords(), radian)
    
    
def getDihedral(coords1, coords2, coords3, coords4, radian=False):
    """Return the dihedral angle in degrees."""
    
    a1 = coords2 - coords1
    a2 = coords3 - coords2
    a3 = coords4 - coords3
    
    v1 = cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5  
    v2 = cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5  
    sign = sign((v1 * a3).sum(-1))
    rad = arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if radian:    
        return sign * rad
    else:
        return sign * rad * RAD2DEG


def calcOmega(residue, radian=False, dist=4.1):
    """Return ω (omega) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0:s} must be a Residue instance')
    next = residue.getNext()
    if not isinstance(next, Residue):
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
    
    return getDihedral(CA._getCoords(), C._getCoords(), _N._getCoords(), 
                       _CA._getCoords(), radian)


def calcPhi(residue, radian=False, dist=4.1):
    """Return φ (phi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0:s} must be a Residue instance')

    C_, N, CA, C = getPhiAtoms(residue, dist=dist)

    return getDihedral(C_._getCoords(), N._getCoords(), CA._getCoords(), 
                       C._getCoords(), radian)


def getPhiAtoms(residue, dist=4.1):
    """Return the four atoms that form the φ (phi) angle of *residue*."""
    
    prev = residue.getPrev()
    if not isinstance(prev, Residue):
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
    if CA_ is None:
        raise ValueError('{0:s} does not have CA atom'.format(str(prev)))

    if dist and dist < calcDistance(CA, CA_):
        raise ValueError('{0:s} and {1:s} does not seem to be connected'
                         .format(str(residue), str(prev)))
    
    return C_, N, CA, C


def calcPsi(residue, radian=False, dist=4.1):
    """Return ψ (psi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0:s} must be a Residue instance')
        
    N, CA, C, _N = getPsiAtoms(residue, dist=dist)
    
    return getDihedral(N._getCoords(), CA._getCoords(), C._getCoords(), 
                       _N._getCoords(), radian)

def getPsiAtoms(residue, dist=4.1):
    """Return the four atoms that form the φ (phi) angle of *residue*."""
    
    next = residue.getNext()
    if not isinstance(next, Residue):
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
                        
    return N, CA, C, _N

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
    
    if weights is not None:
        if not isinstance(weights, ndarray):
            raise TypeError('weights must be a numpy array')
        elif weights.ndim != 1:
            raise ValueError('weights must be a 1 dimensional array')
        elif weights.shape[0] != coords.shape[0]:
            raise ValueError('weights length must be equal to number of atoms')

    return getCenter(coords, weights)


def getCenter(coords, weights=None):
    
    if weights is None:
        return coords.mean(0)
    else:
        return (coords * weights).mean(0) / weights.sum()


def pickCentral(atoms, weights=None):
    """Return :class:`.Atom` that is closest to the center, which is calculated
    using :func:`calcCenter`."""
        
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')
    if isinstance(atoms, Atom):
        return atoms
    if isinstance(atoms, AtomGroup):
        ag = atoms
    else:
        ag = atoms.getAtomGroup()
    if len(atoms) == 1:
        index = atoms._getIndices()[0]
    else:
        index = getCentral(atoms._getCoords(), weights)
        if not isinstance(atoms, AtomGroup):
             index = atoms._getIndices()[index]
    return Atom(ag, index, atoms.getACSIndex())


def getCentral(coords, weights=None):
    """Return index of coordinates closest to the center."""
    
    return ((coords - getCenter(coords, weights))**2).sum(1).argmin()


def calcGyradius(atoms, weights=None):
    """Calculate radius of gyration of *atoms*."""
    
    if not isinstance(atoms, ndarray):
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
        d2sum = array(rgyr)
    return (d2sum / wsum) ** 0.5


def calcDeformVector(from_atoms, to_atoms):
    """Returns deformation from *from_atoms* to *atoms_to* as a 
    :class:`~.Vector` instance."""
    
    name = '{0:s} => {1:s}'.format(repr(from_atoms), repr(to_atoms))
    if len(name) > 30: 
        name = 'Deformation'
    arr = (to_atoms.getCoords() - from_atoms.getCoords()).flatten()
    return prody.dynamics.Vector(arr, name)
            
            
def calcADPAxes(atoms, **kwargs):
    """Return a 3Nx3 array containing principal axes defining anisotropic 
    displacement parameter (ADP, or anisotropic temperature factor) ellipsoids.
    
    :arg atoms: a ProDy object for handling atomic data
    :type atoms: :class:`~.Atomic`

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
    <NMA: ADPs (3 modes; 46 atoms)>
    >>> writeNMD( 'adp_axes.nmd', nma, calphas )
    'adp_axes.nmd'"""
    
    linalg = importLA()
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()

    axes = zeros((n_atoms*3, 3))    
    variances = zeros((n_atoms, 3))
    stddevs = zeros((n_atoms, 3))
    anisou = anisous[0]
    element = zeros((3,3))
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
        vals = vals * sign((vecs * axes[(i-1)*3:(i)*3,:]).sum(0))
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
        torf = tile(torf.reshape((n_atoms,1)), (1,3)).reshape((n_atoms*3, 1))
        axes = axes * torf
    return axes
        
        
def calcADPs(atom):
    """Calculate anisotropic displacement parameters (ADPs) from 
    anisotropic temperature factors (ATFs).
    
    *atom* must have ATF values set for ADP calculation. ADPs are returned
    as a tuple, i.e. (eigenvalues, eigenvectors)."""
    
    linalg = importLA()
    if not isinstance(atom, Atom):
        raise TypeError('atom must be of type Atom, not {0:s}'
                        .format(type(atom)))
    anisou = atom.getAnisou()
    if anisou is None:
        raise ValueError('atom does not have anisotropic temperature factors')
    element = zeros((3,3))
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
    
    >>> from prody import *
    >>> protein = parsePDB('1ejg')  
    >>> calphas = protein.select('calpha')
    >>> adp_matrix = buildADPMatrix( calphas )"""
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be of type Atomic, not {0:s}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()
    n_dof = n_atoms * 3
    adp = zeros((n_dof, n_dof))
    
    for i in range(n_atoms):
        anisou = anisous[i]
        element = zeros((3,3))
        element[0,0] = anisou[0]
        element[1,1] = anisou[1]
        element[2,2] = anisou[2]
        element[0,1] = element[1,0] = anisou[3]
        element[0,2] = element[2,0] = anisou[4]
        element[1,2] = element[2,1] = anisou[5]
        adp[i*3:(i+1)*3, i*3:(i+1)*3] = element
    return adp

