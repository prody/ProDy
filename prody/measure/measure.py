# -*- coding: utf-8 -*-
"""This module defines a class and methods and for comparing coordinate data
and measuring quantities."""
from numbers import Integral, Number

from numpy import ndarray, power, sqrt, array, zeros, arccos, dot
from numpy import sign, tile, concatenate, pi, cross, subtract, var
from numpy import unique, where, divide, power

from prody.atomic import Atomic, Residue, Atom, extendAtomicData
from prody.kdtree import KDTree
from prody.utilities import importLA, solveEig, checkCoords, getDistance, getCoords
from prody.utilities import calcTree, findSubgroups
from prody import LOGGER, PY2K

if PY2K:
    range = xrange

__all__ = ['buildDistMatrix', 'calcDistance', 'calcGyradius',
           'calcCenter', 'calcAngle', 'calcDihedral',
           'getCenter', 'getAngle', 'getDihedral',
           'calcOmega', 'calcPhi', 'calcPsi',
           'calcMSF', 'calcRMSF', 'calcDeformVector',
           'buildADPMatrix', 'calcADPAxes', 'calcADPs',
           'pickCentral', 'pickCentralAtom', 'pickCentralConf', 'getWeights',
           'calcInertiaTensor', 'calcPrincAxes', 'calcDistanceMatrix',
           'assignBlocks']

RAD2DEG = 180 / pi

DISTMAT_FORMATS = set(['mat', 'rcd', 'arr'])


def buildDistMatrix(atoms1, atoms2=None, unitcell=None, format='mat', seqsep=None):
    """Returns distance matrix.  When *atoms2* is given, a distance matrix
    with shape ``(len(atoms1), len(atoms2))`` is built.  When *atoms2* is
    **None**, a symmetric matrix with shape ``(len(atoms1), len(atoms1))``
    is built.  If *unitcell* array is provided, periodic boundary conditions
    will be taken into account.

    :arg atoms1: atom or coordinate data
    :type atoms1: :class:`.Atomic`, :class:`numpy.ndarray`

    :arg atoms2: atom or coordinate data
    :type atoms2: :class:`.Atomic`, :class:`numpy.ndarray`

    :arg unitcell: orthorhombic unitcell dimension array with shape ``(3,)``
    :type unitcell: :class:`numpy.ndarray`

    :arg format: format of the resulting array, one of ``'mat'`` (matrix,
        default), ``'rcd'`` (arrays of row indices, column indices, and
        distances), or ``'arr'`` (only array of distances)
    :type format: bool
    
    :arg seqsep: if provided, distances will only be measured between atoms 
        with resnum differences that are greater than or equal to seqsep. 
    :type seqsep: int
    """

    spacing = 1
    if seqsep is not None:
        spacing = seqsep - 1

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

        if atoms2.ndim == 1:
            atoms2 = atoms2.reshape((1,3))

    if atoms1.shape[-1] != 3 or atoms2.shape[-1] != 3:
        raise ValueError('one and two must have shape ([M,]N,3)')

    if unitcell is not None:
        if not isinstance(unitcell, ndarray):
            raise TypeError('unitcell must be an array')
        elif unitcell.shape != (3,):
            raise ValueError('unitcell.shape must be (3,)')

    dist = zeros((len(atoms1), len(atoms2)))
    if symmetric:
        if format not in DISTMAT_FORMATS:
            raise ValueError('format must be one of mat, rcd, or arr')
        if format == 'mat':
            for i, xyz in enumerate(atoms1[:-1]):
                dist[i, i+spacing:] = dist[i+spacing:, i] = getDistance(xyz, atoms2[i+spacing:],
                                                            unitcell)
        else:
            dist = concatenate([getDistance(xyz, atoms2[i+spacing:], unitcell)
                                for i, xyz in enumerate(atoms1)])
            if format == 'rcd':
                n_atoms = len(atoms1)
                rc = array([(i, j) for i in range(n_atoms)
                            for j in range(i + 1, n_atoms)])
                row, col = rc.T
                dist = (row, col, dist)

    else:
        for i, xyz in enumerate(atoms1):
            dist[i] = getDistance(xyz, atoms2, unitcell)
    return dist


def calcDistance(atoms1, atoms2, unitcell=None):
    """Returns the Euclidean distance between *atoms1* and *atoms2*.  Arguments
    may be :class:`~.Atomic` instances or NumPy arrays.  Shape of numpy arrays
    must be ``([M,]N,3)``, where *M* is number of coordinate sets and *N* is
    the number of atoms.  If *unitcell* array is provided, periodic boundary
    conditions will be taken into account.

    :arg atoms1: atom or coordinate data
    :type atoms1: :class:`.Atomic`, :class:`numpy.ndarray`

    :arg atoms2: atom or coordinate data
    :type atoms2: :class:`.Atomic`, :class:`numpy.ndarray`

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
        raise ValueError('atoms1 and atoms2 must have shape ([M,]N,3)')

    if unitcell is not None:
        if not isinstance(unitcell, ndarray):
            raise TypeError('unitcell must be an array')
        elif unitcell.shape != (3,):
            raise ValueError('unitcell.shape must be (3,)')

    return getDistance(atoms1, atoms2, unitcell)


def calcAngle(atoms1, atoms2, atoms3, radian=False):
    """Returns the angle between atoms in degrees."""

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


def getAngle(coords1, coords2, coords3, radian=False):
    """Returns bond angle in degrees unless ``radian=True``"""

    v1 = coords1 - coords2
    v2 = coords3 - coords2

    rad = arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if radian:
        return rad
    else:
        return rad * RAD2DEG


def calcDihedral(atoms1, atoms2, atoms3, atoms4, radian=False):
    """Returns the dihedral angle between atoms in degrees."""

    if not isinstance(atoms1, Atomic):
        raise TypeError('atoms1 must be an Atomic instance')
    if not isinstance(atoms2, Atomic):
        raise TypeError('atoms2 must be an Atomic instance')
    if not isinstance(atoms3, Atomic):
        raise TypeError('atoms3 must be an Atomic instance')
    if not isinstance(atoms4, Atomic):
        raise TypeError('atoms4 must be an Atomic instance')
    if not (atoms1.numAtoms() == atoms2.numAtoms() ==
            atoms3.numAtoms() == atoms4.numAtoms()):
        raise ValueError('all arguments must have same number of atoms')

    return getDihedral(atoms1._getCoords(), atoms2._getCoords(),
                       atoms3._getCoords(), atoms4._getCoords(), radian)


def getDihedral(coords1, coords2, coords3, coords4, radian=False):
    """Returns the dihedral angle in degrees unless ``radian=True``."""

    a1 = coords2 - coords1
    a2 = coords3 - coords2
    a3 = coords4 - coords3

    v1 = cross(a1, a2)
    v1 = divide(v1, power((v1 * v1).sum(-1), 0.5).reshape(-1,1))
    v2 = cross(a2, a3)
    v2 = divide(v2, power((v2 * v2).sum(-1), 0.5).reshape(-1,1))
    porm = sign((v1 * a3).sum(-1))
    rad = arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not all(porm == 0):
        rad = rad * porm

    if rad.shape[0] == 1:
        rad = rad[0]
        
    if radian:
        return rad
    else:
        return rad * RAD2DEG


def calcOmega(residue, radian=False, dist=4.1):
    """Returns ω (omega) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0} must be a Residue instance')

    next = residue.getNext()
    if not isinstance(next, Residue):
        raise ValueError('{0} is a terminal residue'.format(str(residue)))

    CA = residue['CA']
    if CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0} does not have C atom'.format(str(residue)))
    _N = next['N']
    if _N is None:
        raise ValueError('{0} does not have N atom'.format(str(residue)))
    _CA = next['CA']
    if _CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(next)))

    if dist and dist < calcDistance(CA, _CA):
        raise ValueError('{0} and {1} does not seem to be connected'
                         .format(str(residue), str(next)))

    return getDihedral(CA._getCoords(), C._getCoords(), _N._getCoords(),
                       _CA._getCoords(), radian)


def calcPhi(residue, radian=False, dist=4.1):
    """Returns φ (phi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0} must be a Residue instance')

    C_, N, CA, C = getPhiAtoms(residue, dist=dist)

    return getDihedral(C_._getCoords(), N._getCoords(), CA._getCoords(),
                       C._getCoords(), radian)


def getPhiAtoms(residue, dist=4.1):
    """Returns the four atoms that form the φ (phi) angle of *residue*."""

    prev = residue.getPrev()
    try:
        isaa = prev.isaminoacid
    except AttributeError:
        raise ValueError('{0} is a terminal residue'.format(str(residue)))

    if not isaa:
        raise ValueError('{0} is not an amino acid'.format(str(prev)))

    C_ = prev['C']
    if C_ is None:
        raise ValueError('{0} does not have C atom'.format(str(prev)))
    N = residue['N']
    if N is None:
        raise ValueError('{0} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0} does not have C atom'.format(str(residue)))
    CA_ = prev['CA']
    if CA_ is None:
        raise ValueError('{0} does not have CA atom'.format(str(prev)))

    if dist and dist < calcDistance(CA, CA_):
        raise ValueError('{0} and {1} does not seem to be connected'
                         .format(str(residue), str(prev)))

    return C_, N, CA, C


def calcPsi(residue, radian=False, dist=4.1):
    """Returns ψ (psi) angle of *residue* in degrees.  This function checks
    the distance between Cα atoms of two residues and raises an exception if
    the residues are disconnected.  Set *dist* to **None**, to avoid this."""

    if not isinstance(residue, Residue):
        raise TypeError('{0} must be a Residue instance')

    N, CA, C, _N = getPsiAtoms(residue, dist=dist)

    return getDihedral(N._getCoords(), CA._getCoords(), C._getCoords(),
                       _N._getCoords(), radian)


def getPsiAtoms(residue, dist=4.1):
    """Returns the four atoms that form the φ (phi) angle of *residue*."""

    next = residue.getNext()
    try:
        isaa = next.isaminoacid
    except AttributeError:
        raise ValueError('{0} is a terminal residue'.format(str(residue)))

    if not isaa:
        raise ValueError('{0} is not an amino acid'.format(str(next)))

    N = residue['N']
    if N is None:
        raise ValueError('{0} does not have N atom'.format(str(residue)))
    CA = residue['CA']
    if CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(residue)))
    C = residue['C']
    if C is None:
        raise ValueError('{0} does not have C atom'.format(str(residue)))
    _N = next['N']
    if _N is None:
        raise ValueError('{0} does not have N atom'.format(str(next)))
    _CA = next['CA']
    if _CA is None:
        raise ValueError('{0} does not have CA atom'.format(str(next)))
    if dist and dist < calcDistance(CA, _CA):
        raise ValueError('{0} and {1} does not seem to be connected'
                         .format(str(residue), str(next)))

    return N, CA, C, _N


def calcCenter(atoms, weights=None):
    """Returns geometric center of *atoms*.  If *weights* is given it must
    be a flat array with length equal to number of atoms.  Mass center of
    atoms can be calculated by setting weights equal to atom masses, i.e.
    ``weights=atoms.getMasses()``."""

    try:
        coords = atoms._getCoords()
    except AttributeError:
        try:
            coords = atoms.getCoords()
        except AttributeError:
            coords = atoms
            try:
                checkCoords(coords, csets=True, dtype=None, name='atoms')
            except TypeError:
                raise TypeError('atoms must be an Atomic instance')

    if weights is not None:
        try:
            ndim, shape = weights.ndim, weights.shape
        except AttributeError:
            raise TypeError('weights must be a numpy array')
        else:
            if shape[0] != coords.shape[-2]:
                raise ValueError('weights.shape[0] must be equal to number of '
                                 'atoms')
            if ndim != 2:
                try:
                    weights = weights.reshape((shape[0], 1))
                except ValueError:
                    raise ValueError('weights.shape must be a (n_atoms, 1)')

    return getCenter(coords, weights)


def getCenter(coords, weights=None):

    if weights is None:
        return coords.mean(-2)
    else:
        return (coords * weights).sum(-2) / weights.sum()

def getWeights(pdb):

    weights = zeros(len(pdb))
    for i in range(len(pdb)):
        weights[i]=pdb[i].getMassess()
    return weights


def pickCentral(obj, weights=None):
    """Returns :class:`.Atom` or :class:`.Conformation` that is closest to the
    center of *obj*, which may be an :class:`.Atomic` or :class:`.Ensemble`
    instance.  See also :func:`pickCentralAtom`, and :func:`pickCentralConf`
    functions."""

    try:
        obj.getACSIndex()
    except AttributeError:
        try:
            obj.numConfs()
        except AttributeError:
            raise TypeError('obj must be an Atomic or Ensemble instance')
        else:
            return pickCentralConf(obj, weights)
    else:
        return pickCentralAtom(obj, weights)


def pickCentralAtom(atoms, weights=None):
    """Returns :class:`.Atom` that is closest to the center, which is calculated
    using :func:`calcCenter`."""

    try:
        acsi, coords = atoms.getACSIndex(), atoms._getCoords()
    except AttributeError:
        raise TypeError('atoms must be an Atomic instance')
    else:
        if coords is None:
            raise ValueError('coordinates are not set')

    try:
        ag = atoms.getAtomGroup()
    except AttributeError:
        ag = atoms
        indices = None
    else:
        indices = atoms._getIndices()
        if len(indices) == 1:
            return Atom(ag, indices[0], acsi)

    index = getCentral(atoms._getCoords(), weights)
    if indices is not None:
        index = indices[index]
    return Atom(ag, index, acsi)


def getCentral(coords, weights=None):
    """Returns index of coordinates closest to the center."""

    return ((coords - getCenter(coords, weights))**2).sum(1).argmin()


def pickCentralConf(ens, weights=None):
    """Returns :class:`.Conformation` that is closest to the center of *ens*.
    In addition to :class:`.Ensemble` instances, :class:`.Atomic` instances
    are accepted as *ens* argument. In this case a :class:`.Selection` with
    central coordinate set as active will be returned."""

    try:
        csets = ens._getCoordsets()
    except AttributeError:
        raise TypeError('ens must be an object with multiple '
                        'coordinate sets')
    else:
        if csets is None:
            raise ValueError('coordinates are not set')

    shape = csets.shape
    if shape[0] == 1:
        index = 0
    else:
        csets = csets.reshape((shape[0], shape[1] * shape[2]))
        mean = csets.mean(0)
        index = ((csets - mean)**2).sum(1).argmin()

    try:
        ens.getACSIndex()
    except AttributeError:
        return ens[index]
    else:
        atoms = ens.select('all')
        atoms.setACSIndex(index)
        return atoms


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
            raise ValueError('coords must have shape ([n_csets,]n_atoms,3)')
    if weights is not None:
        weights = weights.flatten()
        weights = weights.reshape(weights.shape[0], 1)
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
            com = (coords * weights).sum(0) / wsum
            d2sum = (((coords - com)**2) * weights).sum()
    else:
        rgyr = []
        for coords in coords:
            if weights is None:
                com = coords.mean(0)
                d2sum = ((coords - com)**2).sum()
                rgyr.append(d2sum)
            else:
                com = (coords * weights).sum(0) / wsum
                d2sum = (((coords - com)**2) * weights).sum()
                rgyr.append(d2sum)
        d2sum = array(rgyr)
    return (d2sum / wsum) ** 0.5

_MSF_DOCSTRING = """  *coordsets* may be an
    instance of :class:`.Ensemble`, :class:`.TrajBase`, or :class:`.Atomic`.
    For trajectory objects, e.g. :class:`.DCDFile`, frames will be considered
    after they are superposed. For other ProDy objects, coordinate sets should
    be aligned prior to MSF calculation.

    Note that using trajectory files that store 32-bit coordinate will result
    in lower precision in calculations.  Over 10,000 frames this may result
    in up to 5% difference from the values calculated using 64-bit arrays.
    To ensure higher-precision calculations for :class:`.DCDFile` instances,
    you may use *astype* argument, i.e. ``astype=float``, to auto recast
    coordinate data to double-precision (64-bit) floating-point format."""


def calcMSF(coordsets):
    """Calculate mean square fluctuation(s) (MSF)."""

    try:
        ncsets = coordsets.numFrames()
    except AttributeError:
        try:
            coordsets = coordsets.getCoordsets()
        except AttributeError:
            pass
        try:
            ndim, shape = coordsets.ndim, coordsets.shape
        except:
            raise TypeError('coordsets must be a Numpy array or a ProDy '
                            'object with `getCoordsets` method')
        if ndim != 3 or shape[0] == 1:
            raise ValueError('coordsets must contain multiple sets')
        msf = var(coordsets, 0).sum(1)
    else:
        nfi = coordsets.nextIndex()
        natoms = coordsets.numSelected()
        total = zeros((natoms, 3))
        sqsum = zeros((natoms, 3))

        LOGGER.progress('Evaluating {0} frames from {1}:'
                        .format(ncsets, str(coordsets)), ncsets,
                        '_prody_calcMSF')
        ncsets = 0
        coordsets.reset()
        for frame in coordsets:
            frame.superpose()
            coords = frame._getCoords()
            total += coords
            sqsum += coords ** 2
            ncsets += 1
            LOGGER.update(ncsets, label='_prody_calcMSF')
        LOGGER.finish()
        msf = (sqsum/ncsets - (total/ncsets)**2).sum(1)
        coordsets.goto(nfi)
    return msf

calcMSF.__doc__ += _MSF_DOCSTRING


def calcRMSF(coordsets):
    """Returns root mean square fluctuation(s) (RMSF)."""

    return calcMSF(coordsets) ** 0.5

calcRMSF.__doc__ += _MSF_DOCSTRING


def calcDeformVector(from_atoms, to_atoms, weights=None):
    """Returns deformation from *from_atoms* to *atoms_to* as a :class:`.Vector`
    instance."""

    from prody.dynamics import Vector

    name = '{0} => {1}'.format(repr(from_atoms), repr(to_atoms))
    if len(name) > 30:
        name = 'Deformation'
    
    from_coords = getCoords(from_atoms)
    to_coords = getCoords(to_atoms)

    arr = to_coords - from_coords
    if weights is not None:
        if weights.ndim > 1:
            weights = weights.flatten()
        if len(weights) != len(arr):
            raise ValueError('weights must have the same length as from_atoms and to_atoms')
        arr = (arr.T * weights).T
    
    arr = arr.flatten()
    
    return Vector(arr, name)


def calcADPAxes(atoms, **kwargs):
    """Returns a 3Nx3 array containing principal axes defining anisotropic
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


    Keyword arguments *fract*, *ratio2*, or *ratio3* can be used to set
    principal axes to 0 for atoms showing relatively lower degree of
    anisotropy.

    3Nx3 axis contains N times 3x3 matrices, one for each given atom. Columns
    of these 3x3 matrices are the principal axes which are weighted by
    square root of their eigenvalues. The first columns correspond to largest
    principal axes.

    The direction of the principal axes for an atom is determined based on the
    correlation of the axes vector with the principal axes vector of the
    previous atom.

    .. ipython:: python

       from prody import *
       protein = parsePDB('1ejg')
       calphas = protein.select('calpha')
       adp_axes = calcADPAxes( calphas )
       adp_axes.shape

    These can be written in NMD format as follows:

    .. ipython:: python

       nma = NMA('ADPs')
       nma.setEigens(adp_axes)
       nma
       writeNMD('adp_axes.nmd', nma, calphas)"""

    linalg = importLA()
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be of type Atomic, not {0}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()

    axes = zeros((n_atoms*3, 3))
    variances = zeros((n_atoms, 3))
    stddevs = zeros((n_atoms, 3))
    anisou = anisous[0]
    element = zeros((3, 3))
    element[0, 0] = anisou[0]
    element[1, 1] = anisou[1]
    element[2, 2] = anisou[2]
    element[0, 1] = element[1, 0] = anisou[3]
    element[0, 2] = element[2, 0] = anisou[4]
    element[1, 2] = element[2, 1] = anisou[5]
    vals, vecs = linalg.eigh(element)
    # If negative eigenvalues are found (when ADP matrix is not positive
    # definite) set them to 0
    vals[vals < 0] = 0
    variances[0] = vals
    vals = vals**0.5
    stddevs[0] = vals
    axes[0:3, :] = vals * vecs

    for i in range(1, n_atoms):
        anisou = anisous[i]
        element[0, 0] = anisou[0]
        element[1, 1] = anisou[1]
        element[2, 2] = anisou[2]
        element[0, 1] = element[1, 0] = anisou[3]
        element[0, 2] = element[2, 0] = anisou[4]
        element[1, 2] = element[2, 1] = anisou[5]
        vals, vecs = linalg.eigh(element)
        # If negative eigenvalues are found (when ADP matrix is not positive
        # definite) set them to 0
        vals[vals < 0] = 0
        variances[i] = vals
        vals = vals**0.5
        stddevs[i] = vals
        # Make sure the direction that correlates with the previous atom
        # is selected
        vals = vals * sign((vecs * axes[(i-1)*3:(i)*3, :]).sum(0))
        axes[i*3:i*3, :] = vals * vecs
    # Resort the columns before returning array
    axes = axes[:, [2, 1, 0]]
    torf = None
    if 'fract' in kwargs:
        fract = float(kwargs['fract'])
        assert 0.33 < fract < 1.0, 'fract must be > 0.33 and < 1.0'
        variances = variances[:, [2, 1, 0]]
        torf = variances[:, 0] / variances.sum(1) > fract
    elif 'ratio' in kwargs or 'ratio3' in kwargs or 'ratio2' in kwargs:
        if 'ratio2' in kwargs:
            ratio = float(kwargs['ratio2'])
            assert 0 < ratio < 1.0, 'ratio2 must be > 0 and < 1.0'
            dim = 1
        else:
            ratio = float(kwargs.get('ratio', kwargs.get('ratio3')))
            assert 0 < ratio < 1.0, 'ratio or ratio3 must be > 0 and < 1.0'
            dim = 2
        variances = variances[:, [2, 1, 0]]
        torf = variances[:, dim] / variances[:, 0] <= ratio
    if torf is not None:
        torf = tile(torf.reshape((n_atoms, 1)), (1, 3)).reshape((n_atoms*3, 1))
        axes = axes * torf
    return axes


def calcADPs(atom):
    """Calculate anisotropic displacement parameters (ADPs) from
    anisotropic temperature factors (ATFs).

    *atom* must have ATF values set for ADP calculation. ADPs are returned
    as a tuple, i.e. (eigenvalues, eigenvectors)."""

    linalg = importLA()
    if not isinstance(atom, Atom):
        raise TypeError('atom must be of type Atom, not {0}'
                        .format(type(atom)))
    anisou = atom.getAnisou()
    if anisou is None:
        raise ValueError('atom does not have anisotropic temperature '
                         'factors')
    element = zeros((3, 3))
    element[0, 0] = anisou[0]
    element[1, 1] = anisou[1]
    element[2, 2] = anisou[2]
    element[0, 1] = element[1, 0] = anisou[3]
    element[0, 2] = element[2, 0] = anisou[4]
    element[1, 2] = element[2, 1] = anisou[5]
    vals, vecs = linalg.eigh(element)
    return vals[[2, 1, 0]], vecs[:, [2, 1, 0]]


def buildADPMatrix(atoms):
    """Returns a 3Nx3N symmetric matrix containing anisotropic displacement
    parameters (ADPs) along the diagonal as 3x3 super elements.

    .. ipython:: python

       from prody import *
       protein = parsePDB('1ejg')
       calphas = protein.select('calpha')
       adp_matrix = buildADPMatrix(calphas)"""

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be of type Atomic, not {0}'
                        .format(type(atoms)))
    anisous = atoms.getAnisous()
    if anisous is None:
        raise ValueError('anisotropic temperature factors are not set')
    n_atoms = atoms.numAtoms()
    n_dof = n_atoms * 3
    adp = zeros((n_dof, n_dof))

    for i in range(n_atoms):
        anisou = anisous[i]
        element = zeros((3, 3))
        element[0, 0] = anisou[0]
        element[1, 1] = anisou[1]
        element[2, 2] = anisou[2]
        element[0, 1] = element[1, 0] = anisou[3]
        element[0, 2] = element[2, 0] = anisou[4]
        element[1, 2] = element[2, 1] = anisou[5]
        adp[i*3:i*3+3, i*3:i*3+3] = element
    return adp


def calcInertiaTensor(coords):
    """"Calculate inertia tensor from coords"""
    coords = getCoords(coords)

    center = calcCenter(coords)
    coords = coords - center
    return dot(coords.transpose(), coords)


def calcPrincAxes(coords, turbo=True):
    """Calculate principal axes from coords"""
    M = calcInertiaTensor(coords)
    _, vectors, _ = solveEig(M, 3, zeros=True, turbo=turbo, reverse=True)
    return vectors.transpose()


def calcDistanceMatrix(coords, cutoff=None):
    """Calculate matrix of distances between coordinates within *cutoff*.
    Other matrix entries are set to maximum of calculated distances.

    :arg coords: a coordinate set or an object with :meth:`getCoords` method.
    :type coords: :class:`~numpy.ndarray`, :class:`.Atomic`

    :arg cutoff: cutoff distance for searching the KDTree.
        Default (**None**) is to use the length of the longest coordinate axis.
    :type cutoff: None, float
    """
    try:
        coords = (coords._getCoords() if hasattr(coords, '_getCoords') else
                coords.getCoords())
    except AttributeError:
        try:
            checkCoords(coords)
        except TypeError:
            raise TypeError('coords must be a Numpy array or an object '
                            'with `getCoords` method')

    n_atoms = coords.shape[0]
    dist_mat = zeros((n_atoms, n_atoms))

    if cutoff is None:
        cutoff = max(coords.max(axis=0) - coords.min(axis=0))

    kdtree = KDTree(coords)
    kdtree.search(cutoff)

    dists = kdtree.getDistances()

    r = 0
    for i, j in kdtree.getIndices():
        dist_mat[i, j] = dist_mat[j, i] = dists[r]
        r += 1

    for i in range(n_atoms):
        for j in range(i, n_atoms):
            if dist_mat[i, j] == 0.:
                dist_mat[i, j] = dist_mat[j, i] = max(dists)

    return dist_mat


def assignBlocks(atoms, res_per_block=None, secstr=False, **kwargs):
    """Assigns blocks to protein from *atoms*
    using a block size of *res_per_block* or 
    secondary structure information if *secstr* is **True**.

    Returns an array of block IDs and an 
    AtomMap corresponding to protein atoms.

    :arg atoms: atoms to be assigned blocks
    :type atoms: :class:`Atomic`

    :arg res_per_block: number of residues per block
        The last block may be smaller or larger than this.
        Default is **None**, allowing *secstr* to be used easily instead.
    :type res_per_block: int

    :arg secstr: use secondary structure information to assign blocks.
        Default is **False**, allowing *res_per_block* to be used easily instead.
        Any set of strings that can be retrieved by :meth:`.getSecstr` is acceptable
        including from PDB header, DSSP or STRIDE.
    :type secstr: bool

    :arg shortest_block: smallest number of residues to be included 
        in a block before merging with the previous block
        Default is **4** as smaller numbers can cause problems for distance matrices.
    :type shortest_block: int

    :arg longest_block: largest number of residues to be included 
        in a block before splitting it in half.
        Default is the length of the protein so it isn't triggered.
    :type longest_block: int

    :arg min_dist_cutoff: minimum distance of a residue from others beyond which 
        it is not included in the same block as them using :meth:`.findSubgroups`.
        Default is 20 A, which was found to work well with *res_per_block*=10.
    :type min_dist_cutoff: Number
    """

    if not isinstance(atoms, Atomic):
        raise TypeError("atoms should be an Atomic object")

    if not atoms.ca:
        raise ValueError("atoms should have Calpha atoms")

    if not isinstance(res_per_block, Integral) and not secstr:
        raise TypeError("res_per_block should be an integer or "
                        "secstr should be set to true")

    if secstr and res_per_block:
        raise ValueError("Either secstr or res_per_block "
                         "should be set, not both")

    if not isinstance(secstr, bool):
        raise TypeError('secstr should be a Boolean')

    shortest_block = kwargs.get('shortest_block', 4)
    shortest_block = kwargs.get('min_size', shortest_block)
    if not isinstance(shortest_block, Integral):
        raise TypeError("shortest_block should be an integer")

    sel_ca = atoms.ca
    n_res = sel_ca.numAtoms()

    longest_block = kwargs.get('max_size', n_res)
    longest_block = kwargs.get('longest_block', longest_block)
    if not isinstance(longest_block, Integral):
        raise TypeError("longest_block should be an integer")

    try:
        min_dist_cutoff = float(kwargs.get('min_dist_cutoff', 20))
    except:
        raise TypeError("min_dist_cutoff should be a number")

    blocks = []

    if res_per_block:
        n_blocks = int(n_res/res_per_block)

        blocks_stack = [[b] * res_per_block for b in range(n_blocks)] + [[n_blocks] * (n_res % res_per_block)]
        for i, block in enumerate(blocks_stack):
            if len(block) < shortest_block:
                # join onto previous block
                block = [blocks_stack[i-1][0] for b in block]

            blocks.extend(block)    
    else:
        secstrs = sel_ca.getSecstrs()
        if secstrs is None:
            raise OSError("Please parse secstr information "
                          "from PDB header or use DSSP or STRIDE "
                          "to use secstr for assigning blocks")
        
        blocks.append(0)
        secstr_prev = secstrs[0]
        i = 0
        for secstr in secstrs[1:-1]:
            if secstr != secstr_prev:
                secstr_prev = secstr
                if len(where(array(blocks) == i)[0]) >= shortest_block:
                    # long enough so next sec str is a new block
                    i += 1
            blocks.append(i)
        
        # include last residue in previous block
        blocks.append(i)

    blocks = array(blocks)

    unique_blocks = unique(blocks)
    for i in unique_blocks:
        block = where(blocks == i)[0]
        dist_mat = calcDistanceMatrix(sel_ca[block])
        dist_tree = calcTree([str(n) for n in block], dist_mat)
        subgroups = [array(list(reversed(n)), dtype=int) for n in
                     reversed(findSubgroups(dist_tree, min_dist_cutoff))]

        if len(subgroups) == 1:
            continue
        elif len(subgroups) == 2:
            # move last subgroup to next block
            blocks[subgroups[-1]] = blocks[subgroups[-1][-1] + 1]
        elif len(subgroups) == 3:
            # move first subgroup to previous block
            blocks[subgroups[0]] = blocks[subgroups[0][0] - 1]
            # move last subgroup to next block
            blocks[subgroups[-1]] = blocks[subgroups[-1][-1] + 1]
        else:
            raise ValueError("Block {0} is getting too separated. "
                             "Please increase min_dist_cutoff".format(i))

    unique_blocks, lengths = unique(blocks, return_counts=True)
    max_length = max(lengths)
    while max_length > longest_block:
        for i in unique_blocks:
            len_block = len(where(blocks == i)[0])

            if len_block > longest_block:
                new_block = max(unique_blocks)+1
                blocks[where(blocks == i)[0][len_block // 2 :]] = new_block

                unique_blocks, lengths = unique(blocks, return_counts=True)
                max_length = max(lengths)
    
    unique_blocks = unique(blocks)
    for i in unique_blocks:
        block = where(blocks == i)[0]
        if len(block) < shortest_block:
            block_im1 = where(blocks == i-1)[0]
            block_ip1 = where(blocks == i+spacing)[0]

            dist_back = calcDistance(atoms[block_im1][-1], atoms[block][0])
            dist_fwd = calcDistance(atoms[block][-1], atoms[block_ip1][0])

            if dist_back < min_dist_cutoff:
                # join onto previous block
                blocks[where(blocks == i)[0]] = i-1
            elif dist_fwd < min_dist_cutoff:
                # join onto next block
                blocks[where(blocks == i)[0]] = i

    blocks, amap = extendAtomicData(blocks, sel_ca, atoms)

    if amap.getHierView().numResidues() < atoms.getHierView().numResidues():
        amap.setTitle("protein from " + amap.getTitle())

    return blocks, amap
