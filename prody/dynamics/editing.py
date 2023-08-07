# -*- coding: utf-8 -*-
"""This module defines functions for editing normal mode data."""

import numpy as np

from prody.atomic import Atomic, AtomGroup, AtomMap, AtomSubset
from prody.atomic import Selection, SELECT, sliceAtoms, extendAtoms
from prody.measure import calcDistance
from prody.utilities import importLA, isListLike, getCoords
from prody import _PY3K, LOGGER

from .nma import NMA
from .mode import VectorBase, Mode, Vector
from .gnm import GNM
from .anm import ANM
from .pca import PCA

if not _PY3K:
    range = xrange

__all__ = ['extendModel', 'extendMode', 'extendVector',
           'sliceMode', 'sliceModel', 'sliceModelByMask', 'sliceVector',
           'reduceModel', 'reduceModelByMask', 'trimModel', 'trimModelByMask',
           'interpolateModel']


def extendModel(model, nodes, atoms, norm=False):
    """Extend a coarse grained *model* built for *nodes* to *atoms*.  *model*
    may be :class:`.ANM`, :class:`.GNM`, :class:`.PCA`, or :class:`.NMA`
    instance.  This function will take part of the normal modes for each node
    (i.e. Cα atoms) and extend it to all other atoms in the same residue.  For
    each atom in *nodes* argument *atoms* argument must contain a corresponding
    residue.  If *norm* is **True**, extended modes are normalized."""

    try:
        evecs = model._getArray()
        evals = model.getEigvals()
    except AttributeError:
        raise ValueError('model must be an NMA instance')

    if model.numAtoms() != nodes.numAtoms():
        raise ValueError('atom numbers must be the same')

    indices, atommap = extendAtoms(nodes, atoms, model.is3d())

    evecs = evecs[indices, :]
    if norm:
        evecs /= np.array([((evecs[:, i]) ** 2).sum() ** 0.5
                           for i in range(evecs.shape[1])])

    if model.is3d():
        extended = NMA('Extended ' + str(model))
    else:
        extended = GNM('Extended ' + str(model))
    extended.setEigens(evecs, evals)
    return extended, atommap


def extendMode(mode, nodes, atoms, norm=False):
    """Extend a coarse grained normal *mode* built for *nodes* to *atoms*.
    This function will take part of the normal modes for each node (i.e. Cα
    atoms) and extend it to all other atoms in the same residue.  For each atom
    in *nodes* argument *atoms* argument must contain a corresponding residue.
    Extended mode is multiplied by the square root of variance of the mode.
    If *norm* is **True**, extended mode is normalized."""

    try:
        vec = mode._getArray()
        std = mode.getVariance() ** 0.5
    except AttributeError:
        raise ValueError('mode must be a normal Mode instance')

    indices, atommap = extendAtoms(nodes, atoms, mode.is3d())
    vec = vec[indices]
    if norm:
        vec /= ((vec) ** 2).sum() ** 0.5
    else:
        vec *= std
    extended = Vector(vec, 'Extended ' + str(mode), mode.is3d())
    return extended, atommap


def extendVector(vector, nodes, atoms):
    """Extend a coarse grained *vector* for *nodes* to *atoms*.  This function
    will take part of the normal modes for each node (i.e. Cα atoms) and extend
    it to all other atoms in the same residue.  For each atom in *nodes*,
    *atoms* argument must contain a corresponding residue."""

    try:
        vec = vector._getArray()
    except AttributeError:
        raise ValueError('vector must be a Vector instance')

    if vector.numAtoms() != nodes.numAtoms():
        raise ValueError('atom numbers must be the same')

    indices, atommap = extendAtoms(nodes, atoms, vector.is3d())
    extended = Vector(vec[indices], 'Extended ' + str(vector), vector.is3d())
    return extended, atommap

def sliceVector(vector, atoms, select):
    """Returns part of the *vector* for *atoms* matching *select*.  Note that
    returned :class:`.Vector` instance is not normalized.

    :arg vector: vector instance to be sliced
    :type vector: :class:`.VectorBase`

    :arg atoms: atoms for which *vector* describes a deformation, motion, etc.
    :type atoms: :class:`.Atomic`

    :arg select: an atom selection or a selection string
    :type select: :class:`.Selection`, str

    :returns: (:class:`.Vector`, :class:`.Selection`)"""

    if not isinstance(vector, VectorBase):
        raise TypeError('vector must be a VectorBase instance, not {0}'
                        .format(type(vector)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0}'
                        .format(type(atoms)))
    if atoms.numAtoms() != vector.numAtoms():
        raise ValueError('number of atoms in model and atoms must be equal')

    which, sel = sliceAtoms(atoms, select)

    if vector.is3d():
        vec = Vector(vector.getArrayNx3()[which, :].flatten(),
                    '{0} slice {1}'.format(str(vector), select),
                    vector.is3d())
    else:
        vec = Vector(vector.getArray()[which].flatten(),
                    '{0} slice {1}'.format(str(vector), select),
                    vector.is3d())
    return (vec, sel)

def trimModel(model, atoms, select):
    """Returns a part of the *model* for *atoms* matching *select*. This method removes 
    columns and rows in the connectivity matrix and fix the diagonal sums. Normal modes 
    need to be calculated again after the trim.

    :arg mode: NMA model instance to be sliced
    :type mode: :class:`.NMA`

    :arg atoms: atoms for which the *model* was built
    :type atoms: :class:`.Atomic`

    :arg select: an atom selection or a selection string
    :type select: :class:`.Selection`, str

    :returns: (:class:`.NMA`, :class:`.Selection`)"""

    if not isinstance(model, NMA):
        raise TypeError('mode must be a NMA instance, not {0}'
                        .format(type(model)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0}'
                        .format(type(atoms)))
    if atoms.numAtoms() != model.numAtoms():
        raise ValueError('number of atoms in model and atoms must be equal')

    which, sel = sliceAtoms(atoms, select)
    nma = trimModelByMask(model, which)

    return (nma, sel)

def trimModelByMask(model, mask):
    """Returns a part of the *model* indicated by *mask*. This method removes 
    columns and rows in the connectivity matrix indicated by *mask* and fix the diagonal sums.
    Normal modes need to be calculated again after the trim.

    :arg mode: NMA model instance to be sliced
    :type mode: :class:`.NMA`

    :arg mask: an Integer array or a Boolean array where ``"True"`` indicates 
        the parts being selected 
    :type mask: list, :class:`~numpy.ndarray`

    :returns: :class:`.NMA`"""

    if not isListLike(mask):
        raise TypeError('mask must be either a list or a numpy.ndarray, not {0}'
                        .format(type(model)))
    
    is_bool = mask.dtype is np.dtype('bool')

    if is_bool:
        if len(mask) != model.numAtoms():
            raise ValueError('number of atoms in model and mask must be equal')
        which = mask
    else:
        if mask.min() < 0 or mask.max() >= model.numAtoms():
            raise ValueError('index in mask exceeds range')
        which = np.zeros(model.numAtoms(), dtype=bool)
        which[mask] = True

    if model.is3d():
        which = np.repeat(which, 3)

    if isinstance(model, GNM):
        matrix = model._kirchhoff
    elif isinstance(model, ANM):
        matrix = model._hessian
    elif isinstance(model, PCA):
        matrix = model._cov
    
    if isinstance(model, PCA):
        ss = matrix[which, :][:, which]
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(ss)
        return eda
    else:
        matrix = matrix[which, :][:, which]

        if isinstance(model, GNM):
            gnm = GNM(model.getTitle() + ' reduced')
            I = np.eye(len(matrix), dtype=bool)
            matrix[I] = - (matrix.sum(axis=0) - np.diag(matrix))
            gnm.setKirchhoff(matrix)
            return gnm
        elif isinstance(model, ANM):
            model_type = type(model)
            anm = model_type(model.getTitle() + ' reduced')
            
            n = len(matrix) // 3
            for i in range(n):
                S = np.zeros((3, 3))
                for j in range(n):
                    if i == j:
                        continue
                    S -= matrix[i*3:i*3+3, j*3:j*3+3]
                matrix[i*3:i*3+3, i*3:i*3+3] = S
            anm.setHessian(matrix)
            if hasattr(anm, 'getMembrane'):
                anm._membrane = model.getMembrane()
                anm._combined = model.getCombined()
            return anm

def sliceMode(mode, atoms, select):
    """Returns part of the *mode* for *atoms* matching *select*.  This works
    slightly different from :func:`.sliceVector`. Mode array (eigenvector) is
    multiplied by square-root of the variance along the mode.  If mode is from
    an elastic network model, variance is defined as the inverse of the
    eigenvalue.  Note that returned :class:`.Vector` instance is not
    normalized.

    :arg mode: mode instance to be sliced
    :type mode: :class:`.Mode`

    :arg atoms: atoms for which *mode* describes a deformation, motion, etc.
    :type atoms: :class:`.Atomic`

    :arg select: an atom selection or a selection string
    :type select: :class:`.Selection`, str

    :returns: (:class:`.Vector`, :class:`.Selection`)"""

    if not isinstance(mode, Mode):
        raise TypeError('mode must be a Mode instance, not {0}'
                        .format(type(mode)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0}'
                        .format(type(atoms)))
    if atoms.numAtoms() != mode.numAtoms():
        raise ValueError('number of atoms in model and atoms must be equal')

    which, sel = sliceAtoms(atoms, select)

    if mode.is3d():
        vec = Vector(mode.getArrayNx3()[which, :].flatten() *
                     mode.getVariance()**0.5,
                    '{0} slice {1}'.format(str(mode), select), mode.is3d())
    else:
        vec = Vector(mode.getArray()[which].flatten() *
                     mode.getVariance()**0.5,
                    '{0} slice {1}'.format(str(mode), select), mode.is3d())
    return (vec, sel)


def sliceModel(model, atoms, select, norm=False):
    """Returns a part of the *model* (modes calculated) for *atoms* matching *select*. 
    Note that normal modes are sliced instead the connectivity matrix. Sliced normal 
    modes (eigenvectors) are not normalized unless *norm* is **True**.

    :arg mode: NMA model instance to be sliced
    :type mode: :class:`.NMA`

    :arg atoms: atoms for which the *model* was built
    :type atoms: :class:`.Atomic`

    :arg select: an atom selection or a selection string
    :type select: :class:`.Selection`, str

    :arg norm: whether to normalize eigenvectors, default **False**
    :type norm: bool

    :returns: (:class:`.NMA`, :class:`.Selection`)"""

    if not isinstance(model, NMA):
        raise TypeError('mode must be a NMA instance, not {0}'
                        .format(type(model)))
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance, not {0}'
                        .format(type(atoms)))
    if atoms.numAtoms() != model.numAtoms():
        raise ValueError('number of atoms in model and atoms must be equal')

    which, sel = sliceAtoms(atoms, select)
    nma = sliceModelByMask(model, which, norm=norm)

    return (nma, sel)

def sliceModelByMask(model, mask, norm=False):
    """Returns a part of the *model* indicated by *mask*.  Note that
    normal modes (eigenvectors) are not normalized unless *norm* is **True**.

    :arg mode: NMA model instance to be sliced
    :type mode: :class:`.NMA`

    :arg mask: an Integer array or a Boolean array where ``"True"`` indicates 
        the parts being selected 
    :type mask: list, :class:`~numpy.ndarray`

    :arg norm: whether to normalize eigenvectors, default **False**
    :type norm: bool

    :returns: :class:`.NMA`"""

    if not isListLike(mask):
        raise TypeError('mask must be either a list or a numpy.ndarray, not {0}'
                        .format(type(model)))
    
    is_bool = mask.dtype is np.dtype('bool')

    if is_bool:
        if len(mask) != model.numAtoms():
            raise ValueError('number of atoms in model and mask must be equal')
        which = mask
    else:
        if mask.min() < 0 or mask.max() >= model.numAtoms():
            raise ValueError('index in mask exceeds range')
        which = np.zeros(model.numAtoms(), dtype=bool)
        which[mask] = True

    array = model._getArray()
    
    nma = type(model)('{0} sliced'
                    .format(model.getTitle()))
    if model.is3d():
        which = np.repeat(which, 3)

    evecs = array[which, :]
    if norm:
        evecs /= np.array([((evecs[:, i]) ** 2).sum() ** 0.5
                           for i in range(evecs.shape[1])])

    nma.setEigens(evecs, model.getEigvals())
    return nma

def reduceModel(model, atoms, select):
    """Returns reduced NMA model.  Reduces a :class:`.NMA` model to a subset of
    *atoms* matching *select*.  This function behaves differently depending on
    the type of the *model* argument.  For :class:`.ANM` and :class:`.GNM` or
    other :class:`.NMA` models, force constant matrix for system of interest
    (specified by the *select*) is derived from the force constant matrix for
    the *model* by assuming that for any given displacement of the system of
    interest, other atoms move along in such a way as to minimize the potential
    energy.  This is based on the formulation in [KH00]_.  For :class:`.PCA`
    models, this function simply takes the sub-covariance matrix for selection.

    .. [KH00] Hinsen K, Petrescu A-J, Dellerue S, Bellissent-Funel M-C, Kneller GR.
       Harmonicity in slow protein dynamics. *Chem Phys* **2000** 261:25-37.

    :arg model: dynamics model
    :type model: :class:`.ANM`, :class:`.GNM`, or :class:`.PCA`

    :arg atoms: atoms that were used to build the model
    :type atoms: :class:`.Atomic`

    :arg select: an atom selection or a selection string
    :type select: :class:`.Selection`, str

    :returns: (:class:`.NMA`, :class:`.Selection`)"""

    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance, not {0}'
                        .format(type(model)))
    if not isinstance(atoms, (AtomGroup, AtomSubset, AtomMap)):
        raise TypeError('atoms type is not valid')
    if len(atoms) <= 1:
        raise TypeError('atoms must contain more than 1 atoms')

    if isinstance(model, GNM):
        matrix = model._kirchhoff
    elif isinstance(model, ANM):
        matrix = model._hessian
    elif isinstance(model, PCA):
        matrix = model._cov
    else:
        raise TypeError('model does not have a valid type derived from NMA')
    if matrix is None:
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not '
                         'built')

    which, select = sliceAtoms(atoms, select)
    nma = reduceModelByMask(model, which)

    return nma, select

def reduceModelByMask(model, mask):
    """Returns NMA model reduced based on *mask*. 

    :arg model: dynamics model
    :type model: :class:`.ANM`, :class:`.GNM`, or :class:`.PCA`

    :arg mask: an Integer array or a Boolean array where ``"True"`` indicates 
        the parts being selected 
    :type mask: list, :class:`~numpy.ndarray`

    :returns: :class:`.NMA`"""

    if not isinstance(model, NMA):
        raise TypeError('model must be an NMA instance, not {0}'
                        .format(type(model)))
    
    if not isListLike(mask):
        raise TypeError('mask must be either a list or a numpy.ndarray, not {0}'
                        .format(type(model)))
    
    is_bool = mask.dtype is np.dtype('bool')

    if is_bool:
        if len(mask) != model.numAtoms():
            raise ValueError('number of atoms in model and mask must be equal')
        system = mask
    else:
        if mask.min() < 0 or mask.max() >= model.numAtoms():
            raise ValueError('index in mask exceeds range')
        system = np.zeros(model.numAtoms(), dtype=bool)
        system[mask] = True

    if isinstance(model, GNM):
        matrix = model._kirchhoff
    elif isinstance(model, ANM):
        matrix = model._hessian
    elif isinstance(model, PCA):
        matrix = model._cov
    else:
        raise TypeError('model does not have a valid type derived from NMA')
    if matrix is None:
        raise ValueError('model matrix (Hessian/Kirchhoff/Covariance) is not '
                         'built')

    if model.is3d():
        system = np.repeat(system, 3)

    if isinstance(model, PCA):
        ss = matrix[system, :][:, system]
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(ss)
        return eda
    else:
        matrix = _reduceModel(matrix, system)

        if isinstance(model, GNM):
            gnm = GNM(model.getTitle() + ' reduced')
            gnm.setKirchhoff(matrix)
            return gnm
        elif isinstance(model, ANM):
            anm = ANM(model.getTitle() + ' reduced')
            anm.setHessian(matrix)
            return anm

def _reduceModel(matrix, system):
    """This is the underlying function that reduces models, which shall 
    remain private. *system* is a Boolean array where **True** indicates 
    system nodes."""
    
    linalg = importLA()

    other = np.invert(system)

    ss = matrix[system, :][:, system]
    so = matrix[system, :][:, other]
    os = matrix[other, :][:, system]
    oo = matrix[other, :][:, other]

    if other.any():
        try:
            invoo = linalg.inv(oo)
        except:
            invoo = linalg.pinv(oo)
        
        matrix = ss - np.dot(so, np.dot(invoo, os))
    else:
        matrix = ss

    return matrix


def interpolateModel(model, nodes, coords, norm=False, **kwargs):
    """Interpolate a coarse grained *model* built for *nodes* to *coords*.  
    
    *model* may be :class:`.ANM`, :class:`.PCA`, or :class:`.NMA`
    instance

    :arg nodes: the coordinate set or object with :meth:`getCoords` method
        that corresponds to the model
    :type nodes: :class:`.Atomic`, :class:`~numpy.ndarray`

    :arg coords: a coordinate set or an object with :meth:`getCoords` method
        onto which the model should be interpolated
    :type coords: :class:`.Atomic`, :class:`~numpy.ndarray`
    
    This function will take the part of the normal modes for each node
    (i.e. Cα atoms) and extend it to nearby atoms. 
    
    If *norm* is **True**, extended modes are normalized.
    
    Adapted from ModeHunter as described in [JS09]_.

    .. [JS09] Stember JN, Wriggers W. Bend-twist-stretch model for coarse 
    elastic network simulation of biomolecular motion. *J Chem Phys* **2009** 131:074112.

    # Legal notice:
    # 
    # This software is copyrighted, (c) 2009-10, by Joseph N. Stember and Willy Wriggers
    # under the following terms:
    #
    # The authors hereby grant permission to use, copy, modify, and re-distribute this
    # software and its documentation for any purpose, provided that existing copyright
    # notices are retained in all copies and that this notice is included verbatim in
    # any distributions. No written agreement, license, or royalty fee is required for
    # any of the authorized uses.
    #
    # IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY FOR DIRECT,
    # INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE
    # OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY DERIVATIVES THEREOF, EVEN IF THE
    # AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    #
    # THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
    # BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
    # PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE IS PROVIDED ON AN "AS
    # IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
    # MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
    #
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """

    quiet = kwargs.pop('quiet', False)

    try:
        eigvecs = model._getArray()
        eigvals = model.getEigvals()
    except AttributeError:
        raise TypeError('model must be an NMA instance')

    if model.numAtoms() != nodes.numAtoms():
        raise ValueError('atom numbers must be the same')

    if not model.is3d():
        raise ValueError('model must be 3D')

    cg_coords = nodes.getCoords()
    len_cg = nodes.numAtoms()

    fg_coords = getCoords(coords)
    len_fg = fg_coords.shape[0]

    n_modes = model.numModes()

    eigvecs_fg = np.zeros((3*len_fg, n_modes))

    if not quiet:
        LOGGER.progress('Interpolating {0} coarse modes to fine modes...'.format(n_modes),
                        n_modes, '_prody_iterp')

    for k in np.arange(n_modes): # loop through cg eigenvecs
        coarse_grained_nms = eigvecs[:, k]
        coarse_grained_nms.shape = (-1,3)

        lmat = np.zeros((len_cg+4,len_cg+4))
        for i in np.arange(len_cg):
            for j in np.arange(len_cg):
                lmat[i,j] = np.linalg.norm(cg_coords[i]-cg_coords[j])

        for i in np.arange(len_cg):
            lmat[i,len_cg] = 1.0
            lmat[i,len_cg+1] = cg_coords[i,0]
            lmat[i,len_cg+2] = cg_coords[i,1]
            lmat[i,len_cg+3] = cg_coords[i,2]
            lmat[len_cg,i] = 1.0
            lmat[len_cg+1,i] = cg_coords[i,0]
            lmat[len_cg+2,i] = cg_coords[i,1]
            lmat[len_cg+3,i] = cg_coords[i,2]

        lmat_inv = np.linalg.inv(lmat)

        vxprime = coarse_grained_nms[:,0]
        for i in np.arange(4):
            vxprime = np.append(vxprime,0.0)

        vyprime = coarse_grained_nms[:,1]
        for i in np.arange(4):
            vyprime = np.append(vyprime,0.0)
            
        vzprime = coarse_grained_nms[:,2]
        for i in np.arange(4):
            vzprime = np.append(vzprime,0.0)

        vx = np.dot(lmat_inv,vxprime)
        vy = np.dot(lmat_inv,vyprime)
        vz = np.dot(lmat_inv,vzprime)

        nmx=np.zeros(len_fg)
        nmy=np.zeros(len_fg)
        nmz=np.zeros(len_fg)

        for j in np.arange(len_fg):
            nmx[j] += vx[len_cg] + fg_coords[j,0]*vx[len_cg+1] + fg_coords[j,1]*vx[len_cg+2] + fg_coords[j,2]*vx[len_cg+3]  
            nmy[j] += vy[len_cg] + fg_coords[j,0]*vy[len_cg+1] + fg_coords[j,1]*vy[len_cg+2] + fg_coords[j,2]*vy[len_cg+3]  
            nmz[j] += vz[len_cg] + fg_coords[j,0]*vz[len_cg+1] + fg_coords[j,1]*vz[len_cg+2] + fg_coords[j,2]*vz[len_cg+3]  

            dist_j = calcDistance(cg_coords, fg_coords[j])
            for i in np.arange(len_cg):
                nmx[j] += vx[i]*dist_j[i]
                nmy[j] += vy[i]*dist_j[i]
                nmz[j] += vz[i]*dist_j[i]

        eigvecs_fg[::3, k] = nmx
        eigvecs_fg[1::3, k] = nmy
        eigvecs_fg[2::3, k] = nmz  

        if not quiet:
            LOGGER.update(k+1, label='_prody_iterp')

    if not quiet:
        LOGGER.finish()

    if norm:
        for k in np.arange(np.size(cg_coords)):
            eigvecs_fg[k] = eigvecs_fg[k]/np.linalg.norm(eigvecs_fg[k])  
            
    interpolated = NMA('Interpolated ' + str(model))
    interpolated.setEigens(eigvecs_fg, eigvals)
    return interpolated, coords
