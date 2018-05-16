# -*- coding: utf-8 -*-
"""This module defines functions for editing normal mode data."""

import numpy as np

from prody.atomic import Atomic, AtomGroup, AtomMap, AtomSubset
from prody.atomic import Selection, SELECT, sliceAtoms, extendAtoms
from prody.utilities import importLA
from prody import _PY3K

from .nma import NMA
from .mode import VectorBase, Mode, Vector
from .gnm import GNM
from .anm import ANM
from .pca import PCA

if not _PY3K:
    range = xrange

__all__ = ['extendModel', 'extendMode', 'extendVector',
           'sliceMode', 'sliceModel', 'sliceVector',
           'reduceModel']


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

    vec = Vector(vector.getArrayNx3()[
                 which, :].flatten(),
                 '{0} slice {1}'.format(str(vector), select),
                 vector.is3d())
    return (vec, sel)


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

    vec = Vector(mode.getArrayNx3()[which, :].flatten() *
                 mode.getVariance()**0.5,
                 '{0} slice {1}'.format(str(mode), select), mode.is3d())
    return (vec, sel)


def sliceModel(model, atoms, select):
    """Returns a part of the *model* for *atoms* matching *select*.  Note that
    normal modes (eigenvectors) are not normalized.

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

    array = model._getArray()

    which, sel = sliceAtoms(atoms, select)

    nma = type(model)('{0} slice {1}'
                      .format(model.getTitle(), select))
    if model.is3d():
        which = [which.reshape((len(which), 1))*3]
        which.append(which[0]+1)
        which.append(which[0]+2)
        which = np.concatenate(which, 1).flatten()
    nma.setEigens(array[which, :], model.getEigvals())
    return (nma, sel)


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

    linalg = importLA()

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
    system = np.zeros(model.numAtoms(), dtype=bool)
    system[which] = True

    other = np.invert(system)

    if model.is3d():
        system = np.repeat(system, 3)
        other = np.repeat(other, 3)
    ss = matrix[system, :][:, system]
    if isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(ss)
        return eda, system
    so = matrix[system, :][:, other]
    os = matrix[other, :][:, system]
    oo = matrix[other, :][:, other]
    try:
        invoo = linalg.inv(oo)
    except:
        invoo = linalg.pinv(oo)
    matrix = ss - np.dot(so, np.dot(invoo, os))

    if isinstance(model, GNM):
        gnm = GNM(model.getTitle() + ' reduced')
        gnm.setKirchhoff(matrix)
        return gnm, select
    elif isinstance(model, ANM):
        anm = ANM(model.getTitle() + ' reduced')
        anm.setHessian(matrix)
        return anm, select
    elif isinstance(model, PCA):
        eda = PCA(model.getTitle() + ' reduced')
        eda.setCovariance(matrix)
        return eda, select
