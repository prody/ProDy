# -*- coding: utf-8 -*-
""" This module defines a class and function for identifying contacts."""

from numpy import array, ndarray

from prody.atomic import Atomic, Atom, AtomGroup, AtomSubset, Selection
from prody.kdtree import KDTree
from prody.utilities import rangeString

__all__ = ['Contacts', 'iterNeighbors', 'findNeighbors']

class Contacts(object):

    """A class for contact identification.
    Contacts are identified using the coordinates of atoms at the time
    of instantiation."""

    def __init__(self, atoms, unitcell=None):
        """*atoms* must be an :class:`.Atomic` instance.  When an orthorhombic
        *unitcell* array is given"""

        try:
            self._acsi = atoms.getACSIndex()
        except AttributeError:
            try:
                self._ag = atoms.getAtoms()
                unitcell = unitcell or atoms.getUnitcell()[:3]
                self._indices = atoms.getSelection()
            except AttributeError:
                try:
                    ndim, shape = atoms.ndim, atoms.shape
                except AttributeError:
                    raise TypeError('atoms must be an Atomic or Frame instance'
                                    ', not a {0}'.format(type(atoms)))
                else:
                    if not (ndim == 2 and shape[1] == 3):
                        raise ValueError('atoms.shape must be (n_atoms, 3) or '
                                         '(3,).')
                    self._ag = None
                    self._indices = None
                    self._kdtree = KDTree(atoms, unitcell=unitcell)
            else:
                if self._ag is not None:
                    self._acsi = self._ag.getACSIndex()
                    if self._indices is not None:
                        self._indices = self._indices.getIndices()
                else:
                    self._acsi = None
                self._kdtree = KDTree(self._atoms._getCoords(),
                                      unitcell=unitcell)
        else:
            try:
                self._ag = atoms.getAtomGroup()
            except AttributeError:
                self._ag = atoms
                self._indices = None
                self._kdtree = KDTree(atoms._getCoords(), unitcell=unitcell)
            else:
                self._indices = atoms._getIndices()
                self._kdtree = KDTree(atoms._getCoords(), unitcell=unitcell)
        self._unitcell = unitcell
        self._atoms = atoms

    def __repr__(self):

        return '<Contacts: {0} (active coordset index: {1})>'.format(
                                                str(self._atoms), self._acsi)

    def __str__(self):

        return 'Contacts ' + str(self._atoms)

    def __call__(self, radius, center):
        """Select atoms radius *radius* (Å) of *center*, which can be point(s)
        in 3-d space (:class:`numpy.ndarray` with shape ``(n_atoms, 3)``) or a
        set of atoms, e.g. :class:`.Selection`."""

        try:
            center = center._getCoords()
        except AttributeError:
            try:
                ndim, shape = center.ndim, center.shape
            except AttributeError:
                raise TypeError('center must be an Atomic instance or a'
                                'coordinate array')
            else:
                if shape == (3,):
                    center = [center]
                elif not ndim == 2 and shape[1] == 3:
                    raise ValueError('center.shape must be (n_atoms, 3) or'
                                     '(3,)')
        else:
            if center is None:
                raise ValueError('center does not have coordinate data')

        search = self._kdtree.search
        get_indices = self._kdtree.getIndices
        get_count = self._kdtree.getCount
        indices = set()
        update = indices.update
        radius = float(radius)
        for xyz in center:
            search(radius, xyz)
            if get_count():
                update(get_indices())
        indices = list(indices)
        if indices:
            indices.sort()
            if self._ag is None:
                return array(indices)
            else:
                if self._indices is not None:
                    indices = self._indices[indices]
                return Selection(self._ag, array(indices), 'index ' +
                                 rangeString(indices), acsi=self._acsi,
                                 unique=True)

    select = __call__

    def getAtoms(self):
        """Returns atoms, or coordinate array, provided at instantiation.."""

        return self._atoms

    def getUnitcell(self):
        """Returns unitcell array, or **None** if one was not provided."""

        return self._unitcell.copy()


def iterNeighbors(atoms, radius, atoms2=None, unitcell=None):
    """Yield pairs of *atoms* that are within *radius* of each other and the
    distance between them.  If *atoms2* is also provided, one atom from *atoms*
    and another from *atoms2* will be yielded.  If one of *atoms* or *atoms2*
    is a coordinate array, pairs of indices and distances will be yielded.
    When orthorhombic *unitcell* dimensions are provided, periodic boundary
    conditions will be taken into account (see :class:`.KDTree` and also
    :func:`wrapAtoms` for details).  If *atoms* is a :class:`.Frame` instance
    and *unitcell* is not provided, unitcell information from frame will be
    if available."""

    radius = float(radius)
    if radius <= 0:
        raise ValueError('radius must be a positive number')

    try:
        acsi = atoms.getACSIndex()
    except AttributeError:
        try:
            ndim, shape = atoms.ndim, atoms.shape
        except AttributeError:
            try:
                uc = atoms.getUnitcell()[:3]
            except AttributeError:
                raise TypeError('atoms must be an Atomic or Frame instance or '
                                'a coordinate array')
            else:
                coords = atoms._getCoords()
                ag = atoms.getAtoms()
                if unitcell is None:
                    unitcell = uc
                if ag is not None:
                    acsi = ag.getACSIndex()
                    _ = atoms.getSelection()
                    if _:
                        indices = _._getIndices()
                        index = lambda i: indices[i]
                    else:
                        index = lambda i: i
        else:
            if ndim > 2:
                raise ValueError('number of dimensions of coordinate array '
                                 'must be 1 or 2')
            coords = atoms
            ag = None
            acsi = None
    else:
        coords = atoms._getCoords()
        try:
            ag = atoms.getAtomGroup()
            indices = atoms.getIndices()
            index = lambda i: indices[i]
        except AttributeError:
            ag = atoms
            index = lambda i: i

    if coords.ndim == 1:
        coords = array([coords])

    if atoms2 is None:
        if len(coords) <= 1:
            raise ValueError('atoms must be more than 1')

        kdtree = KDTree(coords, unitcell=unitcell, none=list)

        _dict = {}
        if ag is None:
            for (i, j), r in zip(*kdtree(radius)):
                yield (i, j, r)
        else:
            for (i, j), r in zip(*kdtree(radius)):
                a1 = _dict.get(i)
                if a1 is None:
                    a1 = Atom(ag, index(i), acsi)
                    _dict[i] = a1
                a2 = _dict.get(j)
                if a2 is None:
                    a2 = Atom(ag, index(j), acsi)
                    _dict[j] = a2
                yield (a1, a2, r)
    else:
        try:
            coords2 = atoms2._getCoords()
        except AttributeError:
            try:
                ndim, shape = atoms2.ndim, atoms2.shape
            except AttributeError:
                raise TypeError('atoms2 must be an Atomic or Frame instance or '
                                'a coordinate array')
            else:
                if ndim > 2:
                    raise ValueError('number of dimensions of second '
                                     'coordinate array must be 1 or 2')
                coords2 = atoms2
                ag2 = None
                acsi2 = None
        else:
            try:
                acsi2 = atoms2.getACSIndex()
            except AttributeError:
                acsi2 = None
                ag2 = None
                index2 = None
            else:
                try:
                    ag2 = atoms2.getAtomGroup()
                    indices2 = atoms2.getIndices()
                    index2 = lambda i: indices2[i]
                except AttributeError:
                    ag2 = atoms2
                    index2 = lambda i: i

        if coords2.ndim == 1:
            coords2 = array([coords2])
        if len(coords) >= len(coords2):
            kdtree = KDTree(coords, unitcell=unitcell, none=list)
            _dict = {}
            if ag is None or ag2 is None:
                for j, xyz in enumerate(coords2):
                    for i, r in zip(*kdtree(radius, xyz)):
                        yield (i, j, r)
            else:
                for a2 in atoms2.iterAtoms():
                    for i, r in zip(*kdtree(radius, a2._getCoords())):
                        a1 = _dict.get(i)
                        if a1 is None:
                            a1 = Atom(ag, index(i), acsi)
                            _dict[i] = a1
                        yield (a1, a2, r)
        else:
            kdtree = KDTree(coords2, unitcell=unitcell, none=list)
            _dict = {}
            if ag is None or ag2 is None:
                for i, xyz in enumerate(coords):
                    for j, r in zip(*kdtree(radius, xyz)):
                        yield (i, j, r)
            else:
                for a1 in atoms.iterAtoms():
                    for i, r in zip(*kdtree(radius, a1._getCoords())):
                        a2 = _dict.get(i)
                        if a2 is None:
                            a2 = Atom(ag2, index2(i), acsi2)
                            _dict[i] = a2
                        yield (a1, a2, r)


def findNeighbors(atoms, radius, atoms2=None, unitcell=None):
    """Returns list of neighbors that are within *radius* of each other and the
    distance between them.  See :func:`iterNeighbors` for more details."""

    return list(iterNeighbors(atoms, radius, atoms2, unitcell))
