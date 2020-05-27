# -*- coding: utf-8 -*-
"""This module defines :class:`Dihedral` for dealing with dihedral information provided
by using :meth:`.AtomGroup.setDihedrals` method."""

from numbers import Integral
import numpy as np

RAD2DEG = 180 / np.pi

__all__ = ['Dihedral']


class Dihedral(object):

    """A pointer class for dihedrald atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns dihedral length, i.e. :meth:`getSize`
    * :func:`iter` yields :class:`~.Atom` instances"""

    __slots__ = ['_ag', '_acsi', '_indices']

    def __init__(self, ag, indices, acsi=None):

        self._ag = ag
        self._indices = np.array(indices)
        if acsi is None:
            self._acsi = ag.getACSIndex()
        else:
            self._acsi = acsi

    def __repr__(self):

        one, two, three, four = self._indices
        names = self._ag._getNames()
        return '<Dihedral: {0}({1})--{2}({3})--{4}({5}--{6}({7})) from {8}>'.format(
            names[one], one, names[two], two,
            names[three], three, names[four], four,
            str(self._ag))

    def __str__(self):

        one, two, three, four = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})--{4}({5})--{6}({7})'.format(
            names[one], one, names[two], two,
            names[three], three, names[four], four)

    def __eq__(self, other):

        return (isinstance(other, Dihedral) and other.getAtomGroup() is self._ag
                and (np.all(other.getIndices() == self._indices) or
                     np.all(other.getIndices() == list(reversed(self._indices)))))

    def __ne__(self, other):

        return not self.__eq__(other)

    def __size__(self):

        return self.getSize()

    def __iter__(self):

        for index in self._indices:
            yield self._ag[index]

    def getAtomGroup(self):
        """Returns atom group."""

        return self._ag

    def getAtoms(self):
        """Returns dihedrald atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]], self._ag[self._indices[2]])

    def getIndices(self):
        """Returns indices of dihedrald atoms."""

        return self._indices.copy()

    def getSize(self, radian=False):
        """Returns dihedral size."""

        a1, a2, a3 = self.getVectors()

        v1 = np.cross(a1, a2)
        v1 = v1 / (v1 * v1).sum(-1)**0.5
        v2 = np.cross(a2, a3)
        v2 = v2 / (v2 * v2).sum(-1)**0.5
        porm = np.sign((v1 * a3).sum(-1))
        rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
        if not porm == 0:
            rad = rad * porm
        if radian:
            return rad
        else:
            return rad * RAD2DEG

    def getVectors(self):
        """Returns bond vectors that originate from the central atom."""

        one, two, three, four = self._indices
        acsi = self.getACSIndex()
        vector1 = self._ag._coords[acsi, two] - self._ag._coords[acsi, one]
        vector2 = self._ag._coords[acsi, three] - self._ag._coords[acsi, two]
        vector3 = self._ag._coords[acsi, four] - self._ag._coords[acsi, three]
        return vector1, vector2, vector3

    def getACSIndex(self):
        """Returns index of the coordinate set."""

        acsi = self._acsi
        if acsi >= self._ag._n_csets:
            raise ValueError('{0} has fewer coordsets than assumed by {1}'
                             .format(str(self._ag), str(self)))
        return acsi

    def setACSIndex(self, index):
        """Set the coordinate set at *index* active."""

        if self._ag._coords is None:
            raise AttributeError('coordinates are not set')

        if not isinstance(index, Integral):
            raise TypeError('index must be an integer')

        n_csets = self._ag._n_csets
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')

        if index < 0:
            index += n_csets

        self._acsi = index


def evalDihedrals(dihedrals, n_atoms):
    """Returns an array mapping atoms to their dihedrald neighbors and an array
    that stores number of dihedrals made by each atom."""

    numdihedrals = np.bincount(dihedrals.reshape((dihedrals.shape[0] * 4)))
    dmap = np.zeros((n_atoms, numdihedrals.max(), 3), int)
    dmap.fill(-1)
    index = np.zeros(n_atoms, int)
    for dihedral in dihedrals:
        a, b, c, d = dihedral
        dmap[a, index[a]] = [b, c, d]
        dmap[b, index[b]] = [a, c, d]
        dmap[c, index[c]] = [a, b, d]
        dmap[d, index[d]] = [a, b, c]
        index[dihedral] += 1
    return dmap, numdihedrals


def trimDihedrals(dihedrals, indices):
    """Returns dihedrals between atoms at given indices."""

    iset = set(indices)
    dihedrals = [dihedral for dihedral in dihedrals if dihedral[0]
                 in iset and dihedral[1] in iset and dihedral[2] in iset and dihedral[3] in iset]
    if dihedrals:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(dihedrals)]
