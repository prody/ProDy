# -*- coding: utf-8 -*-
"""This module defines :class:`Bond` for dealing with bond information provided
by using :meth:`.AtomGroup.setBonds` method."""

import numpy as np

__all__ = ['Bond']

class Bond(object):

    """A pointer class for bonded atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns bond length, i.e. :meth:`getLength`
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

        one, two = self._indices
        names = self._ag._getNames()
        return '<Bond: {0}({1})--{2}({3}) from {4}>'.format(
                            names[one], one, names[two], two, str(self._ag))

    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})'.format(
                                            names[one], one, names[two], two)

    def __eq__(self, other):

        return (isinstance(other, Bond) and other.getAtomGroup() is self._ag
                and (np.all(other.getIndices() == self._indices) or
                 np.all(other.getIndices() == list(reversed(self._indices)))))

    def __ne__(self, other):

        return not self.__eq__(other)

    def __len__(self):

        return self.getLength()

    def __iter__(self):

        for index in self._indices:
            yield self._ag[index]

    def getAtomGroup(self):
        """Returns atom group."""

        return self._ag

    def getAtoms(self):
        """Returns bonded atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Returns indices of bonded atoms."""

        return self._indices.copy()

    def getLength(self):
        """Returns bond length."""

        vector = self.getVector()
        return np.multiply(vector, vector, vector).sum() ** 0.5

    def getVector(self):
        """Returns bond vector that originates from the first atom."""

        one, two = self._indices
        acsi = self.getACSIndex()
        return self._ag._coords[acsi, two] - self._ag._coords[acsi, one]

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

        if not isinstance(index, int):
            raise TypeError('index must be an integer')

        n_csets = self._ag._n_csets
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')

        if index < 0:
            index += n_csets

        self._acsi = index


def evalBonds(bonds, n_atoms):
    """Returns an array mapping atoms to their bonded neighbors and an array
    that stores number of bonds made by each atom."""

    numbonds = np.bincount(bonds.reshape((bonds.shape[0] * 2)))
    bmap = np.zeros((n_atoms, numbonds.max()), int)
    bmap.fill(-1)
    index = np.zeros(n_atoms, int)
    for bond in bonds:
        a, b = bond
        bmap[a, index[a]] = b
        bmap[b, index[b]] = a
        index[bond] += 1
    return bmap, numbonds


def trimBonds(bonds, indices):
    """Returns bonds between atoms at given indices."""

    iset = set(indices)
    bonds = [bond for bond in bonds if bond[0] in iset and bond[1] in iset]
    if bonds:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(bonds)]
