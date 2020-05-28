# -*- coding: utf-8 -*-
"""This module defines :class:`NBExclusion` for dealing with bond information provided
by using :meth:`.AtomGroup.setNBExclusions` method."""

from numbers import Integral
import numpy as np

__all__ = ['NBExclusion']


class NBExclusion(object):

    """A pointer class for nonnonbonded exclusion exclusion atoms.  Following built-in functions are
    customized for this class:

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
        return '<NBExclusion: {0}({1})--{2}({3}) from {4}>'.format(
            names[one], one, names[two], two, str(self._ag))

    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})'.format(
            names[one], one, names[two], two)

    def __eq__(self, other):

        return (isinstance(other, NBExclusion) and other.getAtomGroup() is self._ag
                and (np.all(other.getIndices() == self._indices) or
                     np.all(other.getIndices() == list(reversed(self._indices)))))

    def __ne__(self, other):

        return not self.__eq__(other)

    def __iter__(self):

        for index in self._indices:
            yield self._ag[index]

    def getAtomGroup(self):
        """Returns atom group."""

        return self._ag

    def getAtoms(self):
        """Returns nonbonded exclusion atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Returns indices of nonbonded exclusion atoms."""

        return self._indices.copy()

    def getVector(self):
        """Returns vector that originates from the first atom."""

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

        if not isinstance(index, Integral):
            raise TypeError('index must be an integer')

        n_csets = self._ag._n_csets
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')

        if index < 0:
            index += n_csets

        self._acsi = index


def evalNBExclusions(exclusions, n_atoms):
    """Returns an array mapping atoms to their nonbonded exclusion neighbors and an array
    that stores number of nonbonded exclusions made by each atom."""

    numexclusions = np.bincount(
        exclusions.reshape((exclusions.shape[0] * 2)))
    nbemap = np.zeros((n_atoms, numexclusions.max()), int)
    nbemap.fill(-1)
    index = np.zeros(n_atoms, int)
    for nbexclusion in exclusions:
        a, b = nbexclusion
        nbemap[a, index[a]] = b
        nbemap[b, index[b]] = a
        index[nbexclusion] += 1
    return nbemap, numexclusions


def trimNBExclusions(exclusions, indices):
    """Returns nonbonded exclusions between atoms at given indices."""

    iset = set(indices)
    exclusions = [nbexclusion for nbexclusion in exclusions if nbexclusion[0]
                    in iset and nbexclusion[1] in iset]
    if exclusions:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(exclusions)]
