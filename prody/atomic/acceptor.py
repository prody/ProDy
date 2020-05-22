# -*- coding: utf-8 -*-
"""This module defines :class:`Acceptor` for dealing with acceptor information provided
by using :meth:`.AtomGroup.setAcceptors` method."""

from numbers import Integral
import numpy as np

__all__ = ['Acceptor']

class Acceptor(object):

    """A pointer class for acceptored atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns acceptor length, i.e. :meth:`getLength`
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
        return '<Acceptor: {0}({1})--{2}({3}) from {4}>'.format(
                            names[one], one, names[two], two, str(self._ag))

    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})'.format(
                                            names[one], one, names[two], two)

    def __eq__(self, other):

        return (isinstance(other, Acceptor) and other.getAtomGroup() is self._ag
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
        """Returns acceptored atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Returns indices of acceptored atoms."""

        return self._indices.copy()

    def getLength(self):
        """Returns acceptor length."""

        vector = self.getVector()
        return np.multiply(vector, vector, vector).sum() ** 0.5

    def getVector(self):
        """Returns acceptor vector that originates from the first atom."""

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


def evalAcceptors(acceptors, n_atoms):
    """Returns an array mapping atoms to their common acceptor neighbors and an array
    that stores number of acceptors made by each atom.
    
    This excludes entries corresponding to missing hydrogens (with 0 in them that became -1)."""

    if acceptors.min() < 0:
        # remove entries corresponding to missing hydrogens (with 0 in them that became -1)
        acceptors = acceptors[np.argwhere(acceptors == -1)[-1][0]+1:]

    numacceptors = np.bincount(acceptors.reshape((acceptors.shape[0] * 2)))

    acmap = np.zeros((n_atoms, numacceptors.max()), int)
    index = np.zeros(n_atoms, int)
    
    acmap.fill(-1)
    
    for acceptor in acceptors:
        a, b = acceptor
        acmap[a, index[a]] = b
        acmap[b, index[b]] = a
        index[acceptor] += 1
    return acmap, numacceptors


def trimAcceptors(acceptors, indices):
    """Returns acceptors between atoms at given indices."""

    iset = set(indices)
    acceptors = [acceptor for acceptor in acceptors if acceptor[0] in iset and acceptor[1] in iset]
    if acceptors:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(acceptors)]
