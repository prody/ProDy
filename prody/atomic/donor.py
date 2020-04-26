# -*- coding: utf-8 -*-
"""This module defines :class:`Donor` for dealing with donor information provided
by using :meth:`.AtomGroup.setDonors` method."""

from numbers import Integral
import numpy as np

__all__ = ['Donor']

class Donor(object):

    """A pointer class for donored atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns donor length, i.e. :meth:`getLength`
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
        return '<Donor: {0}({1})--{2}({3}) from {4}>'.format(
                            names[one], one, names[two], two, str(self._ag))

    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})'.format(
                                            names[one], one, names[two], two)

    def __eq__(self, other):

        return (isinstance(other, Donor) and other.getAtomGroup() is self._ag
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
        """Returns donored atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Returns indices of donored atoms."""

        return self._indices.copy()

    def getLength(self):
        """Returns donor length."""

        vector = self.getVector()
        return np.multiply(vector, vector, vector).sum() ** 0.5

    def getVector(self):
        """Returns donor vector that originates from the first atom."""

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


def evalDonors(donors, n_atoms):
    """Returns an array mapping atoms to their common donor neighbors and an array
    that stores number of donors made by each atom.
    
    This excludes entries corresponding to missing hydrogens (with 0 in them that became -1)."""

    if donors.min() < 0:
        # remove entries corresponding to missing hydrogens (with 0 in them that became -1)
        donors = donors[np.argwhere(donors == -1)[-1][0]+1:]

    numdonors = np.bincount(donors.reshape((donors.shape[0] * 2)))        

    domap = np.zeros((n_atoms, numdonors.max()), int)
    index = np.zeros(n_atoms, int)

    domap.fill(-1)
    
    for donor in donors:
        a, b = donor
        domap[a, index[a]] = b
        domap[b, index[b]] = a
        index[donor] += 1
    return domap, numdonors


def trimDonors(donors, indices):
    """Returns donors between atoms at given indices."""

    iset = set(indices)
    donors = [donor for donor in donors if donor[0] in iset and donor[1] in iset]
    if donors:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(donors)]
