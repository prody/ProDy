# -*- coding: utf-8 -*-
"""This module defines :class:`Improper` for dealing with improper information provided
by using :meth:`.AtomGroup.setImpropers` method."""

from numbers import Integral
import numpy as np

RAD2DEG = 180 / np.pi

__all__ = ['Improper']

class Improper(object):

    """A pointer class for improperd atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns improper length, i.e. :meth:`getSize`
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
        return '<Improper: {0}({1})--{2}({3})--{4}({5}--{6}({7})) from {8}>'.format(
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

        return (isinstance(other, Improper) and other.getAtomGroup() is self._ag
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
        """Returns improperd atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]], self._ag[self._indices[2]])

    def getIndices(self):
        """Returns indices of improperd atoms."""

        return self._indices.copy()

    def getSize(self):
        """Returns improper size."""

        one, two, three, four = self._indices
        acsi = self.getACSIndex()
        atoms1 = self._ag._coords[acsi, one]
        atoms2 = self._ag._coords[acsi, two]
        atoms3 = self._ag._coords[acsi, three]
        atoms4 = self._ag._coords[acsi, four]
        return calcImproper(atoms1, atoms2, atoms3, atoms4)

    def getVectors(self):
        """Returns bond vectors that originate from the central atom."""

        one, two, three, four = self._indices
        acsi = self.getACSIndex()
        vector1 = self._ag._coords[acsi, one] - self._ag._coords[acsi, two]
        vector2 = self._ag._coords[acsi, two] - self._ag._coords[acsi, three]
        vector3 = self._ag._coords[acsi, three] - self._ag._coords[acsi, four]
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


def evalImpropers(impropers, n_atoms):
    """Returns an array mapping atoms to their improperd neighbors and an array
    that stores number of impropers made by each atom."""

    numimpropers = np.bincount(impropers.reshape((impropers.shape[0] * 4)))
    imap = np.zeros((n_atoms, numimpropers.max(), 3), int)
    imap.fill(-1)
    index = np.zeros(n_atoms, int)
    for improper in impropers:
        a, b, c, d = improper
        imap[a, index[a]] = [b, c, d]
        imap[b, index[b]] = [a, c, d]
        imap[c, index[c]] = [a, b, d]
        imap[d, index[d]] = [a, b, c]
        index[improper] += 1
    return imap, numimpropers


def trimImpropers(impropers, indices):
    """Returns impropers between atoms at given indices."""

    iset = set(indices)
    impropers = [improper for improper in impropers if improper[0]
                 in iset and improper[1] in iset and improper[2] in iset and improper[3] in iset]
    if impropers:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(impropers)]
