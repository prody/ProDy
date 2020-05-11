# -*- coding: utf-8 -*-
"""This module defines :class:`Angle` for dealing with angle information provided
by using :meth:`.AtomGroup.setAngles` method."""

from numbers import Integral
import numpy as np

RAD2DEG = 180 / np.pi

__all__ = ['Angle']

class Angle(object):

    """A pointer class for angled atoms.  Following built-in functions are
    customized for this class:

    * :func:`len` returns angle length, i.e. :meth:`getSize`
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

        one, two, three = self._indices
        names = self._ag._getNames()
        return '<Angle: {0}({1})--{2}({3})--{4}({5}) from {6}>'.format(
                            names[one], one, names[two], two, 
                            names[three], three, str(self._ag))

    def __str__(self):

        one, two, three = self._indices
        names = self._ag._getNames()
        return '{0}({1})--{2}({3})--{4}({5})'.format(
                                            names[one], one, names[two], two,
                                            names[three], three)

    def __eq__(self, other):

        return (isinstance(other, Angle) and other.getAtomGroup() is self._ag
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
        """Returns angled atoms."""

        return (self._ag[self._indices[0]], self._ag[self._indices[1]], self._ag[self._indices[2]])

    def getIndices(self):
        """Returns indices of angled atoms."""

        return self._indices.copy()

    def getSize(self, radian=False):
        """Returns angle size."""

        v1, v2 = self.getVectors()
        rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
        if radian:
            return rad
        else:
            return rad * RAD2DEG

    def getVectors(self):
        """Returns bond vectors that originate from the central atom."""

        one, two, three = self._indices
        acsi = self.getACSIndex()
        vector1 = self._ag._coords[acsi, one] - self._ag._coords[acsi, two]
        vector2 = self._ag._coords[acsi, three] - self._ag._coords[acsi, two]
        return vector1, vector2

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


def evalAngles(angles, n_atoms):
    """Returns an array mapping atoms to their angled neighbors and an array
    that stores number of angles made by each atom."""

    numangles = np.bincount(angles.reshape((angles.shape[0] * 3)))
    angmap = np.zeros((n_atoms, numangles.max(), 2), int)
    angmap.fill(-1)
    index = np.zeros(n_atoms, int)
    for angle in angles:
        a, b, c = angle
        angmap[a, index[a]] = [b,c]
        angmap[b, index[b]] = [a,c]
        angmap[c, index[c]] = [a,b]
        index[angle] += 1
    return angmap, numangles


def trimAngles(angles, indices):
    """Returns angles between atoms at given indices."""

    iset = set(indices)
    angles = [angle for angle in angles if angle[0] in iset and angle[1] in iset and angle[2] in iset]
    if angles:
        newindices = np.zeros(indices.max()+1, int)
        newindices[indices] = np.arange(len(indices))
        return newindices[np.array(angles)]
