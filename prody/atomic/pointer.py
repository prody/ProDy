# -*- coding: utf-8 -*-
"""This module defines atom pointer base class."""

from numpy import all, array, concatenate, ones, unique

from .atomic import Atomic
from .bond import Bond

from prody import LOGGER

__all__ = ['AtomPointer']


class AtomPointer(Atomic):

    """A base for classes pointing to atoms in :class:`.AtomGroup` instances.
    Derived classes are:

      * :class:`.Atom`
      * :class:`.AtomSubset`
      * :class:`.AtomMap`"""

    __slots__ = ['_ag', '_acsi']

    def __init__(self, ag, acsi):

        if not isinstance(ag, AtomGroup):
            raise TypeError('ag must be an AtomGroup instance, not {0}'
                            .format(type(ag)))
        self._ag = ag
        if acsi is None:
            self._acsi = ag.getACSIndex()
        else:
            self._acsi = int(acsi)

    def __contains__(self, item):

        try:
            ag = item.getAtomGroup()
        except AttributeError:
            return False
        else:
            return (self._ag == ag and len(item) <= len(self) and
                    set(item._getIndices()).issubset(set(self._getIndices())))

    def __eq__(self, other):

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            return False
        else:
            return (self._ag == ag and other.numAtoms() == self.numAtoms() and
                    all(self._getIndices() == other._getIndices()))

    def __ne__(self, other):

        return not self.__eq__(other)

    def __invert__(self):

        torf = ones(self._ag.numAtoms(), bool)
        torf[self._indices] = False
        return Selection(self._ag, torf.nonzero()[0],
                         'not ({0})'.format(self.getSelstr()),
                         self.getACSIndex(), unique=True)

    def __or__(self, other):

        if self is other:
            return self

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('other must be an AtomPointer')

        if self._ag != ag:
            raise ValueError('both selections must be from the same AtomGroup')

        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warn('Active coordinate set indices do not match, it will '
                        'be set to zero.')
            acsi = 0

        indices = unique(concatenate((self._getIndices(),
                                      other._getIndices())))
        if indices[-1] == atommap.DUMMY:
            indices = indices[:-1]
        return Selection(self._ag, indices, '({0}) or ({1})'
                         .format(self.getSelstr(), other.getSelstr()),
                         acsi, unique=True)

    def __and__(self, other):

        if self is other:
            return self

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('other must be an AtomPointer')

        if self._ag != ag:
            raise ValueError('both selections must be from the same AtomGroup')

        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0

        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warn('Active coordinate set indices do not match, it will '
                        'be set to zero.')
            acsi = 0

        indices = set(self._getIndices())

        indices = indices.intersection(other.getIndices())
        if indices:
            indices = unique(indices)
            if indices[-1] == atommap.DUMMY:
                indices = indices[:-1]
            return Selection(self._ag, indices, '({0}) and ({1})'
                             .format(self.getSelstr(), other.getSelstr()),
                             acsi)

    def __add__(self, other):
        """Returns an :class:`.AtomMap` instance. Order of pointed atoms are
        preserved."""

        try:
            ag = other.getAtomGroup()
        except AttributeError:
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        if ag != self._ag:
            raise ValueError('AtomPointer instances must point to the same '
                             'AtomGroup instance')
        acsi = self.getACSIndex()
        if acsi != other.getACSIndex():
            LOGGER.warning('Active coordset indices of atoms are not the same.'
                           ' Result will have ACSI {0}.'.format(acsi))

        title = '({0}) + ({1})'.format(str(self), str(other))
        indices = concatenate([self._getIndices(), other._getIndices()])

        dummies = 0
        try:
            dummies += self.numDummies()
        except AttributeError:
            pass
        try:
            dummies += other.numDummies()
        except AttributeError:
            pass

        return AtomMap(ag, indices, acsi, title=title, intarrays=True,
                       dummies=dummies)

    def _getTimeStamp(self, index=None):

        if index is None:
            index = self.getACSIndex()
        return self._ag._getTimeStamp(index)

    def _getKDTree(self):
        """Return KDTree for the active coordinate set from the atom group."""

        return self._ag._getKDTree(self.getACSIndex())

    def getAtomGroup(self):
        """Return associated atom group."""

        return self._ag

    def numCoordsets(self):
        """Return number of coordinate sets."""

        return self._ag._n_csets

    def getACSIndex(self):
        """Return index of the coordinate set."""

        acsi = self._acsi
        if acsi >= self._ag._n_csets:
            raise ValueError('{0} has fewer coordsets than assumed by {1}'
                             .format(str(self._ag), str(self)))
        return acsi

    def setACSIndex(self, index):
        """Set coordinates at *index* active."""

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

    def getACSLabel(self):
        """Return active coordinate set label."""

        if self._ag._n_csets:
            return self._ag._cslabels[self.getACSIndex()]

    def getCSLabels(self):
        """Return coordinate set labels."""

        return self._ag.getCSLabels()

    def isDataLabel(self, label):
        """Return **True** if data associated with *label* is present."""

        return self._ag.isDataLabel(label)

    def getDataLabels(self, which=None):
        """Return data labels.  For ``which='user'``, return only labels of
        user provided data."""

        return self._ag.getDataLabels(which)

    def getDataType(self, label):
        """Return type of the data (i.e. ``data.dtype``) associated with
        *label*, or **None** label is not used."""

        return self._ag.getDataType(label)

    def getFlagLabels(self, which=None):
        """Return flag labels.  For ``which='user'``,  return labels of user
        or parser (e.g. :term:`hetatm`) provided flags, for ``which='all'``
        return all possible :ref:`flags` labels in addition to those present
        in the instance."""

        return self._ag.getFlagLabels(which)

    def isFlagLabel(self, label):
        """Return **True** if flags associated with *label* are present."""

        return self._ag.isFlagLabel(label)

    def _getFlags(self, label):
        """Return atom flags."""

        flags = self._ag._getFlags(label)
        if flags is not None:
            return flags[self._getIndices()]

    def _getSubset(self, label):

        subset = array(list(set(self._ag._getSubset(label))
                            .intersection(set(self._getIndices()))), int)
        subset.sort()
        return subset

    def _iterBonds(self):
        """Yield pairs of indices for bonded atoms that are within the pointer.
        Use :meth:`setBonds` for setting bonds."""

        if self._ag._bonds is None:
            raise ValueError('bonds are not set, use `AtomGroup.setBonds`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) * 0.5 >= len(self):
            for a, b in self._ag._iterBonds():
                if a in iset and b in iset:
                    yield a, b
        else:
            for a, bmap in zip(indices, self._ag._bmap[indices]):
                for b in bmap:
                    if b > -1 and b in iset:
                        yield a, b
                iset.remove(a)
