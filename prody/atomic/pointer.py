# -*- coding: utf-8 -*-
"""This module defines atom pointer base class."""

from numbers import Integral
from numpy import all, array, concatenate, ones, unique, any

from .atomic import Atomic
from .bond import Bond
from .angle import Angle
from .dihedral import Dihedral
from .improper import Improper
from .donor import Donor
from .acceptor import Acceptor
from .crossterm import Crossterm
from .nbexclusion import NBExclusion

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
        """Returnss an :class:`.AtomMap` instance. Order of pointed atoms are
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
        """Returns KDTree for the active coordinate set from the atom group."""

        return self._ag._getKDTree(self.getACSIndex())

    def getAtomGroup(self):
        """Returns associated atom group."""

        return self._ag

    def numCoordsets(self):
        """Returns number of coordinate sets."""

        return self._ag._n_csets

    def getACSIndex(self):
        """Returns index of the coordinate set."""

        acsi = self._acsi
        if acsi >= self._ag._n_csets:
            raise ValueError('{0} has fewer coordsets than assumed by {1}'
                             .format(str(self._ag), str(self)))
        return acsi

    def setACSIndex(self, index):
        """Set coordinates at *index* active."""

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

    def getACSLabel(self):
        """Returns active coordinate set label."""

        if self._ag._n_csets:
            return self._ag._cslabels[self.getACSIndex()]

    def getCSLabels(self):
        """Returns coordinate set labels."""

        return self._ag.getCSLabels()

    def isDataLabel(self, label):
        """Returns **True** if data associated with *label* is present."""

        return self._ag.isDataLabel(label)

    def getDataLabels(self, which=None):
        """Returns data labels.  For ``which='user'``, return only labels of
        user provided data."""

        return self._ag.getDataLabels(which)

    def getDataType(self, label):
        """Returns type of the data (i.e. ``data.dtype``) associated with
        *label*, or **None** label is not used."""

        return self._ag.getDataType(label)

    def getFlagLabels(self, which=None):
        """Returns flag labels.  For ``which='user'``,  return labels of user
        or parser (e.g. :term:`hetatm`) provided flags, for ``which='all'``
        return all possible :ref:`flags` labels in addition to those present
        in the instance."""

        return self._ag.getFlagLabels(which)

    def isFlagLabel(self, label):
        """Returns **True** if flags associated with *label* are present."""

        return self._ag.isFlagLabel(label)

    def _getFlags(self, label):
        """Returns atom flags."""

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
            LOGGER.warning('bonds are not set, use `setBonds` or `inferBonds`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 2 >= len(self):
            for a, b in self._ag._iterBonds():
                if a in iset and b in iset:
                    yield a, b
        else:
            if any(self._ag._bmap):
                for a, bmap in zip(indices, self._ag._bmap[indices]):
                    for b in bmap:
                        if b > -1 and b in iset:
                            yield a, b
                    iset.remove(a)

    def iterBonds(self):
        """Yield bonds formed by the atom.  Use :meth:`setBonds` or 
        :meth:`inferBonds` for setting bonds."""

        ag = self._ag
        acsi = self.getACSIndex()
        for bond in self._iterBonds():
            yield Bond(ag, bond, acsi)

    def _iterAngles(self):
        """Yield triplets of indices for angled atoms that are within the pointer.
        Use :meth:`setAngles` for setting angles."""

        if self._ag._angles is None:
            LOGGER.warning('angles are not set, use `AtomGroup.setAngles`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 3 >= len(self):
            for a, b, c in self._ag._iterAngles():
                if a in iset and b in iset and c in iset:
                    yield a, b, c
        else:
            if any(self._ag._angmap):
                for a, amap in zip(indices, self._ag._angmap[indices]):
                    for b, c in amap:
                        if b > -1 and b in iset and c > -1 and c in iset:
                            yield a, b, c
                    iset.remove(a)

    def iterAngles(self):
        """Yield angles formed by the atom.  Use :meth:`setAngles` for setting
        angles."""

        ag = self._ag
        acsi = self.getACSIndex()
        for angle in self._iterAngles():
            yield Angle(ag, angle, acsi)

    def _iterDihedrals(self):
        """Yield quadruples of indices for dihedraled atoms that are within the pointer.
        Use :meth:`setDihedrals` for setting dihedrals."""

        if self._ag._dihedrals is None:
            LOGGER.warning('dihedrals are not set, use `AtomGroup.setDihedrals`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 4 >= len(self):
            for a, b, c, d in self._ag._iterDihedrals():
                if a in iset and b in iset and c in iset and d in iset:
                    yield a, b, c, d
        else:
            if any(self._ag._dmap):
                for a, dmap in zip(indices, self._ag._dmap[indices]):
                    for b, c, d in dmap:
                        if b > -1 and b in iset and c > -1 and c in iset \
                        and d > -1 and d in iset:
                            yield a, b, c, d
                    iset.remove(a)

    def iterDihedrals(self):
        """Yield dihedrals formed by the atom.  Use :meth:`setDihedrals` for setting
        dihedrals."""

        ag = self._ag
        acsi = self.getACSIndex()
        for dihedral in self._iterDihedrals():
            yield Dihedral(ag, dihedral, acsi) 

    def _iterImpropers(self):
        """Yield quadruplet of indices for impropered atoms that are within the pointer.
        Use :meth:`setImpropers` for setting impropers."""

        if self._ag._impropers is None:
            LOGGER.warning('impropers are not set, use `AtomGroup.setImpropers`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 4 >= len(self):
            for a, b, c, d in self._ag._iterImpropers():
                if a in iset and b in iset and c in iset and d in iset:
                    yield a, b, c, d
        else:
            if any(self._ag._imap):
                for a, imap in zip(indices, self._ag._imap[indices]):
                    for b, c, d in imap:
                        if b > -1 and b in iset and c > -1 and c in iset \
                        and d > -1 and d in iset:
                            yield a, b, c, d
                    iset.remove(a)

    def iterImpropers(self):
        """Yield impropers formed by the atom.  Use :meth:`setImpropers` for setting
        impropers."""

        ag = self._ag
        acsi = self.getACSIndex()
        for improper in self._iterImpropers():
            yield Improper(ag, improper, acsi) 

    def _iterCrossterms(self):
        """Yield quadruplet of indices for crosstermed atoms that are within the pointer.
        Use :meth:`setCrossterms` for setting crossterms."""

        if self._ag._crossterms is None:
            LOGGER.warning('crossterms are not set, use `AtomGroup.setCrossterms`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 4 >= len(self):
            for a, b, c, d in self._ag._iterCrossterms():
                if a in iset and b in iset and c in iset and d in iset:
                    yield a, b, c, d
        else:
            if any(self._ag._cmap):
                for a, cmap in zip(indices, self._ag._cmap[indices]):
                    for b, c, d in cmap:
                        if b > -1 and b in iset and c > -1 and c in iset \
                        and d > -1 and d in iset:
                            yield a, b, c, d
                    iset.remove(a)

    def iterCrossterms(self):
        """Yield crossterms formed by the atom.  Use :meth:`setCrossterms` for setting
        crossterms."""

        ag = self._ag
        acsi = self.getACSIndex()
        for crossterm in self._iterCrossterms():
            yield Crossterm(ag, crossterm, acsi) 

    def _iterDonors(self):
        """Yield pairs of indices for donored atoms that are within the pointer.
        Use :meth:`setDonors` for setting donors."""

        if self._ag._donors is None:
            LOGGER.warning('donors are not set, use `AtomGroup.setDonors`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 2 >= len(self):
            for a, b in self._ag._iterDonors():
                if a in iset and b in iset:
                    yield a, b
        else:
            if any(self._ag._domap):
                for a, dmap in zip(indices, self._ag._domap[indices]):
                    for b in dmap:
                        if b > -1 and b in iset:
                            yield a, b
                    iset.remove(a)

    def iterDonors(self):
        """Yield donors formed by the atom.  Use :meth:`setDonors` for setting
        donors."""

        ag = self._ag
        acsi = self.getACSIndex()
        for donor in self._iterDonors():
            yield Donor(ag, donor, acsi)

    def _iterAcceptors(self):
        """Yield pairs of indices for acceptored atoms that are within the pointer.
        Use :meth:`setAcceptors` for setting acceptors."""

        if self._ag._acceptors is None:
            LOGGER.warning('acceptors are not set, use `AtomGroup.setAcceptors`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 2 >= len(self):
            for a, b in self._ag._iterAcceptors():
                if a in iset and b in iset:
                    yield a, b
        else:
            if any(self._ag._acmap):
                for a, amap in zip(indices, self._ag._acmap[indices]):
                    for b in amap:
                        if b > -1 and b in iset:
                            yield a, b
                    iset.remove(a)

    def iterAcceptors(self):
        """Yield acceptors formed by the atom.  Use :meth:`setAcceptors` for setting
        acceptors."""

        ag = self._ag
        acsi = self.getACSIndex()
        for acceptor in self._iterAcceptors():
            yield Acceptor(ag, acceptor, acsi)

    def _iterNBExclusions(self):
        """Yield pairs of indices for nbexclusioned atoms that are within the pointer.
        Use :meth:`setNBExclusions` for setting nbexclusions."""

        if self._ag._nbexclusions is None:
            LOGGER.warning('nbexclusions are not set, use `AtomGroup.setNBExclusions`')

        indices = self._getIndices()
        iset = set(indices)
        if len(self._ag) / 2 >= len(self):
            for a, b in self._ag._iterNBExclusions():
                if a in iset and b in iset:
                    yield a, b
        else:
            if any(self._ag._nbemap):
                for a, nbemap in zip(indices, self._ag._nbemap[indices]):
                    for b in nbemap:
                        if b > -1 and b in iset:
                            yield a, b
                    iset.remove(a)

    def iterNBExclusions(self):
        """Yield nbexclusions formed by the atom.  Use :meth:`setNBExclusions` for setting
        nbexclusions."""

        ag = self._ag
        acsi = self.getACSIndex()
        for nbexclusion in self._iterNBExclusions():
            yield NBExclusion(ag, nbexclusion, acsi)
