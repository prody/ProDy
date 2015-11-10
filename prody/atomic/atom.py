# -*- coding: utf-8 -*-
"""This module defines classes to handle individual atoms."""

import numpy as np

from . import flags
from .fields import ATOMIC_FIELDS, READONLY
from .fields import wrapGetMethod, wrapSetMethod
from .pointer import AtomPointer
from .bond import Bond

__all__ = ['Atom']


class Atom(AtomPointer):

    """A class for handling individual atoms in an :class:`.AtomGroup`."""

    __slots__ = ['_ag', '_acsi', '_index']

    def __init__(self, ag, index, acsi):
        AtomPointer.__init__(self, ag, acsi)
        self._index = int(index)

    def __repr__(self):

        n_csets = self._ag.numCoordsets()
        if n_csets == 1:
            return '<Atom: {0} from {1} (index {2})>'.format(
                   self.getName(), self._ag.getTitle(), self._index)
        elif n_csets > 1:
            return ('<Atom: {0} from {1} (index {2}; active #{3} of '
                    '{4} coordsets)>').format(self.getName(),
                     self._ag.getTitle(), self._index, self.getACSIndex(),
                     n_csets)
        else:
            return ('<Atom: {0} from {1} (index {2}; no coordinates)>'
                    ).format(self.getName(), self._ag.getTitle(), self._index)

    def __str__(self):

        return 'Atom {0} (index {1})'.format(self.getName(), self._index)

    def __len__(self):

        return 1

    def __int__(self):

        return self._index

    def numAtoms(self, flag=None):
        """Return number of atoms, or number of atoms with given *flag*."""

        return len(self._getSubset(flag)) if flag else 1

    def getIndex(self):
        """Return index of the atom."""

        return self._index

    def getIndices(self):
        """Return index of the atom in an :class:`numpy.ndarray`."""

        return np.array([self._index])

    _getIndices = getIndices

    def iterAtoms(self):
        """Yield atoms."""

        yield Atom(ag=self._ag, index=self._index, acsi=self.getACSIndex())

    __iter__ = iterAtoms

    def getCoords(self):
        """Return a copy of coordinates of the atom from the active coordinate
        set."""

        if self._ag._coords is not None:
            return self._ag._coords[self.getACSIndex(), self._index].copy()

    def _getCoords(self):
        """Return a view of coordinates of the atom from the active coordinate
        set."""

        if self._ag._coords is not None:
            return self._ag._coords[self.getACSIndex(), self._index]

    def setCoords(self, coords):
        """Set coordinates of the atom in the active coordinate set."""

        acsi = self.getACSIndex()
        self._ag._coords[acsi, self._index] = coords
        self._ag._setTimeStamp(acsi)

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate set(s) at given *indices*."""

        if self._ag._coords is None:
            return None

        if indices is None:
            return self._ag._coords[:, self._index].copy()

        if isinstance(indices, (int, slice)):
            return self._ag._coords[indices, self._index].copy()

        if isinstance(indices, (list, np.ndarray)):
            return self._ag._coords[indices, self._index]

        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')

    def _getCoordsets(self, indices=None):
        """Return a view of coordinate set(s) at given *indices*."""

        if self._ag._coords is None:
            return None

        if indices is None:
            indices = slice(None)

        return self._ag._coords[indices, self._index]

    def iterCoordsets(self):
        """Yield copies of coordinate sets."""

        for i in range(self.numCoordsets()):
            yield self._ag._coords[i, self._index].copy()

    def _iterCoordsets(self):
        """Yield views of coordinate sets."""

        for i in range(self.numCoordsets()):
            yield self._ag._coords[i, self._index]

    def getMassess(self):
        """get the mass atom. """
        mass_dict = {'C':12,'N':14,'S':32,'O':16,'H':1}
        
        return mass_dict[self.getElement()]

    def getData(self, label):
        """Return a copy of data associated with *label*, if it is present."""

        try:
            data = self._ag._getData(label)
        except KeyError:
            pass
        else:
            if data.ndim > 1:
                return data[self._index]
            else:
                return data[self._index].copy()

    _getData = getData

    def setData(self, label, data):
        """Update *data* associated with *label*.

        :raise AttributeError: when *label* is not in use or read-only"""

        if label in READONLY:
            raise AttributeError('{0} is read-only'.format(repr(label)))
        if label in ATOMIC_FIELDS:
            getattr(self, 'set' + ATOMIC_FIELDS[label].meth)(data)
        else:
            try:
                self._ag._data[label][self._index] = data
            except KeyError:
                raise AttributeError('data with label {0} must be set for'
                                       ' AtomGroup first'.format(repr(label)))

    def getFlag(self, label):
        """Return atom flag."""

        return self._ag._getFlags(label)[self._index]

    def setFlag(self, label, value):
        """Update flag associated with *label*.

         :raise AttributeError: when *label* is not in use or read-only"""

        if label in flags.PLANTERS:
            raise AttributeError('flag {0} cannot be changed by user'
                                    .format(repr(label)))
        flags = self._ag._getFlags(label)
        if flags is None:
            raise AttributeError('flags with label {0} must be set for '
                                    'AtomGroup first'.format(repr(label)))
        flags[self._index] = value

    def getSelstr(self):
        """Return selection string that will select this atom."""

        return 'index {0}'.format(self._index)

    def numBonds(self):
        """Return number of bonds formed by this atom.  Bonds must be set first
        using :meth:`.AtomGroup.setBonds`."""

        numbonds = self._ag._data.get('numbonds')
        if numbonds is not None:
            return numbonds[self._index]

    def iterBonds(self):
        """Yield bonds formed by the atom.  Use :meth:`setBonds` for setting
        bonds."""

        ag = self._ag
        acsi = self.getACSIndex()
        for bond in self._iterBonds():
            yield Bond(ag, bond, acsi)

    def _iterBonds(self):
        """Yield pairs of bonded atom indices."""

        ag = self._ag
        if ag._bmap is None:
            raise ValueError('bonds are not set, use `AtomGroup.setBonds`')

        this = self._index
        for other in ag._bmap[this]:
            if other == -1:
                break
            yield this, other

    def iterBonded(self):
        """Yield bonded atoms.  Use :meth:`setBonds` for setting bonds."""

        ag = self._ag
        if ag._bmap is None:
            raise ValueError('bonds are not set, use `AtomGroup.setBonds`')

        acsi = self.getACSIndex()
        this = self._index
        for other in self._ag._bmap[this]:
            if other == -1:
                break
            yield Atom(ag, other, acsi)


for fname, field in ATOMIC_FIELDS.items():

    if field.private:
        continue

    meth = field.meth
    getMeth = 'get' + meth
    setMeth = 'set' + meth
    # Define public method for retrieving a copy of data array
    def getData(self, meth=field.meth_pl, call=field.call):
        data = getattr(self._ag, '_get' + meth)()
        if data is not None:
            return data[self._index]
    getData = wrapGetMethod(getData)
    getData.__name__ = getMeth
    getData.__doc__ = field.getDocstr('get', False)
    setattr(Atom, getMeth, getData)
    setattr(Atom, '_' + getMeth, getData)

    if field.readonly:
        continue

    # Define public method for setting values in data array
    def setData(self, value, var=fname, none=field.none):
        array = self._ag._data[var]
        if array is None:
            raise AttributeError('attribute of the AtomGroup is '
                                 'not set')
        array[self._index] = value
        if none: self._ag._none(none)
    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth
    setData.__doc__ = field.getDocstr('set', False)
    setattr(Atom, setMeth, setData)

del getData
del setData
