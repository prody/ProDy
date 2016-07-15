# -*- coding: utf-8 -*-
import numpy as np

from . import flags
from .atom import Atom
from .fields import ATOMIC_FIELDS, READONLY
from .fields import wrapGetMethod, wrapSetMethod
from .pointer import AtomPointer
from prody import LOGGER

__all__ = ['AtomSubset']


class AtomSubset(AtomPointer):

    """A class for manipulating subset of atoms in an :class:`.AtomGroup`.
    Derived classes are:

      * :class:`.Selection`
      * :class:`.Segment`
      * :class:`.Chain`
      * :class:`.Residue`

    This class stores a reference to an :class:`.AtomGroup` instance, a set of
    atom indices, and active coordinate set index for the atom group."""

    __slots__ = ['_ag', '_indices', '_acsi', '_selstr']

    def __init__(self, ag, indices, acsi, **kwargs):

        AtomPointer.__init__(self, ag, acsi)

        if not isinstance(indices, np.ndarray):
            indices = np.array(indices, int)
        elif not np.issubdtype(indices.dtype, int):
            indices = list(indices[0])
            indices = np.array(indices, int)

        if kwargs.get('unique'):
            self._indices = indices
        else:
            self._indices = np.unique(indices)

        self._selstr = kwargs.get('selstr')

    def __len__(self):

        return len(self._indices)

    def getCoords(self):
        """Returns a copy of coordinates from the active coordinate set."""

        if self._ag._coords is not None:
            # Since this is not slicing, a view is not returned
            return self._ag._coords[self.getACSIndex(), self._indices]

    _getCoords = getCoords

    def setCoords(self, coords):
        """Set coordinates in the active coordinate set."""

        if self._ag._coords is not None:
            self._ag._coords[self.getACSIndex(), self._indices] = coords
            self._ag._setTimeStamp(self.getACSIndex())

    def getCoordsets(self, indices=None):
        """Returns coordinate set(s) at given *indices*, which may be an integer
        or a list/array of integers."""

        if self._ag._coords is None:
            return None
        if indices is None:
            return self._ag._coords[:, self._indices]
        if isinstance(indices, (int, slice)):
            return self._ag._coords[indices, self._indices]
        if isinstance(indices, (list, np.ndarray)):
            return self._ag._coords[indices, self._indices]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')

    _getCoordsets = getCoordsets

    def iterCoordsets(self):
        """Yield copies of coordinate sets."""

        coords = self._ag._getCoordsets()
        if coords is not None:
            indices = self._indices
            for xyz in coords:
                yield xyz[indices]

    _iterCoordsets = iterCoordsets

    def getIndices(self):
        """Returns a copy of the indices of atoms."""

        return self._indices.copy()

    def _getIndices(self):
        """Returns indices of atoms."""

        return self._indices

    def numAtoms(self, flag=None):
        """Returns number of atoms, or number of atoms with given *flag*."""

        return len(self._getSubset(flag)) if flag else len(self._indices)

    def iterAtoms(self):
        """Yield atoms."""

        ag = self._ag
        acsi = self.getACSIndex()
        for index in self._indices:
            yield Atom(ag=ag, index=index, acsi=acsi)

    __iter__ = iterAtoms

    def getData(self, label):
        """Returns a copy of data associated with *label*, if it is present."""

        data = self._ag._getData(label)
        if data is not None:
            return data[self._indices]

    _getData = getData

    def setData(self, label, data):
        """Update *data* associated with *label*.

        :raise AttributeError: when *label* is not in use or read-only"""

        if label in READONLY:
            raise AttributeError('{0} is read-only'.format(repr(label)))
        if label in ATOMIC_FIELDS:
            getattr(self, 'set' + ATOMIC_FIELDS[label].meth_pl)(data)
        else:
            try:
                self._ag._data[label][self._indices] = data
            except KeyError:
                raise AttributeError('data with label {0} must be set for '
                                     'AtomGroup first'.format(repr(label)))

    def getFlags(self, label):
        """Returns a copy of atom flags for given *label*, or **None** when
        flags for *label* is not set."""

        return self._ag._getFlags(label)[self._indices]

    def setFlags(self, label, value):
        """Update flag associated with *label*.

         :raise AttributeError: when *label* is not in use or read-only"""

        if label in flags.PLANTERS:
            raise AttributeError('flag {0} cannot be changed by user'
                                    .format(repr(label)))
        flags = self._ag._getFlags(label)
        if flags is None:
            raise AttributeError('flags with label {0} must be set for '
                                    'AtomGroup first'.format(repr(label)))
        flags[self._indices] = value


for fname, field in ATOMIC_FIELDS.items():

    if field.private:
        continue

    meth = field.meth_pl
    getMeth = 'get' + meth
    setMeth = 'set' + meth
    # Define public method for retrieving a copy of data array
    def getData(self, meth=field.meth_pl, call=field.call):
        data = getattr(self._ag, '_get' + meth)()
        if data is not None:
            return data[self._indices]
    getData = wrapGetMethod(getData)
    getData.__name__ = getMeth
    getData.__doc__ = field.getDocstr('get')
    setattr(AtomSubset, getMeth, getData)
    setattr(AtomSubset, '_' + getMeth, getData)

    if field.readonly:
        continue

    # Define public method for setting values in data array
    def setData(self, value, var=fname, none=field.none):
        array = self._ag._data[var]
        if array is None:
            raise AttributeError(var + ' data is not set')
        array[self._indices] = value
        if none: self._ag._none(none)
    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth
    setData.__doc__ = field.getDocstr('set')
    setattr(AtomSubset, setMeth, setData)

del getData
del setData
