# -*- coding: utf-8 -*-
"""This module defines :class:`AtomMap` class that allows for pointing atoms in
arbitrary order.

.. _atommaps:

How AtomMap's work
===============================================================================

:class:`AtomMap` class adds great flexibility to manipulating atomic data.

First let's see how an instance of :class:`.Selection` (:class:`.Chain`, or
:class:`.Residue`) works.  Below table shows indices for a selection of atoms
in an :class:`~.AtomGroup` and values returned when
:meth:`~.Selection.getNames`, :meth:`~.Selection.getResnames` and
:meth:`~.Selection.getResnums` methods are called.

.. csv-table:: **Atom Subset**
   :header: "Indices", "Names", "Resnames", "Resnums"

   0, N, PHE, 1
   1, CA, PHE, 1
   2, C, PHE, 1
   3, O, PHE, 1
   4, CB, PHE, 1
   5, CG, PHE, 1
   6, CD1, PHE, 1
   7, CD2, PHE, 1
   8, CE1, PHE, 1
   9, CE2, PHE, 1
   10, CZ, PHE, 1

:class:`~.Selection` instances keep indices ordered and do not allow duplicate
values, hence their use is limited. In an :class:`AtomMap`, indices do not need
to be sorted, duplicate indices may exist, even "DUMMY" atoms are allowed.

Let's say we instantiate the following AtomMap::

    amap = AtomMap(atomgroup, indices=[0, 1, 3, 8, 8, 9, 10],
                   mapping=[5, 6, 7, 0, 1, 2, 3])


The size of the :class:`AtomMap` based on this mapping is 8, since the larger
mapping is 7.

Calling the same functions for this AtomMap instance would result in the
following:

.. csv-table:: **Atom Map**
   :header: "Mapping", "Indices", "Names", "Resnames", "Resnums", \
            "MappedFlags", "DummyFlags"

   0, 8, CE1, PHE, 1, 1, 0
   1, 8, CE1, PHE, 1, 1, 0
   2, 9, CE2, PHE, 1, 1, 0
   3, 10, CZ, PHE, 1, 1, 0
   4, , , , 0, 0, 1
   5, 0, N, PHE, 1, 1, 0
   6, 1, CA, PHE, 1, 1, 0
   7, 3, O, PHE, 1, 1, 0

For unmapped atoms, numeric attributes are set to 0, others to empty string,
i.e. ``""``.

.. seealso::
   :class:`AtomMap` are used by :mod:`.proteins` module functions that
   match or map protein chains.  :ref:`pca-xray` and :ref:`pca-dimer`
   examples that make use of these functions and :class:`AtomMap` class.


"""

try:
    from sys import maxint as DUMMY
except ImportError:
    from sys import maxsize as DUMMY

from numpy import arange, array, ndarray, ones, zeros, dtype

from prody.utilities import rangeString

from .atom import Atom
from .fields import ATOMIC_FIELDS
from .fields import wrapGetMethod, wrapSetMethod
from .pointer import AtomPointer

__all__ = ['AtomMap']


class AtomMap(AtomPointer):

    """A class for mapping atomic data."""

    __slots__ = ['_ag', '_indices', '_acsi', '_mapping', '_dummies', '_title',
                 '_len', '_idarray']

    def __init__(self, ag, indices, acsi=None, **kwargs):
        """Instantiate an atom map.

        :arg ag: AtomGroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg acsi: active coordinate set index, defaults is that of *ag*
        :arg mapping: mapping of atom *indices*
        :arg dummies: dummy atom indices
        :arg title: title of the instance, default is 'Unknown'

        *mapping* and *dummies* arrays must be provided together.  Length of
        *mapping* must be equal to length of *indices*.  Elements of *mapping*
        must be an ordered in ascending order.  When dummy atoms are present,
        number of atoms is the sum of lengths of *mapping* and *dummies*.

        Following built-in functions are customized for this class:

        * :func:`len` returns the number of atoms in the instance.
        * :func:`iter` yields :class:`.Atom` instances.
        * Indexing returns an :class:`.Atom` or an :class:`.AtomMap` instance
          depending on the type and value of the index."""

        AtomPointer.__init__(self, ag, acsi)

        self._dummies = self._mapping = None
        mapping = kwargs.get('mapping', None)
        dummies = kwargs.get('dummies', False)

        try:
            len(dummies)
        except TypeError:
            dummy_array = False
        else:
            dummy_array = True

        if mapping is None:
            if not kwargs.get('intarrays'):
                indices = array(indices, int)
            self._len = len(indices)

            if dummy_array:
                raise TypeError('mapping and dummies must be provided '
                                'together')

            if dummies:
                dummies = (indices == DUMMY).nonzero()[0]
                if len(dummies):
                    self._dummies = dummies
                    self._mapping = (indices < DUMMY).nonzero()[0]
                    self._indices = indices[self._mapping]
                    self._idarray = indices
                else:
                    self._indices = self._idarray = indices
            else:
                self._indices = self._idarray = indices
        else:
            if dummies is None:
                raise TypeError('mapping and dummies must be provided '
                                'together')

            if len(mapping) != len(indices):
                raise ValueError('indices and mapping arrays must have the '
                                 'same length')
            if not kwargs.get('intarrays'):
                indices = array(indices, int)
                mapping = array(mapping, int)
                if dummy_array:
                    dummies = array(dummies, int)
            if any(mapping[1:] - mapping[:-1] < 0):
                raise ValueError('mapping must be an ordered array')
            self._len = len(indices)
            if dummy_array:
                self._indices = indices
                self._mapping = mapping
                self._dummies = dummies
                self._len += len(dummies)
                self._idarray = idarray = zeros(self._len, int)
                idarray[self._mapping] = self._indices
                idarray[self._dummies] = DUMMY
            else:
                self._indices = self._idarray = indices[mapping]
        self._title = str(kwargs.get('title', 'Unknown'))

    def __repr__(self):

        rep = '<AtomMap: {0} from {1} ({2} atoms'.format(
                self._title, self._ag.getTitle(), self._len)
        if self.numDummies():
            rep += ', {0} mapped, {1} dummy'.format(self.numMapped(),
                                                        self.numDummies())

        n_csets = self._ag.numCoordsets()
        if n_csets > 1:
            rep += '; active #{0} of {1} coordsets)>'.format(
                    self.getACSIndex(), n_csets)
        elif n_csets == 0:
            rep += '; no coordinates'
        return rep + ')>'

    def __str__(self):

        return 'AtomMap {0}'.format(self._title)

    def __len__(self):

        return self._len

    def __getitem__(self, index):

        indices = self._idarray[index]
        try:
            n_sel = len(indices)
        except TypeError:
            if indices != DUMMY:
                return self._ag[indices]
        else:
            mapping = (indices > -1).nonzero()[0]
            return AtomMap(self._ag, indices, self._acsi,
                       title='({0})[{1}]'.format(self._title, repr(index)),
                       intarrays=True, dummies=self.numDummies())

    def getTitle(self):
        """Return title of the instance."""

        return self._title

    def setTitle(self, title):
        """Set title of the instance."""

        self._title = str(title)

    def numAtoms(self, flag=None):
        """Return number of atoms."""

        return len(self._getSubset(flag)) if flag else self._len

    def iterAtoms(self):
        """Yield atoms, and ``None`` for dummies."""

        ag = self._ag
        acsi = self.getACSIndex()
        for index in self.getIndices():
            yield Atom(ag, index, acsi) if index < DUMMY else None

    __iter__ = iterAtoms

    def getCoords(self):
        """Return a copy of coordinates from the active coordinate set."""

        coords = self._ag._getCoordsets()
        if coords is not None:
            if self._mapping is None:
                xyz = coords[self.getACSIndex(), self._indices]
            else:
                xyz = zeros((self._len, 3), float)
                xyz[self._mapping] = coords[self.getACSIndex(), self._indices]
            return xyz

    _getCoords = getCoords

    def setCoords(self, coords):
        """Set coordinates of atoms in the active coordinate set."""

        coordsets = self._ag._getCoordsets()
        if coordsets is not None:
            if self._mapping is None:
                coordsets[self.getACSIndex(), self._indices] = coords
            elif self._dummies is None:
                coordsets[self.getACSIndex(), self._indices] = coords


    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer
        or a list/array of integers."""

        coords = self._ag._getCoordsets()
        if coords is not None:
            n_csets = self._ag.numCoordsets()
            if indices is None:
                coords = coords[:, self._indices]
            else:
                try:
                    len(indices)
                except TypeError:
                    coords = coords[indices, self._indices]
                else:
                    coords = coords[indices][:, self._indices]

            if self._mapping is None:
                return coords
            else:
                csets = zeros(coords.shape[:-2] + (self._len, 3))
                csets[:, self._mapping] = coords
                return csets

    _getCoordsets = getCoordsets

    def iterCoordsets(self):
        """Yield copies of coordinate sets."""

        coords = self._ag._getCoordsets()
        if coords is not None:
            mapping = self._mapping
            n_atoms = self._len
            indices = self._indices
            for i in range(self._ag.numCoordsets()):
                xyz = zeros((n_atoms, 3), float)
                xyz[mapping] = coords[i, indices]
                yield xyz

    _iterCoordsets = iterCoordsets

    def getData(self, label):
        """Return a copy of data associated with *label*, if it is present."""

        data = self._ag._getData(label)
        if data is not None:
            result = zeros((self._len,) + data.shape[1:], data.dtype)
            result[self._mapping] = data[self._indices]
            return result

    _getData = getData

    def getFlags(self, label):
        """Return a copy of atom flags for given *label*, or **None** when
        flags for *label* is not set."""

        if label == 'dummy':
            flags = zeros(self._len, bool)
            if self._dummies is not None:
                flags[self._dummies] = True
        elif label == 'mapped':
            flags = ones(self._len, bool)
            if self._dummies is not None:
                flags[self._dummies] = False
        else:
            flags = None
            agflags = self._ag._getFlags(label)
            if agflags is not None:
                flags = zeros(self._len, bool)
                flags[self._mapping] = agflags[self._indices]
        return flags

    _getFlags = getFlags

    def _getSubset(self, label):

        return self._idarray[self._getFlags(label)]

    def getIndices(self):
        """Return a copy of indices of atoms, with maximum integer value
        dummies."""

        return self._idarray.copy()

    def _getIndices(self):
        """Return indices of atoms, with maximum integer value dummies."""

        return self._idarray

    def getMapping(self):
        """Return a copy of mapping of indices."""

        mapping = self._mapping
        return arange(self._len) if mapping is None else mapping.copy()

    def _getMapping(self):
        """Return mapping of indices."""

        mapping = self._mapping
        return arange(self._len) if mapping is None else mapping

    def numMapped(self):
        """Return number of mapped atoms."""

        return len(self._indices)

    def numDummies(self):
        """Return number of dummy atoms."""

        return 0 if self._dummies is None else len(self._dummies)

    def getSelstr(self):
        """Return selection string that selects mapped atoms."""

        return 'index ' + rangeString(self._indices)


for fname, field in ATOMIC_FIELDS.items():

    if field.private:
        continue

    meth = field.meth_pl
    getMeth = 'get' + meth

    def getData(self, meth=field.meth_pl, dtype=field.dtype):
        data = getattr(self._ag, '_get' + meth)()
        if data is not None:
            if self._mapping is None:
                return data[self._indices]
            else:
                result = zeros((self._len,) + data.shape[1:], dtype)
                result[self._mapping] = data[self._indices]
                return result

    getData = wrapGetMethod(getData)
    getData.__name__ = getMeth
    getData.__doc__ = (field.getDocstr('get', selex=False) +
                       ' Entries for dummy atoms will be ``{0}``.'
                       .format(repr(dtype(field.dtype).type())))
    setattr(AtomMap, getMeth, getData)
    setattr(AtomMap, '_' + getMeth, getData)

del getData
