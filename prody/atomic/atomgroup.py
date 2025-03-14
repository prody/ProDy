# -*- coding: utf-8 -*-
"""This module defines :class:`AtomGroup` class that stores atomic data and
multiple coordinate sets in :class:`~numpy.ndarray` instances."""

from time import time
from numbers import Integral

import numpy as np

from prody import LOGGER, PY2K
from prody.kdtree import KDTree
from prody.utilities import (checkCoords, checkAnisous, 
                             rangeString, getDistance, copy)

from .atomic import Atomic
from .fields import ATOMIC_FIELDS, READONLY
from .fields import wrapGetMethod, wrapSetMethod
from .flags import PLANTERS as FLAG_PLANTERS
from .flags import ALIASES as FLAG_ALIASES
from .flags import FIELDS as FLAG_FIELDS
from .atom import Atom
from .bond import Bond, evalBonds
from .angle import Angle, evalAngles
from .dihedral import Dihedral, evalDihedrals
from .crossterm import Crossterm, evalCrossterms
from .improper import Improper, evalImpropers
from .donor import Donor, evalDonors
from .acceptor import Acceptor, evalAcceptors
from .nbexclusion import NBExclusion, evalNBExclusions
from .selection import Selection

from . import flags

__all__ = ['AtomGroup']

if PY2K:
    range = xrange


def checkLabel(label):
    """Check suitability of *label* for labeling user data or flags."""

    label = str(label)
    if not label:
        raise ValueError('label cannot be empty string')

    if not label[0].isalpha():
        raise ValueError('label must start with a letter')

    if not (''.join(label.split('_'))).isalnum():
        raise ValueError('label may contain alphanumeric characters and '
                         'underscore, {0} is not valid'.format(label))

    if isReserved(label):
        raise ValueError('{0} is a reserved word and cannot be used '
                         'as a label'.format(repr(label)))

    if label in READONLY:
        raise AttributeError('{0} is read-only'.format(label))

    return label


class AtomGroup(Atomic):

    """A class for storing and accessing atomic data.  The number of atoms of
    the atom group is inferred at the first set method call from the size of
    the data array.

    **Atomic data**

    All atomic data is stored in :class:`~numpy.ndarray` instances.

    **Get and set methods**

    *get* methods, e.g. :meth:`getResnames`, return copies of the data arrays.

    *set* methods, e.g. :meth:`setResnums`, accept data in :func:`list` or
    :class:`~numpy.ndarray` instances.  The length of the list or array must
    match the number of atoms in the atom group.  These methods set attributes
    of all atoms at once.

    **Coordinate sets**

    Atom groups with multiple coordinate sets may have one of these sets as
    the *active coordinate set*.  The active coordinate set may be changed
    using :meth:`setACSIndex()` method.  :meth:`getCoords` returns coordinates
    from the *active set*.

    **Atom subsets**

    To access and modify data associated with a subset of atoms in an atom
    group, :class:`.Selection` instances may be used.  A :class:`.Selection`
    has initially the same coordinate set as the *active coordinate set*, but
    it may be changed using :meth:`.Selection.setACSIndex()` method.

    **Customizations**

    Following built-in functions are customized for this class:

    * :func:`len` returns the number of atoms, i.e. :meth:`numAtoms`
    * :func:`iter` yields :class:`.Atom` instances

    Indexing :class:`AtomGroup` instances by:
         - *int* (:func:`int`), e.g, ``10``, returns an :class:`.Atom`
         - *slice* (:func:`slice`), e.g, ``10:20:2``, returns a
           :class:`.Selection`
         - *segment name* (:func:`str`), e.g. ``'PROT'``, returns a
           a :class:`.Segment`
         - *chain identifier* (:func:`str`), e.g. ``'A'``, returns a
           a :class:`.Chain`
         - *[segment name,] chain identifier, residue number[, insertion code]*
           (:func:`tuple`), e.g. ``'A', 10`` or  ``'A', 10, 'B'`` or
           ``'PROT', 'A', 10, 'B'``, returns a :class:`.Residue`

    *Addition*

    Addition of two :class:`AtomGroup` instances, let's say *A* and *B*,
    results in a new :class:`AtomGroup` instance, say *C*.  *C* stores an
    independent copy of the data of *A* and *B*.  If *A* or *B* is missing
    a certain data type, zero values will be used for that part in *C*.
    If *A* and *B* has same number of coordinate sets, *C* will have a copy
    of all coordinate sets, otherwise *C* will have a single coordinate set,
    which is a copy of of active coordinate sets of *A* and *B*."""

    __slots__ = ['_title', '_n_atoms', '_coords', '_hv', '_sn2i',
                 '_timestamps', '_kdtrees',
                 '_bmap', '_angmap', '_dmap', '_imap',
                 '_domap', '_acmap', '_nbemap', '_cmap',
                 '_bonds', '_bondOrders', '_bondIndex', '_angles',
                 '_dihedrals', '_impropers',
                 '_donors', '_acceptors', '_nbexclusions', '_crossterms',
                 '_cslabels', '_acsi', '_n_csets', '_data',
                 '_fragments', '_flags', '_flagsts', '_subsets',
                 '_msa', '_sequenceMap', '_anisous']

    def __init__(self, title='Unnamed'):

        self._title = str(title)
        self._n_atoms = 0
        self._coords = None
        self._hv = None
        self._sn2i = None
        self._timestamps = None
        self._kdtrees = None
        self._bmap = None
        self._bonds = None
        self._bondOrders = None
        self._bondIndex = None
        self._angmap = None
        self._angles = None
        self._dmap = None
        self._dihedrals = None
        self._imap = None
        self._impropers = None
        self._domap = None
        self._donors = None
        self._acmap = None
        self._acceptors = None
        self._nbemap = None
        self._nbexclusions = None
        self._cmap = None
        self._crossterms = None
        self._fragments = None

        self._cslabels = []
        self._acsi = None
        self._n_csets = 0

        self._data = dict()

        self._flags = None
        self._flagsts = 0
        self._subsets = None
        self._msa = None
        self._sequenceMap = None
        self._anisous = None

    def __repr__(self):

        n_csets = self._n_csets
        if n_csets == 1:
            return '<AtomGroup: {0} ({1} atoms)>'.format(self._title,
                                                         self._n_atoms)
        elif n_csets > 1:
            return ('<AtomGroup: {0} ({1} atoms; active #{2} of {3}'
                    ' coordsets)>').format(self._title, self._n_atoms,
                                           self._acsi, n_csets)
        else:
            return ('<AtomGroup: {0} ({1} atoms; no coordinates)>'
                    ).format(self._title, self._n_atoms)

    def __str__(self):

        return 'AtomGroup ' + self._title

    def __getitem__(self, index):

        acsi = self._acsi

        if isinstance(index, slice):
            start, stop, step = index.indices(self._n_atoms)
            start = start or 0
            index = np.arange(start, stop, step)
            if len(index):
                if start > stop:
                    index = index[::-1]
                selstr = 'index {0}:{1}:{2}'.format(start, stop, step)
                return Selection(self, index, selstr, acsi, unique=True)

        elif isinstance(index, (list, np.ndarray)):
            unique = np.unique(index)
            if unique[0] < 0 or unique[-1] >= self._n_atoms:
                raise IndexError('index out of range')
            return Selection(self, unique, 'index ' + rangeString(index),
                             acsi, unique=True)

        elif isinstance(index, (str, tuple)):
            return self.getHierView()[index]

        else:
            try:
                index = int(index)
                n_atoms = self._n_atoms
                if index >= n_atoms or index < -n_atoms:
                    raise IndexError('index out of bounds')
                return Atom(self, index if index >= 0 else n_atoms + index, acsi)
            except:
                raise TypeError('invalid index')

    def __len__(self):

        return self._n_atoms

    def __add__(self, other, newAG=True):
        """This method adds two atom groups. By default this creates a new atom group.

        e.g.:
           newAG = ag1 + ag2
           len(newAG) == (len(ag1) + len(ag2))
        
        if newAG is set to false, ag1 is extended with the atoms of ag2
           oldLength = len(ag1)
           oldID = id(ag1)
           ag1.__add__(ag2, newAG=False)
           len(ag1) == (oldLength + len(ag2))
           id(ag1) === oldID
        """
        if not isinstance(other, AtomGroup):
            raise TypeError('unsupported operand type(s) for +: {0} and '
                            '{1}'.format(repr(type(self).__name__),
                                         repr(type(other).__name__)))

        oldSize = len(self)
        if newAG:
            new = AtomGroup(self._title + ' + ' + other._title)
        else:
            new = self
        self._n_atoms += other._n_atoms

        if self._n_csets:
            if self._n_csets == other._n_csets:
                new.setCoords(np.concatenate((self._coords, other._coords), 1), overwrite=True)
                this = self._anisous
                that = other._anisous
                if this is not None and that is not None:
                    if (isinstance(this, np.ndarray) and isinstance(that, np.ndarray)
                        and len(this) > 0 and len(that) > 0):
                        new.setAnisous(np.concatenate((self._anisous, other._anisous), 1))
                if self._n_csets > 1:
                    LOGGER.info('All {0} coordinate sets are copied to '
                                '{1}.'.format(self._n_csets, new.getTitle()))
            else:
                new.setCoords(np.concatenate((self._getCoords(),
                                              other._getCoords())))
                new.setAnisous(np.concatenate((self.getAnisous(),
                                               other.getAnisous())))
                LOGGER.info('Active coordinate sets are copied to {0}.'
                            .format(new.getTitle()))
        elif other._n_csets:
            LOGGER.warn('No coordinate sets are copied to {0}'
                        .format(new.getTitle()))

        for key in set(list(self._data) + list(other._data)):
            if key in ATOMIC_FIELDS and ATOMIC_FIELDS[key].readonly:
                continue
            this = self._data.get(key)
            that = other._data.get(key)
            if this is not None or that is not None:
                if this is None:
                    shape = list(that.shape)
                    shape[0] = oldSize
                    this = np.zeros(shape, that.dtype)
                if that is None:
                    shape = list(this.shape)
                    shape[0] = len(other)
                    that = np.zeros(shape, this.dtype)
                new._data[key] = np.concatenate((this, that))

        keys = []
        if self._flags:
            for flag in self._flags:
                if flag not in keys: keys.append(flag)

        if other._flags:
            for flag in other._flags:
                if flag not in keys: keys.append(flag)

        # remove aliases
        skip = []
        uniqueKeys = []
        for k in keys:
            aliases = FLAG_ALIASES.get(k, [k])
            if aliases[0] not in uniqueKeys and aliases[0] not in skip:
                uniqueKeys.append(aliases[0])
                if len(aliases) > 1:
                    skip.extend(list(aliases[1:]))
                
        for key in uniqueKeys:
            this = None
            that = None
            if self._flags:
                this = self._flags.get(key)
            if other._flags:
                that = other._flags.get(key)
            if this is not None or that is not None:
                if this is None:
                    shape = list(that.shape)
                    shape[0] = oldSize
                    this = np.zeros(shape, that.dtype)
                if that is None:
                    shape = list(this.shape)
                    shape[0] = len(other)
                    that = np.zeros(shape, this.dtype)
                new._setFlags(key, np.concatenate((this, that)))
                
        if self._bondOrders is not None:
            if other._bondOrders is not None:
                bo = np.concatenate([self._bondsOrders, other._bondOrders])
            else:
                bo = np.concatenate([self._bondsOrders, np.ones(len(other), np.int8)])
        else:
            if other._bondOrders is not None:
                bo = np.concatenate([np.ones(len(self), np.int8), self._bondsOrders])
            else:
                bo = None

        if self._bonds is not None and other._bonds is not None:
            new.setBonds(np.concatenate([self._bonds,
                                         other._bonds + oldSize]), bo)
        elif self._bonds is not None:
            new.setBonds(self._bonds.copy(), bo)

        elif other._bonds is not None:
            bonds = other._bonds + self._n_atoms
            new.setBonds(bonds, bo)

        if self._angles is not None and other._angles is not None:
            new.setAngles(np.concatenate([self._angles,
                                          other._angles + self._n_atoms]))
        elif self._angles is not None:
            new.setAngles(self._angles.copy())
        elif other._angles is not None:
            new.setAngles(other._angles + self._n_atoms)

        if self._dihedrals is not None and other._dihedrals is not None:
            new.setDihedrals(np.concatenate([self._dihedrals,
                                             other._dihedrals + self._n_atoms]))
        elif self._dihedrals is not None:
            new.setDihedrals(self._dihedrals.copy())
        elif other._dihedrals is not None:
            new.setDihedrals(other._dihedrals + self._n_atoms)

        if self._impropers is not None and other._impropers is not None:
            new.setImpropers(np.concatenate([self._impropers,
                                             other._impropers + self._n_atoms]))
        elif self._impropers is not None:
            new.setImpropers(self._impropers.copy())
        elif other._impropers is not None:
            new.setImpropers(other._impropers + self._n_atoms)

        if self._donors is not None and other._donors is not None:
            new.setDonors(np.concatenate([self._donors,
                                          other._donors + self._n_atoms]))
        elif self._donors is not None:
            new.setDonors(self._donors.copy())
        elif other._donors is not None:
            new.setDonors(other._donors + self._n_atoms)

        if self._acceptors is not None and other._acceptors is not None:
            new.setAcceptors(np.concatenate([self._acceptors,
                                             other._acceptors + self._n_atoms]))
        elif self._acceptors is not None:
            new.setAcceptors(self._acceptors.copy())
        elif other._acceptors is not None:
            new.setAcceptors(other._acceptors + self._n_atoms)

        if self._nbexclusions is not None and other._nbexclusions is not None:
            new.setNBExclusions(np.concatenate([self._nbexclusions,
                                                other._nbexclusions + self._n_atoms]))
        elif self._nbexclusions is not None:
            new.setNBExclusions(self._nbexclusions.copy())
        elif other._nbexclusions is not None:
            new.setNBExclusions(other._nbexclusions + self._n_atoms)

        if self._crossterms is not None and other._crossterms is not None:
            new.setCrossterms(np.concatenate([self._crossterms,
                                              other._crossterms + self._n_atoms]))
        elif self._crossterms is not None:
            new.setCrossterms(self._crossterms.copy())
        elif other._crossterms is not None:
            new.setCrossterms(other._crossterms + self._n_atoms)

        return new

    def __contains__(self, item):

        try:
            item.getACSIndex()
        except AttributeError:
            return False
        else:
            try:
                ag = item.getAtomGroup()
            except AttributeError:
                return item == self
            else:
                return ag == self

    def __eq__(self, other):

        return (isinstance(other, AtomGroup) and
                (self._n_atoms and self._n_atoms == other._n_atoms) and
                set(self._data) == set(other._data) and
                (self._n_csets and self._n_csets == other._n_csets and
                 np.all(self._coords == other._coords)) and
                all(np.all(self._data[key] == other._data[key])
                    for key in self._data))

    def __hash__(self):
        return object.__hash__(self)


    def __iter__(self):
        """Yield atom instances."""

        acsi = self._acsi
        for index in range(self._n_atoms):
            yield Atom(self, index, acsi)

    iterAtoms = __iter__

    def _none(self, attrs):
        """Set *attrs* **None** or remove them from data dictionary."""

        [self.__setattr__(nm,  None) if nm[0] == '_' else
         self._data.pop(nm,  None) for nm in attrs]

    def _getTimeStamp(self, index):
        """Returns time stamp showing when coordinates were last changed."""

        if self._n_cset:
            if index is None:
                return self._timestamps[self._acsi]
            else:
                return self._timestamps[index]
        else:
            return None

    def _setTimeStamp(self, index=None):
        """Set time stamp when :meth:`setCoords` methods of
        atom group or atom pointer instances are called.
        """

        if index is None:
            self._timestamps = np.zeros(self._n_csets)
            self._timestamps.fill(time())
            self._kdtrees = [None] * self._n_csets
        else:
            self._timestamps[index] = time()
            self._kdtrees[index] = None

    def _getKDTree(self, index=None):
        """Returns KDTree for coordinate set at given index."""

        if self._n_csets:
            if index is None:
                index = self._acsi
            kdtree = self._kdtrees[index]
            if kdtree is None:
                kdtree = KDTree(self._coords[index])
                self._kdtrees[index] = kdtree
            return kdtree
        else:
            return None

    def _getSN2I(self):
        """Returns a mapping of serial numbers to indices."""

        if self._sn2i is None and 'serial' in self._data:
            serials = self._data['serial']
            if serials is None:
                raise AttributeError('atom serial numbers are not set')
            unique = np.unique(serials)
            if len(unique) != self._n_atoms:
                raise ValueError('atom serial numbers are not unique')
            if unique[0] < 0:
                raise ValueError('atom serial numbers must be positive')
            sn2i = np.zeros(unique[-1] + 1, int)
            sn2i.fill(-1)
            sn2i[serials] = np.arange(self._n_atoms)
            self._sn2i = sn2i
        return self._sn2i

    def getTitle(self):
        """Returns title of the instance."""

        return self._title

    def setTitle(self, title):
        """Set title of the instance."""

        self._title = str(title)

    def numAtoms(self, flag=None):
        """Returns number of atoms, or number of atoms with given *flag*."""

        return len(self._getSubset(flag)) if flag else self._n_atoms

    def getCoords(self):
        """Returns a copy of coordinates from active coordinate set."""

        if self._coords is not None:
            return self._coords[self._acsi].copy()

    def _getCoords(self):
        """Returns a view of coordinates from active coordinate set."""

        if self._coords is not None:
            return self._coords[self._acsi]

    def setCoords(self, coords, label='', overwrite=False):
        """Set coordinates of atoms.  *coords* may be any array like object
        or an object instance with :meth:`getCoords` method.  If the shape of
        coordinate array is ``(n_csets > 1, n_atoms, 3)``, it will replace all
        coordinate sets and the active coordinate set index  will reset to
        zero.  This situation can be avoided using :meth:`addCoordset`.
        If shape of *coords* is ``(n_atoms, 3)`` or ``(1, n_atoms, 3)``, it
        will replace the active coordinate set.  *label* argument may be used
        to label coordinate set(s).  *label* may be a string or a list of
        strings length equal to the number of coordinate sets. The optional
        argument *overwrite* can be set to True to force resizing the
        coordinates array when the number of atoms in the AtomGroup changed."""

        atoms = coords
        try:
            if self._coords is None and hasattr(atoms, '_getCoords'):
                coords = atoms._getCoords()
            else:
                coords = atoms.getCoords()
        except AttributeError:
            if self._coords is None:
                coords = np.array(coords)
        else:
            if coords is None:
                raise ValueError('coordinates of {0} are not set'
                                 .format(str(atoms)))

        try:
            checkCoords(coords, csets=True, dtype=(float, np.float32))
        except TypeError:
            raise TypeError('coords must be a numpy array or an '
                            'object with `getCoords` method')

        self._setCoords(coords, label=label, overwrite=overwrite)

    def _setCoords(self, coords, label='', overwrite=False):
        """Set coordinates without data type checking.  *coords* must
        be a :class:`~numpy.ndarray`, but may have data type other than
        :class:`~numpy.float64`, e.g. :class:`~numpy.float32`.  *label*
        argument may be used to label coordinate sets.  *label* may be
        a string or a list of strings length equal to the number of
        coordinate sets. The optional argument *overwrite* can be set
        to True to force resizing the coordinates array when the number
        of atoms in the AtomGroup changed."""

        n_atoms = self._n_atoms
        if n_atoms:
            if coords.shape[-2] != n_atoms:
                raise ValueError('coords array has incorrect number of atoms')
        else:
            self._n_atoms = n_atoms = coords.shape[-2]

        ndim = coords.ndim
        shape = coords.shape
        if self._coords is None or overwrite or (ndim == 3 and shape[0] > 1):
            if ndim == 2:
                self._coords = coords.reshape((1, n_atoms, 3))
                self._cslabels = [str(label)]
                self._n_csets = n_csets = 1

            else:
                self._coords = coords
                self._n_csets = n_csets = shape[0]

                if isinstance(label, list):
                    if len(label) == n_csets:
                        self._cslabels = list(label)

                    else:
                        self._cslabels = [''] * n_csets
                        LOGGER.warn('Number of labels does not match number '
                                    'of coordinate sets.')
                else:
                    self._cslabels = [str(label)] * n_csets
            self._acsi = 0
            self._setTimeStamp()

        else:
            acsi = self._acsi
            if ndim == 2:
                self._coords[acsi] = coords
            else:
                self._coords[acsi] = coords[0]
            self._setTimeStamp(acsi)
            self._cslabels[acsi] = str(label)

    def getAnisous(self):
        """Returns a copy of anisotropic temperature factors from active coordinate set."""

        if self._anisous is not None:
            return self._anisous[self._acsi].copy()

    def _getAnisous(self):
        """Returns a view of anisotropic temperature factors from active coordinate set."""

        if self._anisous is not None:
            return self._anisous[self._acsi]

    def setAnisous(self, anisous, label=''):
        """Set anisotropic temperature factors of atoms. *anisous* may be any array like object
        or an object instance with :meth:`getAnisous` method.  If the shape of
        anisou array is ``(n_csets > 1, n_atoms, 3)``, it will replace all
        coordinate sets and the active coordinate set index  will reset to
        zero.  This situation can be avoided using :meth:`addCoordset`.
        If shape of *coords* is ``(n_atoms, 3)`` or ``(1, n_atoms, 3)``, it
        will replace the active coordinate set.  *label* argument may be used
        to label coordinate set(s).  *label* may be a string or a list of
        strings length equal to the number of coordinate sets."""

        atoms = anisous
        try:
            if self._anisous is None and hasattr(atoms, '_getAnisous'):
                anisous = atoms._getAnisous()
            else:
                anisous = atoms.getAnisous()
        except AttributeError:
            if self._anisous is None:
                anisous = np.array(anisous)
        else:
            if anisous is None:
                raise ValueError('anisous of {0} are not set'
                                 .format(str(atoms)))

        try:
            checkAnisous(anisous, csets=True, dtype=(float, np.float32))
        except TypeError:
            raise TypeError('anisous must be a numpy array or an '
                            'object with `getAnisous` method')

        self._setAnisous(anisous, label=label)

    def _setAnisous(self, anisous, label='', overwrite=False):
        """Set anisotropic temperature factors without data type checking.
        *anisous* must be a :class:`~numpy.ndarray`, but may have data type
        other than :class:`~numpy.float64`, e.g. :class:`~numpy.float32`.
        *label* argument may be used to label coordinate sets.  *label* may be
        a string or a list of strings length equal to the number of
        coordinate sets."""

        n_atoms = self._n_atoms
        if n_atoms:
            if anisous.shape[-2] != n_atoms:
                raise ValueError('anisous array has incorrect number of atoms')
        else:
            self._n_atoms = n_atoms = anisous.shape[-2]

        ndim = anisous.ndim
        shape = anisous.shape
        if self._anisous is None or overwrite or (ndim == 6 and shape[0] > 1):
            if ndim == 2:
                self._anisous = anisous.reshape((1, n_atoms, 6))
                self._cslabels = [str(label)]
                self._n_csets = n_csets = 1

            else:
                self._anisous = anisous
                self._n_csets = n_csets = shape[0]

                if isinstance(label, list):
                    if len(label) == n_csets:
                        self._cslabels = list(label)

                    else:
                        self._cslabels = [''] * n_csets
                        LOGGER.warn('Number of labels does not match number '
                                    'of coordinate sets.')
                else:
                    self._cslabels = [str(label)] * n_csets
            self._acsi = 0
            self._setTimeStamp()

        else:
            acsi = self._acsi
            if ndim == 2:
                self._anisous[acsi] = anisous
            else:
                self._anisous[acsi] = anisous[0]
            self._setTimeStamp(acsi)
            self._cslabels[acsi] = str(label)

    def addCoordset(self, coords, label=None, anisous=None):
        """Add a coordinate set.  *coords* argument may be an object with
        :meth:`getCoordsets` method."""

        if self._coords is None:
            return self.setCoords(coords)

        n_atoms = self._n_atoms
        atoms = coords
        try:
            coords = (atoms._getCoordsets()
                      if hasattr(coords, '_getCoordsets') else
                      atoms.getCoordsets())
        except AttributeError:
            pass
        else:
            if coords is None:
                raise ValueError('coordinates of {0} are not set'
                                 .format(str(atoms)))

        try:
            checkCoords(coords, csets=True, natoms=n_atoms, dtype=None)
        except TypeError:
            raise TypeError('coords must be a numpy array or an '
                            'object with `getCoords` method')

        if coords.ndim == 2:
            coords = coords.reshape((1, n_atoms, 3))

        if anisous is not None and anisous.ndim == 2:
            anisous = anisous.reshape((1, n_atoms, 6))

        diff = coords.shape[0]
        self._coords = np.concatenate((self._coords, coords), axis=0)
        if anisous is not None and self._anisous is not None:
            self._anisous = np.concatenate((self._anisous, anisous/10000), axis=0)
        self._n_csets = self._coords.shape[0]
        timestamps = self._timestamps
        self._timestamps = np.zeros(self._n_csets)
        self._timestamps[:len(timestamps)] = timestamps
        self._timestamps[len(timestamps):] = time()
        self._kdtrees.extend([None] * diff)
        if label is None or isinstance(label, str):
            self._cslabels.extend([label] * diff)
        elif isinstance(label, (list, tuple)):
            if len(label) == diff:
                self._cslabels.extend([str(lbl) for lbl in label])
            else:
                LOGGER.warn('Number of labels does not match number '
                            'of coordinate sets.')
        else:
            LOGGER.warn('Wrong type for `label` argument.')

    def delCoordset(self, index):
        """Delete a coordinate set from the atom group."""

        n_csets = self._n_csets
        if not n_csets:
            raise AttributeError('coordinates are not set')

        which = np.ones(n_csets, bool)
        which[index] = False
        which = which.nonzero()[0]
        if len(which) == 0:
            self._coords = None
            self._n_csets = 0
            self._acsi = None
            self._cslabels = None
            self._kdtrees = None
        else:
            self._coords = self._coords[which]
            self._n_csets = self._coords.shape[0]
            self._acsi = 0
            self._cslabels = [self._cslabels[i] for i in which]
            self._kdtrees = [self._kdtrees[i] for i in which]
        self._timestamps = self._timestamps[which]

    def getCoordsets(self, indices=None):
        """Returns a copy of coordinate set(s) at given *indices*.  *indices*
        may  be an integer, a list of integers, or **None** meaning all
        coordinate sets."""

        if self._coords is None:
            return None
        if indices is None:
            return self._coords.copy()
        if isinstance(indices, (Integral, slice)):
            return self._coords[indices].copy()

        # following fancy indexing makes a copy, so .copy() is not needed
        if isinstance(indices, (list, np.ndarray)):
            return self._coords[indices]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')

    def _getCoordsets(self, indices=None):
        """Returns a view of coordinate set(s) at given *indices*."""

        if self._coords is None:
            return None

        try:
            return self._coords if indices is None else self._coords[indices]
        except:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')

    def numBytes(self, all=False):
        """Returns number of bytes used by atomic data arrays, such as
        coordinate, flag, and attribute arrays.  If *all* is **True**,
        internal arrays for indexing hierarchical views, bonds, and
        fragments will also be included.  Note that memory usage of
        Python objects is not taken into account and that this may
        change in the future."""

        arrays = {}
        getbase = lambda arr: arr if arr.base is None else getbase(arr.base)
        getpair = lambda arr: (id(arr), arr)
        getboth = lambda arr: getpair(getbase(arr))

        if self._coords is not None:
            arrays[id(self._coords)] = self._coords
        arrays.update(getboth(val)
                      for key, val in self._data.items() if val is not None)
        if self._bonds is not None:
            arrays[id(self._bonds)] = self._bonds
        if self._angles is not None:
            arrays[id(self._angles)] = self._angles
        if self._dihedrals is not None:
            arrays[id(self._dihedrals)] = self._dihedrals
        if self._impropers is not None:
            arrays[id(self._impropers)] = self._impropers
        if self._flags:
            arrays.update(getboth(val)
                          for key, val in self._flags.items() if val is not None)
        if all:
            if self._subsets:
                arrays.update(getboth(val)
                              for key, val in self._subsets.items() if val is not None)
            if self._fragments:
                for val in self._fragments:
                    val = getbase(val)
                    arrays[id(val)] = val
            if self._bmap is not None:
                arrays[id(self._bonds)] = self._bmap
            if self._hv is not None:
                arrays.update(getboth(val) if hasattr(val, 'base') else
                              getboth(val._indices) for val in self._hv._residues)
                arrays.update(getboth(val) if hasattr(val, 'base') else
                              getboth(val._indices) for val in self._hv._chains)
                arrays.update(getboth(val) if hasattr(val, 'base') else
                              getboth(val._indices) for val in self._hv._segments)

        return sum(getbase(arr).nbytes for arr in arrays.values())

    def numCoordsets(self):
        """Returns number of coordinate sets."""

        return self._n_csets

    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate
        set."""

        for i in range(self._n_csets):
            yield self._coords[i].copy()

    def _iterCoordsets(self):
        """Iterate over coordinate sets by returning a view of each coordinate
        set."""

        for i in range(self._n_csets):
            yield self._coords[i]

    def getACSIndex(self):
        """Returns index of the coordinate set."""

        return self._acsi

    def setACSIndex(self, index):
        """Set the coordinate set at *index* active."""

        n_csets = self._n_csets
        if n_csets == 0:
            self._acsi = 0
        if not isinstance(index, Integral):
            raise TypeError('index must be an integer')
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += n_csets
        self._acsi = index

    def getHierView(self, **kwargs):
        """Returns a hierarchical view of the atom group."""

        if self._hv is None:
            self._hv = HierView(self, **kwargs)
        else:
            self._hv.update(**kwargs)

        return self._hv

    def numSegments(self):
        """Returns number of segments."""

        return self.getHierView().numSegments()

    def numChains(self):
        """Returns number of chains."""

        return self.getHierView().numChains()

    def numResidues(self):
        """Returns number of residues."""

        return self.getHierView().numResidues()

    def iterSegments(self):
        """Iterate over segments."""

        return self.getHierView().iterSegments()

    def iterChains(self):
        """Iterate over chains."""

        return self.getHierView().iterChains()

    def iterResidues(self):
        """Iterate over residues."""

        return self.getHierView().iterResidues()

    def setData(self, label, data):
        """Store atomic *data* under *label*, which must:

            * start with a letter
            * contain only alphanumeric characters and underscore
            * not be a reserved word (see :func:`.listReservedWords`)

        *data* must be a :func:`list` or a :class:`~numpy.ndarray` and its
        length must be equal to the number of atoms.  If the dimension of the
        *data* array is 1, i.e. ``data.ndim==1``, *label* may be used to make
        atom selections, e.g. ``"label 1 to 10"`` or ``"label C1 C2"``.  Note
        that, if data with *label* is present, it will be overwritten."""

        if label in ATOMIC_FIELDS:
            getattr(self, 'set' + ATOMIC_FIELDS[label].meth_pl)(data)
        else:
            label = checkLabel(label)

            if np.isscalar(data):
                data = [data] * self._n_atoms

            if label in self._data.keys():
                data = np.asarray(data, dtype=self._data[label].dtype)
            else:
                data = np.asarray(data)

            ndim, dtype, shape = data.ndim, data.dtype, data.shape

            if ndim == 1 and dtype == bool:
                raise TypeError('1 dimensional boolean arrays are not '
                                'accepted, use `setFlags` instead')

            if len(data) != self._n_atoms:
                raise ValueError('len(data) must match number of atoms')

            self._data[label] = data

    def delData(self, label):
        """Returns data associated with *label* and remove from the instance.
        If data associated with *label* is not found, return **None**."""

        return self._data.pop(label, None)

    def getData(self, label):
        """Returns a copy of the data array associated with *label*, or **None**
        if such data is not present."""

        data = self._getData(label)
        if data is not None:
            return data.copy()

    def _getData(self, label):
        """Returns data array associated with *label*, or **None** if such data
        is not present."""

        try:
            return self._data[label]
        except KeyError:
            try:
                field = ATOMIC_FIELDS[label]
            except KeyError:
                return None
            else:
                return getattr(self, '_get' + field.meth_pl)()

    def isDataLabel(self, label):
        """Returns **True** if data associated with *label* is present."""

        if label in self._data:
            return True
        else:
            try:
                return self._getData(label) is not None
            except:
                return False

    def getDataLabels(self, which=None):
        """Returns data labels.  For ``which='user'``, return only labels of
        user provided data."""

        if str(which).startswith('u'):  # user
            labels = [key for key in (self._data or {})
                      if not key in ATOMIC_FIELDS]
        else:
            labels = list(self._data or [])
        labels.sort()
        return labels

    def getDataType(self, label):
        """Returns type of the data (i.e. ``data.dtype``) associated with
        *label*, or **None** label is not used."""

        try:
            return self._data[label].dtype
        except KeyError:
            return None

    def isFlagLabel(self, label):
        """Returns **True** if flags associated with *label* are present."""

        return label in FLAG_PLANTERS or label in (self._flags or {})

    def getFlags(self, label):
        """Returns a copy of atom flags for given *label*, or **None** when
        flags for *label* is not set."""

        flags = self._getFlags(label)
        if flags is not None:
            return flags.copy()

    def _getFlags(self, label):
        """Returns atom flag values for given *label*, or **None** when
        flags for *label* is not set."""

        if self._flags is None:
            self._flags = {}
            self._subsets = {}
        elif flags.TIMESTAMP != self._flagsts:
            self._resetFlags()
        self._flagsts = flags.TIMESTAMP

        try:
            return self._flags[label]
        except KeyError:
            try:
                return FLAG_PLANTERS[label](self, label)
            except KeyError:
                pass

    def setFlags(self, label, flags):
        """Set atom *flags* for *label*."""

        label = checkLabel(label)
        if np.isscalar(flags):
            flags = np.array([flags] * self._n_atoms)

        try:
            ndim, dtype = flags.ndim, flags.dtype
        except AttributeError:
            flags = np.array(flags)
            ndim, dtype = flags.ndim, flags.dtype
        if ndim != 1:
            raise ValueError('flags.ndim must be 1')
        if dtype != bool:
            raise ValueError('flags.dtype must be bool')
        if len(flags) != self._n_atoms:
            raise ValueError('len(flags) must be equal to number of atoms')
        self._setFlags(label, flags)

    def _setFlags(self, label, flags):
        """Set atom flags."""

        if self._flags is None:
            self._flags = {}
            self._subsets = {}
        for label in FLAG_ALIASES.get(label, [label]):
            self._flags[label] = flags

    def delFlags(self, label):
        """Returns flags associated with *label* and remove from the instance.
        If flags associated with *label* is not found, return **None**."""

        return self._flags.pop(label, None)

    def _setSubset(self, label, indices):
        """Set indices of a subset of atoms."""

        for label in FLAG_ALIASES.get(label, [label]):
            self._subsets[label] = indices

    def _getSubset(self, label):
        """Returns indices of atoms."""

        if self._flags is None:
            self._flags = {}
            self._subsets = {}
        elif flags.TIMESTAMP != self._flagsts:
            self._resetFlags()
        self._flagsts = flags.TIMESTAMP

        try:
            return self._subsets[label]
        except KeyError:
            flgs = self._getFlags(label)
            try:
                return self._subsets[label]
            except KeyError:
                indices = flgs.nonzero()[0]
                self._setSubset(label, indices)
                return indices.copy()

    def getFlagLabels(self, which=None):
        """Returns flag labels.  For ``which='user'``,  return labels of user
        or parser (e.g. :term:`hetatm`) provided flags, for ``which='all'``
        return all possible :ref:`flags` labels in addition to those present
        in the instance."""

        which = str(which)
        if which.startswith('a'):  # all possible
            labels = set(self._flags or [])
            labels.update(FLAG_PLANTERS)
            labels = list(labels)
        elif which.startswith('u'):  # user
            labels = [key for key in (self._flags or {})
                      if not key in FLAG_PLANTERS]
        else:
            labels = list(self._flags or [])
        labels.sort()
        return labels

    def _resetFlags(self, field=None):
        """Reset flags and subsets associated with *field*."""

        flags = self._flags
        if flags is None:
            return
        if field:
            labels = FLAG_FIELDS[field]
        else:
            labels = list(FLAG_PLANTERS)
        subsets = self._subsets
        for label in labels:
            flags.pop(label, None)
            subsets.pop(label, None)

    def getBySerial(self, serial, stop=None, step=None):
        """Get an atom(s) by *serial* number (range).  *serial* must be zero or
        a positive integer. *stop* may be **None**, or an integer greater than
        *serial*.  ``getBySerial(i, j)`` will return atoms whose serial numbers
        are i+1, i+2, ..., j-1.  Atom whose serial number is *stop* will be
        excluded as it would be in indexing a Python :class:`list`.  *step*
        (default is 1) specifies increment.  If atoms with matching serial
        numbers are not found, **None** will be returned."""

        if not isinstance(serial, Integral):
            raise TypeError('serial must be an integer')
        if serial < 0:
            raise ValueError('serial must be greater than or equal to zero')
        sn2i = self._getSN2I()
        if sn2i is None:
            raise ValueError('serial numbers are not set')
        if stop is None:
            if serial < len(sn2i):
                index = sn2i[serial]
                if index != -1:
                    return Atom(self, index)
        else:
            if not isinstance(stop, Integral):
                raise TypeError('stop must be an integer')
            if stop <= serial:
                raise ValueError('stop must be greater than serial')

            if step is None:
                step = 1
            else:
                if not isinstance(step, Integral):
                    raise TypeError('step must be an integer')
                if step < 1:
                    raise ValueError('step must be greater than zero')

            indices = sn2i[serial:stop:step]
            indices = indices[indices > -1]
            return Selection(self, indices, 'serial {0}:{1}:{2}'
                                            .format(serial, stop, step))

    def getACSLabel(self):
        """Returns active coordinate set label."""

        if self._n_csets:
            return self._cslabels[self._acsi]

    def setACSLabel(self, label):
        """Set active coordinate set label."""

        if self._n_csets:
            if label is None or isinstance(label, str):
                self._cslabels[self._acsi] = label
            else:
                raise TypeError('label must be a string')

    def getCSLabels(self):
        """Returns coordinate set labels."""

        if self._n_csets:
            return list(self._cslabels)

    def setCSLabels(self, labels):
        """Set coordinate set labels. *labels* must be a list of strings."""

        if isinstance(labels, list):
            if len(labels) == self._n_csets:
                if all((lbl is None or isinstance(lbl, str))
                       for lbl in labels):
                    self._cslabels = list(labels)
                else:
                    raise ValueError('all items of labels must be strings')
            else:
                raise ValueError('length of labels must be equal to the '
                                 'number of coordinate sets')
        else:
            raise TypeError('labels must be a list')

    def setBondOrders(self, bondOrders):
        """Set covalent bond order. *bondOrders* must be a list or an array
        of integers and provide a value for each bond. Possible values are
        1:single, 2:double, 3:triple, 4:aromatic, 5:amide. The bond order of
        all bonds must be set at once. This method must be called after
        the setBonds() has been called. The bond order is stored in the
        *_bondOrders* array."""

        if bondOrders is None:
          self._bondOrders = bondOrders
          return
      
        if len(bondOrders)!=len(self._bonds):
            raise ValueError('invalid bond order list, bond and bond order length mismatch')
        if min(bondOrders)<1 or max(bondOrders)>5:
            raise ValueError('invalid bond order value, values must range from 1 to 5')

        self._bondOrders = np.array(bondOrders, np.int8)

    def setBonds(self, bonds, bondOrders=None):
        """Set covalent bonds between atoms.  *bonds* must be a list or an
        array of pairs of indices. Optionally *bondOrders* can be specified
        (see :meth:`setBondOrder` method).  All bonds must be
        set at once.  Bonding information can be used to make atom selections,
        e.g. ``"bonded to index 1"``.  See :mod:`.select` module documentation
        for details. Also, a data array with number of bonds will be generated
        and stored with label *numbonds*.  This can be used in atom selections,
        e.g. ``'numbonds 0'`` can be used to select ions in a system. The keys
        in the *_bondIndex* dictionary is a string representation of the bond's
        atom indices i.e. '%d %d'%(i, j) with i<j. If  *bonds*  is empty or
        **None**, then all bonds will be removed for this
        :class:`.AtomGroup`. """

        if bonds is None or len(bonds) == 0:
            self._bmap = None
            self._bonds = None
            self._fragments = None
            return

        if isinstance(bonds, list):
            bonds = np.array(bonds, int)
        if bonds.ndim != 2:
            raise ValueError('bonds.ndim must be 2')
        if bonds.shape[1] != 2:
            raise ValueError('bonds.shape must be (n_bonds, 2)')
        if bonds.min() < 0:
            raise ValueError('negative atom indices are not valid')

        n_atoms = self._n_atoms
        if bonds.max() >= n_atoms:
            raise ValueError('atom indices are out of range')

        bonds.sort(1)
        bonds = bonds[bonds[:, 1].argsort(), ]
        bonds = bonds[bonds[:, 0].argsort(), ]
        bonds = np.unique(bonds, axis=0)

        d = {}
        for n, b in enumerate(bonds):
            key = '%d %d'%(b[0], b[1])
            d[key] = n
        self._bondIndex = d

        self._bmap, self._data['numbonds'] = evalBonds(bonds, n_atoms)
        self._bonds = bonds
        self._fragments = None

        self.setBondOrders(bondOrders)
            
    def numBonds(self):
        """Returns number of bonds.  Use :meth:`setBonds` or 
        :meth:`inferBonds` for setting bonds."""

        if self._bonds is not None:
            return self._bonds.shape[0]
        return 0

    def getBonds(self):
        """Returns bonds.  Use :meth:`setBonds` or 
        :meth:`inferBonds` for setting bonds."""

        if self._bonds is not None:
            acsi = self._acsi
            return np.array([Bond(self, bond, acsi) for bond in self._bonds])
        return None

    def inferBonds(self, max_bond=1.6, min_bond=0, set_bonds=True):
        """Returns bonds based on distances **max_bond** and **min_bond**."""

        bonds = []
        for atom_i in self.iterAtoms():
            sele = self.select('index > {0} and exwithin {1} of index {0}'
                               .format(atom_i.getIndex(), max_bond))
            if sele is not None:
                for atom_j in sele.iterAtoms():
                    distance = getDistance(atom_i.getCoords(),
                                           atom_j.getCoords())
                    if distance > min_bond:
                        bonds.append([atom_i._index, atom_j._index])

        if set_bonds:
            self.setBonds(bonds)

        acsi = self._acsi
        return np.array([Bond(self, bond, acsi) for bond in bonds])

    def iterBonds(self):
        """Yield bonds.  Use :meth:`setBonds` or `inferBonds` for setting bonds."""

        if self._bonds is not None:
            acsi = self._acsi
            for bond in self._bonds:
                yield Bond(self, bond, acsi)

    def _iterBonds(self):
        """Yield pairs of bonded atom indices. Use :meth:`setBonds` 
        or `inferBonds` for setting bonds."""

        if self._bonds is not None:
            for a, b in self._bonds:
                yield a, b

    def setAngles(self, angles):
        """Set covalent angles between atoms.  *angles* must be a list or an
        array of triplets of indices.  All angles must be set at once.  Angle
        information can be used to make atom selections, e.g. ``"angle to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of angles will be generated and stored
        with label *numangles*.  This can be used in atom selections, e.g.
        ``'numangles 0'`` can be used to select ions in a system."""

        if isinstance(angles, list):
            angles = np.array(angles, int)
        if angles.ndim != 2:
            raise ValueError('angles.ndim must be 2')
        if angles.shape[1] != 3:
            raise ValueError('angles.shape must be (n_angles, 3)')
        if angles.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if angles.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        angles.sort(1)
        angles = angles[angles[:, 2].argsort(), ]
        angles = angles[angles[:, 1].argsort(), ]
        angles = angles[angles[:, 0].argsort(), ]

        self._angmap, self._data['numangles'] = evalAngles(angles, n_atoms)
        self._angles = angles

    def numAngles(self):
        """Returns number of angles.  Use :meth:`setAngles` for setting angles."""

        if self._angles is not None:
            return self._angles.shape[0]
        return 0

    def getAngles(self):
        """Returns angles.  Use :meth:`setAngles` for setting angles."""

        if self._angles is not None:
            acsi = self._acsi
            return np.array([Angle(self, angle, acsi) for angle in self._angles])
        return None

    def iterAngles(self):
        """Yield angles.  Use :meth:`setAngles` for setting angles."""

        if self._angles is not None:
            acsi = self._acsi
            for angle in self._angles:
                yield Angle(self, angle, acsi)

    def _iterAngles(self):
        """Yield triplets of angled atom indices. Use :meth:`setAngles` for setting
        angles."""

        if self._angles is not None:
            for a, b, c in self._angles:
                yield a, b, c

    def setDihedrals(self, dihedrals):
        """Set covalent dihedrals between atoms.  *dihedrals* must be a list or an
        array of triplets of indices.  All dihedrals must be set at once.  Dihedral
        information can be used to make atom selections, e.g. ``"dihedral to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of dihedrals will be generated and stored
        with label *numdihedrals*.  This can be used in atom selections, e.g.
        ``'numdihedrals 0'`` can be used to select ions in a system."""

        if isinstance(dihedrals, list):
            dihedrals = np.array(dihedrals, int)
        if dihedrals.ndim != 2:
            raise ValueError('dihedrals.ndim must be 2')
        if dihedrals.shape[1] != 4:
            raise ValueError('dihedrals.shape must be (n_dihedrals, 4)')
        if dihedrals.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if dihedrals.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        dihedrals.sort(1)
        dihedrals = dihedrals[dihedrals[:, 3].argsort(), ]
        dihedrals = dihedrals[dihedrals[:, 2].argsort(), ]
        dihedrals = dihedrals[dihedrals[:, 1].argsort(), ]
        dihedrals = dihedrals[dihedrals[:, 0].argsort(), ]

        self._dmap, self._data['numdihedrals'] = evalDihedrals(
            dihedrals, n_atoms)
        self._dihedrals = dihedrals

    def numDihedrals(self):
        """Returns number of dihedrals.  Use :meth:`setDihedrals` for setting dihedrals."""

        if self._dihedrals is not None:
            return self._dihedrals.shape[0]
        return 0

    def getDihedrals(self):
        """Returns dihedrals.  Use :meth:`setDihedrals` for setting dihedrals."""

        if self._dihedrals is not None:
            acsi = self._acsi
            return np.array([Dihedral(self, dihedral, acsi) for dihedral in self._dihedrals])
        return None

    def iterDihedrals(self):
        """Yield dihedrals.  Use :meth:`setDihedrals` for setting dihedrals."""

        if self._dihedrals is not None:
            acsi = self._acsi
            for dihedral in self._dihedrals:
                yield Dihedral(self, dihedral, acsi)

    def _iterDihedrals(self):
        """Yield quadruplets of dihedraled atom indices. Use :meth:`setDihedrals` for setting
        dihedrals."""

        if self._dihedrals is not None:
            for a, b, c, d in self._dihedrals:
                yield a, b, c, d

    def setImpropers(self, impropers):
        """Set covalent impropers between atoms.  *impropers* must be a list or an
        array of triplets of indices.  All impropers must be set at once.  Improper
        information can be used to make atom selections, e.g. ``"improper to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of impropers will be generated and stored
        with label *numimpropers*.  This can be used in atom selections, e.g.
        ``'numimpropers 0'`` can be used to select ions in a system."""

        if isinstance(impropers, list):
            impropers = np.array(impropers, int)
        if impropers.ndim != 2:
            raise ValueError('impropers.ndim must be 2')
        if impropers.shape[1] != 4:
            raise ValueError('impropers.shape must be (n_impropers, 4)')
        if impropers.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if impropers.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        impropers.sort(1)
        impropers = impropers[impropers[:, 3].argsort(), ]
        impropers = impropers[impropers[:, 2].argsort(), ]
        impropers = impropers[impropers[:, 1].argsort(), ]
        impropers = impropers[impropers[:, 0].argsort(), ]

        self._imap, self._data['numimpropers'] = evalImpropers(
            impropers, n_atoms)
        self._impropers = impropers

    def numImpropers(self):
        """Returns number of impropers.  Use :meth:`setImpropers` for setting impropers."""

        if self._impropers is not None:
            return self._impropers.shape[0]
        return 0

    def getImpropers(self):
        """Returns impropers.  Use :meth:`setImpropers` for setting impropers."""

        if self._impropers is not None:
            acsi = self._acsi
            return np.array([Improper(self, improper, acsi) for improper in self._impropers])
        return None

    def iterImpropers(self):
        """Yield impropers.  Use :meth:`setImpropers` for setting impropers."""

        if self._impropers is not None:
            acsi = self._acsi
            for improper in self._impropers:
                yield Improper(self, improper, acsi)

    def _iterImpropers(self):
        """Yield quadruplets of impropered atom indices. Use :meth:`setImpropers` for setting
        impropers."""

        if self._impropers is not None:
            for a, b, c, d in self._impropers:
                yield a, b, c, d

    def setDonors(self, donors):
        """Set covalent donors between atoms.  *donors* must be a list or an
        array of pairs of indices.  All donors must be set at once.  Donoring
        information can be used to make atom selections, e.g. ``"donored to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of donors will be generated and stored
        with label *numdonors*.  This can be used in atom selections, e.g.
        ``'numdonors 0'`` can be used to select ions in a system. If *donors* 
        is empty or **None**, then all donors will be removed for this 
        :class:`.AtomGroup`. """

        if donors is None or len(donors) == 0:
            self._domap = None
            self._donors = None
            return

        if isinstance(donors, list):
            donors = np.array(donors, int)
        if donors.ndim != 2:
            raise ValueError('donors.ndim must be 2')
        if donors.shape[1] != 2:
            raise ValueError('donors.shape must be (n_donors, 2)')
        if donors.min() < -1:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if donors.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        donors.sort(1)
        donors = donors[donors[:, 1].argsort(), ]
        donors = donors[donors[:, 0].argsort(), ]
        donors = np.unique(donors, axis=0)

        self._domap, self._data['numdonors'] = evalDonors(donors, n_atoms)
        self._donors = donors

    def numDonors(self):
        """Returns number of donors.  Use :meth:`setDonors` for setting donors."""

        if self._donors is not None:
            return self._donors.shape[0]
        return 0

    def getDonors(self):
        """Returns donors.  Use :meth:`setDonors` for setting donors."""

        if self._donors is not None:
            acsi = self._acsi
            return np.array([Donor(self, donor, acsi) for donor in self._donors])
        return None

    def iterDonors(self):
        """Yield donors.  Use :meth:`setDonors` for setting donors."""

        if self._donors is not None:
            acsi = self._acsi
            for donor in self._donors:
                yield Donor(self, donor, acsi)

    def _iterDonors(self):
        """Yield pairs of donored atom indices. Use :meth:`setDonors` for setting
        donors."""

        if self._donors is not None:
            for a, b in self._donors:
                yield a, b

    def setAcceptors(self, acceptors):
        """Set covalent acceptors between atoms.  *acceptors* must be a list or an
        array of pairs of indices.  All acceptors must be set at once.  Acceptoring
        information can be used to make atom selections, e.g. ``"acceptored to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of acceptors will be generated and stored
        with label *numacceptors*.  This can be used in atom selections, e.g.
        ``'numacceptors 0'`` can be used to select ions in a system. If *acceptors* 
        is empty or **None**, then all acceptors will be removed for this 
        :class:`.AtomGroup`. """

        if acceptors is None or len(acceptors) == 0:
            self._acmap = None
            self._acceptors = None
            return

        if isinstance(acceptors, list):
            acceptors = np.array(acceptors, int)
        if acceptors.ndim != 2:
            raise ValueError('acceptors.ndim must be 2')
        if acceptors.shape[1] != 2:
            raise ValueError('acceptors.shape must be (n_acceptors, 2)')
        if acceptors.min() < -1:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if acceptors.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        acceptors.sort(1)
        acceptors = acceptors[acceptors[:, 1].argsort(), ]
        acceptors = acceptors[acceptors[:, 0].argsort(), ]
        acceptors = np.unique(acceptors, axis=0)

        self._acmap, self._data['numacceptors'] = evalAcceptors(acceptors, n_atoms)
        self._acceptors = acceptors

    def numAcceptors(self):
        """Returns number of acceptors.  Use :meth:`setAcceptors` for setting acceptors."""

        if self._acceptors is not None:
            return self._acceptors.shape[0]
        return 0

    def getAcceptors(self):
        """Returns acceptors.  Use :meth:`setAcceptors` for setting acceptors."""

        if self._acceptors is not None:
            acsi = self._acsi
            return np.array([Acceptor(self, acceptor, acsi) for acceptor in self._acceptors])
        return None

    def iterAcceptors(self):
        """Yield acceptors.  Use :meth:`setAcceptors` for setting acceptors."""

        if self._acceptors is not None:
            acsi = self._acsi
            for acceptor in self._acceptors:
                yield Acceptor(self, acceptor, acsi)

    def _iterAcceptors(self):
        """Yield pairs of acceptored atom indices. Use :meth:`setAcceptors` for setting
        acceptors."""

        if self._acceptors is not None:
            for a, b in self._acceptors:
                yield a, b

    def setNBExclusions(self, nbexclusions):
        """Set nbexclusions between atoms.  *nbexclusions* must be a list or an
        array of pairs of indices.  All nbexclusions must be set at once.  Acceptoring
        information can be used to make atom selections, e.g. ``"nbexclusioned to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of nbexclusions will be generated and stored
        with label *numnbexclusions*.  This can be used in atom selections, e.g.
        ``'numnbexclusions 0'`` can be used to select ions in a system. If *nbexclusions* 
        is empty or **None**, then all nbexclusions will be removed for this 
        :class:`.AtomGroup`. """

        if nbexclusions is None or len(nbexclusions) == 0:
            self._nbemap = None
            self._nbexclusions = None
            return

        if isinstance(nbexclusions, list):
            nbexclusions = np.array(nbexclusions, int)
        if nbexclusions.ndim != 2:
            raise ValueError('nbexclusions.ndim must be 2')
        if nbexclusions.shape[1] != 2:
            raise ValueError('nbexclusions.shape must be (n_nbexclusions, 2)')
        if nbexclusions.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if nbexclusions.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        nbexclusions.sort(1)
        nbexclusions = nbexclusions[nbexclusions[:, 1].argsort(), ]
        nbexclusions = nbexclusions[nbexclusions[:, 0].argsort(), ]
        nbexclusions = np.unique(nbexclusions, axis=0)

        self._nbemap, self._data['numnbexclusions'] = evalNBExclusions(
            nbexclusions, n_atoms)
        self._nbexclusions = nbexclusions

    def numNBExclusions(self):
        """Returns number of nbexclusions.  Use :meth:`setNBExclusions` for setting nbexclusions."""

        if self._nbexclusions is not None:
            return self._nbexclusions.shape[0]
        return 0

    def getNBExclusions(self):
        """Returns nbexclusions.  Use :meth:`setNBExclusions` for setting nbexclusions."""

        if self._nbexclusions is not None:
            acsi = self._acsi
            return np.array([Acceptor(self, nbexclusion, acsi) for nbexclusion in self._nbexclusions])
        return None

    def iterNBExclusions(self):
        """Yield nbexclusions.  Use :meth:`setNBExclusions` for setting nbexclusions."""

        if self._nbexclusions is not None:
            acsi = self._acsi
            for nbexclusion in self._nbexclusions:
                yield Acceptor(self, nbexclusion, acsi)

    def _iterNBExclusions(self):
        """Yield pairs of nbexclusioned atom indices. Use :meth:`setNBExclusions` for setting
        nbexclusions."""

        if self._nbexclusions is not None:
            for a, b in self._nbexclusions:
                yield a, b

    def setCrossterms(self, crossterms):
        """Set covalent crossterms between atoms.  *crossterms* must be a list or an
        array of triplets of indices.  All crossterms must be set at once.  Crossterm
        information can be used to make atom selections, e.g. ``"crossterm to
        index 1"``.  See :mod:`.select` module documentation for details.
        Also, a data array with number of crossterms will be generated and stored
        with label *numcrossterms*.  This can be used in atom selections, e.g.
        ``'numcrossterms 0'`` can be used to select ions in a system."""

        if isinstance(crossterms, list):
            crossterms = np.array(crossterms, int)
        if crossterms.ndim != 2:
            raise ValueError('crossterms.ndim must be 2')
        if crossterms.shape[1] != 4:
            raise ValueError('crossterms.shape must be (n_crossterms, 4)')
        if crossterms.min() < 0:
            raise ValueError('negative atom indices are not valid')
        n_atoms = self._n_atoms
        if crossterms.max() >= n_atoms:
            raise ValueError('atom indices are out of range')
        crossterms.sort(1)
        crossterms = crossterms[crossterms[:, 3].argsort(), ]
        crossterms = crossterms[crossterms[:, 2].argsort(), ]
        crossterms = crossterms[crossterms[:, 1].argsort(), ]
        crossterms = crossterms[crossterms[:, 0].argsort(), ]

        self._cmap, self._data['numcrossterms'] = evalCrossterms(
            crossterms, n_atoms)
        self._crossterms = crossterms

    def numCrossterms(self):
        """Returns number of crossterms.  Use :meth:`setCrossterms` for setting crossterms."""

        if self._crossterms is not None:
            return self._crossterms.shape[0]
        return 0

    def getCrossterms(self):
        """Returns crossterms.  Use :meth:`setCrossterms` for setting crossterms."""

        if self._crossterms is not None:
            acsi = self._acsi
            return np.array([Crossterm(self, crossterm, acsi) for crossterm in self._crossterms])
        return None

    def iterCrossterms(self):
        """Yield crossterms.  Use :meth:`setCrossterms` for setting crossterms."""

        if self._crossterms is not None:
            acsi = self._acsi
            for crossterm in self._crossterms:
                yield Crossterm(self, crossterm, acsi)

    def _iterCrossterms(self):
        """Yield quadruplets of crosstermed atom indices. Use :meth:`setCrossterms` for setting
        crossterms."""

        if self._crossterms is not None:
            for a, b, c, d in self._crossterms:
                yield a, b, c, d

    def numFragments(self):
        """Returns number of connected atom subsets."""

        self._fragment()
        if self._fragments is None:
            return 0
        return self._data['fragindex'].max() + 1

    def iterFragments(self):
        """Yield connected atom subsets as :class:`.Selection` instances."""

        if self._bmap is not None:

            acsi = self._acsi
            if self._fragments is None:
                self._fragment()
                if self._fragments is None:
                    return

            for i, frag in enumerate(self._fragments):
                try:
                    frag.getAtomGroup()
                except AttributeError:
                    frag = Selection(self, frag, 'fragment ' + str(i),
                                     acsi=acsi, unique=True)
                finally:
                    self._fragments[i] = frag
                yield frag

    def _fragment(self):
        """Set unique fragment indices to connected atom subsets using bond
        information."""

        if self._bmap is None:
            LOGGER.warn('bonds must be set for fragment determination, '
                        'use `setBonds` or `inferBonds` to set them')
            self._data['fragindex'] = None
            self._fragments = None
            return

        fids = np.zeros(self._n_atoms, int)
        fdict = {}
        c = 0
        for a, b in self._bonds:
            af = fids[a]
            bf = fids[b]
            if af and bf:
                if af != bf:
                    frag = fdict[af]
                    temp = fdict[bf]
                    fids[temp] = af
                    frag.extend(temp)
                    fdict.pop(bf)
            elif af:
                fdict[af].append(b)
                fids[b] = af
            elif bf:
                fdict[bf].append(a)
                fids[a] = bf
            else:
                c += 1
                fdict[c] = [a, b]
                fids[a] = fids[b] = c
        fragindices = np.zeros(self._n_atoms, int)
        fragments = []
        append = fragments.append
        fidset = set()
        c = 0
        for i, fid in enumerate(fids):
            if fid in fidset:
                continue
            elif fid:
                fidset.add(fid)
                indices = fdict[fid]
                indices.sort()
                append(indices)
                fragindices[indices] = c
                c += 1
            else:
                # these are non-bonded atoms, e.g. ions
                fragindices[i] = c
                append([i])
                c += 1
        self._data['fragindex'] = fragindices
        self._fragments = fragments


for fname, field in ATOMIC_FIELDS.items():

    meth = field.meth_pl
    getMeth = 'get' + meth
    setMeth = 'set' + meth
    if field.call:
        # Define public method for retrieving a copy of data array
        if not field.private:
            def getData(self, var=fname, call=field.call):
                try:
                    return copy(self._data[var])
                except KeyError:
                    [getattr(self, meth)() for meth in call]
                    return copy(self._data[var])

        # Define private method for retrieving actual data array
        def _getData(self, var=fname, call=field.call):
            try:
                return self._data[var]
            except KeyError:
                [getattr(self, meth)() for meth in call]
                return copy(self._data[var])
    else:
        if not field.private:
            def getData(self, var=fname):
                try:
                    return copy(self._data[var])
                except KeyError:
                    pass

        def _getData(self, var=fname):
            return self._data.get(var)

    if not field.private:
        getData = wrapGetMethod(getData)
        getData.__name__ = getMeth
        getData.__doc__ = field.getDocstr('get')
        setattr(AtomGroup, getMeth, getData)

    _getData = wrapGetMethod(_getData)
    _getData.__name__ = '_' + getMeth
    _getData.__doc__ = field.getDocstr('_get')
    setattr(AtomGroup, '_' + getMeth, _getData)

    if field.readonly or field.private:
        continue

    # Define public method for setting values in data array
    def setData(self, array, var=fname, dtype=field.dtype,
                ndim=field.ndim, none=field.none, flags=field.flags):
        if array is None:
            self._data.pop(var, None)
        else:
            if np.isscalar(array):
                self._data[var][:] = array
            else:
                if self._n_atoms == 0:
                    self._n_atoms = len(array)
                elif len(array) != self._n_atoms:
                    raise ValueError('length of array must match number '
                                     'of atoms')

                if not np.isscalar(array):
                    if var == 'chain':
                        max_len = 0
                        for val in array:
                            if len(val) > max_len:
                                max_len = len(val)

                        if max_len > int(dtype[1:]):
                            dtype = dtype[0] + str(max_len)

                    array = np.asarray(array, dtype)
                else:
                    raise TypeError('array must be an ndarray or a list')

                if array.ndim != ndim:
                    raise ValueError('array must be {0} '
                                     'dimensional'.format(ndim))
                elif array.dtype != dtype:
                    try:
                        array = array.astype(dtype)
                    except ValueError:
                        raise ValueError('array cannot be assigned type '
                                         '{0}'.format(dtype))
                self._data[var] = array
                if none:
                    self._none(none)
                if flags and self._flags:
                    self._resetFlags(var)

    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth
    setData.__doc__ = field.getDocstr('set')
    setattr(AtomGroup, setMeth, setData)

del getData
del _getData
del setData
