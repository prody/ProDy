# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines :class:`AtomGroup` class that stores atomic data and 
multiple coordinate sets in :class:`numpy.ndarray` instances.

Below example shows how to build an :class:`AtomGroup` from scratch.  We start 
by importing everything from the ProDy package and the NumPy package:

>>> from prody import *
>>> import numpy as np

.. _atomgroup:

Atom Group
===============================================================================

>>> wtr1 = AtomGroup('Water')
>>> wtr1
<AtomGroup: Water (0 atoms; no coordinates)>


Coordinates
-------------------------------------------------------------------------------

The best way to start constructing an :class:`AtomGroup` is by setting the
coordinates first. Number of atoms will be automatically set according to
the size of the coordinate data array:

>>> coords = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=float)
>>> print( coords )
[[ 1.  0.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  1.]]
>>> wtr1.setCoords( coords )
>>> wtr1
<AtomGroup: Water (3 atoms)>

Attributes
-------------------------------------------------------------------------------

Attributes must be passed in a list or an array whose size is the same
as the number of atoms.

>>> wtr1.setNames( ['H', 'O', 'H'] )
>>> wtr1.setResnums( [1, 1, 1] )
>>> wtr1.setResnames( ['WAT', 'WAT', 'WAT'] )

Accessing data will return a copy of the data:

>>> print( wtr1.getNames() )
['H' 'O' 'H']

Atoms
===============================================================================

Atoms are represented by instance of :class:`.Atom`.

Iteration
-------------------------------------------------------------------------------

Atoms in an :class:`AtomGroup` can be iterated over

>>> for a in wtr1: a
... 
<Atom: H from Water (index 0)>
<Atom: O from Water (index 1)>
<Atom: H from Water (index 2)>

Indexing
-------------------------------------------------------------------------------

Atoms in an atom group can be accessed via indexing:

>>> a = wtr1[0]
>>> a
<Atom: H from Water (index 0)>
>>> print( a.getCoords() )
[ 1.  0.  0.]

Coordinate sets
===============================================================================

Let's add another coordinate set to the atom group:

>>> wtr1.addCoordset(np.array([[0, 1, 0], [0, 0, 0], [0, 0, 1]], dtype=float))
>>> wtr1
<AtomGroup: Water (3 atoms; active #0 of 2 coordsets)>

Note that number of coordinate sets is now 2, but active coordinate set index
is still 0. Active coordinate set incex can be changed for :class:`AtomGroup`

>>> a.setACSIndex(1)
>>> a
<Atom: H from Water (index 0; active #1 of 2 coordsets)>

Changing active coordinate set for an atom group, does not affect the active 
coordinate set of the atom group:

>>> wtr1
<AtomGroup: Water (3 atoms; active #0 of 2 coordsets)>

Coordinates for the atom group will be returned from the active coordinate set

>>> print( a.getCoords() )
[ 0.  1.  0.]

**Iterations**

Coordinate sets can also be iterated over for :class:`.Atom` and
:class:`AtomGroup` instances:

>>> for xyz in a.iterCoordsets(): print( xyz )
... 
[ 1.  0.  0.]
[ 0.  1.  0.]

Copying
===============================================================================

Now let's make another copy of this water:

>>> wtr2 = wtr1.copy()
>>> wtr2
<AtomGroup: Water (3 atoms; active #0 of 2 coordsets)>

**Translate copy**

Let's translate the coordinates of wtr2 so that it does not overlap with wtr1

>>> wtr2.setCoords( wtr2.getCoords() + 2 )
>>> print( wtr2.getCoords() )
[[ 3.  2.  2.]
 [ 2.  2.  2.]
 [ 2.  2.  3.]]

Above operation only translated the coordinate set at index 0

>>> wtr2.setACSIndex(1)
>>> print( wtr2.getCoords() )
[[ 0.  1.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  1.]]
>>> wtr2.setCoords( wtr2.getCoords() + 2 ) # translate the 2nd coordset as well

**Change copy attributes**

Before we merge wtr1 and wtr2, let's change resid's of wtr2:

>>> wtr2.setResnums( [2, 2, 2] )
>>> print( wtr2.getResnums() )
[2 2 2]

We can do this in an alternate way too:

>>> wtr2.select('all').setResnums(2)
>>> print( wtr2.getResnums() )
[2 2 2]

Note that the following won't work:

>>> wtr2.setResnums(2)
Traceback (most recent call last):
  File "/usr/lib/python2.6/doctest.py", line 1248, in __run
    compileflags, 1) in test.globs
  File "<doctest __main__[29]>", line 1, in <module>
    wtr2.resids = 2
  File "....", line 424, in set_resnums
    if len(resids) != self._n_atoms:
TypeError: object of type 'int' has no len()

Merging
===============================================================================

Let's merge two water atom groups:

>>> wtrs = wtr1 + wtr2
>>> wtrs
<AtomGroup: Water + Water (6 atoms; active #0 of 2 coordsets)>
>>> print( wtrs.getCoords() )
[[ 1.  0.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  1.]
 [ 3.  2.  2.]
 [ 2.  2.  2.]
 [ 2.  2.  3.]]
>>> print( wtrs.getNames() )
['H' 'O' 'H' 'H' 'O' 'H']
>>> print( wtrs.getResnums() )
[1 1 1 2 2 2]

Hierarchical view
===============================================================================

Hierarchical views of atom groups are represented by :class:`.HierView`.

Residues (and also chains) in an atom group can also be iterated over

>>> for res in wtrs.getHierView().iterResidues(): res
<Residue: WAT 1 from Water + Water (3 atoms; active #0 of 2 coordsets)>
<Residue: WAT 2 from Water + Water (3 atoms; active #0 of 2 coordsets)>

Finally, it's is possible to change the name of *wtrs* from 
"Water + Water" to something shorter:

>>> wtrs.setTitle('2Waters')
>>> wtrs
<AtomGroup: 2Waters (6 atoms; active #0 of 2 coordsets)>

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from time import time
from types import NoneType

import numpy as np

from prody import LOGGER, PY2K
from prody.kdtree import KDTree
from prody.utilities import checkCoords, rangeString

from .atomic import Atomic
from .fields import ATOMIC_FIELDS, READONLY
from .fields import wrapGetMethod, wrapSetMethod
from .flags import PLANTERS as FLAG_PLANTERS
from .flags import ALIASES as FLAG_ALIASES
from .flags import FIELDS as FLAG_FIELDS
from .atom import Atom
from .bond import Bond, evalBonds
from .selection import Selection

from . import flags

__all__ = ['AtomGroup']

if PY2K: range = xrange

def checkLabel(label):
    """Check suitability of *label* for labeling user data or flags."""
    
    label = str(label)
    if not label:
        raise ValueError('label cannot be empty string')
        
    label = str(label)

    if not label:
        raise ValueError('label cannot be empty string')

    if not label[0].isalpha():
        raise ValueError('label must start with a letter')

    if not (''.join(label.split('_'))).isalnum():
        raise ValueError('label may contain alphanumeric characters and '
                         'underscore, {0:s} is not valid'.format(label))

    if isReserved(label):
        raise ValueError('{0:s} is a reserved word and cannot be used '
                           'as a label'.format(repr(label)))

    if label in READONLY:
        raise AttributeError('{0:s} is read-only'.format(label))

    return label


class AtomGroup(Atomic):
    
    """A class for storing and accessing atomic data.  The number of atoms of 
    the atom group is inferred at the first set method call from the size of 
    the data array. 

    **Atomic data**
    
    All atomic data is stored in :class:`numpy.ndarray` instances.

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
                 '_timestamps', '_kdtrees', '_bmap', '_bonds', '_cslabels',
                 '_acsi', '_n_csets', '_data', '_fragments',
                 '_flags', '_flagsts', '_subsets']
    
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
        self._fragments = None
        
        self._cslabels = []
        self._acsi = None
        self._n_csets = 0
        
        self._data = dict()

        self._flags = None
        self._flagsts = 0
        self._subsets = None 

    def __repr__(self):

        n_csets = self._n_csets
        if n_csets == 1:
            return '<AtomGroup: {0:s} ({1:d} atoms)>'.format(
                    self._title, self._n_atoms)
        elif n_csets > 1:
            return ('<AtomGroup: {0:s} ({1:d} atoms; active #{2:d} of {3:d}' 
                    ' coordsets)>').format(self._title, self._n_atoms, 
                                           self._acsi, n_csets)
        else:
            return '<AtomGroup: {0:s} ({1:d} atoms; no coordinates)>'.format(
                    self._title, self._n_atoms)
        
    def __str__(self):
        
        return 'AtomGroup ' + self._title

    def __getitem__(self, index):
        
        acsi = self._acsi
        
        if isinstance(index, int):
            n_atoms = self._n_atoms
            if index >= n_atoms or index < -n_atoms:
                raise IndexError('index out of bounds')
            return Atom(self, index if index >= 0 else n_atoms + index, acsi)
        
        elif isinstance(index, slice):
            start, stop, step = index.indices(self._n_atoms)
            start = start or 0
            index = np.arange(start, stop, step)
            if len(index):
                if start > stop:
                    index = index[::-1]
                selstr = 'index {0:d}:{1:d}:{2:d}'.format(start, stop, step)
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
            raise TypeError('invalid index') 
    
    def __len__(self):
    
        return self._n_atoms
    
    def __add__(self, other):
        
        if not isinstance(other, AtomGroup):
            raise TypeError('unsupported operand type(s) for +: {0:s} and '
                            '{1:s}'.format(repr(type(self).__name__), 
                                           repr(type(other).__name__)))
                                                       
        new = AtomGroup(self._title + ' + ' + other._title)
        if self._n_csets:
            if self._n_csets == other._n_csets:
                new.setCoords(np.concatenate((self._coords, other._coords), 1))
                if self._n_csets > 1:
                    LOGGER.info('All {0:d} coordinate sets are copied to '
                                '{1:s}.'.format(self._n_csets, new.getTitle()))
            else:
                new.setCoords(np.concatenate((self._getCoords(), 
                                              other._getCoords())))
                LOGGER.info('Active coordinate sets are copied to {0:s}.'
                            .format(new.getTitle()))
        elif other._n_csets:
            LOGGER.warn('No coordinate sets are copied to {0:s}'
                        .format(new.getTitle()))
        
        for key in set(list(self._data) + list(other._data)):
            if key in ATOMIC_FIELDS and ATOMIC_FIELDS[key].readonly:
                continue
            this = self._data.get(key)
            that = other._data.get(key)
            if this is not None or that is not None:
                if this is None:
                    this = np.zeros(that.shape, that.dtype)
                if that is None:
                    that = np.zeros(this.shape, this.dtype)
                new._data[key] = np.concatenate((this, that))

        if self._bonds is not None and other._bonds is not None:
            new.setBonds(np.concatenate([self._bonds, 
                                         other._bonds + self._n_atoms]))
        elif self._bonds is not None:
            new.setBonds(self._bonds.copy())
        elif other._bonds is not None:
            new.setBonds(other._bonds + self._n_atoms)
        
        return new        

    def __contains__(self, item):

        try:        
            acsi = item.getACSIndex()
        except AttributeError:
            return False
        else:
            try:
                ag = item.getAtomGroup()
            except AttributeError:
                return item == self
            else:
                return ag == self
    
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
        """Return time stamp showing when coordinates were last changed."""

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
        """Return KDTree for coordinate set at given index."""

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
        """Return a mapping of serial numbers to indices."""
        
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
        """Return title of the instance."""
        
        return self._title
    
    def setTitle(self, title):
        """Set title of the instance."""
        
        self._title = str(title)
    
    def numAtoms(self, flag=None):
        """Return number of atoms, or number of atoms with given *flag*."""
        
        return len(self._getSubset(flag)) if flag else self._n_atoms
    
    def getCoords(self):
        """Return a copy of coordinates from active coordinate set."""
        
        if self._coords is not None:
            return self._coords[self._acsi].copy()
    
    def _getCoords(self): 
        """Return a view of coordinates from active coordinate set."""
        
        if self._coords is not None:
            return self._coords[self._acsi]

    def setCoords(self, coords, label=None):
        """Set coordinates of atoms.  *coords* may be any array like object
        or an object instance with :meth:`getCoords` method.  If the shape of 
        coordinate array is ``(n_csets > 1, n_atoms, 3)``, it will replace all 
        coordinate sets and the active coordinate set index  will reset to 
        zero.  This situation can be avoided using :meth:`addCoordset`.  
        If shape of *coords* is ``(n_atoms, 3)`` or ``(1, n_atoms, 3)``, it 
        will replace the active coordinate set.  *label* argument may be used 
        to label coordinate set(s).  *label* may be a string or a list of 
        strings length equal to the number of coordinate sets."""

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
                raise ValueError('coordinates of {0:s} are not set'
                                 .format(str(atoms)))
        
        try:
            checkCoords(coords, csets=True, dtype=(float, np.float32))
        except TypeError:
            raise TypeError('coords must be a numpy array or an '
                            'object with `getCoords` method')

        self._setCoords(coords, label=label)
        
    def _setCoords(self, coords, label=None, overwrite=False):
        """Set coordinates without data type checking.  *coords* must 
        be a :class:`~numpy.ndarray`, but may have data type other than 
        :class:`numpy.float64`, e.g. :class:`numpy.float32`.  *label* 
        argument may be used to label coordinate sets.  *label* may be 
        a string or a list of strings length equal to the number of 
        coordinate sets."""
        
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
                if label is None:
                    self._cslabels = [None]
                else:
                    self._cslabels = [str(label)]
                self._n_csets = n_csets = 1

            else:
                self._coords = coords
                self._n_csets = n_csets = shape[0]

                if isinstance(label, (NoneType, str)):
                    self._cslabels = [label] * n_csets
                
                elif isinstance(label, (list, tuple)):
                    if len(label) == n_csets:
                        self._cslabels = list(label)

                    else:
                        self._cslabels = [None] * n_csets
                        LOGGER.warn('Number of labels does not match number '
                                    'of coordinate sets.')
                else:
                    LOGGER.warn('Wrong type for `label` argument.')
            self._acsi = 0
            self._setTimeStamp()

        else:
            acsi = self._acsi
            if ndim == 2:
                self._coords[acsi] = coords
            else:
                self._coords[acsi] = coords[0]
            self._setTimeStamp(acsi)
            if label is not None:
                self._cslabels[acsi] = str(label)
    
    def addCoordset(self, coords, label=None):
        """Add a coordinate set.  *coords* argument may be an object with 
        :meth:`getCoordsets` method."""
        
        if self._coords is None:
            return self.setCoords(coords)

        n_atoms = self._n_atoms
        try:
            coords = (coords._getCoordsets() 
                      if hasattr(coords, '_getCoordsets') else
                      coords.getCoordsets())
        except AttributeError:
            pass
        else:
            if coords is None:
                raise ValueError('coordinates of {0:s} are not set'
                                 .format(str(atoms)))

        try:
            checkCoords(coords, csets=True, natoms=n_atoms, dtype=None)
        except TypeError:
            raise TypeError('coords must be a numpy array or an '
                            'object with `getCoords` method')
                            
        if coords.ndim == 2:
            coords = coords.reshape((1, n_atoms, 3))
            
        diff = coords.shape[0]
        self._coords = np.concatenate((self._coords, coords), axis=0)
        self._n_csets = self._coords.shape[0]
        timestamps = self._timestamps
        self._timestamps = np.zeros(self._n_csets)
        self._timestamps[:len(timestamps)] = timestamps
        self._timestamps[len(timestamps):] = time()
        self._kdtrees.extend([None] * diff)
        if isinstance(label, (str, NoneType)):
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
        """Return a copy of coordinate set(s) at given *indices*.  *indices* 
        may  be an integer, a list of integers, or **None** meaning all 
        coordinate sets."""
        
        if self._coords is None:
            return None
        if indices is None:
            return self._coords.copy()
        if isinstance(indices, (int, slice)):
            return self._coords[indices].copy()
    
        # following fancy indexing makes a copy, so .copy() is not needed
        if isinstance(indices, (list, np.ndarray)): 
            return self._coords[indices]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')
        
    def _getCoordsets(self, indices=None):
        """Return a view of coordinate set(s) at given *indices*."""
        
        if self._coords is None:
            return None

        try:
            return self._coords if indices is None else self._coords[indices] 
        except:
            raise IndexError('indices must be an integer, a list/array of '
                               'integers, a slice, or None')

    def numCoordsets(self):
        """Return number of coordinate sets."""
        
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
        """Return index of the coordinate set."""
        
        return self._acsi 
            
    def setACSIndex(self, index):
        """Set the coordinate set at *index* active."""
        
        n_csets = self._n_csets
        if n_csets == 0:
            self._acsi = 0
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if n_csets <= index or n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += n_csets 
        self._acsi = index
        
    def getHierView(self):
        """Return a hierarchical view of the atom group."""
        
        if self._hv is None:
            self._hv = HierView(self)
        return self._hv
    
    def numSegments(self):
        """Return number of segments."""
        
        return self.getHierView().numSegments()

    def numChains(self):
        """Return number of chains."""
        
        return self.getHierView().numChains()
    
    def numResidues(self):
        """Return number of residues."""
        
        return self.getHierView().numResidues()

    def iterSegments(self):
        """Iterate over chains."""
        
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

            try:
                ndim, dtype, shape = data.ndim, data.dtype, data.shape
            except AttributeError:
                data = np.array(data)
                ndim, dtype, shape = data.ndim, data.dtype, data.shape

            if ndim == 1 and dtype == bool:
                raise TypeError('1 dimensional boolean arrays are not '
                                  'accepted, use `setFlags` instead')

            if len(data) != self._n_atoms:
                raise ValueError('len(data) must match number of atoms')

            self._data[label] = data
    
    def delData(self, label):
        """Return data associated with *label* and remove from the instance.
        If data associated with *label* is not found, return **None**."""
        
        return self._data.pop(label, None)
    
    def getData(self, label):
        """Return a copy of the data array associated with *label*, or **None** 
        if such data is not present."""
        
        data = self._getData(label)
        if data is not None:
            return data.copy()

    def _getData(self, label):
        """Return data array associated with *label*, or **None** if such data 
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

    def isData(self, label):
        """Deprecated, use :meth:`isDataLabel` instead."""
        
        from prody import deprecate
        deprecate('isData', 'isDataLabel', 'v1.3')
        return self.isDataLabel(label)
        
    def isDataLabel(self, label):
        """Return **True** if data associated with *label* is present."""
        
        if label in self._data:
            return True
        else:
            try:
                return self._getData(label) is not None
            except:
                return False
  
    def getDataLabels(self, which=None):
        """Return data labels.  For ``which='user'``, return only labels of 
        user provided data."""

        if str(which).startswith('u'): # user
            labels = [key for key in (self._data or {}) 
                      if not key in ATOMIC_FIELDS]
        else:
            labels = list(self._data or [])
        labels.sort()
        return labels
        
    def getDataType(self, label):
        """Return type of the data (i.e. ``data.dtype``) associated with 
        *label*, or **None** label is not used."""

        try:
            return self._data[label].dtype
        except KeyError:
            return None
    
    def isFlagLabel(self, label):
        """Return **True** if flags associated with *label* are present."""
        
        return label in FLAG_PLANTERS or label in (self._flags or {})
    
    def getFlags(self, label):
        """Return a copy of atom flags for given *label*, or **None** when 
        flags for *label* is not set."""
        
        flags = self._getFlags(label)
        if flags is not None:
            return flags.copy()
        
    def _getFlags(self, label):
        """Return atom flag values for given *label*, or **None** when 
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
        """Return flags associated with *label* and remove from the instance.
        If flags associated with *label* is not found, return **None**."""
        
        return self._flags.pop(label, None)
    
    def _setSubset(self, label, indices):
        """Set indices of a subset of atoms."""

        for label in FLAG_ALIASES.get(label, [label]):
            self._subsets[label] = indices
            
    def _getSubset(self, label):
        """Return indices of atoms."""
        
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
        """Return flag labels.  For ``which='user'``,  return labels of user 
        or parser (e.g. :term:`hetatm`) provided flags, for ``which='all'`` 
        return all possible :ref:`flags` labels in addition to those present 
        in the instance."""

        which = str(which)
        if which.startswith('a'): # all possible
            labels = set(self._flags or [])
            labels.update(FLAG_PLANTERS)
            labels = list(labels)
        elif which.startswith('u'): # user
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

        if not isinstance(serial, int):
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
            if not isinstance(stop, int):
                raise TypeError('stop must be an integer')
            if stop <= serial:
                raise ValueError('stop must be greater than serial')
                
            if step is None:
                step = 1
            else:
                if not isinstance(step, int):
                    raise TypeError('step must be an integer')
                if step < 1:
                    raise ValueError('step must be greater than zero')
            
            indices = sn2i[serial:stop:step]
            indices = indices[indices > -1]
            return Selection(self, indices, 'serial {0:d}:{1:d}:{2:d}'
                                            .format(serial, stop, step))

    def getACSLabel(self):
        """Return active coordinate set label."""
        
        if self._n_csets:
            return self._cslabels[self._acsi]

    def setACSLabel(self, label):
        """Set active coordinate set label."""

        if self._n_csets:
            if isinstance(label, (str, NoneType)):
                self._cslabels[self._acsi] = label 
            else:
                raise TypeError('label must be a string')
    
    def getCSLabels(self):
        """Return coordinate set labels."""
        
        if self._n_csets:
            return list(self._cslabels)

    def setCSLabels(self, labels):
        """Set coordinate set labels. *labels* must be a list of strings."""
        
        if isinstance(labels, list):
            if len(labels) == self._n_csets:
                if all(isinstance(lbl, (str, NoneType)) for lbl in labels):
                    self._cslabels = list(labels)
                else:
                    raise ValueError('all items of labels must be strings')
            else:
                raise ValueError('length of labels must be equal to the '
                                 'number of coordinate sets')
        else:
            raise TypeError('labels must be a list')    
            
    def setBonds(self, bonds):
        """Set covalent bonds between atoms.  *bonds* must be a list or an
        array of pairs of indices.  All bonds must be set at once.  Bonding
        information can be used to make atom selections, e.g. ``"bonded to 
        index 1"``.  See :mod:`.select` module documentation for details.   
        Also, a data array with number of bonds will be generated and stored 
        with label *numbonds*.  This can be used in atom selections, e.g. 
        ``'numbonds 0'`` can be used to select ions in a system."""
        
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
        bonds = bonds[bonds[:,1].argsort(),]
        bonds = bonds[bonds[:,0].argsort(),]
        
        self._bmap, self._data['numbonds'] = evalBonds(bonds, n_atoms)
        self._bonds = bonds
        self._fragments = None

    def numBonds(self):
        """Return number of bonds.  Use :meth:`setBonds` for setting bonds."""
        
        if self._bonds is not None:
            return self._bonds.shape[0]
        return 0

    def iterBonds(self):
        """Yield bonds.  Use :meth:`setBonds` for setting bonds."""
        
        if self._bonds is not None:
            acsi = self._acsi
            for bond in self._bonds:
                yield Bond(self, bond, acsi)
    
    def _iterBonds(self):
        """Yield pairs of bonded atom indices. Use :meth:`setBonds` for setting
        bonds."""
        
        if self._bonds is not None:
            for a, b in self._bonds:     
                yield a, b

    def numFragments(self):
        """Return number of connected atom subsets."""
        
        self._fragment()
        return self._data['fragindex'].max() + 1
        
    def iterFragments(self):
        """Yield connected atom subsets as :class:`.Selection` instances."""
        
        if self._bmap is not None:

            frags = self._data.get('fragindex')
            acsi = self._acsi
            if self._fragments is None:
                self._fragment()
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
            raise ValueError('bonds must be set for fragment determination, '
                             'use `setBonds`')
        
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
                    return self._data[var].copy()
                except KeyError:
                    [getattr(self, meth)() for meth in call]
                    return self._data[var].copy()
        # Define private method for retrieving actual data array
        def _getData(self, var=fname, call=field.call):
            try:
                return self._data[var]
            except KeyError:
                [getattr(self, meth)() for meth in call]
                return self._data[var].copy()
    else:
        if not field.private:
            def getData(self, var=fname):
                try:
                    return self._data[var].copy()
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
            if self._n_atoms == 0:
                self._n_atoms = len(array)
            elif len(array) != self._n_atoms:
                raise ValueError('length of array must match number '
                                   'of atoms')
                
            if isinstance(array, list):
                array = np.array(array, dtype)
            elif not isinstance(array, np.ndarray):
                raise TypeError('array must be an ndarray or a list')
            elif array.ndim != ndim:
                    raise ValueError('array must be {0:d} '
                                       'dimensional'.format(ndim))
            elif array.dtype != dtype:
                try:
                    array = array.astype(dtype)
                except ValueError:
                    raise ValueError('array cannot be assigned type '
                                     '{0:s}'.format(dtype))
            self._data[var] = array
            if none: self._none(none)
            if flags and self._flags:
                self._resetFlags(var)
            
    setData = wrapSetMethod(setData)
    setData.__name__ = setMeth 
    setData.__doc__ = field.getDocstr('set')
    setattr(AtomGroup, setMeth, setData)

del getData
del _getData
del setData
