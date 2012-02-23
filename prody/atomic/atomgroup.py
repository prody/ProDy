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

.. _atomgroup:

Constructing an :class:`AtomGroup`
===============================================================================

We start by importing everything from the ProDy package and the NumPy package:

>>> from prody import *
>>> import numpy as np

Instantiate an AtomGroup
-------------------------------------------------------------------------------

>>> wtr1 = AtomGroup('Water')
>>> wtr1
<AtomGroup: Water (0 atoms; no coordinates)>


Set coordinates
-------------------------------------------------------------------------------

The best way to start constructing an :class:`AtomGroup` is by setting the
coordinates first. Number of atoms will be automatically set according to
the size of the coordinate data array:

>>> coords = np.array( [ [1, 0, 0], [0, 0, 0], [0, 0, 1] ] )
>>> print( coords )
[[1 0 0]
 [0 0 0]
 [0 0 1]]
>>> wtr1.setCoords( coords )
>>> wtr1
<AtomGroup: Water (3 atoms)>

Set attributes
-------------------------------------------------------------------------------

Attributes must be passed in a list or an array whose size is the same
as the number of atoms.

>>> wtr1.setNames( ['H', 'O', 'H'] )
>>> wtr1.setResnums( [1, 1, 1] )
>>> wtr1.setResnames( ['WAT', 'WAT', 'WAT'] )

Accessing data will return a copy of the data:

>>> print( wtr1.getNames() )
['H' 'O' 'H']

Individual atoms
-------------------------------------------------------------------------------

Atoms are represented by instance of :class:`~.Atom`.

**Iteration**

Atoms in an :class:`AtomGroup` can be iterated over

>>> for a in wtr1: a
... 
<Atom: H from Water (index 0)>
<Atom: O from Water (index 1)>
<Atom: H from Water (index 2)>

**Indexing**

Atoms in an atom group can be accessed via indexing:

>>> a = wtr1[0]
>>> a
<Atom: H from Water (index 0)>
>>> print( a.getCoords() )
[ 1.  0.  0.]

Coordinate sets
-------------------------------------------------------------------------------

Let's add another coordinate set to the atom group:

>>> wtr1.addCoordset( np.array( [ [0, 1, 0], [0, 0, 0], [0, 0, 1] ] ) )
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

Coordinate sets can also be iterated over for :class:`~.Atom` and 
:class:`AtomGroup` instances:

>>> for xyz in a.iterCoordsets(): print( xyz )
... 
[ 1.  0.  0.]
[ 0.  1.  0.]

Copy atom groups
-------------------------------------------------------------------------------

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
>>> wtr2.setCoords( wtr2.getCoords() + 2 ) # translate the second coordinate set as well

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

Merge atom groups
-------------------------------------------------------------------------------

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
-------------------------------------------------------------------------------

Hierarchical views of atom groups are represented by :class:`~.HierView`.

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

from prody.tools import checkCoords
from prody.KDTree import getKDTree

from atomic import Atomic
from fields import ATOMIC_ATTRIBUTES, ATOMIC_DATA_FIELDS, READONLY
from fields import wrapGetMethod, wrapSetMethod
from atom import Atom
from bond import Bond, evalBonds, trimBonds
from atommap import AtomMap
from subset import AtomSubset
from selection import Selection



__all__ = ['AtomGroup']

pkg = __import__(__package__)
LOGGER = pkg.LOGGER

SELECT = None

class AtomGroupMeta(type):

    def __init__(cls, name, bases, dict):
    
        for field in ATOMIC_DATA_FIELDS.values():

            meth = field.meth_pl
            getMeth = 'get' + meth
            setMeth = 'set' + meth
            # Define public method for retrieving a copy of data array
            if field.call:
                def getData(self, var=field.var, call=field.call):
                    for meth in call:
                        getattr(self, meth)()
                    array = self._data[var]
                    return array.copy()                 
            else:
                def getData(self, var=field.var):
                    array = self._data[var]
                    if array is None:
                        return None
                    return array.copy() 
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            getData.__doc__ = field.getDocstr('get')
            setattr(cls, getMeth, getData)
            
            # Define private method for retrieving actual data array
            def _getData(self, var=field.var):
                return self._data[var]
            _getData = wrapGetMethod(_getData)
            _getData.__name__ = '_' + getMeth
            _getData.__doc__ = field.getDocstr('_get')
            setattr(cls, '_' + getMeth, _getData)
            
            if field.readonly:
                continue
            
            # Define public method for setting values in data array
            def setData(self, array, var=field.var, dtype=field.dtype, 
                        ndim=field.ndim, none=field.none):
                if self._n_atoms == 0:
                    self._n_atoms = len(array)
                elif len(array) != self._n_atoms:
                    raise ValueError('length of array must match numAtoms')
                    
                if isinstance(array, list):
                    array = np.array(array, dtype)
                elif not isinstance(array, np.ndarray):
                    raise TypeError('array must be an ndarray or a list')
                elif array.ndim != ndim:
                        raise ValueError('array must be {0:d} dimensional'
                                         .format(ndim))
                elif array.dtype != dtype:
                    try:
                        array = array.astype(dtype)
                    except ValueError:
                        raise ValueError('array cannot be assigned type '
                                         '{0:s}'.format(dtype))
                self._data[var] = array
                if none:
                    self.__setattr__('_'+none,  None)
            setData = wrapSetMethod(setData)
            setData.__name__ = setMeth 
            setData.__doc__ = field.getDocstr('set')
            setattr(cls, setMeth, setData)


class AtomGroup(Atomic):
    
    """A class for storing and accessing atomic data.  The number of atoms of 
    the atom group is inferred at the first set method call from the size of 
    the data array. 

    **Atomic Data**
    
    All atomic data is stored in :class:`numpy.ndarray` instances.

    **Get and Set Methods**
    
    *get* methods, e.g. :meth:`getResnames`, return copies of the data arrays. 
    
    *set* methods, e.g. :meth:`setResnums`, accept data in :class:`list` or 
    :class:`~numpy.ndarray` instances.  The length of the list or array must 
    match the number of atoms in the atom group.  These methods set attributes 
    of all atoms at once.
    
    Atom groups with multiple coordinate sets may have one of these sets as 
    the *active coordinate set*.  The active coordinate set may be changed 
    using :meth:`setACSIndex()` method.  :meth:`getCoords` returns coordinates
    from the *active set*.
    
    To access and modify data associated with a subset of atoms in an atom 
    group, :class:`~.Selection` instances may be used.  A selection from an 
    atom group has initially the same coordinate set as the *active coordinate 
    set*.
    
    Some :class:`object` methods are customized as follows:
    
    * :func:`len` returns the number of atoms, i.e. :meth:`numAtoms`
    * :func:`iter` yields :class:`~.Atom` instances
    * indexing by:
         - *int* (:func:`int`), e.g, ``10``, returns an :class:`~.Atom`
         - *slice* (:func:`slice`), e.g, ``10:20:2``, returns a 
           :class:`~.Selection`
         - *segment name* (:func:`str`), e.g. ``"PROT"``, returns a 
           a :class:`~.Segment` 
         - *chain identifier* (:func:`str`), e.g. ``"A"``, returns a 
           a :class:`~.Chain`
         - *[segment name,] chain identifier, residue number[, insertion code]* 
           (:func:`tuple`), e.g. ``"A", 10`` or  ``"A", 10, "B"`` or
           ``"PROT", "A", 10, "B"``, returns a :class:`~.Residue`
    """
    
    __metaclass__ = AtomGroupMeta
    
    __slots__ = ['_title', '_n_atoms', '_coords', '_hv', '_sn2i', 
                 '_timestamps', '_kdtrees', '_bmap', '_bonds', '_cslabels',
                 '_acsi', '_n_csets', '_data']
    
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
        
        self._cslabels = []
        self._acsi = None
        self._n_csets = 0
        
        self._data = dict([(field.var, None) 
                           for field in ATOMIC_DATA_FIELDS.itervalues()])

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
            if index < 0:
                index = n_atoms + index
            return Atom(self, index, acsi)
        
        elif isinstance(index, slice):
            start, stop, step = index.indices(self._n_atoms)
            if start is None:
                start = 0
            if step is None:
                step = 1
            index = np.arange(start,stop,step)
            if len(index) > 0:
                selstr = 'index {0:d}:{1:d}:{2:d}'.format(start, stop, step)
                return Selection(self, index, selstr, acsi)
        
        elif isinstance(index, (list, np.ndarray)):
            unique = np.unique(index)
            if unique[0] < 0 or unique[-1] >= self._n_atoms:
                raise IndexError('index out of range')
            return Selection(self, unique,  
                             'index ' + ' '.join(np.array(index, '|S')), 
                             acsi, unique=True)
        
        elif isinstance(index, (str, tuple)):
            return self.getHierView()[index]
        
        else:
            raise TypeError('invalid index') 
    
    def __len__(self):
    
        return self._n_atoms
    
    def __add__(self, other):
        
        if not isinstance(other, AtomGroup):
            raise TypeError('can only concatenate two AtomGroup`s or can '
                            'deform AtomGroup along a Vector/Mode')
            
        new = AtomGroup(self._title + ' + ' + other._title)
        n_csets = self._n_csets
        if n_csets != other._n_csets:
            LOGGER.warning('AtomGroups {0:s} and {1:s} do not have same '
                           'number of coordinate sets.  First from both '
                           'AtomGroups will be merged.'
              .format(str(self._title), str(other._title), n_csets))
            n_csets = 1
        coordset_range = range(n_csets)
        new.setCoords(np.concatenate((self._coords[coordset_range],
                                      other._coords[coordset_range]), 1))
        
        for key in set(self._data.keys() + other._data.keys()):
            if key in ATOMIC_ATTRIBUTES and \
                ATOMIC_ATTRIBUTES[key].readonly:
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
        
        torf = False
        if isinstance(item, Atomic):
            torf = self == item or self == item.getAtomGroup()
        return torf
    
    def __iter__(self):
        """Yield atom instances."""
        
        acsi = self._acsi
        for index in xrange(self._n_atoms):
            yield Atom(self, index, acsi)

    iterAtoms = __iter__

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
                kdtree = getKDTree(self._coords[index])
                self._kdtrees[index] = kdtree
            return kdtree
        else:
            return None

    def _getSN2I(self):
        """Return a mapping of serial numbers to indices."""
        
        if self._sn2i is None:
            serials = self._serials  
            if serials is None:
                raise AttributeError('atom serial numbers are not set')
            unique = np.unique(serials) 
            if len(unique) != self._n_atoms:
                raise ValueError('atom serial numbers must be unique')
            if unique[0] < 0:
                raise ValueError('atoms must not have negative serial numbers')
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
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
    
    def getCoords(self):
        """Return a copy of coordinates from active coordinate set."""
        
        if self._coords is not None:
            return self._coords[self._acsi].copy()
    
    def _getCoords(self): 
        """Return a view of coordinates from active coordinate set."""
        
        if self._coords is not None:
            return self._coords[self._acsi]

    def setCoords(self, coords, label=None):
        """Set coordinates.  *coords* may be a :class:`numpy.ndarray` instance 
        or an object instance with ``getCoordsets`` method.  If the shape of 
        the coordinates array is (n_csets,n_atoms,3), the given array will 
        replace all coordinate sets.  To avoid it, :meth:`addCoordset` may be 
        used.  If the shape of the *coords* array is (n_atoms,3) or 
        (1,n_atoms,3), the coordinate set will replace the coordinates of the 
        currently active set.  *label* argument may be used to label coordinate
        sets.  *label* may be a string or a list of strings length equal to 
        the number of coordinate sets."""

        if not isinstance(coords, np.ndarray):
            coords = coords.getCoordsets()

        coords = checkCoords(coords, 'coords',
                                  cset=True, n_atoms=self._n_atoms,
                                  reshape=True)
        if self._n_atoms == 0:
            self._n_atoms = coords.shape[-2] 
            
        acsi = None
        if self._coords is None:
            self._coords = coords
            self._n_csets = coords.shape[0]
            self._acsi = 0
            self._setTimeStamp()
            if isinstance(label, (NoneType, str)):
                self._cslabels = [label] * self._n_csets
            elif isinstance(label, (list, tuple)):
                if len(label) == self._n_csets:
                    self._cslabels = label
                else:
                    self._cslabels = [None] * self._n_csets
                    LOGGER.warning('Length of `label` does not match number '
                                   'of coordinate sets.')
        else:
            if coords.shape[0] == 1:
                acsi = self._acsi
                self._coords[acsi] = coords[0]
                self._setTimeStamp(acsi)
                if isinstance(label, str):
                    self._cslabels[self._acsi] = label
            else:
                self._coords = coords
                self._n_csets = coords.shape[0]
                self._acsi = min(self._n_csets - 1, self._acsi)
                self._setTimeStamp()
                
        if acsi is None:
            if isinstance(label, (str, NoneType)):
                self._cslabels = [label] * self._n_csets
            elif isinstance(label, (list, tuple)):
                if len(label) == self._n_csets:
                    if all([isinstance(lbl, str) for lbl in label]):
                        self._cslabels += label
                    else:
                        LOGGER.warning('all items of label must be strings')
                else:
                    LOGGER.warning('label must have same length as the '
                                   'coords array')
            else:
                LOGGER.warning('label must be a string or list of strings')
                
        elif label is not None:
            if isinstance(label, str):
                self._cslabels[acsi] = label
            elif isinstance(label, (list, tuple)):
                if len(label) == 1:
                    if isinstance(label[0], str):
                        self._cslabels[acsi] = label
                    else:
                        LOGGER.warning('all items of label must be strings')
                else:
                    LOGGER.warning('length of label must be one')
            else:
                LOGGER.warning('label must be a string or list of strings')    
    
    def _setCoords(self, coords, label=None):
        """Set coordinates without date type checking.  *coords* must 
        be a :class:`numpy.ndarray`, but may have data type other than 
        :class:`numpy.float64`, e.g. :class:`numpy.float32`.  *label* 
        argument may be used to label coordinate sets.  *label* may be 
        a string or a list of strings length equal to the number of 
        coordinate sets."""

        n_atoms = self._n_atoms
        if self._n_atoms: 
            if coords.shape[-2] != self._n_atoms:
                raise ValueError('coords have incorrect number of atoms')
        else:
            n_atoms = coords.shape[-2]
            self._n_atoms = n_atoms
        
        if coords.ndim == 2:
            coords = coords.reshape((1, n_atoms, 3))
            self._cslabels = [label]
            self._n_csets = 1
        else:
            self._n_csets = coords.shape[0]
            if isinstance(label, (NoneType, str)):
                self._cslabels = [label] * self._n_csets
            elif isinstance(label, (list, tuple)):
                if len(label) == self._n_csets:
                    self._cslabels = label
                else:
                    self._cslabels = [None] * self._n_csets
                    LOGGER.warning('Length of `label` does not match number '
                                   'of coordinate sets.')
        self._coords = coords
        self._acsi = 0
        self._setTimeStamp()
    
    def addCoordset(self, coords, label=None):
        """Add a coordinate set.  *coords* argument may be an object instance 
        with ``getCoordsets`` method."""
        
        if not isinstance(coords, np.ndarray):
            coords = coords.getCoordsets()

        if self._coords is None:
            self.setCoords(coords)
            return

        coords = checkCoords(coords, 'coords', cset=True, 
                             n_atoms=self._n_atoms, reshape=True)
        diff = coords.shape[0]
        self._coords = np.concatenate((self._coords, coords), axis=0)
        self._n_csets = self._coords.shape[0]
        timestamps = self._timestamps
        self._timestamps = np.zeros(self._n_csets)
        self._timestamps[:len(timestamps)] = timestamps
        self._timestamps[len(timestamps):] = time()
        self._kdtrees.extend([None] * diff)
        if isinstance(label, (str, NoneType)):
            self._cslabels += [label] * diff
        elif isinstance(label, (list, tuple)):
            if len(label) == diff:
                if all([isinstance(lbl, str) for lbl in label]):
                    self._cslabels += label
                else:
                    LOGGER.warning('all items of `label` must be strings')
            else:
                LOGGER.warning('`label` list must have same length as the '
                               '`coords` array')
        else:
            LOGGER.warning('`label` must be a string or list of strings')
        
    def delCoordset(self, index):
        """Delete a coordinate set from the atom group."""
        
        if self._n_csets == 0:
            raise AttributeError('coordinates are not set')

        which = np.ones(self._n_csets, bool)
        which[index] = False
        n_csets = self._n_csets
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
        may  be an integer, a list of integers, or ``None`` meaning all 
        coordinate sets."""
        
        if self._coords is None:
            return None
        if indices is None:
            return self._coords.copy()
        if isinstance(indices, (int, slice)):
            return self._coords[indices].copy()
        if isinstance(indices, (list, np.ndarray)):
            return self._coords[indices]
        raise IndexError('indices must be an integer, a list/array of '
                         'integers, a slice, or None')
        
    def _getCoordsets(self, indices=None):
        """Return a view of coordinate set(s) at given *indices*."""
        
        if self._coords is None:
            return None
        if indices is None:
            return self._coords
        if isinstance(indices, (int, slice, list, np.ndarray)):
            return self._coords[indices]
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
        
    def copy(self, which=None):
        """Return a copy of atoms indicated *which* as a new AtomGroup 
        instance.
        
        *which* may be:
            * ``None``, make a copy of the AtomGroup
            * a Selection, Residue, Chain, or Atom instance
            * a list or an array of indices
            * a selection string"""
        
        title = self._title
        if which is None:
            indices = None
            newmol = AtomGroup('{0:s}'.format(title))
            newmol.setCoords(self._coords.copy())
            
        elif isinstance(which, int):
            indices = [which]
            newmol = AtomGroup('{0:s} index {1:d}'.format(title, which))
            
        elif isinstance(which, str):
            indices = SELECT.getIndices(self, which)
            if len(indices) == 0:
                return None
            newmol = AtomGroup('{0:s} selection "{1:s}"'.format(title, which))
            
        elif isinstance(which, (list, np.ndarray)):
            if isinstance(which, list):
                indices = np.array(which)
            elif which.ndim != 1:
                raise ValueError('which must be a 1d array')
            else:
                indices = which
            newmol = AtomGroup('{0:s} subset'.format(title))
            
        else:
            if isinstance(which, Atom):
                indices = [which.getIndex()]
            elif isinstance(which, (AtomSubset, AtomMap)):
                indices = which.getIndices()
            else:
                raise TypeError('{0:s} is not a valid type'.format(
                                                                type(which)))            
            newmol = AtomGroup('{0:s} selection "{1:s}"'.format(title, 
                                                                str(which)))
        if indices is not None:
            newmol.setCoords(self._coords[:, indices])
        for key, array in self._data.iteritems():
            if key == 'numbonds':
                continue
            if array is not None:
                if indices is None:
                    newmol._data[key] = array.copy()
                else:
                    newmol._data[key] = array[indices]
        
        newmol._cslabels = list(self._cslabels)
        bonds = self._bonds
        bmap = self._bmap
        if bonds is not None and bmap is not None:
            if indices is None:
                newmol._bonds = bonds.copy()
                newmol._bmap = bmap.copy()
                newmol._data['numbonds'] = self._data['numbonds'].copy()
            else:
                bonds = trimBonds(bonds, indices)
                if bonds is not None:
                    newmol.setBonds(bonds)
        return newmol
    
    __copy__ = copy
    
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
        """Store atomic *data* under *label*.
        
        *label* must:
            
            * start with a letter
            * contain only alphanumeric characters and underscore
            * not be a reserved word 
              (see :func:`~.getReservedWords`)

        *data* must be a :func:`list` or a :class:`numpy.ndarray`, its length 
        must be equal to the number of atoms, and the type of data array must 
        be one of:
            
            * :func:`bool`
            * :func:`float`
            * :func:`int`
            * :func:`str`
        
        If a :class:`list` is given, its type must match one of the above after 
        it is converted to an :class:`numpy.ndarray`.  If the dimension of the 
        *data* array is 1 (i.e. ``data.ndim==1``), *label* can be used to make
        atom selections, e.g. ``"label 1 to 10"`` or ``"label C1 C2"``.  Note 
        that, if data with *label* is present, it will be overridden."""
        
        if not isinstance(label, str):
            raise TypeError('label must be a string')
        if label == '':
            raise ValueError('label cannot be empty string')
        if not label[0].isalpha():
            raise ValueError('label must start with a letter')
        if label in READONLY:
            raise AttributeError("{0:s} is read-only".format(label))
        if not (''.join((''.join(label.split('_'))).split())).isalnum():
            raise ValueError('label may contain alphanumeric characters and '
                             'underscore, {0:s} is not valid'.format(label))
            
        if isReserved(label):
            raise ValueError('label cannot be a reserved word or a selection '
                             'keyword, "{0:s}" is invalid'.format(label))
        if len(data) != self._n_atoms:
            raise ValueError('length of data array must match number of atoms')
        if isinstance(data, list):
            data = np.array(data)
        elif not isinstance(data, np.ndarray):
            raise TypeError('data must be a numpy.ndarray instance')
        if not data.dtype in (np.float, np.int, np.bool) and \
              data.dtype.type != np.string_:
            raise TypeError('type of data array must be float, int, or '
                            'string_, {0:s} is not valid'.format(
                            str(data.dtype)))
            
        self._data[label] = data
    
    def delData(self, label):
        """Return data associated with *label* and remove it from the atom 
        group.  If data associated with *label* is not found, ``None`` will 
        be returned."""
        
        if not isinstance(label, str):
            raise TypeError('label must be a string')
        return self._data.pop(label, None)
    
    def getData(self, label):
        """Return a copy of the data array associated with *label*, or ``None`` 
        if such data is not present."""
        
        data = self._data.get(label, None)
        if data is None:
            return None
        else:
            return data.copy()

    def _getData(self, label):
        """Return data array associated with *label*, or ``None`` if such data 
        is not present."""
        
        data = self._data.get(label, None)
        if data is None:
            return None
        else:
            return data

    def isData(self, label):
        """Return **True** if data with *label* is set by user."""
        
        return label in self._data and self._data[label] is not None
  
    def getDataLabels(self):
        """Return list of user data labels."""
        
        return [key for key, data in self._data.iteritems() 
                    if data is not None]
        
    def getDataType(self, label):
        """Return type of the user data (i.e. data.dtype) associated with
        *label*, or ``None`` label is not used."""
        
        try:
            return self._data[label].dtype
        except KeyError:
            return None
    
    def getBySerial(self, serial, stop=None, step=None):
        """Get an atom(s) by *serial* number (range).  *serial* must be zero or 
        a positive integer. *stop* may be ``None``, or an integer greater than 
        *serial*.  ``getBySerial(i, j)`` will return atoms whose serial numbers
        are i+1, i+2, ..., j-1.  Atom whose serial number is *stop* will be 
        excluded as it would be in indexing a Python :class:`list`.  *step* 
        (default is 1) specifies increment.  If atoms with matching serial 
        numbers are not found, ``None`` will be returned."""

        if not isinstance(serial, int):
            raise TypeError('serial must be an integer')
        if serial < 0:
            raise ValueError('serial must be greater than or equal to zero')
        sn2i = self._getSN2I()
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
        array of pairs of indices.  All bonds must be set at once.  An array
        with number of bonds will be generated and stored as *numbonds*.
        This can be used in atom selections, e.g. ``ag.select('numbonds 0')``
        can be used to select ions in a system."""
        
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

    def numBonds(self):
        """Return number of bonds.  Bonds must be set using :meth:`setBonds`.
        """
        
        if self._bonds is not None:
            return self._bonds.shape[0]

    def iterBonds(self):
        """Yield bonds.  Bonds must be set using :meth:`setBonds`."""
        
        if self._bonds is not None:
            acsi = self._acsi
            for bond in self._bonds:
                yield Bond(self, bond, acsi)
