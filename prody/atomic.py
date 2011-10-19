# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

"""This module defines classes for handling atomic data.

.. _atomic:
    
Atomic data
===============================================================================

ProDy stores atomic data in instances of :class:`AtomGroup` class. This
class is designed to be efficient and responsive, i.e. facilitates user
to access atomic data quickly for any subset of atoms. An :class:`AtomGroup`
instance can be obtained by parsing a PDB file as follows: 
    
>>> from prody import *
>>> ag = parsePDB('1aar')

To read this page in a Python session, type:
    
>>> # help(atomic)


:class:`AtomGroup` instances can store multiple coordinate sets, which may
be models from an NMR structure, snapshots from an MD simulation.


ProDy stores all atomic data in :class:`AtomGroup` instances and comes
with other classes acting as pointers to provide convenient read/write access 
to such data. These classes are:

* :class:`Atom` - Points to a single atom in an :class:`AtomGroup` instance.                          

* :class:`Selection` - Points to an arbitrary subset of atoms. See 
  :ref:`selections` and :ref:`selection-operations` for usage examples.

* :class:`Chain` - Points to atoms that have the same chain identifier.

* :class:`Residue` - Points to atoms that have the same chain identifier, 
  residue number and insertion code.
                      
* :class:`AtomMap` - Points to arbitrary subsets of atoms while allowing for 
  duplicates and missing atoms. Indices of atoms are stored in the order 
  provided by the user.
    

Atom selections
-------------------------------------------------------------------------------

Flexible and powerful atom selections is one of the most important features 
of ProDy. The details of the selection grammar is described in 
:ref:`selections`. 

.. versionadded:: 0.7.1

Using the flexibility of Python, atom selections are made much easier by
overriding the ``.`` operator i.e. the :meth:`__getattribute__` 
method of :class:`Atomic` class. So the following will be interpreted
as atom selections:
    
>>> ag.chain_A # selects chain A
<Selection: "chain A" from 1aar (608 atoms; 1 coordinate sets, active set index: 0)>
>>> ag.calpha # selects alpha carbons
<Selection: "calpha" from 1aar (152 atoms; 1 coordinate sets, active set index: 0)>
>>> ag.resname_ALA # selects alanine residues
<Selection: "resname ALA" from 1aar (20 atoms; 1 coordinate sets, active set index: 0)>

It is also possible to combine selections with ``and`` and ``or`` operators:

>>> ag.chain_A_and_backbone
<Selection: "chain A and backbone" from 1aar (304 atoms; 1 coordinate sets, active set index: 0)>
>>> ag.acidic_or_basic
<Selection: "acidic or basic" from 1aar (422 atoms; 1 coordinate sets, active set index: 0)>


Using dot operator will behave like the logical ``and`` operator:
    
>>> ag.chain_A.backbone
<Selection: "(backbone) and (chain A)" from 1aar (304 atoms; 1 coordinate sets, active set index: 0)>
  
For this to work, the first word following the dot operator must be a selection
keyword, e.g. ``resname``, ``name``, ``apolar``, ``protein``, etc. 
Underscores will be interpreted as white space, as obvious from the
previous examples. The limitation of this is that parentheses, special 
characters cannot be used.     

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

from collections import defaultdict
import time

import numpy as np

import prody
LOGGER = prody.LOGGER
checkCoordsArray = prody.checkCoordsArray 

__all__ = ['Atomic', 'AtomGroup', 'AtomPointer', 'Atom', 'AtomSubset', 
           'Selection', 'Chain',
           'Residue', 'AtomMap', 'HierView', 'ATOMIC_DATA_FIELDS',
           'loadAtoms', 'saveAtoms',]

def isBooleanKeyword(name):
    return True

class Field(object):
    __slots__ = ('_name', '_var', '_dtype',  '_doc', '_doc_pl', 
                 '_meth', '_meth_pl', '_ndim', '_none')
    def __init__(self, name, var, dtype, doc, meth, 
                 doc_pl=None, meth_pl=None, ndim=1, none=None):
        self._name = name
        self._var = var
        self._dtype = dtype
        self._doc = doc
        self._ndim = ndim
        if doc_pl is None:
            self._doc_pl = doc + 's'
        else:
            self._doc_pl = doc_pl
        self._meth = meth
        if meth_pl is None:
            self._meth_pl = meth + 's'
        else:
            self._meth_pl = meth_pl
        self._none = none
        
    def name(self):
        return self._name
    name = property(name, 
        doc='Atomic data field name.')
    def var(self):
        return self._var
    var = property(var, 
        doc='Atomic data variable name.')
    def dtype(self):
        return self._dtype
    dtype = property(dtype, 
        doc='Atomic data type (NumPy types).')
    def doc(self):
        return self._doc
    doc = property(doc, 
        doc='Atomic data field name, used in documentation.')
    def doc_pl(self):
        return self._doc_pl
    doc_pl = property(doc_pl, 
        doc='Atomic data field name in plural form, used in documentation.')
    def meth(self):
        return self._meth
    meth = property(meth, 
        doc='Atomic data field name, used in get/set methods.')
    def meth_pl(self):
        return self._meth_pl
    meth_pl = property(meth_pl, 
        doc='Atomic data field name in plural form.')
    def ndim(self):
        return self._ndim
    ndim = property(ndim, 
        doc='Dimensionality of the NumPy array storing atomic data.')
    def none(self):
        return self._none
    none = property(none, 
        doc='Set the value of the variable to None.')

ATOMIC_DATA_FIELDS = {
    'name':      Field('name',      'names',       '|S6',      'atom name',                      'AtomName'),
    'altloc':    Field('altloc',    'altlocs',     '|S1',      'alternate location indicator',   'AltLocIndicator'),
    'anisou':    Field('anisou',    'anisou',      np.float64, 'anisotropic temperature factor', 'AnisoTempFactor', ndim=2),
    'chain':     Field('chain',     'chids',       '|S1',      'chain identifier',               'ChainIdentifier', none='hv'),
    'element':   Field('element',   'elements',    '|S2',      'element symbol',                 'ElementSymbol'),
    'hetero':    Field('hetero',    'hetero',      np.bool,    'hetero flag',                    'HeteroFlag'),
    'occupancy': Field('occupancy', 'occupancies', np.float64, 'occupancy value',                'Occupancy',      meth_pl='Occupancies'),
    'resname':   Field('resname',   'resnames',    '|S6',      'residue name',                   'ResidueName'),
    'resnum':    Field('resnum',    'resnums',     np.int64,   'residue number',                 'ResidueNumber', none='hv'),
    'secondary': Field('secondary', 'secondary',   '|S1',      'secondary structure assignment', 'SecondaryStr'),
    'segment':   Field('segment',   'segments',    '|S6',      'segment name',                   'SegmentName'),
    'siguij':    Field('siguij',    'siguij',      np.float64, 'standard deviations for the anisotropic temperature factor',                   
                                                                                                 'AnisoStdDev', ndim=2),
    'serial':    Field('serial',    'serials',     np.int64,   'atom serial number (from file)', 'SerialNumber',  none='sn2i'),
    'beta':      Field('beta',      'bfactors',    np.float64, 'temperature (B) factor',         'TempFactor'),
    'icode':     Field('icode',     'icodes',      '|S1',      'insertion code',                 'InsertionCode', none='hv'),
    'type':      Field('type',      'types',       '|S6',      'atom type',                      'AtomType'),
    'charge':    Field('charge',    'charges',     np.float64, 'atomic partial charge',          'Charge'),
    'mass':      Field('mass',      'masses',      np.float64, 'atomic mass',                    'Mass', 'atomic masses', 'Masses'),
    'radius':    Field('radius',    'radii',       np.float64, 'atomic radius',                  'Radius', 'atomic radii', 'Radii'),
}

ATOMIC_ATTRIBUTES = {}
for field in ATOMIC_DATA_FIELDS.values():
    ATOMIC_ATTRIBUTES[field.var] = field

def wrapGetMethod(fn):
    def wrapped(self):
        return fn(self)
    return wrapped
def wrapSetMethod(fn):
    def wrapped(self, data):
        return fn(self, data)
    return wrapped

__doc__ += """

Common methods
-------------------------------------------------------------------------------

Atomic data contained in a PDB file can be accessed and changed using 
``get`` and ``set`` methods defined for :class:`Atomic` classes. To provide 
a coherent interface, these methods are defined for :class:`AtomGroup`, 
:class:`Atom`, :class:`Selection`, :class:`Chain`, :class:`Residue`, and 
:class:`AtomMap` classes, with the following exceptions: 

* Names of methods of the :class:`Atom` class are in singular form.
* ``set`` methods are not defined for the :class:`AtomMap` class.

The list of methods are below (they link to the documentation of the 
:class:`AtomGroup` methods):
 
===========================  ==================================================
Get/set method               Description
===========================  ==================================================
"""

keys = ATOMIC_DATA_FIELDS.keys()
keys.sort()

for key in keys:
    field = ATOMIC_DATA_FIELDS[key]
    __doc__ += '``get/set{0:18s}  get/set {1:s}\n'.format(field.meth_pl+'``', field.doc_pl)

__doc__ += """
===========================  ==================================================

.. note:: Note that ``get`` methods return a copy of the data. Changes in the 
   array obtained by calling one of the above methods will not be saved in the
   :class:`AtomGroup` instance. To change the data stored in :class:`AtomGroup`
   instance, use ``set`` methods.

Other functions common to all atomic classes is given below:

===========================  ==================================================
Method name                  Description
===========================  ==================================================
``copy``                     returns a deep copy atomic data.
``select``                   selects a subset of atoms (see :ref:`selections`).
``getCoordinates``           changes the index of the active coordinate set.
``getNumOfAtoms``            returns number of atoms.
``getNumOfCoordset``         returns number of coordinate sets.
``getCoordsets``             returns specified coordinate sets.
``getActiveCoordsetIndex``   returns the index of the active coordinate set.
``setActiveCoordsetIndex``   changes the index of the active coordinate set.
``iterCoordsets``            iterate over coordinate sets.
``isAttribute``              checks whether a user set attribute exists.
``getAttribute``             returns user set attribute data.
``setAttribute``             changes user set attribute data.
===========================  ==================================================



Special methods
-------------------------------------------------------------------------------

Atomic classes also have the following class specific methods: 
    
========================  =====================================================
Method                    Description
========================  =====================================================
:class:`AtomGroup`  
* ``getName``             returns name of the instance.
* ``setName``             changes name of the instance.
* ``delAttribute``        deletes a user set attribute from the instance.
* ``addCoordset``
* ``getNumOfChains``      returns the number of chains.
* ``iterChains``          iterates over chains.
* ``getNumOfResidues``    returns the total number of residues from all chains.
* ``iterResidues``        iterates over all residues.

                      
:class:`Atom`              
* ``getIndex``            returns atom index.
* ``getName``             return atom name.
* ``getSelectionString``  returns string that selects the atom.
                    
:class:`Selection`         
* ``getIndices``          returns indices of atoms.
* ``getSelectionString``  returns selection string of the instance.

:class:`Chain`
* ``getIdentifier``       returns chain identifier.
* ``setIdentifier``       changes chain identifier.
* ``getResidue``          returns residue with given number.
* ``iterResidues``        iterates over residues.
* ``getNumOfResidues``    returns the number of residues in the instance.
* ``getSequence``         returns single letter amino acid sequence. 
* ``getSelectionString``  returns a string that selects chain atoms.
                      
:class:`Residue`
* ``getIndices``          returns indices of atoms.
* ``getAtom``             returns :class:`Atom` with given name.
* ``getChain``            returns :class:`Chain` of residue instance.
* ``getChainIdentifier``  returns chain identifier.
* ``getInsertionCode``    returns insertion code.
* ``setInsertionCode``    changes insertion code.
* ``getName``             returns residue name.
* ``setName``             changes residue name.
* ``getNumber``           returns residue number.
* ``setNumber``           changes residue number.
* ``getSelectionString``  returns a string that selects residue atoms.

:class:`AtomMap`
* ``getIndices``          returns indices of atoms.
* ``getName``             returns name of the instance.
* ``setName``             changes name of the instance.
* ``getNumOfMapped``      returns number of mapped atoms.
* ``getNumOfUnmapped``    returns number of unmapped atoms.
* ``getMapping``          returns mapping of indices.
* ``getMappedFlags``      returns an boolean array indicating mapped atoms.
* ``getUnmappedFlags``    returns an boolean array indicating unmapped atoms.
========================  =====================================================

Functions common to :class:`Atom`, :class:`Selection`, :class:`Chain`,
:class:`Residue`, and :class:`AtomMap` include: 
    
========================  =====================================================
Method                    Description
========================  =====================================================  
* ``getAtomGroup``        returns the associated :class:`AtomGroup`.
* ``getIndices``          returns the indices of atoms.
========================  =====================================================


Behavioral differences
-------------------------------------------------------------------------------

Atomic classes behave differently to indexing and to calls of certain built-in 
functions. These differences are:

=========  ====================================================================
Class               Properties and differences
=========  ====================================================================
AtomGroup  * :func:`len` returns the number of atoms.
           * :func:`iter` yields :class:`Atom` instances.
           * Indexing by:
               
             - *atom index* (:func:`int`), e.g, ``10`` returns an 
               :class:`Atom`.
             - *slice* (:func:`slice`), e.g, ``10:20:2`` returns a 
               :class:`Selection`.
             - *chain identifier* (:func:`str`), e.g. ``"A"`` return 
               a :class:`Chain`.
             - *chain identifier, residue number [, insertion code]* 
               (:func:`tuple`), e.g. ``"A", 10`` or  ``"A", 10, "B"`` 
               returns a :class:`Residue`.
                       
Atom       * :func:`len` returns 1.
           * :func:`iter` is not applicable.
           * Indexing is not applicable.
                      
Selection  * :func:`len` returns the number of selected atoms.
           * :func:`iter` yields :class:`Atom` instances.
           * Indexing is not available.

Chain      * :func:`len` returns the number of residues in the chain.
           * :func:`iter` yields :class:`Residue` instances.
           * Indexing by:
                
             - *residue number [, insertion code]* (:func:`tuple`), 
               e.g. ``10`` or  ``10, "B"`` returns a :class:`Residue`.
                    
Residue    * :func:`len` returns the number of atoms in the instance.
           * :func:`iter` yields :class:`Atom` instances.
           * Indexing by:
              
             - *atom name* (:func:`str`), e.g. ``"CA"`` returns 
               an :class:`Atom`.

AtomMap    * :func:`len` returns the number of atoms in the instance.
           * :func:`iter` yields :class:`Atom` instances.
           * Indexing is not available.
=========  ====================================================================


Hierarchical views
-------------------------------------------------------------------------------

:class:`HierView` instances can be built for :class:`AtomGroup` and 
:class:`Selection` instances.

Some overridden functions are:

* :func:`len` return the number of chains.
* :func:`iter()` iterates over chains.
* Indexing:
    
  - *chain identifier* (:func:`str`), e.g. ``"A"`` returns a :class:`Chain`.
  - *chain identifier, residue number [, insertion code]* 
    (:func:`tuple`), e.g. ``"A", 10`` or  ``"A", 10, "B"`` 
    returns a :class:`Residue`


"""

__doc__ += """

:mod:`prody.atomic`
===============================================================================

Classes
-------

    * :class:`AtomGroup`
    * :class:`Atom`
    * :class:`Chain`
    * :class:`Residue`
    * :class:`Selection`
    * :class:`AtomMap`
    * :class:`HierView`
    
Base Classes
------------

    * :class:`Atomic`
    * :class:`AtomPointer`
    * :class:`AtomSubset`


Functions
---------

    * :func:`saveAtoms`
    * :func:`loadAtoms`

Inheritance Diagram
-------------------

.. inheritance-diagram:: prody.atomic
   :parts: 1
    
"""

class Atomic(object):
    """Base class for all atomic classes.
    
    Derived classes are:
        
      * :class:`AtomGroup`
      * :class:`AtomPointer`"""
      
    __slots__ = ('_acsi')
    
    def __contains__(self, item):
        """.. versionadded:: 0.5.3"""
        
        if isinstance(item, Atomic):
            if isinstance(item, AtomGroup) and self == item: 
                return True
            elif isinstance(self, AtomGroup) and self == item.getAtomGroup():
                return True
            elif len(item) <= len(self):
                if set(item.getIndices()).issubset(set(self.getIndices())):
                    return True
        return False        
      
    def __eq__(self, other):
        """
        .. versionadded:: 0.5.3
        
        .. versionchanged:: 0.8.1
           A :class:`Selection` (:class:`AtomPointer`) of all atoms is 
           considered not equal to the :class:`AtomGroup` anymore as 
           this causes problems in :mod:`select` module. 
        
        """
        
        if isinstance(other, Atomic):
            # AtomMaps may need special handling
            if self is other:
                return True
            elif isinstance(self, AtomPointer) and \
                isinstance(other, AtomPointer):
                self_indices = self._indices
                if len(self_indices) == len(other):
                    other_indices = other.getIndices()
                    if np.all(self_indices == other_indices):
                        return True
        return False
    
    def __ne__(self, other):
        """.. versionadded:: 0.5.3"""
        
        return not self.__eq__(other)
      
    def __getattribute__(self, name):
        """.. versionadded:: 0.7.1"""
        
        try:
            return object.__getattribute__(self, name)
        except AttributeError:
            selstr = name
            items = name.split('_')
            if prody.select.isKeyword(items[0]) or items[0] == 'not':
                selstr = ' '.join(items)
                return prody.ProDyAtomSelect.select(self, selstr)
        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))
        
    def getActiveCoordsetIndex(self):
        """Return index of the active coordinate set."""
        
        return self._acsi
    
    def select(self, selstr, **kwargs):
        """Return atoms matching the criteria in *selstr*.
        
        .. seealso:: :meth:`~prody.select.Select.select()` for more usage 
           details.
           
        """
        
        return prody.ProDyAtomSelect.select(self, selstr, **kwargs)

class AtomGroupMeta(type):
    def __init__(cls, name, bases, dict):
        for field in ATOMIC_DATA_FIELDS.values():
            
            # Define public method for retrieving a copy of data array
            def getData(self, var=field.var):
                array = self._data[var]
                if array is None:
                    return None
                return array.copy() 
            getData = wrapGetMethod(getData)
            getData.__name__ = field.meth_pl
            getData.__doc__ = 'Return a copy of {0:s} array.'.format(
                                                                field.doc_pl)
            setattr(cls, 'get'+field.meth_pl, getData)
            
            # Define private method for retrieving actual data array
            def _getData(self, var=field.var):
                return self._data[var]
            _getData = wrapGetMethod(_getData)
            _getData.__name__ = field.meth_pl
            _getData.__doc__ = ('Return {0:s} array.').format(field.doc_pl)
            setattr(cls, '_get'+field.meth_pl, _getData)
            
            # Define public method for setting values in data array
            if field.none:
                def setData(self, array, var=field.var, dtype=field.dtype, 
                            ndim=field.ndim, none=field.none):
                    if self._n_atoms == 0:
                        self._n_atoms = len(array)
                    elif len(array) != self._n_atoms:
                        raise ValueError('length of array must match n_atoms')
                        
                    if isinstance(array, list):
                        array = np.array(array, dtype)
                    elif not isinstance(array, np.ndarray):
                        raise TypeError('array must be a NumPy array or a list')
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
                    self.__setattr__('_'+none,  None)
            else:
                def setData(self, array, var=field.var, dtype=field.dtype, 
                            ndim=field.ndim):
                    if self._n_atoms == 0:
                        self._n_atoms = len(array)
                    elif len(array) != self._n_atoms:
                        raise ValueError('length of array must match n_atoms')
                        
                    if isinstance(array, list):
                        array = np.array(array, dtype)
                    elif not isinstance(array, np.ndarray):
                        raise TypeError('array must be a NumPy array or a list')
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
            setData = wrapSetMethod(setData)
            setData.__name__ = field.meth_pl 
            setData.__doc__ = 'Set {0:s}.'.format(field.doc_pl)  
            setattr(cls, 'set'+field.meth_pl, setData)


class AtomGroup(Atomic):
    
    """A class for storing and accessing atomic data.
    
    The number of atoms of the atom group is inferred at the first set method
    call from the size of the data array. 

    **Atomic Data**
    
    All atomic data is stored in :class:`numpy.ndarray` instances.

    **Get and Set Methods** 
    
    ``get()`` methods return copies of the data arrays. 
    
    ``set()`` methods accepts data in :func:`list` or :class:`~numpy.ndarray` 
    instances. The length of the list or array must match the number of atoms 
    in the atom group. Set method sets attributes of all atoms at once.
    
    Atom groups with multiple coordinate sets may have one of these sets as 
    the active coordinate set. The active coordinate set may be changed using
    :meth:`setActiveCoordsetIndex()` method. :meth:`getCoordinates` returns
    coordinates from the active set.
    
    To access and modify data associated with a subset of atoms in an atom 
    group, :class:`Selection` instances may be used. A selection from an atom 
    group has initially the same coordinate set as the active coordinate set.
    
    User can iterate over atoms and coordinate sets in an atom group. To 
    iterate over residues and chains, get a hierarchical view of the atom 
    group by calling :meth:`getHierView()`.
    
    """
    __metaclass__ = AtomGroupMeta
    
    __slots__ = ('_acsi', '_name', '_n_atoms', '_coordinates', '_n_csets', 
                 '_hv', '_userdata', '_sn2i',  
                 '_trajectory', '_frameindex', '_tcsi', '_timestamps',
                 '_data')
    
    def __init__(self, name='Unnamed'):
        """Instantiate an AtomGroup with a *name*."""
        self._name = str(name)
        self._n_atoms = 0
        self._coordinates = None
        self._acsi = None                  # Active Coordinate Set Index
        self._n_csets = 0
        self._hv = None
        self._userdata = {}
        self._sn2i = None
        self._trajectory = None
        self._frameindex = None
        self._tcsi = None # Trajectory Coordinate Set Index
        self._timestamps = None
        self._data = dict()
        
        for field in ATOMIC_DATA_FIELDS.values():
            self._data[field.var] = None

    def _getTimeStamp(self, index):
        """Return time stamp showing when coordinates were last changed."""

        if self._n_csets > 0:
            if index is None:
                return self._timestamps[self._acsi]
            else:
                return self._timestamps[index]
        else:
            return None
    
    def _setTimeStamp(self, index=None):
        """Set time stamp when:
           
            * :meth:`setCoordinates` method of :class:`AtomGroup` or 
              :class:`AtomPointer` instances are called.
              
            * one of :meth:`nextFrame`, :meth:`skipFrame`, or :meth:`gotoFrame` 
              methods is called.
        """
        
        if index is None:
            self._timestamps = np.zeros(self._n_csets)
            self._timestamps[:] = time.time()
        else:
            self._timestamps[index] = time.time()

    def __repr__(self):
        if self._trajectory is None:
            if self._n_csets > 0:
                return ('<AtomGroup: {0:s} ({1:d} atoms; {2:d} coordinate '
                        'sets, active set index: {3:d})>').format(self._name, 
                                    self._n_atoms, self._n_csets, self._acsi)
            else:
                return ('<AtomGroup: {0:s} ({1:d} atoms; {2:d} coordinate '
                        'sets)>').format(self._name,  self._n_atoms, 
                        self._n_csets)
        else:
            return ('<AtomGroup: {0:s} ({1:d} atoms; trajectory {2:s}, '
                    'frame index {3:d})>').format(self._name, 
                    self._n_atoms, self._trajectory.getName(), self._tcsi, 
                    len(self._trajectory))
        
    def __str__(self):
        return ('AtomGroup {0:s}').format(self._name)
        return ('{0:s} ({1:d} atoms; {2:d} coordinate sets, active '
               'set index: {3:d})').format(self._name, 
              self._n_atoms, self._n_csets, self._acsi)

    def __getitem__(self, indices):
        acsi = self._acsi
        if isinstance(indices, int):
            if indices < 0:
                indices = self._n_atoms + indices
            return Atom(self, indices, acsi)
        elif isinstance(indices, slice):
            start, stop, step = indices.indices(self._n_atoms)
            if start is None:
                start = 0
            if step is None:
                step = 1
            selstr = 'index {0:d}:{1:d}:{2:d}'.format(start, stop, step)
            return Selection(self, np.arange(start,stop,step), selstr, acsi)
        elif isinstance(indices, (list, np.ndarray)):
            return Selection(self, np.array(indices), 'Some atoms', 
                 'index {0:s}'.format(' '.join(np.array(indices, '|S'))), acsi)
        elif isinstance(indices, (str, tuple)):
            hv = self.getHierView()
            return hv[indices]
        else:
            raise TypeError('invalid index') 
    
    def __iter__(self):
        """Iterate over atoms in the atom group."""
        
        acsi = self._acsi
        for index in xrange(self._n_atoms):
            yield Atom(self, index, acsi)

    def __len__(self):
        return self._n_atoms
    
    def __add__(self, other):
        """.. versionadded:: 0.5"""
        
        if isinstance(other, AtomGroup):
            if self == other:
                raise ValueError('an atom group cannot be added to itself')
            
            new = AtomGroup(self._name + ' + ' + other._name)
            n_coordsets = self._n_csets
            if n_coordsets != other._n_csets:
                LOGGER.warning('AtomGroups {0:s} and {1:s} do not have same number'
                  ' of coordinate sets. First from both AtomGroups will be merged.'
                  .format(str(self._name), str(other._name), n_coordsets))
                n_coordsets = 1
            coordset_range = range(n_coordsets)
            new.setCoordinates(np.concatenate((self._coordinates[coordset_range],
                                        other._coordinates[coordset_range]), 1))
            for field in ATOMIC_DATA_FIELDS.values():
                var = field.var
                this = self._data[var]
                that = other._data[var]
                if this is not None and that is not None:
                    new._data[var] = np.concatenate((this, that))
            return new
        elif isinstance(other, prody.VectorBase):
            if self._n_atoms != other.getNumOfAtoms(): 
                raise ValueError('Vector/Mode must have same number of atoms '
                                 'as the AtomGroup')
            self.addCoordset(self._coordinates[self._acsi] + 
                             other._getArrayNx3())
            self.setActiveCoordsetIndex(self._n_csets - 1)
        else:
            raise TypeError('can only concatenate two AtomGroup`s or can '
                            'deform AtomGroup along a Vector/Mode')

    def _buildSN2I(self):
        """Builds a mapping from serial numbers to atom indices."""
        
        serials = self._serials  
        if serials is None:
            raise AttributeError('atom serial numbers are not set')
        if len(np.unique(serials)) != self._n_atoms:
            raise ValueError('atom serial numbers must be unique')
        if serials.min() < 0:
            raise ValueError('atoms must not have negative serial numbers')
        sn2i = -np.ones(serials.max() + 1)
        sn2i[serials] = np.arange(self._n_atoms)
        self._sn2i = sn2i

    def _getSN2I(self):
        if self._sn2i is None:
            self._buildSN2I()
        return self._sn2i

    def getName(self):
        """Return name of the atom group instance."""
        
        return self._name
    
    def setName(self, name):
        """Set name of the atom group instance."""
        
        self._name = str(name)
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
    
    def getCoordinates(self): 
        """Return a copy of coordinates from active coordinate set."""
        
        if self._coordinates is None:
            return None
        return self._coordinates[self._acsi].copy()
    
    def _getCoordinates(self): 
        """Return a view of coordinates from active coordinate set."""
        
        if self._coordinates is None:
            return None
        return self._coordinates[self._acsi]


    def setCoordinates(self, coordinates):
        """Set coordinates.  Coordinates must be a :class:`numpy.ndarray` 
        instance.  If the shape of the coordinates array is 
        (n_coordsets,n_atoms,3), the given array will replace all coordinate 
        sets. To avoid it, :meth:`addCoordset` may be used.
        
        If the shape of the coordinates array is (n_atoms,3) or (1,n_atoms,3),
        the coordinate set will replace the coordinates of the currently active
        coordinate set."""

        coordinates = checkCoordsArray(coordinates, 'coordinates',
                                       cset=True, n_atoms=self._n_atoms,
                                       reshape=True)
        if self._n_atoms == 0:
            self._n_atoms = coordinates.shape[-2] 
        
        if self._coordinates is None:
            self._coordinates = coordinates
            self._n_csets = coordinates.shape[0]
            self._acsi = 0
            self._setTimeStamp()
        else:
            if coordinates.shape[0] == 1:
                self._coordinates[self._acsi] = coordinates[0]
                self._setTimeStamp(self._acsi)
            else:
                self._coordinates = coordinates
                self._n_csets = coordinates.shape[0]
                self._acsi = min(self._n_csets - 1, self._acsi)
                self._setTimeStamp()
            
    def addCoordset(self, coords):
        """Add a coordinate set to the atom group.
        
        .. versionchanged:: 0.6.2
            :class:`~prody.ensemble.Ensemble` and :class:`Atomic` instances are 
            accepted as *coords* argument.
        """
        
        if self._trajectory is not None:
            raise AttributeError('AtomGroup is locked for coordinate set '
                                 'addition/deletion when its associated with '
                                 'a trajectory')
        if isinstance(coords, (prody.ensemble.Ensemble, Atomic)):
            if self._n_atoms != coords.getNumOfAtoms(): 
                raise ValueError('coords must have same number of atoms')
            coords = coords.getCoordsets()

        if self._coordinates is None:
            self.setCoordinates(coords)
            return

        coords = checkCoordsArray(coords, 'coords', cset=True, 
                                  n_atoms=self._n_atoms, reshape=True)

        self._coordinates = np.concatenate((self._coordinates, coords), axis=0)
        self._n_csets = self._coordinates.shape[0]
        timestamps = self._timestamps
        self._timestamps = np.zeros(self._n_csets)
        self._timestamps[:len(timestamps)] = timestamps
        self._timestamps[len(timestamps):] = time.time()
        
    def delCoordset(self, index):
        """Delete a coordinate set from the atom group."""
        
        if self._n_csets == 0:
            raise AttributeError('coordinates are not set')
        if self._trajectory is not None:
            raise AttributeError('AtomGroup is locked for coordinate set '
                                 'addition/deletion when its associated with '
                                 'a trajectory')
        which = np.ones(self._n_csets, np.bool)
        which[index] = False
        n_csets = self._n_csets
        if which.sum() == 0:
            self._coordinates = None
            self._n_csets = 0
            self._acsi = None
        else:
            self._coordinates = self._coordinates[which]
            self._n_csets = self._coordinates.shape[0]
            self._acsi = 0
        self._timestamps = self._timestamps[which]        

    def getCoordsets(self, indices=None):
        """Return a copy of coordinate sets at given indices.
        
        *indices* may be an integer, a list of integers or ``None``. ``None``
        returns all coordinate sets."""
        
        if self._coordinates is None:
            return None
        if indices is None:
            return self._coordinates.copy()
        elif isinstance(indices, (int, slice)):
            return self._coordinates[indices].copy()
        elif isinstance(indices, (list, np.ndarray)):
            return self._coordinates[indices]
        else:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')
        
    def _getCoordsets(self, indices=None):
        """Return a view of coordinate sets at given indices."""
        
        if self._coordinates is None:
            return None
        if indices is None:
            return self._coordinates
        elif isinstance(indices, (int, slice)):
            return self._coordinates[indices]
        elif isinstance(indices, (list, np.ndarray)):
            return self._coordinates[indices]
        else:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')

    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        
        return self._n_csets
    
    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate
        set."""
        
        for i in range(self._n_csets):
            yield self._coordinates[i].copy()
    
    def _iterCoordsets(self):
        """Iterate over coordinate sets by returning a view of each coordinate
        set."""
        
        for i in range(self._n_csets):
            yield self._coordinates[i]

    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set."""
        
        if self._n_csets == 0:
            self._acsi = 0
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._n_csets <= index or self._n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._n_csets 
        self._acsi = index

    def copy(self, which=None):
        """Return a copy of atoms indicated *which* as a new AtomGroup 
        instance.
        
        *which* may be:
            * ``None``, make a copy of the AtomGroup         else:
            * a Selection, Residue, Chain, or Atom instance
            * a list or an array of indices
            * a selection string
            
        .. versionchanged:: 0.7.1
           If selection string does not select any atoms, ``None`` is returned.
        
        .. versionchanged:: 0.8
           User data stored in the atom group is also copied.
        
        Note that association of an atom group with a trajectory is not copied.
        
        """
        
        name = self._name
        if which is None:
            indices = None
            newmol = AtomGroup('Copy of {0:s}'.format(name))
            newmol.setCoordinates(self._coordinates.copy())
            for field in ATOMIC_DATA_FIELDS.values():
                var = field.var
                array = self._data[var]
                if array is not None:
                    newmol._data[var] = array.copy()
        elif isinstance(which, int):
            indices = [which]
            newmol = AtomGroup('Copy of {0:s} index {1:d}'.format(name, which))
        elif isinstance(which, str):
            indices = prody.ProDyAtomSelect.getIndices(self, which)
            if len(indices) == 0:
                return None
            newmol = AtomGroup('Copy of {0:s} selection "{1:s}"'
                               .format(name, which))
        elif isinstance(which, (list, np.ndarray)):
            if isinstance(which, list):
                indices = np.array(which)
            elif which.ndim != 1:
                raise ValueError('which must be a 1d array')
            else:
                indices = which
            newmol = AtomGroup('Copy of a {0:s} subset'.format(name))
        else:
            if isinstance(which, Atom):
                indices = [which.getIndex()]
            elif isinstance(which, (AtomSubset, AtomMap)):
                indices = which.getIndices()
            else:
                raise TypeError('{0:s} is not a valid type'.format(type(which)))            
            newmol = AtomGroup('Copy of {0:s} selection "{1:s}"'.format(name, 
                                                                str(which)))
        if indices is not None:
            newmol.setCoordinates(self._coordinates[:, indices])
        for field in ATOMIC_DATA_FIELDS.values():
            var = field.var
            array = self._data[var]
            if array is not None:
                if indices is None:
                    newmol._data[var] = array.copy()
                else:
                    newmol._data[var] = array[indices]
        for name, data in self._userdata.iteritems():
            if indices is None:
                newmol._userdata[name] = data.copy()
            else:
                newmol._userdata[name] = data[indices]
        return newmol
    
    def getHierView(self):
        """Return a hierarchical view of the atom group."""
        hv = self._hv
        if hv is None:
            hv = HierView(self)
            self._hv = hv
        return hv
    
    def getNumOfChains(self):
        """Return number of chains.
        
        .. versionadded:: 0.7.1
        
        """
        
        return self.getHierView().getNumOfChains()
    
    def iterChains(self):
        """Iterate over chains.
        
        .. versionadded:: 0.7.1
        
        """
        
        return self.getHierView().iterChains()
    
    def getNumOfResidues(self):
        """Return number of chains.
        
        .. versionadded:: 0.7.1
        
        """
        
        return self.getHierView().getNumOfResidues()

    def iterResidues(self):
        """Iterate over residues.
        
        .. versionadded:: 0.7.1
        
        """
        
        return self.getHierView().iterResidues()

    def getAttrNames(self):
        """Return a list of user attribute names.
        
        .. versionadded:: 0.8"""
        
        return self._userdata.keys()
        

    def setAttribute(self, name, data):
        """Set a new attribute called *name* storing atomic *data*.
        
        .. versionadded:: 0.7.1
        
        *name* must:
            
            * start with a letter
            * contain alphanumeric characters and underscore
            * not be a selection keyword or one of the reserved names 
              listed below
        
        *data* must be a :func:`list` or a :class:`numpy.ndarray`, its length 
        must be equal to the number of atoms, and the type of data array must 
        be one of the following:
            
            * :class:`numpy.bool_`
            * :class:`numpy.float64`
            * :class:`numpy.int64`
            * :class:`numpy.string_`
        
        If a :func:`list` is given, its type must match one of the above 
        after it is converted to an array.  
        
        If the dimension of the *data* array is 1, name can be used
        as a selection string for making selections
        
        Note that, if an attribute with given *name* exists, it will be 
        overridden.

        :arg name: name of the attribute
        :type name: str
        
        :arg data: atomic data
        :type data: :class:`numpy.ndarray`
        
        """
        
        if not isinstance(name, str):
            raise TypeError('name must be a string')
        elif name == '':
            raise ValueError('name cannot be empty string')
        elif not name[0].isalpha():
            raise ValueError('name must start with a letter')
        elif not (''.join((''.join(name.split('_'))).split())).isalnum():
            raise ValueError('name may contain alphanumeric characters and '
                             'underscore, {0:s} is not valid'.format(name))
            
        if name in ATOMIC_DATA_FIELDS or prody.select.isReserved(name):
            raise ValueError('name cannot be a reserved name or a selection '
                             'keyword, "{0:s}" is invalid'.format(name))
        if isinstance(data, list):
            data = np.array(data)
        if not isinstance(data, np.ndarray):
            raise TypeError('data must be a numpy.ndarray instance')
        elif len(data) != self._n_atoms:
            raise ValueError('length of data array must match number of atoms')
        elif not data.dtype.type in (np.float64, np.int64, np.string_, np.bool_):
            raise TypeError('type of data array must be float64, int64, or '
                            'string_, {0:s} is not valid'.format(
                            str(data.dtype.type)))
            
        self._userdata[name] = data
    
    setAttribute.__doc__ += """Name of the attribute cannot be one of the 
        following reserved names:
            
          * {0:s}
    """.format('\n          * '.join(ATOMIC_DATA_FIELDS.keys()))
    
    def delAttribute(self, name):
        """Delete the attribute with given *name* and return the stored data.
        
        .. versionadded:: 0.7.1
        
        """
        
        if not isinstance(name, str):
            raise TypeError('name must be a string')
        return self._userdata.pop(name, None)
    
    def getAttribute(self, name):
        """Return a copy of the attribute *name* data array, if it exists.
        
        .. versionadded:: 0.7.1
        
        """
        
        data = self._userdata.get(name, None)
        if data is not None:
            return data.copy()
        return None


    def _getAttribute(self, name):
        """Return the attribute *name* data array, if it exists (for internal 
        use)."""
        
        data = self._userdata.get(name, None)
        if data is not None:
            return data
        return None


    def isAttribute(self, name):    
        """Return ``True`` if *name* is a user set attribute.
        
        .. versionadded:: 0.7.1
        
        """
        
        return name in self._userdata
  
    def getBySerial(self, serial):
        """Get an atom by its *serial* number.
        
        .. versionadded:: 0.8
        
        """

        sn2i = self._getSN2I()
        if not isinstance(serial, int):
            raise TypeError('serial must be an integer')
        if serial < 0:
            raise ValueError('serial must be greater than or equal to zero')
        if serial < len(sn2i):
            index = sn2i[serial]
            if index != -1:
                return Atom(self, index)

    def getBySerialRange(self, start, stop, step=None):
        """Get atoms by a *serial* number range.
        
        .. versionadded:: 0.8
        
        *start* must be zero or a positive integer and *stop* must be an 
        integer greater than *start*. ``getBySerialRange(i, j)`` will return
        atoms whose serial numbers are i+1, i+2, ..., j-1. Atom with serial
        number *stop* will be excluded as in indexing a Python :func:`list`.
        When *step* (default is 1) is given, it specified increment.  
        
        """
        
        sn2i = self._getSN2I()
        if not isinstance(start, int):
            raise TypeError('start must be an integer')
        if start < 0:
            raise ValueError('start must be greater than or equal to zero')
    
        if not isinstance(stop, int):
            raise TypeError('stop must be an integer')
        if stop <= start:
            raise ValueError('stop must be greater than start')
            
        if step is None:
            step = 1
        else:
            if not isinstance(step, int):
                raise TypeError('step must be an integer')
            if step < 1:
                raise ValueError('step must be greater than zero')
        
        indices = sn2i[start:stop:step]
        indices = indices[ indices > -1 ]
        return Selection(self, indices, 'serial {0:d}:{1:d}:{2:d}'
                                        .format(start, stop, step))

    def setTrajectory(self, trajectory):              
        """Associates atom group with a *trajectory*. *trajectory* may be
        a filename or a :class:`~prody.ensemble.Trajectory` instance.
        Number of atoms in the atom group and the trajectory must match. 
        At association a new coordinate set will be added to the atom group.
        :meth:`nextFrame`, :meth:`skipFrame`, and :meth:`gotoFrame` methods
        can be used to read coordinate sets from the trajectory.  To remove 
        association with a trajectory, pass ``None`` as trajectory argument.
        When atom group is associated with a trajectory, it will be locked
        for coordinate set addition/deletion operations.
        
        .. versionadded:: 0.8
        
        """
        
        if trajectory is None:
            self._tcsi = None
            self._trajectory = None
            self.delCoordset(self._acsi)
        else:
            if isinstance(trajectory, str):
                trajectory = prody.Trajectory(trajectory)
            elif not isinstance(trajectory, prody.TrajectoryBase):
                raise TypeError('trajectory must be a file name or a '
                                'TrajectoryBase instance')
            if self._n_atoms != trajectory.getNumOfAtoms():
                raise ValueError('trajectory must have same number of atoms')
                
            self._tcsi = trajectory.getNextFrameIndex()
            self.addCoordset(trajectory.nextCoordset())
            self._acsi = self._n_csets - 1
            self._trajectory = trajectory
        
    def getTrajectory(self):
        """Return trajectory associated with the atom group."""
        
        return self._trajectory
    
    def nextFrame(self):
        """Read the next frame from the trajectory and update coordinates.
        
        .. versionadded:: 0.8
        """
        nfi = self._trajectory.getNextFrameIndex()
        if nfi - self._tcsi == 1:
            self._tcsi = nfi
            self._coordinates[self._acsi] = self._trajectory.nextCoordset()
            self._setTimeStamp(self._acsi)
        else:
            self._gotoFrame(self._tcsi+1)
                
    def skipFrame(self, n=1): 
        """Read the frame after skipping *n* frames and update coordinates.
        
        .. versionadded:: 0.8
        """
        nfi = self._trajectory.getNextFrameIndex()
        self._trajectory.skip(n)
        if nfi - self._tcsi == 1:         
            self._tcsi = self._trajectory.getNextFrameIndex()
            self._coordinates[self._acsi] = self._trajectory.nextCoordset()
            self._setTimeStamp(self._acsi)
        else:
            self.gotoFrame(self._tcsi+1+n)
    
    def gotoFrame(self, n):
        """Read frame *n* from the trajectory and update coordinates.
        
        .. versionadded:: 0.8
        """
        
        self._trajectory.goto(n)
        self._tcsi = self._trajectory.getNextFrameIndex()
        self._coordinates[self._acsi] = self._trajectory.nextCoordset()
        self._setTimeStamp(self._acsi)
    
    def getFrameIndex(self):
        """Return current frame index.
        
        .. versionadded:: 0.8
        """
        
        return self._tcsi
    
    
class AtomPointer(Atomic):
    """Base class for classes pointing to atom(s) in :class:`AtomGroup` 
    instances.
    
    Derived classes are:
        
      * :class:`Atom`
      * :class:`AtomSubset`
      * :class:`AtomMap`
      
    """
    
    def __init__(self, atomgroup, acsi=None):
        if not isinstance(atomgroup, AtomGroup):
            raise TypeError('atomgroup must be AtomGroup, not {0:s}'
                            .format(type(atomgroup)))
        self._ag = atomgroup
        if acsi is None:
            self._acsi = atomgroup.getActiveCoordsetIndex()
        else: 
            self._acsi = int(acsi)

    def __add__(self, other):
        """Returns an :class:`AtomMap` instance. Order of pointed atoms are
        preserved.
        
        .. versionadded:: 0.5
        
        """
        
        if not isinstance(other, AtomPointer):
            raise TypeError('an AtomPointer instance cannot be added to a '
                            '{0:s} instance'.format(type(other)))
        ag = self._ag
        if ag != other._ag:
            raise ValueError('AtomPointer instances must point to same '
                             'AtomGroup instance')
        acsi = self._acsi
        if self._acsi != other._acsi:
            LOGGER.warning('Active coordinate set indices of operands are not '
                           'the same. Result will have {0:d}'.format(acsi))
        
        name = '({0:s}) + ({1:s})'.format(str(self), str(other))
        indices = np.concatenate([self.getIndices(), other.getIndices()])
        length = len(self)
        if isinstance(self, AtomMap):
            mapping = [self.getMapping()]
            unmapped = [self._unmapped]
        else:
            mapping = [np.arange(length)]
            unmapped = [np.array([])]
        
        if isinstance(other, AtomMap):
            mapping.append(other.getMapping() + length)
            unmapped.append(other._unmapped + length) 
        else:
            mapping.append(np.arange(length, length+len(other)))
            unmapped.append(np.array([]))
        return AtomMap(ag, indices, np.concatenate(mapping), 
                       np.concatenate(unmapped), name, acsi)
    
    def _getTimeStamp(self, index=None):
        
        if index is None:
            return self._ag._getTimeStamp(self._acsi)
        else:
            return self._ag._getTimeStamp(index)
    
    def isAttribute(self, name):    
        """Return ``True`` if *name* is a user set attribute.
        
        .. versionadded:: 0.7.1
        
        """
        
        return self._ag.isAttribute(name)

    
    def getAtomGroup(self):
        """Return associated atom group."""
        
        return self._ag
    
    def getNumOfCoordsets(self):
        """Return number of coordinate sets."""
        
        return self._ag._n_csets

    def setActiveCoordsetIndex(self, index):
        """Set the index of the active coordinate set."""
        
        if self._ag._coordinates is None:
            raise AttributeError('coordinates are not set')
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._ag._n_csets <= index or \
           self._ag._n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._ag._n_csets
        self._acsi = index
        
    def copy(self, selstr=None):
        """Make a copy of atoms."""
        
        if selstr is None:
            return self._ag.copy(self)
        elif isinstance(selstr, str):
            return self._ag.copy(self.select(selstr))
        raise TypeError('selstr must be a string')
        
    def nextFrame(self):
        """Read the next frame from the trajectory and update coordinates.
        
        .. versionadded:: 0.8
        """
        
        self._ag.nextFrame()
                
    def skipFrame(self, n=1): 
        """Read the frame after skipping *n* frames and update coordinates.
        
        .. versionadded:: 0.8
        """
        
        self._ag.skipFrame(n)
    
    def gotoFrame(self, n):
        """Read frame *n* from the trajectory and update coordinates.
        
        .. versionadded:: 0.8
        """
        
        self._ag.gotoFrame(n)
    
    def getFrameIndex(self):
        """Return current frame index.
        
        .. versionadded:: 0.8
        """
        
        return self._ag.getFrameIndex()
            

class AtomMeta(type):
    def __init__(cls, name, bases, dict):
        for field in ATOMIC_DATA_FIELDS.values():
            
            # Define public method for retrieving a copy of data array
            def getData(self, var=field.var):
                array = self._ag._data[var]
                if array is None:
                    return None
                return array[self._index] 
            getData = wrapGetMethod(getData)
            getData.__name__ = field.meth
            getData.__doc__ = 'Return {0:s} of the atom.'.format(field.doc)
            setattr(cls, 'get'+field.meth, getData)
            setattr(cls, '_get'+field.meth, getData)
            
            # Define public method for setting values in data array
            if field.none:
                def setData(self, value, var=field.var, none=field.none):
                    array = self._ag._data[var]
                    if array is None:
                        raise AttributeError('attribute of the AtomGroup is not '
                                             'set')
                    array[self._index] = value
                    self._ag.__setattr__('_'+none,  None)
            else:
                def setData(self, value, var=field.var):
                    array = self._ag._data[var]
                    if array is None:
                        raise AttributeError('attribute of the AtomGroup is not '
                                             'set')
                    array[self._index] = value
            setData = wrapSetMethod(setData)
            setData.__name__ = field.meth 
            setData.__doc__ = 'Set {0:s} of the atom.'.format(field.doc)  
            setattr(cls, 'set'+field.meth, setData)
        setattr(cls, 'getName', getattr(cls, 'getAtomName'))
        setattr(cls, 'setName', getattr(cls, 'setAtomName'))


class Atom(AtomPointer):
    
    """A class for accessing and manipulating attributes of an atom 
    in a :class:`AtomGroup` instance.
    
    """
    
    __metaclass__ = AtomMeta
    __slots__ = ('_ag', '_index', '_acsi')
    
    def __init__(self, atomgroup, index, acsi=None):
        AtomPointer.__init__(self, atomgroup, acsi)
        self._index = int(index)
        
    def __repr__(self):
        n_csets = self._ag.getNumOfCoordsets()
        if n_csets > 0:
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets, active set index: {4:d})>').format(
                    self.getAtomName(), self._ag.getName(), self._index,  
                    n_csets, self._acsi)
        else:
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets)>').format(self.getAtomName(), 
                    self._ag.getName(), self._index, n_csets)
                    
        sn = self.getSerialNumber()
        if sn is None: 
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets, active set index: {4:d})>').format(
                    self.getAtomName(), self._ag.getName(), self._index,  
                    self._ag.getNumOfCoordsets(), self._acsi)
        return ('<Atom: {0:s} from {1:s} (index {2:d}; sn {5:d}; {3:d} '
                'coordinate sets, active set index: {4:d})>').format(
                self.getAtomName(), self._ag.getName(), self._index,  
                self._ag.getNumOfCoordsets(), self._acsi, sn)

    def __str__(self):
        return ('Atom {0:s} (index {1:d})').format(self.getAtomName(), 
                                                   self._index)
        sn = self.getSerialNumber()
        if sn is None: 
            return ('Atom {0:s} (index {1:d})').format(self.getAtomName(), 
                                                       self._index)
        return ('Atom {0:s} (index {1:d}; sn {2:d})').format(
                                    self.getAtomName(), self._index, sn)

    def __len__(self):
        return 1
    
    def getIndex(self):
        """Return index of the atom."""
        
        return self._index
    
    def getAttribute(self, name):
        """Return the attribute *name*, if it exists.
        
        .. versionadded:: 0.7.1

        """
        
        if self._ag.isAttribute(name):
            return self._ag._userdata[name][self._index]
    
    _getAttribute = getAttribute
    
    def setAttribute(self, name, data):
        """Set data for the attribute *name*.
        
        .. versionadded:: 0.7.1
        
        :raise AttributeError: when attribute *name* does not exist in the
            the :class:`AtomGroup` instance.
        
        """
        
        if self._ag.isAttribute(name):
            self._ag._userdata[name][self._index] = data 
        else:
            raise AttributeError("AtomGroup '{0:s}' has no attribute '{1:s}'"
                                 .format(self._ag.getName(), name))

    def getIndices(self):
        """Return index of the atom in a :class:`numpy.ndarray`."""
        
        return np.array([self._index])
    
    def getCoordinates(self):
        """Return a copy of coordinates of the atom from the active coordinate 
        set."""
        
        if self._ag._coordinates is None:
            return None
        return self._ag._coordinates[self._acsi, self._index].copy()
    
    def _getCoordinates(self):
        """Return a view of coordinates of the atom from the active coordinate 
        set."""
        
        if self._ag._coordinates is None:
            return None
        return self._ag._coordinates[self._acsi, self._index]
    
    def setCoordinates(self, coordinates):
        """Set coordinates of the atom in the active coordinate set."""
        
        self._ag._coordinates[self._acsi, self._index] = coordinates
        self._ag._setTimeStamp(self._acsi)
        
    def getCoordsets(self, indices=None):
        """Return a copy of coordinate sets at given indices. 
        *indices* may be an integer or a list/array of integers."""
        
        if self._ag._coordinates is None:
            return None
        if indices is None:
            return self._ag._coordinates[:, self._index].copy()
        elif isinstance(indices, (int, slice)):
            return self._ag._coordinates[indices, self._index].copy()
        elif isinstance(indices, (list, np.ndarray)):
            return self._ag._coordinates[indices, self._index]
        else:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')
       
    def _getCoordsets(self, indices=None): 
        """Return a view of coordinate sets at given indices."""
        
        if self._ag._coordinates is None:
            return None
        if indices is None:
            return self._ag._coordinates[:, self._index]
        elif isinstance(indices, (int, slice)):
            return self._ag._coordinates[indices, self._index]
        elif isinstance(indices, (list, np.ndarray)):
            return self._ag._coordinates[indices, self._index]
        else:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')

    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each 
        coordinate set."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coordinates[i, self._index].copy()


    def _iterCoordsets(self):
        """Iterate over coordinate sets by returning a view of each coordinate
        set."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coordinates[i, self._index]

    def getSelectionString(self):
        """Return selection string that will select this atom."""
        
        return 'index {0:d}'.format(self._index)


class AtomSubsetMeta(type):
    def __init__(cls, name, bases, dict):
        for field in ATOMIC_DATA_FIELDS.values():
            
            # Define public method for retrieving a copy of data array
            def getData(self, var=field.var):
                array = self._ag._data[var]
                if array is None:
                    return None
                return array[self._indices] 
            getData = wrapGetMethod(getData)
            getData.__name__ = field.meth_pl
            getData.__doc__ = 'Return a copy of {0:s} of the atoms.'.format(
                                                                field.doc_pl)
            setattr(cls, 'get'+field.meth_pl, getData)
            setattr(cls, '_get'+field.meth_pl, getData)
            
            # Define public method for setting values in data array
            if field.none:
                def setData(self, value, var=field.var, none=field.none):
                    array = self._ag._data[var]
                    if array is None:
                        raise AttributeError('attribute of the AtomGroup is '
                                             'not set')
                    array[self._indices] = value
                    self._ag.__setattr__('_'+none,  None)
            else:
                def setData(self, value, var=field.var):
                    array = self._ag._data[var]
                    if array is None:
                        raise AttributeError('attribute of the AtomGroup is '
                                             'not set')
                    array[self._indices] = value
            setData = wrapSetMethod(setData)
            setData.__name__ = field.meth_pl 
            setData.__doc__ = 'Set {0:s} of the atoms.'.format(field.doc_pl)  
            setattr(cls, 'set'+field.meth_pl, setData)
        
class AtomSubset(AtomPointer):
    
    """A class for manipulating subset of atomic data in an :class:`AtomGroup`.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, and active coordinate set index for the atom group.
    
    """
    
    __metaclass__ = AtomSubsetMeta    
    __slots__ = ['_ag', '_indices', '_acsi']
    def __init__(self, atomgroup, indices, acsi=None):
        """Instantiate atom group base class. 
        
        :arg atomgroup: an atom group
        :type atomgroup: :class:`AtomGroup`
        
        :arg indices: list of indices of atoms in the subset
        :type indices: list of integers
        
        :arg acsi: active coordinate set index
        :type acsi: integer
        
        """
        
        AtomPointer.__init__(self, atomgroup, acsi)
        if not isinstance(indices, np.ndarray):
            indices = np.array(indices, np.int64)
        elif not indices.dtype == np.int64:
            indices = indices.astype(np.int64)
        self._indices = np.unique(indices)
    
    def __iter__(self):
        """Iterate over atoms."""
        acsi = self._acsi
        ag = self._ag 
        for index in self._indices:
            yield Atom(ag, index, acsi)
    
    def __len__(self):
        return len(self._indices)
    
    def __invert__(self):
        arange = range(self._ag.getNumOfAtoms())
        indices = list(self._indices)
        while indices:
            arange.pop(indices.pop())
        sel = Selection(self._ag, arange, "not ({0:s}) ".format(
                                self.getSelectionString()), self._acsi)        
        return sel
    
    def __or__(self, other):
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
        if self is other:
            return self
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
        if isinstance(other, Atom):
            other_indices = np.array([other._index])
        else:
            other_indices = other._indices
        indices = np.unique(np.concatenate((self._indices, other_indices)))
        return Selection(self._ag, indices, '({0:s}) or ({1:s})'.format(
                self.getSelectionString(), other.getSelectionString()), acsi)

    def __and__(self, other):
        """
        .. versionchanged:: 0.7.1
           If intersection of selections does not contain any atoms, ``None``
           is returned.
        
        """
        
        if not isinstance(other, AtomSubset):
            raise TypeError('other must be an AtomSubset')
        if self._ag != other._ag:
            raise ValueError('both selections must be from the same AtomGroup')
        if self is other:
            return self
        acsi = self._acsi
        if acsi != other._acsi:
            LOGGER.warning('active coordinate set indices do not match, '
                           'so it will be set to zero in the union.')
            acsi = 0
        indices = set(self._indices)
        if isinstance(other, Atom):
            other_indices = set([other._index])
        else:
            other_indices = set(other._indices)
        indices = indices.intersection(other_indices)
        if indices:
            indices = np.unique(indices)
            return Selection(self._ag, indices, '({0:s}) and ({1:s})'.format(
               self.getSelectionString(), other.getSelectionString()), acsi)
               
    def getAttribute(self, name):
        """Return a copy of the attribute *name*, if it exists.
        
        .. versionadded:: 0.7.1
        
        """
        
        if self._ag.isAttribute(name):
            return self._ag._userdata[name][self._indices].copy()
    
    _getAttribute = getAttribute
    
    def setAttribute(self, name, data):
        """Set data for the attribute *name*.
        
        .. versionadded:: 0.7.1
        
        :raise AttributeError: when attribute *name* does not exist in the
            the :class:`AtomGroup` instance.
        
        """
        
        if self._ag.isAttribute(name):
            self._ag._userdata[name][self._indices] = data 
        else:
            raise AttributeError("AtomGroup '{0:s}' has no attribute '{1:s}'"
                                 .format(self._ag.getName(), name))
    
    def getIndices(self):
        """Return a copy of the indices of atoms."""
        
        return self._indices.copy()
    
    def getNumOfAtoms(self):
        """Return number of atoms."""
        
        return self._indices.__len__()

    def getCoordinates(self):
        """Return a copy of coordinates from the active coordinate set."""
        
        if self._ag._coordinates is None:
            return None
        # Since this is not slicing, a view is not returned
        return self._ag._coordinates[self._acsi, self._indices]
    
    _getCoordinates = getCoordinates
    
    def setCoordinates(self, coordinates):
        """Set coordinates in the active coordinate set."""
        
        self._ag._coordinates[self._acsi, self._indices] = coordinates
        self._ag._setTimeStamp(self._acsi)
        
    def getCoordsets(self, indices=None):
        """Return coordinate sets at given *indices*.
        *indices* may be an integer or a list/array of integers."""
        
        if self._ag._coordinates is None:
            return None
        if indices is None:
            return self._ag._coordinates[:, self._indices]
        elif isinstance(indices, (int, slice)):
            return self._ag._coordinates[indices, self._indices]
        elif isinstance(indices, (list, np.ndarray)):
            return self._ag._coordinates[indices, self._indices]
        else:
            raise IndexError('indices must be an integer, a list/array of '
                             'integers, a slice, or None')
    _getCoordsets = getCoordsets

    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate 
        set.
        
        """
        
        for i in range(self._ag._n_csets):
            yield self._ag._coordinates[i, self._indices]

    _iterCoordsets = iterCoordsets

class Chain(AtomSubset):
    
    """Instances are generated by :class:`HierView` class.
    
    Indexing a :class:`Chain` instance by residue number returns 
    :class:`Residue` instances.
    
    >>> from prody import *
    >>> pdb = parsePDB('1p38')
    >>> hv = pdb.getHierView()
    >>> chA = hv['A']
    >>> chA[4]
    <Residue: GLU 4 from Chain A from 1p38 (9 atoms; 1 coordinate sets, active set index: 0)>
    >>> print chA[3] # Residue 3 does not exist in chain A
    None
    
    Iterating over a chain yields residue instances:
        
    >>> for res in chA: print res
    GLU 4
    ARG 5
    PRO 6
    THR 7
    ...
        
    """
        
    __slots__ = AtomSubset.__slots__ + ['_seq', '_dict']
    
    def __init__(self, atomgroup, indices, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._seq = None
        self._dict = dict()
        
    def __len__(self):
        return len(self._dict)
    
    def __repr__(self):
        n_csets = self._ag.getNumOfCoordsets()
        if n_csets > 0:
            return ('<Chain: {0:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets, active set index: {4:d})>').format(
                    self.getIdentifier(), self._ag.getName(), 
                    self.getNumOfAtoms(), n_csets, self._acsi)
        else:
            return ('<Chain: {0:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets)>').format(self.getIdentifier(), 
                    self._ag.getName(), self.getNumOfAtoms(), n_csets)

    def __str__(self):
        return ('Chain {0:s}').format(self.getIdentifier())

    def __iter__(self):
        return self.iterResidues()
    
    def __getitem__(self, number):
        """Returns the residue with given number, if it exists.
        
        .. versionchanged:: 6.2
           Tuples composed of chain identifier, residue number, and residue
           insertion code is accepted.
        
        """
        
        if isinstance(number, tuple): 
            if len(number) == 2:
                return self.getResidue(number[0], number[1]) 
            else:
                return self.getResidue(number[0])
        return self.getResidue(number)
    
    def getResidue(self, number, insertcode=''):
        """Return residue with given number."""
        
        return self._dict.get((number, insertcode), None)

    def iterResidues(self):
        """Iterate residues in the chain."""
        
        keys = self._dict.keys()
        keys.sort()
        for key in keys:
            yield self._dict[key]
    
    def getNumOfResidues(self):
        """Return number of residues."""
        
        return len(self._dict)

    def getIdentifier(self):
        """Return chain identifier."""
        
        return self._ag._data['chids'][self._indices[0]]
    
    def setIdentifier(self, identifier):
        """Set chain identifier."""
        
        self.setChainIdentifiers(identifier)
    
    def getSequence(self):
        """Return sequence, if chain is a polypeptide."""
        
        if self._seq:
            return self._seq
        CAs = self.select('name CA').select('protein')
        if len(CAs) > 0:
            self._seq = prody.compare.getSequence(CAs.getResidueNames())
        else:
            self._seq = ''
        return self._seq

    def getSelectionString(self):
        """Return selection string that selects this chain."""
        
        return 'chain {0:s}'.format(self.getIdentifier())


class Residue(AtomSubset):
    """Instances are generated by :class:`HierView` class.
    
    Indexing a :class:`Residue` by atom name returns :class:`Atom` instances.
    
    >>> from prody import *
    >>> pdb = parsePDB('1p38')
    >>> hv = pdb.getHierView()
    >>> chA = hv['A']
    >>> res = chA[4]
    >>> res['CA']
    <Atom: CA from 1p38 (index 1; 1 coordinate sets, active set index: 0)>
    >>> res['CB']
    <Atom: CB from 1p38 (index 4; 1 coordinate sets, active set index: 0)>
    >>> print res['H'] # X-ray structure 1p38 does not contain H atoms
    None
    
    """
     
    __slots__ = AtomSubset.__slots__ + ['_chain']
    
    def __init__(self, atomgroup, indices, chain, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._chain = chain

    def __repr__(self):
        n_csets = self._ag.getNumOfCoordsets()
        if n_csets > 0:
            return ('<Residue: {0:s} {1:d}{2:s} from Chain {3:s} from {4:s} '
                    '({5:d} atoms; {6:d} coordinate sets, active set index: '
                    '{7:d})>').format(self.getName(), self.getNumber(), 
                                      self.getInsertionCode(), 
                                      self.getChain().getIdentifier(), 
                                      self._ag.getName(), len(self), 
                                      n_csets, self._acsi)
        else:        
            return ('<Residue: {0:s} {1:d}{2:s} from Chain {3:s} from {4:s} '
                    '({5:d} atoms; {6:d} coordinate sets)>').format(
                        self.getName(), self.getNumber(), 
                        self.getInsertionCode(), 
                        self.getChain().getIdentifier(), 
                        self._ag.getName(), len(self), n_csets)
            
    def __str__(self):
        return '{0:s} {1:d}{2:s}'.format(self.getName(), self.getNumber(), 
                                         self.getInsertionCode())

    def __getitem__(self, name):
        return self.getAtom(name)
    
    def getAtom(self, name):
        """Return atom with given *name*, ``None`` if not found.
        
        Assumes that atom names in a residue are unique. If more than one atoms 
        with the given *name* exists, the one with the smaller index will be 
        returned.
        
        """
        if isinstance(name, str):
            nz = (self.getAtomNames() == name).nonzero()[0]
            if len(nz) > 0:
                return Atom(self._ag, self._indices[nz[0]], self._acsi)
    
    def getChain(self):
        """Return the chain that the residue belongs to."""
        
        return self._chain
    
    def getNumber(self):
        """Return residue number."""
        
        return int(self._ag._data['resnums'][self._indices[0]])
    
    def setNumber(self, number):
        """Set residue number."""
        
        self.setResidueNumbers(number)
    
    def getName(self):
        """Return residue name."""
        
        return self._ag._data['resnames'][self._indices[0]]
    
    def setName(self, name):
        """Set residue name."""
        
        self.setResidueNames(name)

    def getInsertionCode(self):
        """Return residue insertion code."""
        
        return self._ag._data['icodes'][self._indices[0]]
        
    def setInsertionCode(self, icode):
        """Set residue insertion code."""
        
        self.setInsertionCodes(icode)
    
    def getChainIdentifier(self):
        """Return chain identifier."""
        
        return self._chain.getIdentifier()
    
    def getSelectionString(self):
        """Return selection string that will select this residue."""
        
        return 'chain {0:s} and resnum {1:d}{2:s}'.format(
          self.getChainIdentifier(), self.getNumber(), self.getInsertionCode())


class Selection(AtomSubset):
    """A class for accessing and manipulating attributes of select of atoms 
    in an :class:`AtomGroup` instance.
    
    """
    
    __slots__ = AtomSubset.__slots__ + ['_selstr']
    
    def __init__(self, atomgroup, indices, selstr, acsi=None):
        AtomSubset.__init__(self, atomgroup, indices, acsi)
        self._selstr = str(selstr)
        
    def __repr__(self):
        n_csets = self._ag.getNumOfCoordsets()
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        if n_csets > 0:
            return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets, active set index: {4:d})>').format(
                    selstr, self._ag.getName(), len(self), n_csets, self._acsi)
        else:
            return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets)>').format(
                    selstr, self._ag.getName(), len(self), n_csets)

    def __str__(self):
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        return 'Selection "{0:s}" from {1:s}'.format(selstr, self._ag.getName())
    
    def getSelectionString(self):
        """Return selection string that selects this atom subset."""
        
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom subset."""
        
        return HierView(self)


class AtomMapMeta(type):
    
    def __init__(cls, name, bases, dict):

        for field in ATOMIC_DATA_FIELDS.values():
            def getData(self, name=field.name, var=field.var):
                array = self._ag._data[var]
                if array is None:
                    return None
                data = self._ag._data[var][self._indices]
                result = np.zeros((self._len,) + data.shape[1:], 
                                 ATOMIC_DATA_FIELDS[name].dtype)
                result[self._mapping] = data
                return result 
            getData = wrapGetMethod(getData)
            getData.__name__ = field.meth_pl
            getData.__doc__ = ('Return {0:s} of the atoms. Unmapped atoms '
                               'will have ``0`` or ``""`` as entries.'
                               .format(field.doc_pl))
            setattr(cls, 'get'+field.meth_pl, getData)
            setattr(cls, '_get'+field.meth_pl, getData)


class AtomMap(AtomPointer):
    """A class for mapping atomic data.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, active coordinate set index, mapping for indices, and
    indices of unmapped atoms.
    
    """
    
    __metaclass__ = AtomMapMeta
    __slots__ = ['_ag', '_indices', '_acsi', '_name', '_mapping', '_unmapped', 
                 '_len']
    
    def __init__(self, atomgroup, indices, mapping, unmapped, name='Unnamed', 
                 acsi=None):
        """Instantiate with an AtomMap with following arguments:        
        
        :arg atomgroup: the atomgroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms as a list of indices
        :arg unmapped: list of indices for unmapped atoms
        :arg name: name of the AtomMap instance
        :arg acsi: active coordinate set index, if ``None`` defaults to that 
            of *atomgrup*
        
        Length of *mapping* must be equal to length of *indices*. Number of 
        atoms (including unmapped dummy atoms) are determined from the 
        sum of lengths of *mapping* and *unmapped* arrays.
                 
        """
        
        AtomPointer.__init__(self, atomgroup, acsi)
        
        if not isinstance(indices, np.ndarray):
            self._indices = np.array(indices, np.int64)
        elif not indices.dtype == np.int64:
            self._indices = indices.astype(np.int64)
        else:
            self._indices = indices

        if not isinstance(mapping, np.ndarray):
            self._mapping = np.array(mapping, np.int64)
        elif not mapping.dtype == np.int64:
            self._mapping = mapping.astype(np.int64)
        else:
            self._mapping = mapping

        if not isinstance(unmapped, np.ndarray):
            self._unmapped = np.array(unmapped, np.int64)
        elif not unmapped.dtype == np.int64:
            self._unmapped = unmapped.astype(np.int64)
        else:
            self._unmapped = unmapped
        
        self._name = str(name)
        self._len = len(self._unmapped) + len(self._mapping)
        
    def __repr__(self):
        n_csets = self._ag.getNumOfCoordsets()
        if n_csets > 0:
            return ('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                    '{3:d} mapped; {4:d} unmapped; {5:d} coordinate sets, '
                    'active set index: {6:d})>').format(self._name,
                    self._ag.getName(), self._len, len(self._mapping), 
                    len(self._unmapped), n_csets, self._acsi)
        else:
            return (('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                    '{3:d} mapped; {4:d} unmapped; {5:d} coordinate sets)>')
                    .format(self._name, self._ag.getName(), self._len, 
                    len(self._mapping), len(self._unmapped), n_csets))
            
    def __str__(self):
        return 'AtomMap {0:s}'.format(self._name)
    
    def __iter__(self):
        indices = np.zeros(self._len, np.int64)
        indices[self._unmapped] = -1
        indices[self._mapping] = self._indices
        ag = self._ag
        acsi = self._acsi
        for index in indices:
            if index > -1:
                yield Atom(ag, index, acsi)
            else:
                yield None
    
    def __len__(self):
        return self._len
    
    def getAttribute(self, name):
        """Return a copy of the attribute *name*, if it exists.
        
        .. versionadded:: 0.7.1
        
        """
        
        if self._ag.isAttribute(name):
            data = self._ag._userdata[name][self._indices]
            result = np.zeros((self._len,) + data.shape[1:], data.dtype)
            result[self._mapping] = data
            return result

    _getAttribute = getAttribute

    def getName(self):
        """Return name of the atom map instance."""
        
        return self._name
    
    def setName(self, name):
        """Set name of the atom map instance."""
        
        self._name = str(name)

    def getNumOfAtoms(self):
        """Return number of mapped atoms."""
        
        return self._len

    def getNumOfUnmapped(self):
        """Return number of unmapped atoms."""
        
        return len(self._unmapped)

    def getNumOfMapped(self):
        """Return number of mapped atoms."""
        
        return len(self._mapping)

    def getIndices(self):
        """Return indices of mapped atoms."""
        
        return self._indices.copy()

    def getMapping(self):
        """Return mapping of indices."""
        
        return self._mapping.copy()

    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each coordinate 
        set.
        
        """
        
        for i in range(self._ag._n_csets):
            coordinates = np.zeros((self._len, 3), np.float64)
            coordinates[self._mapping] = self._ag._coordinates[i, 
                                                               self._indices] 
            yield coordinates
    
    _iterCoordsets = iterCoordsets

    def getCoordinates(self):
        """Return coordinates from the active coordinate set."""
        
        if self._ag._coordinates is None:
            return None
        coordinates = np.zeros((self._len, 3), np.float64)
        coordinates[self._mapping] = self._ag._coordinates[self._acsi, 
                                                           self._indices] 
        return coordinates
    
    _getCoordinates = getCoordinates
    
    def setCoordinates(self, coordinates):
        """Set coordinates in the active coordinate set.
        
        Length of the *coordinates* array must match the number of mapped 
        atoms.
        
        """
        
        self._ag._coordinates[self._acsi, self._indices] = coordinates
    

    def getCoordsets(self, indices=None):
        """Return coordinate sets at given indices.
        *indices* may be an integer or a list/array of integers."""
        
        if self._ag._coordinates is None:
            return None
        if indices is None:
            indices = np.arange(self._ag._n_csets)
        elif isinstance(indices, (int, long)):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(indices.indices(self._ag._n_csets))
        try:
            coordsets = np.zeros((len(indices), self._len, 3))
            coordsets[:, self._mapping] = self._ag._coordinates[indices][:, 
                                                                self._indices]  
            return coordsets
        except IndexError:
            raise IndexError('indices may be an integer or a list/array '
                             'of integers')

    _getCoordsets = getCoordsets

    def getUnmappedFlags(self):
        """Return an array with 1s for unmapped atoms."""
        flags = np.zeros(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 1
        return flags
    
    def getMappedFlags(self):
        """Return an array with 1s for mapped atoms."""
        flags = np.ones(self._len)
        if len(self._unmapped):
            flags[self._unmapped] = 0
        return flags

class HierView(object):
    
    """Hierarchical views can be generated for :class:`AtomGroup` 
    and :class:`Selection` instances.
    
    Indexing a :class:`HierView` instance returns a :class:`Chain` instance.
    
    >>> from prody import *
    >>> pdb = parsePDB('1p38')
    >>> hv = pdb.getHierView()
    >>> chA = hv['A']
    >>> chA
    <Chain: A from 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>
    >>> print hv['B'] # Chain B does not exist in 1p38
    None
    
    """
    
    __slots__ = ['_atoms', '_chains']

    def __init__(self, atoms):
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be an atomic instance')
        self._atoms = atoms
        self._chains = dict()
        self.update()

    def getAtoms(self):
        """Return atoms for which the hierarchical view is built.
        
        .. versionadded:: 0.6.2
        
        """
        
        return self._atoms
    
    def build(self):
        """Calls :meth:`update` method. This method will is deprecated and 
        will be removed in v0.8."""
        
        LOGGER.warning('HierView.build() method is deprecated. '
                       'Use HierView.update() method instead.')
        self.update()
    
    def update(self):
        """Rebuild hierarchical view of atoms.
        
        This method is called at instantiation, but can be used to rebuild
        the hierarchical view when attributes of atoms change.
        
        """
        what = 'built'
        if self._chains:
            what = 'updated'
        #start = time.time()
        acsi = self._atoms.getActiveCoordsetIndex()
        atoms = self._atoms
        if isinstance(atoms, AtomGroup):
            atomgroup = atoms
            _indices = np.arange(atomgroup._n_atoms)
            chainids = atomgroup.getChainIdentifiers() 
            if chainids is None:
                chainids = np.zeros(atomgroup._n_atoms, 
                                    dtype=ATOMIC_DATA_FIELDS['chain'].dtype)
                atomgroup.setChainIdentifiers(chainids)
        else:
            atomgroup = atoms._ag
            _indices = atoms._indices
            chainids = atomgroup.getChainIdentifiers() 
            if chainids is None:
                chainids = np.zeros(atomgroup._n_atoms, 
                                    dtype=ATOMIC_DATA_FIELDS['chain'].dtype)
                atomgroup.setChainIdentifiers(chainids)
            chainids = chainids[_indices]


        for chid in np.unique(chainids):
            ch = Chain(atomgroup, _indices[chainids == chid], acsi)
            self._chains[chid] = ch
        
        if atomgroup.getResidueNumbers() is None:
            atomgroup.setResidueNumbers(np.zeros(atomgroup._n_atoms, 
                                     dtype=ATOMIC_DATA_FIELDS['resnum'].dtype))
        if atomgroup.getResidueNames() is None:
            atomgroup.setResidueNames(np.zeros(atomgroup._n_atoms, 
                                    dtype=ATOMIC_DATA_FIELDS['resname'].dtype))
        if atomgroup.getInsertionCodes() is None:
            atomgroup.setInsertionCodes(np.zeros(atomgroup._n_atoms, 
                                    dtype=ATOMIC_DATA_FIELDS['icode'].dtype))

        icodes = atomgroup.getInsertionCodes()

        for chain in self.iterChains():
            chid = chain.getIdentifier()
            rd = defaultdict(list)
            indices = chain.getIndices()
            resnums = chain.getResidueNumbers()
            for i in xrange(len(resnums)):
                rd[resnums[i]].append(indices[i])
            resnums = rd.keys()
            resnums.sort()
            for resnum in resnums:
                resindices = np.array(rd[resnum])
                res_icodes = icodes[resindices]
                
                for ic in np.unique(res_icodes): 
                    subindices = resindices[res_icodes == ic]
                    temp = subindices[0]
                    res = Residue(atomgroup, subindices, chain, acsi)   
                    chain._dict[(resnum, ic)] = res
        #LOGGER.debug('Hierarchical view was {0:s} in {1:.2f}s.'
        #             .format(what, time.time()-start))
        
    def __repr__(self):
        return '<HierView: {0:s}>'.format(str(self._atoms))
    
    def __str__(self):
        return 'HierView of {0:s}'.format(str(self._atoms))
    
    def __iter__(self):
        """Iterate over chains."""
        return self.iterChains()
    
    def __len__(self):
        return len(self._chains)
    
    def __getitem__(self, chid):
        """
        .. versionchanged:: 6.2
           Tuples composed of chain identifier, residue number, and residue
           insertion code is accepted.
        
        """
        
        if isinstance(chid, str):
            return self._chains.get(chid, None)
        elif isinstance(chid, tuple):
            ch = self._chains.get(chid[0], None)
            if ch is not None:
                return ch[chid[1:]]

    def iterResidues(self):
        """Iterate over residues."""
        
        chids = self._chains.keys()
        chids.sort()
        for chid in chids:
            chain = self._chains[chid]
            for res in chain.iterResidues():
                yield res
                
    def getResidue(self, chid, resnum, icode=''):
        """Return residue with number *resnum* and insertion code *icode* from 
        the chain with identifier *chid*, if it exists.
        
        """
        
        ch = self._chains.get(chid, None)
        if ch is not None:
            return ch.getResidue(resnum, icode)
        return None

    def getNumOfResidues(self):
        """Returns number of residues."""
        
        return sum([ch.getNumOfResidues() for ch in self._chains.itervalues()])    

    def iterChains(self):
        """Iterate over chains."""
        
        chids = self._chains.keys()
        chids.sort()
        for chid in chids:
            yield self._chains[chid]
    
    def getChain(self, chid):
        """Return chain with identifier *chid*, if it exists."""
        
        return self._chains.get(chid, None)

    def getNumOfChains(self):
        """Return number of chains."""
        
        return len(self._chains)


def saveAtoms(atoms, filename=None):
    """Save *atoms* in ProDy internal format.
    
    All classes derived from :class:`Atomic` are accepted as *atoms* argument.
    
    .. versionadded:: 0.7.1
    
    This function saves custom atomic attributes as well.
    
    Note that name of the :class:`AtomGroup` instance is used as the filename
    when *atoms* is not an :class:`AtomGroup`. This is because names for
    selections and atom maps may be too long and may contain special 
    characters. To avoid overwriting an existing file with the same name, 
    specify a *filename*.
    
    """
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if isinstance(atoms, AtomGroup):
        ag = atoms
        name = ag.getName()
    else:
        ag = atoms.getAtomGroup()
        name = str(atoms)
    singular = False
    if isinstance(atoms, Atom):
        singular = True
    
    
    if filename is None:
        filename = ag.getName().replace(' ', '_')
    filename += '.ag.npz'
    attr_dict = {'_name': name}
    attr_dict['_coordinates'] = atoms._getCoordsets()
    for name, field in ATOMIC_DATA_FIELDS.items():
        if singular:
            data = atoms.__getattribute__('_get'+field.meth)()
        else:
            data = atoms.__getattribute__('_get'+field.meth_pl)()
        if data is not None:
            attr_dict[field.var] = data 
    for attr in ag._userdata.keys():
        attr_dict[attr] = atoms._getAttribute(attr)
    np.savez(filename, **attr_dict)
    return filename


def loadAtoms(filename):
    """Return :class:`AtomGroup` instance loaded from *filename*.
    
    .. versionadded:: 0.7.1
    
    .. seealso: :func:`saveAtoms`
    
    This function makes use of :func:`numpy.load` function.
    
    """
    
    start = time.time()
    attr_dict = np.load(filename)
    if not '_coordinates' in attr_dict.files:
        raise ValueError("'{0:s}' is not a valid atomic data file"
                         .format(filename))
    ag = AtomGroup(str(attr_dict['_name']))
    for attr in attr_dict.files:
        if attr == '_name':
            continue  
        elif attr == '_coordinates':
            data = attr_dict[attr]
            if data.ndim > 0:
                ag.setCoordinates(data)
        elif attr in ATOMIC_ATTRIBUTES: 
            field = ATOMIC_ATTRIBUTES[attr]
            data = attr_dict[attr]
            if data.ndim > 0:
               ag.__getattribute__('set' + field.meth_pl)(data)
            else:
                ag.__getattribute__('set' + field.meth_pl)([data])
        else:            
            data = attr_dict[attr]
            if data.ndim > 0:
                ag.setAttribute(attr, data)
            else:
                ag.setAttribute(attr, [data])
            
    LOGGER.debug('Atoms were loaded in {0:.2f}s.'.format(time.time() - start))
    return ag

if __name__ == '__main__':
    from prody import *
    p = parsePDB('1aar')
    saveAtoms(p)
