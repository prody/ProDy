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

"""This module defines classes for handling atomic data.

.. _atomic:
    
Atomic data
===============================================================================

ProDy stores atomic data in instances of :class:`AtomGroup` class.  This
class is designed to be efficient and responsive, i.e. facilitates user
to access atomic data quickly for any subset of atoms.  An :class:`AtomGroup`
instance can be obtained by parsing a PDB file as follows: 
    
>>> from prody import *
>>> ag = parsePDB('1aar')

To read this page in a Python session, type::
    
  help(atomic)

:class:`AtomGroup` instances can store multiple coordinate sets, which may
be models from an NMR structure, snapshots from an MD simulation.


ProDy stores all atomic data in :class:`AtomGroup` instances and comes
with other classes acting as pointers to provide convenient read/write access 
to such data.  These classes are:

* :class:`Atom` - Points to a single atom in an :class:`AtomGroup` instance.                          

* :class:`Selection` - Points to an arbitrary subset of atoms. See 
  :ref:`selections` and :ref:`selection-operations` for usage examples.

* :class:`Segment` - Points to atoms that have the same segment name.

* :class:`Chain` - Points to atoms in a segment that have the same chain 
  identifier.

* :class:`Residue` - Points to atoms in a chain that have the same residue 
  number and insertion code.
                      
* :class:`AtomMap` - Points to arbitrary subsets of atoms while allowing for 
  duplicates and missing atoms.  Indices of atoms are stored in the order 
  provided by the user.
    
Atom selections
-------------------------------------------------------------------------------

Flexible and powerful atom selections is one of the most important features 
of ProDy.  The details of the selection grammar is described in 
:ref:`selections`. 

.. versionadded:: 0.7.1

Using the flexibility of Python, atom selections are made much easier by
overriding the ``.`` operator i.e. the :meth:`__getattribute__` 
method of :class:`Atomic` class.  So the following will be interpreted
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
previous examples.  The limitation of this is that parentheses, special 
characters cannot be used.     

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from collections import defaultdict
from types import NoneType
import sys
import time

if sys.version_info[:2] < (2,7):
    from ordereddict import OrderedDict
else:
    from collections import OrderedDict

import numpy as np

from tools import *
import prody
LOGGER = prody.LOGGER

__all__ = ['Atomic', 'AtomGroup', 'AtomPointer', 'Atom', 'AtomSubset', 
           'Selection', 'Segment', 'Chain', 'Residue', 
           'AtomMap', 'HierView', 'Bond', 'ATOMIC_DATA_FIELDS',
           'loadAtoms', 'saveAtoms',]

READONLY = set(['numbonds', 'resindex', 'chindex', 'segindex'])

class Field(object):
    __slots__ = ['_name', '_var', '_dtype',  '_doc', '_doc_pl', 
                 '_meth', '_meth_pl', '_ndim', '_none', '_selstr',
                 '_depr', '_depr_pl', '_synonym', '_readonly', '_call']
    def __init__(self, name, dtype, **kwargs):
        self._name = name
        self._dtype = dtype
        self._var = kwargs.get('var', name + 's')
        self._doc = kwargs.get('doc', name)
        self._ndim = kwargs.get('ndim', 1)
        self._meth = kwargs.get('meth', name.capitalize())
        self._doc_pl = kwargs.get('doc_pl', self._doc + 's')
        self._meth_pl = kwargs.get('meth_pl', self._meth + 's')
        self._none = kwargs.get('none')
        self._selstr = kwargs.get('selstr')
        self._depr = kwargs.get('depr')
        if self._depr is None:
            self._depr_pl = None
        else:
            self._depr_pl = kwargs.get('depr_pl', self._depr + 's')
        self._synonym = kwargs.get('synonym')
        self._readonly = kwargs.get('readonly', False)
        self._call = kwargs.get('call', None)
        
    def name(self):
        return self._name
    name = property(name, doc='Data field name used in atom selections.')
    def var(self):
        return self._var
    var = property(var, doc='Internal variable name.')
    def dtype(self):
        return self._dtype
    dtype = property(dtype, doc='Data type (primitive Python types).')
    def doc(self):
        return self._doc
    doc = property(doc, doc='Data field name, as used in documentation.')
    def doc_pl(self):
        return self._doc_pl
    doc_pl = property(doc_pl, doc='Plural form for documentation.')
    def meth(self):
        return self._meth
    meth = property(meth, doc='Atomic get/set method name.')
    def meth_pl(self):
        return self._meth_pl
    meth_pl = property(meth_pl, doc='get/set method name in plural form.')
    def ndim(self):
        return self._ndim
    ndim = property(ndim, doc='Expected number of data array dimensions.')
    def none(self):
        return self._none
    none = property(none, doc='When to set the value of the variable to None.')
    def selstr(self):
        return self._selstr
    selstr = property(selstr, doc='Selection string examples.')
    def synonym(self):
        return self._synonym
    synonym = property(synonym, doc='Synonym used in atom selections.')    
    def readonly(self):
        return self._readonly
    readonly = property(readonly, 
                        doc='Read-only attribute without a set method.')    
    def call(self):
        return self._call
    call = property(call, doc='List of AtomGroup methods to call.')    
    def depr(self):
        return self._depr
    depr = property(depr, doc='Deprecated method name.')
    def depr_pl(self):
        return self._depr_pl
    depr_pl = property(depr_pl, doc='Deprecated method name in plural form.')
    def getDocstr(self, meth, plural=True, selex=True):
        """Return documentation string for the field."""
        
        assert meth in ('set', 'get', '_get'), "meth must be 'set' or 'get'"
        assert isinstance(plural, bool), 'plural must be a boolean'
        assert isinstance(selex, bool), 'selex must be a boolean'
        
        if meth == 'get':
            if plural:
                docstr = 'Return a copy of {0:s}.'.format(self.doc_pl)
            else:
                docstr = 'Return {0:s} of the atom.'.format(self.doc)
        elif meth == 'set':
            if plural:
                docstr = 'Set {0:s}.'.format(self.doc_pl)
            else:
                docstr = 'Set {0:s} of the atom.'.format(self.doc)
        else:
            selex = False
            if plural:
                docstr = 'Return {0:s} array.'.format(self.doc_pl) 
            
        selstr = self.selstr
        if selex and selstr:
            if plural:
                doc = self.doc_pl
            else:
                doc = self.doc
            if '(' in doc:
                doc = doc[:doc.index('(')]
            selex = "'``, ``'".join(selstr)
            selex = ("  {0:s} can be used in atom selections, e.g. "
                     "``'{1:s}'``.").format(doc.capitalize(), selex)
            if self.synonym is not None:
                selex = selex + ('  Note that *{0:s}* is a synonym for '
                    '*{1:s}*.').format(self.synonym, self.name)
            return docstr + selex
        else:
            return docstr

ATOMIC_DATA_FIELDS = {
    'name':      Field('name', '|S6', selstr=('name CA CB',), depr='AtomName'),
    'altloc':    Field('altloc', '|S1', doc='alternate location indicator', 
                       selstr=('altloc A B', 'altloc _'), 
                       depr='AltLocIndicator'),
    'anisou':    Field('anisou', float, doc='anisotropic temperature factor', 
                       ndim=2, depr='AnisoTempFactor'),
    'chain':     Field('chain', '|S1', var='chids', doc='chain identifier', 
                       meth='Chid', none='hv', synonym='chid', 
                       selstr=('chain A', 'chid A B C', 'chain _'), 
                       depr='ChainIdentifier'),
    'element':   Field('element', '|S2', doc='element symbol', 
                       selstr=('element C O N',), depr='ElementSymbol'),
    'hetero':    Field('hetero', bool, doc='hetero flag', 
                       selstr=('hetero', 'hetero and not water'), 
                       depr='HeteroFlag'),
    'occupancy': Field('occupancy', float, var='occupancies', 
                       doc='occupancy value', meth_pl='Occupancies',
                       selstr=('occupancy 1', 'occupancy > 0')),
    'resname':   Field('resname', '|S6', doc='residue name', 
                       selstr=('resname ALA GLY',), depr='ResidueName'),
    'resnum':    Field('resnum', int, doc='residue number', none='hv',
                       selstr=('resnum 1 2 3', 'resnum 120A 120B', 
                               'resnum 10 to 20', 'resnum 10:20:2', 
                               'resnum < 10'), synonym='resid',
                       depr='ResidueNumber'),
    'secondary': Field('secondary', '|S1', var='secondaries', 
                       doc='secondary structure assignment', 
                       meth='Secstr', synonym='secstr',
                       selstr=('secondary H E', 'secstr H E'),  
                       depr='SecondaryStr'),
    'segment':   Field('segment', '|S6', doc='segment name', meth='Segname',
                       selstr=('segment PROT', 'segname PROT'), 
                       synonym='segname', depr='SegmentName'),
    'siguij':    Field('siguij', float, doc='standard deviations for '
                       'anisotropic temperature factor', meth='Anistd', ndim=2, 
                       depr='AnisoStdDev'),
    'serial':    Field('serial', int, doc='serial number (from file)', 
                       doc_pl='serial numbers (from file)', none='sn2i', 
                       selstr=('serial 1 2 3', 'serial 1 to 10', 
                       'serial 1:10:2', 'serial < 10'), depr='SerialNumber'),
    'beta':      Field('beta', float, doc='β-value (temperature factor)', 
                       doc_pl='β-values (or temperature factors)', 
                       selstr=('beta 555.55', 'beta 0 to 500', 'beta 0:500', 
                       'beta < 500'), depr='TempFactor'),
    'icode':     Field('icode', '|S1', doc='insertion code', none='hv', 
                       selstr=('icode A', 'icode _'), depr='InsertionCode'),
    'type':      Field('type', '|S6', selstr=('type CT1 CT2 CT3',), 
                       depr='AtomType'),
    'charge':    Field('charge', float, doc='partial charge',  
                       selstr=('charge 1', 'abs(charge) == 1', 'charge < 0')),
    'mass':      Field('mass', float, var='masses', doc_pl='masses', 
                       meth_pl='Masses', selstr=('12 <= mass <= 13.5',)),
    'radius':    Field('radius', float, var='radii', doc='radius',  
                       doc_pl='radii', meth_pl='Radii', 
                       selstr=('radii < 1.5', 'radii ** 2 < 2.3')),
    'resindex':  Field('resindex', int, var='resindices', doc='residue index',  
                       doc_pl='residue indices', meth_pl='Resindices',
                       selstr=('resindex 0'), readonly=True, 
                       call=['getHierView']),
    'chindex':   Field('chindex', int, var='chindices', doc='chain index',  
                       doc_pl='chain indices', meth_pl='Chindices',
                       selstr=('chindex 0'), readonly=True, 
                       call=['getHierView']),
    'segindex':  Field('segindex', int, var='segindices', doc='segment index',  
                       doc_pl='segment indices', meth_pl='Segindices',
                       selstr=('segindex 0'), readonly=True, 
                       call=['getHierView']),
}

ATOMIC_ATTRIBUTES = {}
for field in ATOMIC_DATA_FIELDS.values():
    ATOMIC_ATTRIBUTES[field.var] = field

def wrapGetMethod(fn):
    def getMethod(self):
        return fn(self)
    return getMethod
def wrapSetMethod(fn):
    def setMethod(self, data):
        return fn(self, data)
    return setMethod

__doc__ += """

Common methods
-------------------------------------------------------------------------------

Atomic data contained in a PDB file can be accessed and changed using ``get`` 
and ``set`` methods defined for :class:`Atomic` classes.  To provide a coherent
interface, these methods are defined for :class:`AtomGroup`, :class:`Atom`, 
:class:`Selection`, :class:`Chain`, :class:`Residue`, and :class:`AtomMap` 
classes, with the following exceptions: 

* Names of methods of the :class:`Atom` class are in singular form.
* ``set`` methods are not defined for the :class:`AtomMap` class.

The list of methods are below (they link to the documentation of the 
:class:`AtomGroup` methods):
 
======================  =======================================================
Get/set method          Description
======================  =======================================================
``get/setCoords``       get/set coordinates of atoms
"""

keys = ATOMIC_DATA_FIELDS.keys()
keys.sort()

for key in keys:
    field = ATOMIC_DATA_FIELDS[key]
    __doc__ += '``get/set{0:13s}  get/set {1:s}\n'.format(field.meth_pl+'``', 
                                                          field.doc_pl)

__doc__ += """
======================  =======================================================

.. note:: Note that ``get`` methods return a copy of the data. Changes in the 
   array obtained by calling one of the above methods will not be saved in the
   :class:`AtomGroup` instance. To change the data stored in :class:`AtomGroup`
   instance, use ``set`` methods.

Other functions common to all atomic classes is given below:

=================  ==========================================================
Method name        Description
=================  ==========================================================
``copy``           returns a deep copy of atomic data
``select``         selects a subset of atoms (see :ref:`selections`)
``numAtoms``       returns number of atoms
``numCoordsets``   returns number of coordinate sets
``getCoordsets``   returns specified coordinate sets
``getACSIndex``    returns the index of the active coordinate set
``setACSIndex``    changes the index of the active coordinate set
``getACSLabel``    returns the label of the active coordinate set
``setACSLabel``    changes the label of the active coordinate set
``iterCoordsets``  iterate over coordinate sets
``isData``         checks whether a user set attribute exists
``getData``        returns user set attribute data
``setData``        changes user set attribute data
=================  ==========================================================


Special methods
-------------------------------------------------------------------------------

Atomic classes also have the following class specific methods: 
    
======================  =======================================================
Method                  Description
======================  =======================================================
:class:`AtomGroup`  
* ``getTitle``          returns title of the atom group
* ``setTitle``          changes title of the atom group
* ``delData``           deletes a user data from the atom group
* ``addCoordset``       add a coordinate set to the atom group
* ``numChains``         returns the number of chains
* ``numResidues``       returns the total number of residues from all chains
* ``iterChains``        iterates over chains
* ``iterResidues``      iterates over all residues

                      
:class:`Atom`              
* ``getIndex``          returns atom index
* ``getName``           return atom name
* ``getSelstr``         returns string that selects the atom
                    
:class:`Selection`         
* ``getIndices``        returns indices of atoms
* ``getSelstr``         returns selection string that reproduces the selection

:class:`Segment`
* ``getSegname``        returns segment name
* ``setSegname``        changes segment name
* ``getChain``          returns chain with given identifier
* ``iterChains``        iterates over chains
* ``numChains``         returns the number of chains in the instance
* ``getSelstr``         returns a string that selects segment atoms

:class:`Chain`
* ``getChid``           returns chain identifier
* ``setChid``           changes chain identifier
* ``getResidue``        returns residue with given number
* ``iterResidues``      iterates over residues
* ``numResidues``       returns the number of residues in the instance
* ``getSequence``       returns single letter amino acid sequence
* ``getSelstr``         returns a string that selects chain atoms
                      
:class:`Residue`
* ``getIndices``        returns indices of atoms
* ``getAtom``           returns :class:`Atom` with given name
* ``getChain``          returns :class:`Chain` of the residue
* ``getChid``           returns chain identifier
* ``getIcode``          returns residue insertion code
* ``setIcode``          changes residue insertion code 
* ``getResname``        returns residue name
* ``setResname``        changes residue name
* ``getResnum``         returns residue number
* ``setResnum``         changes residue number
* ``getSelstr``         returns a string that selects residue atoms

:class:`AtomMap`
* ``getIndices``        returns indices of atoms
* ``getTitle``          returns name of the atom map
* ``setTitle``          changes name of the atom map
* ``numMapped``         returns number of mapped atoms
* ``numUnmapped``       returns number of unmapped atoms
* ``getMapping``        returns mapping of indices
* ``getMappedFlags``    returns an boolean array indicating mapped atoms
* ``getUnmappedFlags``  returns an boolean array indicating unmapped atoms
======================  =======================================================

Functions common to :class:`Atom`, :class:`Selection`, :class:`Chain`,
:class:`Residue`, and :class:`AtomMap` include: 
    
======================  =======================================================
Method                  Description
======================  =======================================================
* ``getAtomGroup``      returns the associated :class:`AtomGroup`
* ``getIndices``        returns the indices of atoms
======================  =======================================================


Behavioral differences
-------------------------------------------------------------------------------

Atomic classes behave differently to indexing and to calls of certain built-in 
functions.  These differences are:

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

Segment    * :func:`len` returns the number of chains in the segment.
           * :func:`iter` yields :class:`Chain` instances.
           * Indexing by:
                
             - *chain identifier* (:func:`str`), 
               e.g. ``A`` returns a :class:`Chain`.

Chain      * :func:`len` returns the number of residues in the chain.
           * :func:`iter` yields :class:`Residue` instances.
           * Indexing by:
                
             - *residue number [, insertion code]* (:func:`tuple`), 
               e.g. ``10`` or  ``10, "B"`` returns a :class:`Residue`.
             - *slice* (:func:`slice`), e.g, ``10:20`` returns a list of  
               :class:`Residue` instances.
                    
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
  - *segment name, chain identifier, residue number [, insertion code]* 
    (:func:`tuple`), e.g. ``"PROT", "A", 10`` or  ``"PROT", "A", 10, "B"`` 
    returns a :class:`Residue`
    

"""

__doc__ += """

:mod:`prody.atomic`
===============================================================================

Classes
-------

    * :class:`AtomGroup`
    * :class:`Atom`
    * :class:`Segment`
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
      
    __slots__ = ['_acsi']
    
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
           this causes problems in :mod:`select` module."""
        
        if isinstance(other, Atomic):
            # AtomMaps may need special handling
            if self is other:
                return True
            elif isinstance(self, AtomPointer) and \
                isinstance(other, AtomPointer):
                self_indices = self._indices
                if len(self_indices) == len(other):
                    if np.all(self_indices == other._getIndices()):
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
            if prody.select.isKeyword(items[0]) or items[0] == 'not' or \
               items[0] in prody.select.MACROS:
                selstr = ' '.join(items)
                return prody.ProDyAtomSelect.select(self, selstr)
        raise AttributeError("'{0:s}' object has no attribute '{1:s}' "
                             "and '{2:s}' is not a valid selection string"
                             .format(self.__class__.__name__, name, selstr))
    
    def getActiveCoordsetIndex(self):
        """Deprecated, use :meth:`getACSIndex`."""
        
        prody.deprecate('getActiveCoordsetIndex', 'getACSIndex')
        return self.getACSIndex()
    
    def getACSIndex(self):
        """Return index of the active coordinate set."""
        
        return self._acsi
    
    def select(self, selstr, **kwargs):
        """Return atoms matching the criteria in *selstr*.
        
        .. seealso:: :meth:`~prody.select.Select.select()` for more usage 
           details."""
        
        return prody.ProDyAtomSelect.select(self, selstr, **kwargs)


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
            
            
            # DEPRECATIONS
            if field.depr:
                depr = field.depr_pl
                getDepr = 'get' + depr
                setDepr = 'set' + depr
                # Define public method for retrieving a copy of data array
                def getData(self, old=getDepr, new=getMeth):
                    prody.deprecate(old, new, 4)
                    return self.__getattribute__(new)() 
                getData = wrapGetMethod(getData)
                getData.__name__ = getDepr
                getData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(getMeth)
                setattr(cls, getDepr, getData)
                
                # Define public method for setting values in data array
                def setData(self, value, old=setDepr, new=setMeth):
                    prody.deprecate(old, new, 4)
                    self.__getattribute__(new)(value) 
                setData = wrapSetMethod(setData)
                setData.__name__ = setDepr 
                setData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(setMeth)
                setattr(cls, setDepr, setData)


class AtomGroup(Atomic):
    
    """A class for storing and accessing atomic data.
    
    The number of atoms of the atom group is inferred at the first set method
    call from the size of the data array. 

    **Atomic Data**
    
    All atomic data is stored in :class:`numpy.ndarray` instances.

    **Get and Set Methods**
    
    :meth:`get` methods return copies of the data arrays. 
    
    :meth:`set` methods accept data in :class:`list` or :class:`~numpy.ndarray` 
    instances. The length of the list or array must match the number of atoms 
    in the atom group. Set method sets attributes of all atoms at once.
    
    Atom groups with multiple coordinate sets may have one of these sets as 
    the active coordinate set. The active coordinate set may be changed using
    :meth:`setACSIndex()` method.  :meth:`getCoors` returns coordinates from 
    the active set.
    
    To access and modify data associated with a subset of atoms in an atom 
    group, :class:`Selection` instances may be used. A selection from an atom 
    group has initially the same coordinate set as the active coordinate set.
    
    User can iterate over atoms and coordinate sets in an atom group. To 
    iterate over residues and chains, get a hierarchical view of the atom 
    group by calling :meth:`getHierView()`.
    
    """
    
    __metaclass__ = AtomGroupMeta
    
    __slots__ = ['_acsi', '_title', '_n_atoms', '_coords', '_n_csets',
                 '_cslabels', 
                 '_hv', '_sn2i',
                 '_traj', '_tcsi', 
                 '_timestamps', '_kdtrees',
                 '_data', '_bonds', '_bmap']
    
    def __init__(self, title='Unnamed'):
        """Instantiate an AtomGroup with a *title*."""
        
        self._title = str(title)
        self._n_atoms = 0
        self._coords = None
        self._cslabels = []
        self._acsi = None                   # Active Coordinate Set Index
        self._n_csets = 0
        self._hv = None
        self._sn2i = None
        self._traj = None
        self._tcsi = None # Trajectory Coordinate Set Index
        self._timestamps = None
        self._kdtrees = None
        self._data = dict()
        self._bmap = None
        self._bonds = None
        
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
              
            * one of :meth:`nextFrame`, or :meth:`gotoFrame` methods is called.
        """
        
        if index is None:
            self._timestamps = np.zeros(self._n_csets)
            self._timestamps.fill(time.time())
            self._kdtrees = [None] * self._n_csets
        else:
            self._timestamps[index] = time.time()
            self._kdtrees[index] = None

    def _getKDTree(self, index):
        """Return KDTree for coordinate set at given index."""

        if self._n_csets > 0:
            if index is None:
                index = self._acsi
            kdtree = self._kdtrees[index]
            if kdtree is None:
                kdtree = prody.measure.getKDTree(self._coords[index])
                self._kdtrees[index] = kdtree
            return kdtree
        else:
            return None

    def __repr__(self):
        if self._traj is None:
            if self._n_csets:
                return ('<AtomGroup: {0:s} ({1:d} atoms; {2:d} coordinate '
                        'sets, active set index: {3:d})>').format(self._title, 
                                    self._n_atoms, self._n_csets, self._acsi)
            else:
                return ('<AtomGroup: {0:s} ({1:d} atoms; {2:d} coordinate '
                        'sets)>').format(self._title,  self._n_atoms, 
                        self._n_csets)
        else:
            return ('<AtomGroup: {0:s} ({1:d} atoms; trajectory {2:s}, '
                    'frame index {3:d})>').format(self._title, 
                    self._n_atoms, self._traj.getTitle(), self._tcsi, 
                    len(self._traj))
        
    def __str__(self):
        return ('AtomGroup {0:s}').format(self._title)
        return ('{0:s} ({1:d} atoms; {2:d} coordinate sets, active '
               'set index: {3:d})').format(self._title, 
              self._n_atoms, self._n_csets, self._acsi)

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
    
    def __iter__(self):
        """Yield atom instances."""
        
        acsi = self._acsi
        for index in xrange(self._n_atoms):
            yield Atom(self, index, acsi)

    def __len__(self):
        return self._n_atoms
    
    def __add__(self, other):
        """.. versionadded:: 0.5"""
        
        if isinstance(other, AtomGroup):
            
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
        
        elif isinstance(other, prody.VectorBase):
            if self._n_atoms != other.numAtoms(): 
                raise ValueError('Vector/Mode must have same number of atoms '
                                 'as the AtomGroup')
            self.addCoordset(self._coords[self._acsi] + 
                             other._getArrayNx3())
            self.setACSIndex(self._n_csets - 1)
        else:
            raise TypeError('can only concatenate two AtomGroup`s or can '
                            'deform AtomGroup along a Vector/Mode')

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

    def getName(self):
        """Deprecated, use :meth:`getTitle`."""

        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
        
    def getTitle(self):
        """Return title of the atom group instance."""
        
        return self._title
    
    def setName(self, name):
        """Deprecated, use :meth:`setTitle`."""

        prody.deprecate('setName', 'setTitle')
        return self.setTitle(name)
        
    def setTitle(self, title):
        """Set title of the atom group instance."""
        
        self._title = str(title)
    
    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""

        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._n_atoms
    
    iterAtoms = __iter__
    
    def getCoordinates(self):
        """Deprecated, use :meth:`getCoords`."""
        
        prody.deprecate('getCoordinates', 'getCoords')
        return self.getCoords()
        
    def getCoords(self):
        """Return a copy of coordinates from active coordinate set."""
        
        if self._coords is None:
            return None
        return self._coords[self._acsi].copy()
    
    def _getCoords(self): 
        """Return a view of coordinates from active coordinate set."""
        
        if self._coords is None:
            return None
        return self._coords[self._acsi]

    def setCoordinates(self, coordinates):
        """Deprecated, use :meth:`setCoords`."""
        
        prody.deprecate('setCoordinates', 'setCoords')
        return self.setCoords(coordinates)
        
    def setCoords(self, coords, label=None):
        """Set coordinates.  *coords* must be a :class:`numpy.ndarray` 
        instance.  If the shape of the coordinates array is 
        (n_csets,n_atoms,3), the given array will replace all coordinate sets. 
        To avoid it, :meth:`addCoordset` may be used.  If the shape of the 
        coordinates array is (n_atoms,3) or (1,n_atoms,3), the coordinate set 
        will replace the coordinates of the currently active coordinate set.
        
        .. versionadded:: 0.9.3
           *label* argument is added to allow labeling coordinate sets.  
           *label* may be a string or a list of strings length equal to the
           number of coordinate sets."""

        coordinates = checkCoords(coords, 'coords',
                                  cset=True, n_atoms=self._n_atoms,
                                  reshape=True)
        if self._n_atoms == 0:
            self._n_atoms = coordinates.shape[-2] 
        acsi = None
        if self._coords is None:
            self._coords = coordinates
            self._n_csets = coordinates.shape[0]
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
            if coordinates.shape[0] == 1:
                acsi = self._acsi
                self._coords[acsi] = coordinates[0]
                self._setTimeStamp(acsi)
                if isinstance(label, str):
                    self._cslabels[self._acsi] = label
            else:
                self._coords = coordinates
                self._n_csets = coordinates.shape[0]
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
                        LOGGER.warning('all items of `label` must be strings')
                else:
                    LOGGER.warning('`label` must have same length as the '
                                   '`coords` array')
            else:
                LOGGER.warning('`label` must be a string or list of strings')
        elif label is not None:
            if isinstance(label, str):
                self._cslabels[acsi] = label
            elif isinstance(label, (list, tuple)):
                if len(label) == 1:
                    if isinstance(label[0], str):
                        self._cslabels[acsi] = label
                    else:
                        LOGGER.warning('all items of `label` must be strings')
                else:
                    LOGGER.warning('length of `label` must be one')
            else:
                LOGGER.warning('`label` must be a string or list of strings')
                    
            
    def addCoordset(self, coords, label=None):
        """Add a coordinate set to the atom group.
        
        .. versionchanged:: 0.6.2
            :class:`~prody.ensemble.Ensemble` and :class:`Atomic` instances are 
            accepted as *coords* argument."""
        
        if self._traj is not None:
            raise AttributeError('AtomGroup is locked for coordinate set '
                                 'addition/deletion when its associated with '
                                 'a trajectory')
        if isinstance(coords, (prody.ensemble.Ensemble, Atomic)):
            if self._n_atoms != coords.numAtoms(): 
                raise ValueError('coords must have same number of atoms')
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
        self._timestamps[len(timestamps):] = time.time()
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
        if self._traj is not None:
            raise AttributeError('AtomGroup is locked for coordinate set '
                                 'addition/deletion when its associated with '
                                 'a trajectory')
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

    def getNumOfCoordsets(self):
        """Deprecated, use :meth:`numCoordsets`."""
        
        prody.deprecate('getNumOfCoordsets', 'numCoordsets')
        return self.numCoordsets()
        
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

    def setActiveCoordsetIndex(self, index):
        """Deprecated, use :meth:`setACSIndex`."""
        
        prody.deprecate('setActiveCoordsetIndex', 'setACSIndex')
        return self.setACSIndex(index)
        
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
            * a selection string
            
        .. versionchanged:: 0.7.1
           If selection string does not select any atoms, ``None`` is returned.
        
        .. versionchanged:: 0.8
           User data stored in the atom group is also copied.
           
        .. versionchanged:: 0.9.2
           Copy AtomGroup title does not start with 'Copy of'.
        
        Note that association of an atom group with a trajectory is not copied.
        """
        
        title = self._title
        if which is None:
            indices = None
            newmol = AtomGroup('{0:s}'.format(title))
            newmol.setCoords(self._coords.copy())
        elif isinstance(which, int):
            indices = [which]
            newmol = AtomGroup('{0:s} index {1:d}'.format(title, which))
        elif isinstance(which, str):
            indices = prody.ProDyAtomSelect.getIndices(self, which)
            if len(indices) == 0:
                return None
            newmol = AtomGroup('{0:s} selection "{1:s}"'
                               .format(title, which))
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
                bonds = trimBonds(bonds, bmap, indices)
                if bonds is not None:
                    newmol.setBonds(bonds)
        return newmol
    
    __copy__ = copy
    
    def getHierView(self):
        """Return a hierarchical view of the atom group."""
        
        if self._hv is None:
            self._hv = HierView(self)
        return self._hv
    
    def getNumOfChains(self):
        """Deprecated, use :meth:`numChains`."""
        
        prody.deprecate('getNumOfChains', 'numChains')
        return self.numChains()
    
    def numChains(self):
        """Return number of chains.
        
        .. versionadded:: 0.7.1"""
        
        return self.getHierView().numChains()
    
    def iterChains(self):
        """Iterate over chains.
        
        .. versionadded:: 0.7.1"""
        
        return self.getHierView().iterChains()

    def numSegments(self):
        """Return number of segments.
        
        .. versionadded:: 0.9.3"""
        
        return self.getHierView().numSegments()
    
    def iterSegments(self):
        """Iterate over chains.
        
        .. versionadded:: 0.9.3"""
        
        return self.getHierView().iterSegments()


    def getNumOfResidues(self):
        """Deprecated, use :meth:`numResidues`."""
        
        prody.deprecate('getNumOfResidues', 'numResidues')
        return self.numResidues()
        
    def numResidues(self):
        """Return number of residues.
        
        .. versionadded:: 0.7.1"""
        
        return self.getHierView().numResidues()

    def iterResidues(self):
        """Iterate over residues.
        
        .. versionadded:: 0.7.1"""
        
        return self.getHierView().iterResidues()


    def getAttrNames(self):
        """Deprecated, use :meth:`getDataLabels`."""
        
        prody.deprecate('getAttrNames', 'getDataLabels')
        return self.getDataLabels()

    def getDataLabels(self):
        """Return list of user data labels.
        
        .. versionadded:: 0.8"""
        
        return [key for key, data in self._data.iteritems() 
                    if data is not None]
        
    def getAttrType(self, name):
        """Deprecated, use :meth:`getDataType`."""
        
        prody.deprecate('getAttrType', 'getDataType')
        return self.getDataType(name)
        
    def getDataType(self, label):
        """Return type of the user data (i.e. data.dtype) associated with
        *label*, or ``None`` label is not used.
        
        .. versionadded:: 0.9"""
        
        try:
            return self._data[label].dtype
        except KeyError:
            return None

    def setAttribute(self, name, data):
        """Deprecated, use :meth:`setData`."""
        
        prody.deprecate('setAttribute', 'setData')
        return self.setData(name, data)
        
    def setData(self, label, data):
        """Store atomic *data* under *label*.
        
        .. versionadded:: 0.7.1
        
        *label* must:
            
            * start with a letter
            * contain only alphanumeric characters and underscore
            * not be a reserved word 
              (see :func:`~prody.select.getReservedWords`)

        *data* must be a :func:`list` or a :class:`numpy.ndarray`, its length 
        must be equal to the number of atoms, and the type of data array must 
        be one of:
            
            * :class:`bool`
            * :class:`float`
            * :class:`int`
            * :class:`string`
        
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
            
        if prody.select.isReserved(label):
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
    
    def delAttribute(self, name):
        """Deprecated, use :meth:`delData`."""
        
        prody.deprecate('delAttribute', 'delData')
        return self.delData(name)
        
    def delData(self, label):
        """Return data associated with *label* and remove it from the atom 
        group.  If data associated with *label* is not found, ``None`` will 
        be returned.
        
        .. versionadded:: 0.7.1"""
        
        if not isinstance(label, str):
            raise TypeError('label must be a string')
        return self._data.pop(label, None)
    
    def getAttribute(self, name):
        """Deprecated, use :meth:`getData`."""
        
        prody.deprecate('getAttribute', 'getData')
        return self.getData(name)
        
    def getData(self, label):
        """Return a copy of the data array associated with *label*, or ``None`` 
        if such data is not present.
        
        .. versionadded:: 0.7.1"""
        
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

    def isAttribute(self, name):
        """Deprecated, use :meth:`isData`."""
        
        prody.deprecate('isAttribute', 'isData')
        return self.isData(name)

    def isData(self, label):
        """Return **True** if *label* is user data.
        
        .. versionadded:: 0.7.1"""
        
        return label in self._data and self._data[label] is not None
  
    def getBySerial(self, serial, stop=None, step=None):
        """Get an atom(s) by *serial* number (range).  *serial* must be zero or 
        a positive integer. *stop* may be ``None``, or an integer greater than 
        *serial*.  ``getBySerial(i, j)`` will return atoms whose serial numbers
        are i+1, i+2, ..., j-1.  Atom whose serial number is *stop* will be 
        excluded as it would be in indexing a Python :class:`list`.  *step* 
        (default is 1) specifies increment.  If atoms with matching serial 
        numbers are not found, ``None`` will be returned. 
        
        .. versionadded:: 0.8"""

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

    def getBySerialRange(self, start, stop, step=None):
        """Deprecated, use :meth:`getBySerial`."""
        
        prody.deprecate('getBySerialRange', 'getBySerial')
        return self.getBySerial(start, stop, step)


    def setTrajectory(self, trajectory):              
        """Associates atom group with a *trajectory*.  *trajectory* may be a 
        filename or a :class:`~prody.ensemble.Trajectory` instance.  Number of 
        atoms in the atom group and the trajectory must match.  At association
        a new coordinate set will be added to the atom group.  
        :meth:`nextFrame`, and :meth:`gotoFrame` methods can be used to read 
        coordinate sets from the trajectory.  To remove association with a 
        trajectory, pass ``None`` as trajectory argument.  When atom group is 
        associated with a trajectory, it will be locked for coordinate set 
        addition/deletion operations.
        
        .. versionadded:: 0.8"""
        
        if trajectory is None:
            self._tcsi = None
            self._traj = None
            self.delCoordset(self._acsi)
            self._cslabels.pop()
        else:
            if isinstance(trajectory, str):
                trajectory = prody.Trajectory(trajectory)
            elif not isinstance(trajectory, prody.TrajectoryBase):
                raise TypeError('trajectory must be a file name or a '
                                'TrajectoryBase instance')
            if self._n_atoms != trajectory.numAtoms():
                raise ValueError('trajectory must have same number of atoms')
            self._tcsi = trajectory.getNextIndex()
            self._cslabels.append(trajectory.getTitle())
            self.addCoordset(trajectory.nextCoordset())
            self._acsi = self._n_csets - 1
            self._traj = trajectory
        
    def getTrajectory(self):
        """Return trajectory associated with the atom group."""
        
        return self._traj
    
    def nextFrame(self, step=1):
        """Read the next frame from the trajectory and update coordinates.
        *step* can be incremented to skip frames.
        
        .. versionadded:: 0.8"""
        
        if not isinstance(step, int) or step < 1:
            raise TypeError('step must be a positive integer')
        nfi = self._traj.getNextIndex()
        if step > 1:
            self._traj.skip(step - 1)
        if nfi - self._tcsi == 1:
            self._tcsi = nfi
            self._coords[self._acsi] = self._traj.nextCoordset()
            self._setTimeStamp(self._acsi)
        else:
            self._gotoFrame(self._tcsi + step)
                
    def skipFrame(self, n=1): 
        """Deprecated, use :meth:`nextFrame`."""
        
        prody.deprecate('skipFrame', 'nextFrame')
        return self.nextFrame(n+1)
    
    def gotoFrame(self, n):
        """Read frame *n* from the trajectory and update coordinates.
        
        .. versionadded:: 0.8"""
        
        self._traj.goto(n)
        self._tcsi = self._traj.getNextIndex()
        self._coords[self._acsi] = self._traj.nextCoordset()
        self._setTimeStamp(self._acsi)
    
    def getFrameIndex(self):
        """Return current trajectory frame index, ``None`` if atoms are not
        associated with a trajectory.
        
        .. versionadded:: 0.8"""
        
        return self._tcsi
    
    def getACSLabel(self):
        """Return active coordinate set label.
        
        .. versionadded:: 0.9.3"""
        
        if self._n_csets:
            return self._cslabels[self._acsi]

    def setACSLabel(self, label):
        """Set active coordinate set label.
        
        .. versionadded:: 0.9.3"""

        if self._n_csets:
            if isinstance(label, (str, NoneType)):
                self._cslabels[self._acsi] = label 
            else:
                raise TypeError('`label` must be a string')
    
    def getCSLabels(self):
        """Return coordinate set labels.
        
        .. versionadded:: 0.9.3"""
        
        if self._n_csets:
            return list(self._cslabels)

    def setCSLabels(self, labels):
        """Set coordinate set labels. *labels* must be a list of strings.
        
        .. versionadded:: 0.9.3"""
        
        if isinstance(labels, list):
            if len(labels) == self._n_csets:
                if all(isinstance(lbl, (str, NoneType)) for lbl in labels):
                    self._cslabels = list(labels)
                else:
                    raise ValueError('all items of labels must be strings')
            else:
                raise ValueError('length of labels must be equal to number of '
                                 'coordinate sets')
        else:
            raise TypeError('labels must be a list')                

    def setBonds(self, bonds):
        """Set covalent bonds between atoms.  *bonds* must be a list or an
        array of pairs of indices.  All bonds must be set at once.  An array
        with number of bonds will be generated and stored as *numbonds*.
        This can be used in atom selections, e.g. ``ag.select('numbonds 0')``
        can be used to select ions in a system.
        
        .. versionadded:: 0.9.3"""
        
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

    def iterBonds(self):
        """Yield bonds in the atom group.  Bonds must be set first using 
        :meth:`setBonds`.
        
        .. versionadded:: 0.9.3"""
        
        if self._bonds is None:
            return
        acsi = self._acsi
        for bond in self._bonds:
            yield Bond(self, bond, acsi)
        
    def numBonds(self):
        """Return number of bonds.  Bonds must be set first using 
        :meth:`setBonds`.
        
        .. versionadded:: 0.9.3"""
        
        if self._bonds is not None:
            return self._bonds.shape[0]

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
            self._acsi = atomgroup.getACSIndex()
        else: 
            self._acsi = int(acsi)

    def __add__(self, other):
        """Returns an :class:`AtomMap` instance. Order of pointed atoms are
        preserved.
        
        .. versionadded:: 0.5"""
        
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
                           'the same.  Result will have {0:d}'.format(acsi))
        
        title = '({0:s}) + ({1:s})'.format(str(self), str(other))
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
                       np.concatenate(unmapped), title, acsi)
    
    def _getTimeStamp(self, index=None):
        
        if index is None:
            return self._ag._getTimeStamp(self._acsi)
        else:
            return self._ag._getTimeStamp(index)
    
    def iterAtoms(self):
        """Yield atoms."""
        
        ag = self._ag
        acsi = self._acsi
        for i in self._getIndices():
            yield Atom(ag, i, acsi)
    
    def isAttribute(self, name):    
        """Deprecated, use :meth:`isData`."""
        
        prody.deprecate('isAttribute', 'isData')
        return self.isData(name)
        
    def isData(self, label):
        """Return ``True`` if *label* is a user data.
        
        .. versionadded:: 0.7.1"""
        
        return self._ag.isData(label)


    def getAttrType(self, name):
        """Deprecated, use :meth:`getDataType`."""
        
        prody.deprecate('getAttrType', 'getDataType')
        return self.getDataType(name)
        
    def getDataType(self, label):
        """Return type of the user data, ``None`` if data label is not present.
        
        .. versionadded:: 0.9"""
        
        return self._ag.getDataType(label)
    
    def getAtomGroup(self):
        """Return associated atom group."""
        
        return self._ag
    
    def getNumOfCoordsets(self):
        """Deprecated, use :meth:`numCoordsets`."""
        
        prody.deprecate('getNumOfCoordsets', 'numCoordsets')
        return self.numCoordsets()
        
    def numCoordsets(self):
        """Return number of coordinate sets."""
        
        return self._ag._n_csets

    def setActiveCoordsetIndex(self, index):
        """Deprecated, use :meth:`setACSIndex`."""
        
        prody.deprecate('setActiveCoordsetIndex', 'setACSIndex')
        self.setACSIndex(index)
        
    def setACSIndex(self, index):
        """Set the coordinate set at *index* active."""
        
        if self._ag._coords is None:
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
    
    __copy__ = copy
     
    def nextFrame(self):
        """Read the next frame from the trajectory and update coordinates.
        
        .. versionadded:: 0.8"""
        
        self._ag.nextFrame()
                
    def skipFrame(self, n=1): 
        """Deprecated, use :meth:`nextFrame`"""
        
        prody.deprecate('skipFrame', 'nextFrame') 
        self._ag.nextFrame(n+1)
    
    def gotoFrame(self, n):
        """Read frame *n* from the trajectory and update coordinates.
        
        .. versionadded:: 0.8"""
        
        self._ag.gotoFrame(n)
    
    def getFrameIndex(self):
        """Return current frame index.
        
        .. versionadded:: 0.8"""
        
        return self._ag.getFrameIndex()
            
    def getACSLabel(self):
        """Return active coordinate set label.
        
        .. versionadded:: 0.9.3"""
        
        if self._ag._n_csets:
            return self._ag._cslabels[self._acsi]

class AtomMeta(type):

    def __init__(cls, name, bases, dict):
        
        for field in ATOMIC_DATA_FIELDS.values():
            
            meth = field.meth
            getMeth = 'get' + meth
            setMeth = 'set' + meth
            # Define public method for retrieving a copy of data array
            if field.call:
                def getData(self, var=field.var, call=field.call):
                    for meth in call:
                        getattr(self._ag, meth)()
                    array = self._ag._data[var]
                    return array[self._index] 
            else:
                def getData(self, var=field.var):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    return array[self._index] 
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            getData.__doc__ = field.getDocstr('set', False)
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)
            
            if field.readonly:
                continue
            
            # Define public method for setting values in data array
            def setData(self, value, var=field.var, none=field.none):
                array = self._ag._data[var]
                if array is None:
                    raise AttributeError('attribute of the AtomGroup is '
                                         'not set')
                array[self._index] = value
                if None:
                    self._ag.__setattr__('_' + none,  None)
            setData = wrapSetMethod(setData)
            setData.__name__ = setMeth 
            setData.__doc__ = field.getDocstr('set', False)
            setattr(cls, setMeth, setData)
            
            if field.depr:
                depr = field.depr
                getDepr = 'get' + depr
                setDepr = 'set' + depr
                
                # Define public method for retrieving a copy of data array
                def getData(self, old=getDepr, new=getMeth):
                    prody.deprecate(old, new, 4)
                    return self.__getattribute__(new)() 
                getData = wrapGetMethod(getData)
                getData.__name__ = getDepr
                getData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(getMeth)
                setattr(cls, getDepr, getData)
                
                # Define public method for setting values in data array
                def setData(self, value, old=setDepr, new=setMeth):
                    prody.deprecate(old, new, 4)
                    self.__getattribute__(new)(value)
                setData = wrapSetMethod(setData)
                setData.__name__ = setDepr 
                setData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(setMeth)
                setattr(cls, setDepr, setData)
                

class Atom(AtomPointer):
    
    """A class for accessing and manipulating attributes of an atom 
    in a :class:`AtomGroup` instance.
    
    """
    
    __metaclass__ = AtomMeta
    __slots__ = ['_ag', '_index', '_acsi']
    
    def __init__(self, atomgroup, index, acsi=None):
        AtomPointer.__init__(self, atomgroup, acsi)
        self._index = int(index)
        
    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        if n_csets > 0:
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets, active set index: {4:d})>').format(
                    self.getName(), self._ag.getTitle(), self._index,  
                    n_csets, self._acsi)
        else:
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets)>').format(self.getName(), 
                    self._ag.getTitle(), self._index, n_csets)
                    
        sn = self.getSerial()
        if sn is None: 
            return ('<Atom: {0:s} from {1:s} (index {2:d}; {3:d} '
                    'coordinate sets, active set index: {4:d})>').format(
                    self.getName(), self._ag.getTitle(), self._index,  
                    self._ag.numCoordsets(), self._acsi)
        return ('<Atom: {0:s} from {1:s} (index {2:d}; sn {5:d}; {3:d} '
                'coordinate sets, active set index: {4:d})>').format(
                self.getName(), self._ag.getTitle(), self._index,  
                self._ag.numCoordsets(), self._acsi, sn)

    def __str__(self):
        return 'Atom {0:s} (index {1:d})'.format(self.getName(), self._index)
        sn = self.getSerial()
        if sn is None: 
            return 'Atom {0:s} (index {1:d})'.format(self.getName(), 
                                                     self._index)
        return 'Atom {0:s} (index {1:d}; sn {2:d})'.format(self.getName(), 
                                                           self._index, sn)

    def __len__(self):
        return 1
    
    def numAtoms(self):
        """Return number of atoms."""
        
        return 1
    
    def getIndex(self):
        """Return index of the atom."""
        
        return self._index
    
    def getAttribute(self, name):
        """Deprecated, use :meth:`getData`."""
        
        prody.deprecate('getAttribute', 'getData')
        return self.getData(name)
        
    def getData(self, label):
        """Return data *label*, if it exists.
        
        .. versionadded:: 0.7.1"""
        
        if self._ag.isData(label):
            return self._ag._data[label][self._index]
    
    _getData = getData
    
    def setAttribute(self, name, data):
        """Deprecated, use :meth:`setData`."""
        
        prody.deprecate('setAttribute', 'setData')
        return self.setData(name, data)
        
    def setData(self, label, data):
        """Update *data* with *label* for the atom.
        
        .. versionadded:: 0.7.1
        
        :raise AttributeError: when data *label* is not present"""
        
        if self._ag.isData(label):
            if label in READONLY:
                raise AttributeError("{0:s} is read-only".format(label))
            self._ag._data[label][self._index] = data 
        else:
            raise AttributeError("AtomGroup '{0:s}' has no data associated "
                      "with label '{1:s}'".format(self._ag.getTitle(), label))

    def getIndices(self):
        """Return index of the atom in an :class:`numpy.ndarray`."""
        
        return np.array([self._index])
    
    _getIndices = getIndices
    
    def getCoordinates(self):
        """Deprecated, use :meth:`getCoords`."""
        
        prody.deprecate('getCoordinates', 'getCoords')
        return self.getCoords()
        
    def getCoords(self):
        """Return a copy of coordinates of the atom from the active coordinate 
        set."""
        
        if self._ag._coords is None:
            return None
        return self._ag._coords[self._acsi, self._index].copy()
    
    def _getCoords(self):
        """Return a view of coordinates of the atom from the active coordinate 
        set."""
        
        if self._ag._coords is None:
            return None
        return self._ag._coords[self._acsi, self._index]
    
    def setCoordinates(self, coordinates):
        """Deprecated, use :meth:`setCoords`."""
        
        prody.deprecate('setCoordinates', 'setCoords')
        return self.setCoords(coordinates)
        
    def setCoords(self, coords):
        """Set coordinates of the atom in the active coordinate set."""
        
        self._ag._coords[self._acsi, self._index] = coords
        self._ag._setTimeStamp(self._acsi)
        
    def getCoordsets(self, indices=None):
        """Return a copy of coordinate set(s) at given *indices*, which may be 
        an integer or a list/array of integers."""
        
        if self._ag._coords is None:
            return None
        if indices is None:
            return self._ag._coords[:, self._index].copy()
        if isinstance(indices, (int, slice)):
            return self._ag._coords[indices, self._index].copy()
        if isinstance(indices, (list, np.ndarray)):
            return self._ag._coords[indices, self._index]
        raise IndexError('indices must be an integer, a list/array of integers, '
                         'a slice, or None')
       
    def _getCoordsets(self, indices=None): 
        """Return a view of coordinate set(s) at given *indices*."""
        
        if self._ag._coords is None:
            return None
        if indices is None:
            return self._ag._coords[:, self._index]
        if isinstance(indices, (int, slice)):
            return self._ag._coords[indices, self._index]
        if isinstance(indices, (list, np.ndarray)):
            return self._ag._coords[indices, self._index]
        raise IndexError('indices must be an integer, a list/array of integers, '
                         'a slice, or None')

    def iterCoordsets(self):
        """Iterate over coordinate sets by returning a copy of each 
        coordinate set."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coords[i, self._index].copy()


    def _iterCoordsets(self):
        """Iterate over coordinate sets by returning a view of each coordinate
        set."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coords[i, self._index]

    def getSelectionString(self):
        """Deprecated, use :meth:`getSelstr`."""
        
        prody.deprecate('getSelectionString', 'getSelstr')
        return self.getSelstr()
        
    def getSelstr(self):
        """Return selection string that will select this atom."""
        
        return 'index {0:d}'.format(self._index)

    def iterBonds(self):
        """Yield bonds formed by this atom.  Bonds must be set first using 
        :meth:`~AtomGroup.setBonds`.
        
        .. versionadded:: 0.9.3"""
        
        ag = self._ag
        if ag._bmap is not None:
            acsi = self._acsi
            this = self._index
            for other in self._ag._bmap[this]:
                if other == -1:
                    break
                yield Bond(ag, [this, other], acsi) 
                    
    def numBonds(self):
        """Return number of bonds formed by this atom.  Bonds must be set first
        using :meth:`~AtomGroup.setBonds`.
        
        .. versionadded:: 0.9.3"""
        
        numbonds = self._ag._data.get('numbonds')
        if numbonds is not None:
            return numbonds[self._index]
    
    def iterBonded(self):
        """Yield bonded atoms.  Bonds must be set first using 
        :meth:`~AtomGroup.setBonds`.
        
        .. versionadded:: 0.9.3"""
        
        ag = self._ag
        if ag._bmap is not None:
            acsi = self._acsi
            this = self._index
            for other in self._ag._bmap[this]:
                if other == -1:
                    break
                yield Atom(ag, other, acsi) 

class AtomSubsetMeta(type):

    def __init__(cls, name, bases, dict):

        for field in ATOMIC_DATA_FIELDS.values():
            meth = field.meth_pl
            getMeth = 'get' + meth
            setMeth = 'set' + meth
            # Define public method for retrieving a copy of data array
            if field.call:
                def getData(self, var=field.var, call=field.call):
                    for meth in call:
                        getattr(self._ag, meth)()
                    array = self._ag._data[var]
                    return array[self._indices]
            else:
                def getData(self, var=field.var):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    return array[self._indices] 
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            getData.__doc__ = field.getDocstr('get')
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)
            
            if field.readonly:
                continue
            
            # Define public method for setting values in data array
            def setData(self, value, var=field.var, none=field.none):
                array = self._ag._data[var]
                if array is None:
                    raise AttributeError(var + ' data is not set')
                array[self._indices] = value
                if none:
                    self._ag.__setattr__('_'+none,  None)
            setData = wrapSetMethod(setData)
            setData.__name__ = setMeth 
            setData.__doc__ = field.getDocstr('set')  
            setattr(cls, setMeth, setData)

            # DEPRECATIONS
            if field.depr:
                depr = field.depr_pl
                getDepr = 'get' + depr
                setDepr = 'set' + depr
                # Define public method for retrieving a copy of data array
                def getData(self, old=getDepr, new=getMeth):
                    prody.deprecate(old, new, 4)
                    return self.__getattribute__(new)() 
                getData = wrapGetMethod(getData)
                getData.__name__ = getDepr
                getData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(getMeth)
                setattr(cls, getDepr, getData)
                setattr(cls, '_' + getDepr, getData)
                
                # Define public method for setting values in data array
                def setData(self, value, old=setDepr, new=setMeth):
                    prody.deprecate(old, new, 4)
                    self.__getattribute__(new)(value)
                setData = wrapSetMethod(setData)
                setData.__name__ = setDepr 
                setData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(setMeth)
                setattr(cls, setDepr, setData)
                        
class AtomSubset(AtomPointer):
    
    """A class for manipulating subset of atomic data in an :class:`AtomGroup`.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, and active coordinate set index for the atom group.
    
    """
    
    __metaclass__ = AtomSubsetMeta    
    __slots__ = ['_ag', '_indices', '_acsi', '_selstr']
    
    def __init__(self, atomgroup, indices, acsi=None, **kwargs):
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
            indices = np.array(indices, int)
        elif not indices.dtype == int:
            indices = indices.astype(int)
        if kwargs.get('unique'):
            self._indices = indices
        else:
            self._indices = np.unique(indices)
        self._selstr = kwargs.get('selstr')
    
    def __iter__(self):
        """Iterate over atoms."""
        
        acsi = self._acsi
        ag = self._ag 
        for index in self._indices:
            yield Atom(ag, index, acsi)
    
    def __len__(self):
        return len(self._indices)
    
    def __invert__(self):
        arange = range(self._ag.numAtoms())
        indices = list(self._indices)
        while indices:
            arange.pop(indices.pop())
        sel = Selection(self._ag, arange, "not ({0:s}) ".format(
                                                 self.getSelstr()), self._acsi)        
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
                                    self.getSelstr(), other.getSelstr()), acsi)

    def __and__(self, other):
        """
        .. versionchanged:: 0.7.1
           If intersection of selections does not contain any atoms, ``None``
           is returned."""
        
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
                                    self.getSelstr(), other.getSelstr()), acsi)
               
    def getAttribute(self, name):
        """Deprecated, use :meth:`getData`."""
        
        prody.deprecate('getAttribute', 'getData')
        return self.getData(name)
        
    def getData(self, label):
        """Return a copy of the data associated with *label*, if it exists.
        
        .. versionadded:: 0.7.1"""
        
        if self._ag.isData(label):
            return self._ag._data[label][self._indices]
    
    _getData = getData
    
    def setAttribute(self, name, data):
        """Deprecated, use :meth:`setData`."""
        
        prody.deprecate('setAttribute', 'setData')
        return self.setData(name, data)
        
    def setData(self, label, data):
        """Update *data* with label *label* for the atom subset.
        
        .. versionadded:: 0.7.1
        
        :raise AttributeError: when data associated with *label* is not present
        """
        
        if self._ag.isData(label):
            if label in READONLY:
                raise AttributeError("{0:s} is read-only".format(label))
            self._ag._data[label][self._indices] = data 
        else:
            raise AttributeError("AtomGroup '{0:s}' has no data with label "
                            "'{1:s}'".format(self._ag.getTitle(), label))
    
    def getIndices(self):
        """Return a copy of the indices of atoms."""
        
        return self._indices.copy()
    
    def _getIndices(self):
        """Return indices of atoms."""
        
        return self._indices
    
    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""
        
        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of atoms."""
        
        return self._indices.__len__()

    def getCoordinates(self):
        """Deprecated, use :meth:`getCoords`."""
        
        prody.deprecate('getCoordinates', 'getCoords')
        return self.getCoords()
        
    def getCoords(self):
        """Return a copy of coordinates from the active coordinate set."""
        
        if self._ag._coords is None:
            return None
        # Since this is not slicing, a view is not returned
        return self._ag._coords[self._acsi, self._indices]
    
    _getCoords = getCoords
    
    def setCoordinates(self, coordinates):
        """Deprecated, use :meth:`setCoords`."""
        
        prody.deprecate('setCoordinates', 'setCoords')
        return self.setCoords(coordinates)
        
    def setCoords(self, coords):
        """Set coordinates in the active coordinate set."""
        
        self._ag._coords[self._acsi, self._indices] = coords
        self._ag._setTimeStamp(self._acsi)
        
    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer 
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
        """Iterate over coordinate sets by returning a copy of each coordinate 
        set."""
        
        for i in range(self._ag._n_csets):
            yield self._ag._coords[i, self._indices]

    _iterCoordsets = iterCoordsets

class Segment(AtomSubset):
    
    """Instances are generated by :class:`HierView` class.
    
    Indexing a :class:`Segment` instance by chain identifier returns 
    :class:`Chain` instances."""
    
    __slots__ = AtomSubset.__slots__ + ['_dict']

    def __init__(self, atomgroup, indices, acsi=None, **kwargs):
        AtomSubset.__init__(self, atomgroup, indices, acsi, **kwargs)
        self._dict = OrderedDict()
    
    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        if n_csets > 0:
            return ('<Segment: {0:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets, active set index: {4:d})>').format(
                    self.getSegname(), self._ag.getTitle(), 
                    self.numAtoms(), n_csets, self._acsi)
        else:
            return ('<Segment: {0:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets)>').format(self.getSegname(), 
                    self._ag.getTitle(), self.numAtoms(), n_csets)

    def __getitem__(self, chid):
        
        return self._dict.get(chid)
    
    def __str__(self):
        return 'Segment {0:s}'.format(self.getSegname())

    def __len__(self):
        return len(self._dict)
    
    def __iter__(self):
        
        return self._dict.itervalues()

    def getSegname(self):
        """Return segment name."""
        
        return self._ag._data['segments'][self._indices[0]]
    
    def setSegname(self, segname):
        """Set segment name."""
        
        self.setSegnames(segname)
    
    def numChains(self):
        """Return number of residues."""
        
        return len(self._dict)
    
    def getChain(self, chid):
        """Return chain with identifier *chid*."""
        
        return self._dict.get(chid)
    
    def iterChains(self):
        """Iterate chains in the segment."""
        
        return self._dict.itervalues()

    def getSelstr(self):
        """Return selection string that selects atoms in this segment."""
        
        if self._selstr:
            return 'segname {0:s} and ({1:s})'.format(self.getSegname(), 
                                                        self._selstr)
        else:
            return 'segname {0:s}'.format(self.getSegname())

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
        
    __slots__ = AtomSubset.__slots__ + ['_seq', '_dict', '_list', '_segment']
    
    def __init__(self, atomgroup, indices, segment, acsi=None, **kwargs):
        AtomSubset.__init__(self, atomgroup, indices, acsi, **kwargs)
        self._segment = segment
        self._seq = None
        self._dict = dict()
        self._list = list()
        
    def __len__(self):
        return len(self._list)
    
    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        if self._segment:
            segment = ' from Segment {0:s}'.format(self.getSegname())
        else:
            segment = ''
        if n_csets > 0:
            return ('<Chain: {0:s}{5:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets, active set index: {4:d})>').format(
                    self.getChid(), self._ag.getTitle(), 
                    self.numAtoms(), n_csets, self._acsi, segment)
        else:
            return ('<Chain: {0:s}{4:s} from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets)>').format(self.getChid(), 
                    self._ag.getTitle(), self.numAtoms(), n_csets, segment)

    def __str__(self):
        return ('Chain {0:s}').format(self.getChid())

    def __iter__(self):
        return self.iterResidues()
    
    def __getitem__(self, key):
        
        if isinstance(key, tuple): 
            return self.getResidue(*key) 
        elif isinstance(key, slice):
            resnums = self._getResnums()
            resnums = set(np.arange(*key.indices(resnums.max()+1)))
            _list = self._list
            return [_list[i] for (rn, ic), i in self._dict.iteritems() 
                    if rn in resnums]
        else:
            return self.getResidue(key)
    
    def getSegment(self):
        """Return segment that this chain belongs to."""
        
        return self._segment
    
    def getSegname(self):
        """Return name of the segment that this chain belongs to."""
        
        if self._segment:
            return self._segment.getSegname()
    
    def getResidue(self, number, insertcode=None):
        """Return residue with given number."""
        
        i = self._dict.get((number, insertcode or None))
        if i is not None:
            return self._list[i]

    def iterResidues(self):
        """Iterate residues in the chain."""
        
        return self._list.__iter__()
    
    def getNumOfResidues(self):
        """Deprecated, use :meth:`numResidues`."""
        
        prody.deprecate('getNumOfResidues', 'numResidues')
        return self.numResidues()
        
    def numResidues(self):
        """Return number of residues."""
        
        return len(self._list)

    def getIdentifier(self):
        """Deprecated, use :meth:`getChid`."""
        
        prody.deprecate('getIdentifier', 'getChid')
        return self.getChid()
        
    def getChid(self):
        """Return chain identifier."""
        
        return self._ag._data['chids'][self._indices[0]]
    
    def setIdentifier(self, chid):
        """Deprecated, use :meth:`setChid`."""
        
        prody.deprecate('setIdentifier', 'setChid')
        self.setChid(chid)
        
    def setChid(self, chid):
        """Set chain identifier."""
        
        self.setChids(chid)
    
    def getSequence(self):
        """Return sequence, if chain is a polypeptide."""
        
        if self._seq:
            return self._seq
        CAs = self.select('name CA').select('protein')
        if len(CAs) > 0:
            self._seq = prody.compare.getSequence(CAs.getResnames())
        else:
            self._seq = ''
        return self._seq

    def getSelectionString(self):
        """Deprecated, use :meth:`getSelstr`."""
        
        prody.deprecate('getSelectionString', 'getSelstr')
        return self.getSelstr()
        
    def getSelstr(self):
        """Return selection string that selects atoms in this chain."""

        if self._segment is None:        
            if self._selstr:
                return 'chain {0:s} and ({1:s})'.format(self.getChid(),
                                                        self._selstr)
            else:
                return 'chain {0:s}'.format(self.getChid())
        else:
            selstr = self._segment.getSelstr()
            return 'chain {0:s} and ({1:s})'.format(self.getChid(),
                                                    selstr)


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
    
    def __init__(self, atomgroup, indices, chain, acsi=None, **kwargs):
        AtomSubset.__init__(self, atomgroup, indices, acsi, **kwargs)
        self._chain = chain

    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        if self._chain:
            chain = ' from Chain {0:s}'.format(self.getChid())
        else:
            chain = ''
        
        if n_csets > 0:
            return ('<Residue: {0:s} {1:d}{2:s}{3:s} from {4:s} '
                    '({5:d} atoms; {6:d} coordinate sets, active set index: '
                    '{7:d})>').format(self.getResname(), self.getResnum(), 
                                      self.getIcode() or '', chain, 
                                      self._ag.getTitle(), len(self), 
                                      n_csets, self._acsi)
        else:        
            return ('<Residue: {0:s} {1:d}{2:s}{3:s} from {4:s} '
                    '({5:d} atoms; {6:d} coordinate sets)>').format(
                        self.getResname(), self.getResnum(), 
                        self.getIcode() or '',  chain, 
                        self._ag.getTitle(), len(self), n_csets)
            
    def __str__(self):
        return '{0:s} {1:d}{2:s}'.format(self.getResname(), self.getResnum(), 
                                         self.getIcode())

    def __getitem__(self, name):
        return self.getAtom(name)
    
    __iter__ = AtomPointer.iterAtoms
    
    def getAtom(self, name):
        """Return atom with given *name*, ``None`` if not found.  Assumes that 
        atom names in the residue are unique.  If more than one atoms with the 
        given *name* exists, the one with the smaller index will be returned.
        """
        
        if isinstance(name, str):
            nz = (self.getNames() == name).nonzero()[0]
            if len(nz) > 0:
                return Atom(self._ag, self._indices[nz[0]], self._acsi)
    
    def getChain(self):
        """Return the chain that the residue belongs to."""
        
        return self._chain
    
    def getNumber(self):
        """Deprecated, use :meth:`getResnum`."""
        
        prody.deprecate('getNumber', 'getResnum')
        return self.getResnum()
    
    def getResnum(self):
        """Return residue number."""
        
        return int(self._ag._data['resnums'][self._indices[0]])
    
    def setNumber(self, number):
        """Deprecated, use :meth:`setResnum`."""
        
        prody.deprecate('setNumber', 'setResnum')
        return self.setResnum(number)
    
    def setResnum(self, number):
        """Set residue number."""
        
        self.setResnums(number)
    
    def getName(self):
        """Deprecated, use :meth:`getName`."""
        
        prody.deprecate('getName', 'getResname')
        return self.getResname()
    
    def getResname(self):
        """Return residue name."""
        
        data = self._ag._data['resnames']
        if data is not None:
            return data[self._indices[0]]
    
    def setName(self, name):
        """Deprecated, use :meth:`setResname`."""
        
        prody.deprecate('setName', 'setResname')
        return self.setResname(name)
        
    def setResname(self, name):
        """Set residue name."""
        
        self.setResnames(name)

    def getInsertionCode(self):
        """Deprecated, use :meth:`getIcode`."""
        
        prody.deprecate('getInsertionCode', 'getIcode')
        return self.getIcode()
        
    def getIcode(self):
        """Return residue insertion code."""
        
        data = self._ag._data['icodes']
        if data is not None:
            return data[self._indices[0]]
        
    def setInsertionCode(self, icode):
        """Deprecated, use :meth:`setIcode`."""
        
        prody.deprecate('setInsertionCode', 'setIcode')
        return self.setIcode(icode)
        
    def setIcode(self, icode):
        """Set residue insertion code."""
        
        self.setIcodes(icode)
    
    def getChainIdentifier(self):
        """Deprecated, use :meth:`getChid`."""
        
        prody.deprecate('getChainIdentifier', 'getChid')
        return self.getChid()
        
    def getChid(self):
        """Return chain identifier."""
        
        if self._chain:
            return self._chain.getChid()
    
    def getSelectionString(self):
        """Deprecated, use :meth:`getSelstr`."""
        
        prody.deprecate('getSelectionString', 'getSelstr')
        return self.getSelstr()
        
    def getSelstr(self):
        """Return selection string that will select this residue."""
        
        icode = self.getIcode() or ''
        if self._chain is None:        
            if self._selstr:
                return 'resnum {0:s}{1:s} and ({1:s})'.format(
                            self.getResnum(), icode, self._selstr)
            else:
                return 'resnum {0:s}{1:s}'.format(self.getResnum(), icode)
        else:
            selstr = self._chain.getSelstr()
            return 'resnum {0:d}{1:s} and ({2:s})'.format(
                                self.getResnum(), icode, selstr)

    def getPrev(self):
        """Return preceding residue in the chain."""
        
        i = self._chain._dict.get((self.getResnum(), self.getIcode() or None))
        if i is not None and i > 0:
            return self._chain._list[i-1]
        
    def getNext(self):
        """Return following residue in the chain."""

        i = self._chain._dict.get((self.getResnum(), self.getIcode() or None))
        if i is not None and i + 1 < len(self._chain):
            return self._chain._list[i+1]


class Selection(AtomSubset):
    
    """A class for accessing and manipulating attributes of select of atoms 
    in an :class:`AtomGroup` instance.
    
    """
    
    __slots__ = AtomSubset.__slots__
    
    def __init__(self, atomgroup, indices, selstr, acsi=None, **kwargs):
        kwargs['selstr'] = selstr
        AtomSubset.__init__(self, atomgroup, indices, acsi, **kwargs)
        
    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        if n_csets > 0:
            return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets, active set index: {4:d})>').format(
                    selstr, self._ag.getTitle(), len(self), n_csets, 
                    self._acsi)
        else:
            return ('<Selection: "{0:s}" from {1:s} ({2:d} atoms; '
                    '{3:d} coordinate sets)>').format(
                    selstr, self._ag.getTitle(), len(self), n_csets)

    def __str__(self):
        selstr = self._selstr
        if len(selstr) > 33:
            selstr = selstr[:15] + '...' + selstr[-15:]  
        return 'Selection "{0:s}" from {1:s}'.format(selstr, 
                                                     self._ag.getTitle())
    
    __iter__ = AtomPointer.iterAtoms
    
    def getSelectionString(self):
        """Deprecated, use :meth:`getSelstr`."""
        
        prody.deprecate('getSelectionString', 'getSelstr')
        return self.getSelstr()
        
    def getSelstr(self):
        """Return selection string that selects this atom subset."""
        
        return self._selstr

    def getHierView(self):
        """Return a hierarchical view of the atom selection."""
        
        return HierView(self)

    def update(self):    
        """Update selection.
        
        .. versionadded:: 0.9.3"""
        
        self._indices = prody.ProDyAtomSelect.getIndices(self._ag, 
                                                         self._selstr)


class AtomMapMeta(type):
    
    def __init__(cls, name, bases, dict):
        for field in ATOMIC_DATA_FIELDS.values():
            meth = field.meth_pl
            getMeth = 'get' + meth
            if field.call:
                def getData(self, var=field.var, call=field.call):
                    for meth in call:
                        getattr(self._ag, meth)()
                    data = self._ag._data[var][self._indices]
                    result = np.zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result 
            else:
                def getData(self, var=field.var, dtype=field.dtype):
                    array = self._ag._data[var]
                    if array is None:
                        return None
                    data = self._ag._data[var][self._indices]
                    result = np.zeros((self._len,) + data.shape[1:], dtype)
                    result[self._mapping] = data
                    return result 
            getData = wrapGetMethod(getData)
            getData.__name__ = getMeth
            if field.dtype in (int, float):
                zero = '0'
            elif field.dtype == bool:
                zero = 'True'
            else:
                zero = '""'
            getData.__doc__ = field.getDocstr('get', selex=False) + \
                   'Entries for unmapped atoms will be ``{0:s}``.'.format(zero) 
            setattr(cls, getMeth, getData)
            setattr(cls, '_' + getMeth, getData)
        
            if field.depr:
                depr = field.depr_pl
                getDepr = 'get' + depr
                def getData(self, old=getDepr, new=getMeth):
                    prody.deprecate(old, new, 4)
                    return self.__getattribute__(new)()
                getData = wrapGetMethod(getData)
                getData.__name__ = getDepr
                getData.__doc__ = 'Deprecated, use :meth:`{0:s}`'.format(getMeth) 
                setattr(cls, getDepr, getData)
             


class AtomMap(AtomPointer):
    
    """A class for mapping atomic data.
    
    This class stores a reference to an :class:`AtomGroup` instance, a set of 
    atom indices, active coordinate set index, mapping for indices, and
    indices of unmapped atoms.
    
    """
    
    __metaclass__ = AtomMapMeta
    __slots__ = ['_ag', '_indices', '_acsi', '_title', '_mapping', '_unmapped', 
                 '_len']
    
    def __init__(self, atomgroup, indices, mapping, unmapped, title='Unnamed', 
                 acsi=None):
        """Instantiate with an AtomMap with following arguments:        
        
        :arg atomgroup: the atomgroup instance from which atoms are mapped
        :arg indices: indices of mapped atoms
        :arg mapping: mapping of the atoms as a list of indices
        :arg unmapped: list of indices for unmapped atoms
        :arg title: title of the AtomMap instance
        :arg acsi: active coordinate set index, if ``None`` defaults to that 
            of *atomgrup*
        
        Length of *mapping* must be equal to length of *indices*.  Number of 
        atoms (including unmapped dummy atoms) are determined from the sum of 
        lengths of *mapping* and *unmapped* arrays."""
        
        AtomPointer.__init__(self, atomgroup, acsi)
        
        if not isinstance(indices, np.ndarray):
            self._indices = np.array(indices, int)
        elif not indices.dtype == int:
            self._indices = indices.astype(int)
        else:
            self._indices = indices

        if not isinstance(mapping, np.ndarray):
            self._mapping = np.array(mapping, int)
        elif not mapping.dtype == int:
            self._mapping = mapping.astype(int)
        else:
            self._mapping = mapping

        if not isinstance(unmapped, np.ndarray):
            self._unmapped = np.array(unmapped, int)
        elif not unmapped.dtype == int:
            self._unmapped = unmapped.astype(int)
        else:
            self._unmapped = unmapped
        
        self._title = str(title)
        self._len = len(self._unmapped) + len(self._mapping)
        
    def __repr__(self):
        n_csets = self._ag.numCoordsets()
        if n_csets > 0:
            return ('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                    '{3:d} mapped; {4:d} unmapped; {5:d} coordinate sets, '
                    'active set index: {6:d})>').format(self._title,
                    self._ag.getTitle(), self._len, len(self._mapping), 
                    len(self._unmapped), n_csets, self._acsi)
        else:
            return (('<AtomMap: {0:s} (from {1:s}; {2:d} atoms; '
                    '{3:d} mapped; {4:d} unmapped; {5:d} coordinate sets)>')
                    .format(self._title, self._ag.getTitle(), self._len, 
                    len(self._mapping), len(self._unmapped), n_csets))
            
    def __str__(self):
        return 'AtomMap {0:s}'.format(self._title)
    
    def __iter__(self):
        indices = np.zeros(self._len, int)
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
        """Deprecated, use :meth:`getData`."""
        
        prody.deprecate('getAttribute', 'getData')
        return self.getData(name)
        
    def getData(self, label):
        """Return a copy of data associated with *label*, if it exists.
        
        .. versionadded:: 0.7.1"""
        
        if self._ag.isData(label):
            data = self._ag._data[label][self._indices]
            result = np.zeros((self._len,) + data.shape[1:], data.dtype)
            result[self._mapping] = data
            return result

    _getData = getData

    def getName(self):
        """Deprecated, use :meth:`getTitle`."""

        prody.deprecate('getName', 'getTitle')
        return self.getTitle()
        
    def getTitle(self):
        """Return title of the atom map instance."""
        
        return self._title
    
    def setName(self, name):
        """Deprecated, use :meth:`setTitle`."""

        prody.deprecate('setName', 'setTitle')
        return self.setTitle(name)
        
    def setTitle(self, title):
        """Set title of the atom map instance."""
        
        self._title = str(title)

    def getNumOfAtoms(self):
        """Deprecated, use :meth:`numAtoms`."""

        prody.deprecate('getNumOfAtoms', 'numAtoms')
        return self.numAtoms()
        
    def numAtoms(self):
        """Return number of mapped atoms."""
        
        return self._len

    def getNumOfUnmapped(self):
        """Deprecated, use :meth:`numUnmapped`."""

        prody.deprecate('getNumOfUnmapped', 'numUnmapped')
        return self.numUnmapped()
        
    def numUnmapped(self):
        """Return number of unmapped atoms."""
        
        return len(self._unmapped)

    def getNumOfMapped(self):
        """Deprecated, use :meth:`numMapped`."""

        prody.deprecate('getNumOfMapped', 'numMapped')
        return self.numMapped()
        
    def numMapped(self):
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
        set."""
        
        for i in range(self._ag._n_csets):
            coordinates = np.zeros((self._len, 3), float)
            coordinates[self._mapping] = self._ag._coords[i, self._indices] 
            yield coordinates
    
    _iterCoordsets = iterCoordsets

    def getCoordinates(self):
        """Deprecated, use :meth:`getCoords`."""
        
        prody.deprecate('getCoordinates', 'getCoords')
        return self.getCoords()
        
    def getCoords(self):
        """Return coordinates from the active coordinate set."""
        
        if self._ag._coords is None:
            return None
        coordinates = np.zeros((self._len, 3), float)
        coordinates[self._mapping] = self._ag._coords[self._acsi, 
                                                           self._indices] 
        return coordinates
    
    _getCoords = getCoords
    
    def setCoordinates(self, coordinates):
        """Deprecated, use :meth:`setCoords`."""
        
        prody.deprecate('setCoordinates', 'setCoords')
        return self.setCoords(coordinates)
        
    def setCoords(self, coords):
        """Set coordinates in the active coordinate set.  Length of the 
        *coordinates* array must match the number of mapped atoms."""
        
        self._ag._coords[self._acsi, self._indices] = coords
    

    def getCoordsets(self, indices=None):
        """Return coordinate set(s) at given *indices*, which may be an integer 
        or a list/array of integers."""
        
        if self._ag._coords is None:
            return None
        if indices is None:
            indices = np.arange(self._ag._n_csets)
        elif isinstance(indices, (int, long)):
            indices = np.array([indices])
        elif isinstance(indices, slice):
            indices = np.arange(indices.indices(self._ag._n_csets))
        try:
            coordsets = np.zeros((len(indices), self._len, 3))
            coordsets[:, self._mapping] = self._ag._coords[indices][:, 
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
    
    """Hierarchical views can be generated for :class:`AtomGroup` and 
    :class:`Selection` instances.  Indexing a :class:`HierView` instance 
    returns a :class:`Chain` instance.
    
    >>> from prody import *
    >>> pdb = parsePDB('1p38')
    >>> hv = pdb.getHierView()
    >>> chA = hv['A']
    >>> chA
    <Chain: A from 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>
    >>> print hv['B'] # Chain B does not exist in 1p38
    None
    
    Note that is the atom group instance have distinct segments, they will
    be considered when building the hierarchical view.  A :class:`Segment`
    instance will be generated for each distinct segment name.  Then,
    for each segment chains and residues will be evaluated.  Having
    segments in the structure will not change most behaviors of this class,
    except indexing.  For example, when indexing a hierarchical view for 
    chain P in segment PROT needs to be indexed as ``hv['PROT', 'P']``."""
    
    __slots__ = ['_atoms', '_dict', '_segments', '_chains', '_residues']

    def __init__(self, atoms, **kwargs):
        """Build hierarchical view for *atoms*."""
        
        
        if not isinstance(atoms, (AtomGroup, Selection)):
            raise TypeError('atoms must be an AtomGroup or Selection instance')
        self._atoms = atoms
        self.update(**kwargs)

        
    def __repr__(self):
        
        if self._segments:
            return ('<HierView: {0:s} ({1:d} segments, {2:d} chains, {3:d} '
                    'residues)>').format(str(self._atoms), self.numSegments(),
                                         self.numChains(), self.numResidues())
        else:
            return ('<HierView: {0:s} ({1:d} chains, {2:d} residues)>'
               ).format(str(self._atoms), self.numChains(), self.numResidues())
    
    def __str__(self):
        return 'HierView of {0:s}'.format(str(self._atoms))
    
    def __iter__(self):
        """Iterate over chains."""
        
        return self._chains.__iter__()
    
    def __len__(self):
        return len(self._chains)
    
    def __getitem__(self, key):
        
        if isinstance(key, str):
            return self._dict.get(key, self._dict.get((None, key)))
        elif isinstance(key, tuple):
            length = len(key)
            if length == 1:
                return self.__getitem__(key[0])
            elif length == 2:
                return self._dict.get(key,
                        self._dict.get((None, None, key[0], key[1] or None), 
                        self._dict.get((None, key[0], key[1], None))))
            elif length == 3:
                return self._dict.get((None, key[0] or None, 
                                       key[1], key[2] or None), 
                          self._dict.get((key[0] or None, key[1] or None, 
                                          key[2], None)))
            elif length == 4:
                return self._dict.get((key[0] or None, key[1] or None, 
                                       key[2], key[3] or None))
        elif isinstance(key, int):
            return self._dict.get((None, None, key, None))

    def getAtoms(self):
        """Return atoms for which the hierarchical view was built.
        
        .. versionadded:: 0.6.2"""
        
        return self._atoms
    
    def update(self, **kwargs):
        """Rebuild hierarchical view of atoms.  This method is called at 
        instantiation, but can be used to rebuild the hierarchical view 
        when attributes of atoms change."""
        
        array = np.array
        acsi = self._atoms.getACSIndex()
        atoms = self._atoms
        if isinstance(atoms, AtomGroup):
            ag = atoms
            _indices = np.arange(ag._n_atoms)
            selstr = False
        else:
            ag = atoms._ag
            _indices = atoms._indices
            selstr = atoms.getSelstr()
        n_atoms = len(ag)    
        self._dict = _dict = dict()
        self._chains = _chains = list()
        self._residues = _residues = list()
        self._segments = _segments = list()

        segindex = -1
        segindices = np.zeros(n_atoms, int)
        chindex = -1
        chindices = np.zeros(n_atoms, int)
        resindex = -1
        resindices = np.zeros(n_atoms, int)
        
        sgnms = ag._getSegnames()
        if sgnms is None:
            _segments = None
        else:
            if selstr:
                sgnms = sgnms[_indices]
            unique = np.unique(sgnms)
            s = sgnms[0]
            if len(unique) == 1:
                if  s != '':
                    segment = Segment(ag, _indices, acsi, unique=True, 
                                      selstr=selstr)
                    _dict[s] = segment
                    _segments.append(segment)
                    LOGGER.info('Hierarchical view contains segments.')
                else: 
                    _segments = None
            else:
                ps = None
                for i, s in enumerate(sgnms):
                    if s == ps or s in _dict:
                        continue
                    ps = s
                    segindex += 1
                    idx = _indices[i:][sgnms[i:] == s]
                    segment = Segment(ag, idx, acsi, unique=True, 
                                      selstr=selstr)
                    segindices[idx] = segindex
                    _dict[s] = segment
                    _segments.append(segment)
                LOGGER.info('Hierarchical view contains segments.')

        chids = ag._getChids()
        if chids is None:
            _chains = None
        else:
            if selstr:
                chids = chids[_indices]
            if _segments is None:
                if len(np.unique(chids)) == 1:
                    chain = Chain(ag, _indices, None, acsi, unique=True)
                    _dict[(None, chids[0] or None)] = chain
                    _chains.append(chain)
                else:
                    pc = None
                    for i, c in enumerate(chids):
                        if c == pc or (None, c) in _dict:
                            continue
                        pc = c
                        chindex += 1
                        idx = _indices[i:][chids[i:] == c]
                        chain = Chain(ag, idx, None, acsi, unique=True)
                        chindices[idx] = chindex
                        _dict[(None, c)] = chain
                        _chains.append(chain)
            else:
                pc = chids[0]
                ps = sgnms[0]
                _i = 0
                for i, c in enumerate(chids):
                    s = sgnms[i]
                    if c == pc and s == ps:
                        continue
                    s_c = (ps, pc or None)
                    chain = _dict.get(s_c)
                    if chain is None:
                        segment = _dict[ps]
                        chindex += 1
                        idx = _indices[_i:i]
                        chain = Chain(ag, idx, segment, acsi, unique=True)
                        chindices[idx] = chindex
                        _dict[s_c] = chain
                        segment._dict[pc] = chain
                        _chains.append(chain)
                    else:
                        idx = _indices[_i:i]
                        chindices[idx] = chain._indices[0]
                        chain._indices = np.concatenate((chain._indices, idx))
                    pc = c
                    ps = s
                    _i = i
                s_c = (ps, pc or None)
                chain = _dict.get(s_c)
                idx = _indices[_i:]
                if chain is None:
                    segment = _dict[ps]
                    chindex += 1
                    chindices[idx] = chindex
                    chain = Chain(ag, idx, segment, acsi, unique=True)
                    _dict[s_c] = chain
                    segment._dict[pc] = chain
                    _chains.append(chain)
                else:
                    chindices[idx] = chain._indices[0]
                    chain._indices = np.concatenate((chain._indices, idx)) 
        
        if kwargs.get('chain') == True:
            return
    
        rnums = ag._getResnums()
        if rnums is None:
            raise ValueError('resnums are not set')
        if selstr:
            rnums = rnums[_indices]
        nones = None
        if _segments is None:
            if nones is None:
                nones = [None] * len(rnums)
            sgnms = nones
        if _chains is None:
            if nones is None:
                nones = [None] * len(rnums)
            chids = nones
        icods = ag._getIcodes()
        if icods is None:
            if nones is None:
                nones = [None] * len(rnums)
            icods = nones
        elif selstr:
                icods = icods[_indices]

        pr = rnums[0]
        pi = icods[0] or None          
        pc = chids[0]
        ps = sgnms[0] or None
        _j = 0
        for j, r in enumerate(rnums):
            i = icods[j] or None
            c = chids[j] or None
            s = sgnms[j]
            if r != pr or i != pi or c != pc or s != ps:
                s_c_r_i = (ps, pc, pr, pi)
                res = _dict.get(s_c_r_i)
                if res is None:
                    chain = _dict.get((ps, pc))
                    resindex += 1
                    idx = _indices[_j:j]
                    res = Residue(ag, idx, chain, acsi, unique=True, 
                                  selstr=selstr)
                    resindices[idx] = resindex
                    if chain is not None:
                        chain._dict[(pr, pi)] = len(chain._list)
                        chain._list.append(res)
                    _residues.append(res)
                    _dict[s_c_r_i] = res
                else:
                    res._indices = np.concatenate((res._indices, 
                                                   _indices[_j:j]))
                ps = s
                pc = c
                pr = r
                pi = i
                _j = j 
        s_c_r_i = (ps, pc, pr, pi)
        res = _dict.get(s_c_r_i)
        idx = _indices[_j:]
        if res is None:
            chain = _dict.get((ps, pc))
            resindex += 1
            res = Residue(ag, idx, chain, acsi, unique=True, selstr=selstr)
            resindices[idx] = resindex
            if chain is not None:
                chain._dict[(pr, pi)] = len(chain._list)
                chain._list.append(res)
            _residues.append(res)
            _dict[s_c_r_i] = res
        else:
            resindices[idx] = res._indices[0]
            res._indices = np.concatenate((res._indices, idx))
        
        ag._data['segindices'] = segindices
        ag._data['chindices'] = chindices
        ag._data['resindices'] = resindices

    def getResidue(self, chid, resnum, icode=None, segname=None):
        """Return residue with number *resnum* and insertion code *icode* from 
        the chain with identifier *chid* in segment with name *segname*."""
        
        return self._dict.get((segname or None, chid or None, 
                               resnum, icode or None))

    def numResidues(self):
        """Return number of residues."""
        
        return len(self._residues)    

    def iterResidues(self):
        """Iterate over residues."""
        
        return self._residues.__iter__()
                
    def getChain(self, chid, segname=None):
        """Return chain with identifier *chid*, if it exists."""
        
        return self._dict.get((segname or None, chid or None))

    def iterChains(self):
        """Iterate over chains."""

        return self._chains.__iter__()
    
    def numChains(self):
        """Return number of chains."""
        
        return len(self._chains)

    def getSegment(self, segname):
        """Return segment with name *segname*, if it exists."""
        
        return self._segments.get(segname or None)

    def numSegments(self):
        """Return number of chains."""
        
        return len(self._segments)

    def iterSegments(self):
        """Iterate over segments."""

        return self._segments.__iter__()

    def getNumOfResidues(self):
        """Deprecated, use :meth:`numResidues`."""
        
        prody.deprecate('getNumOfResidues', 'numResidues')
        return self.numResidues()
        
    def getNumOfChains(self):
        """Deprecated, use :meth:`numChains`."""
        
        prody.deprecate('getNumOfChains', 'numChains')
        return self.numChains()

class Bond(object):
    
    """A pointer class for bonded atoms."""
    
    __slots__ = ['_ag', '_acsi', '_indices']
    
    def __init__(self, ag, indices, acsi):
        
        self._ag = ag
        self._indices = indices
        self._acsi = acsi

    def __repr__(self):
        
        one, two = self._indices
        names = self._ag._getNames()
        return '<Bond: {0:s}({1:d})--{2:s}({3:d}) from {4:s}>'.format(
                            names[one], one, names[two], two, str(self._ag))
    
    def __str__(self):

        one, two = self._indices
        names = self._ag._getNames()
        return '{0:s}({1:d})--{2:s}({3:d})'.format(
                                            names[one], one, names[two], two)

    def getACSIndex(self):
        """Set active coordinate set index."""

        return self._acsi
        
    def setACSIndex(self):
        """Set the coordinate set at *index* active."""

        if self._ag._coords is None:
            raise AttributeError('coordinates are not set')
        if not isinstance(index, int):
            raise TypeError('index must be an integer')
        if self._ag._n_csets <= index or \
           self._ag._n_csets < abs(index):
            raise IndexError('coordinate set index is out of range')
        if index < 0:
            index += self._ag._n_csets
        self._acsi = index
    
    def getAtomGroup(self):
        """Return atom group."""
        
        return self._ag

    def getAtoms(self):
        """Return bonded atoms."""
        
        return (self._ag[self._indices[0]], self._ag[self._indices[1]])

    def getIndices(self):
        """Return indices of bonded atoms."""
        
        return self._indices.copy()
    
    def getLength(self):
        
        vector = self.getVector()
        return np.multiply(vector, vector, vector).sum() ** 0.5

    def getVector(self):
        """Return bond vector that originates from the first atom."""
        
        one, two = self._indices
        acsi = self._acsi
        return self._ag._coords[acsi, two] - self._ag._coords[acsi, one] 

def evalBonds(bonds, n_atoms):
    numbonds = np.bincount(bonds.reshape((bonds.shape[0] * 2)))
    bmap = np.zeros((n_atoms, numbonds.max()), int)
    bmap.fill(-1)
    index = np.zeros(n_atoms, int)
    for bond in bonds:
        a, b = bond
        bmap[a, index[a]] = b
        bmap[b, index[b]] = a
        index[bond] += 1
    return bmap, numbonds

def trimBonds(bonds, bmap, indices):
    
    newindices = np.zeros(indices.max()+1, int)
    newindices[indices] = np.arange(len(indices))
    iset = set(indices)
    newbonds = []
    for i, bond in enumerate(bonds):
        if bond[0] in iset and bond[1] in iset:
            newbonds.append(newindices[bond])

    newbonds = np.array(newbonds)
    if len(newbonds) > 0:
        return newbonds
    

def saveAtoms(atoms, filename=None, **kwargs):
    """Save *atoms* in ProDy internal format.  All classes derived from 
    :class:`Atomic` are accepted as *atoms* argument.
    
    .. versionadded:: 0.7.1
    
    This function saves user set atomic attributes as well.  Note that name of 
    the :class:`AtomGroup` instance is used as the filename when *atoms* is not
    an :class:`AtomGroup`.  This is because names for selections and atom maps
    may be too long and may contain special characters.  To avoid overwriting 
    an existing file with the same name, specify a *filename*."""
    
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be Atomic instance, not {0:s}'
                        .format(type(atoms)))
    if isinstance(atoms, AtomGroup):
        ag = atoms
        title = ag.getTitle()
    else:
        ag = atoms.getAtomGroup()
        title = str(atoms)
    
    if filename is None:
        filename = ag.getTitle().replace(' ', '_')
    filename += '.ag.npz'
    attr_dict = {'title': title}
    attr_dict['n_atoms'] = atoms.numAtoms()
    attr_dict['n_csets'] = atoms.numCoordsets()
    attr_dict['cslabels'] = atoms.getCSLabels()
    coords = atoms._getCoordsets()
    if coords is not None:
        attr_dict['coordinates'] = coords
    bonds = ag._bonds
    bmap = ag._bmap
    if bonds is not None and bmap is not None:
        if isinstance(atoms, AtomGroup):
            attr_dict['bonds'] = bonds
            attr_dict['bmap'] = bmap
            attr_dict['numbonds'] = ag._data['numbonds']
        else:
            bonds = trimBonds(bonds, bmap, atoms._getIndices())
            attr_dict['bonds'] = bonds
            attr_dict['bmap'], attr_dict['numbonds'] = evalBonds(bonds)
        
    for key, data in ag._data.iteritems():
        if key == 'numbonds':
            continue
        if data is not None:
            attr_dict[key] = data 
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename

SKIP = set(['_name', '_title', 'title', 'n_atoms', 'n_csets', 'bonds', 'bmap',
            'coordinates', '_coordinates', 'cslabels', 'numbonds'])

def loadAtoms(filename):
    """Return :class:`AtomGroup` instance from *filename*.  This function makes
    use of :func:`numpy.load` function.  See also :func:`saveAtoms`.
    
    .. versionadded:: 0.7.1"""
    
    LOGGER.timeit()
    attr_dict = np.load(filename)
    files = set(attr_dict.files)
    # REMOVE support for _coordinates IN v1.0
    if not '_coordinates' in files and not 'n_atoms' in files:
        raise ValueError("'{0:s}' is not a valid atomic data file"
                         .format(filename))
    if '_coordinates' in files:
        ag = AtomGroup(str(attr_dict['_name']))
        for attr in attr_dict.files:
            if attr == '_name':
                continue
            elif attr == '_coordinates':
                data = attr_dict[attr]
                if data.ndim > 0:
                    ag.setCoords(data)
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
                    ag.setData(attr, data)
                else:
                    ag.setData(attr, [data])
    else:        
        ag = AtomGroup(str(attr_dict['title']))
        if 'coordinates' in files:
            ag._coords = attr_dict['coordinates']
        ag._n_atoms = int(attr_dict['n_atoms'])
        ag._n_csets = int(attr_dict['n_csets'])
        ag._setTimeStamp()
        if 'bonds' in files and 'bmap' in files and 'numbonds' in files:
            ag._bonds = attr_dict['bonds']
            ag._bmap = attr_dict['bmap']
            ag._data['numbonds'] = attr_dict['numbonds']
        for key, data in attr_dict.iteritems():
            if key in SKIP:
                continue
            if key in ATOMIC_ATTRIBUTES:
                ag._data[key] = data
            else:
                ag.setData(key, data)
        if ag.numCoordsets() > 0:
            ag._acsi = 0
        if 'cslabels' in files:
            ag.setCSLabels(list(attr_dict['cslabels']))
    LOGGER.timing('Atom group was loaded in %.2fs.')
    return ag

if __name__ == '__main__':
    from prody import *
    p = parsePDB('1aar')
    saveAtoms(p)
