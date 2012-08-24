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

"""This module defines classes for handling atomic data.  Read this page using
``help(atomic)``.

.. _atomic:

Atomic classes
===============================================================================

ProDy stores atomic data in instances of :class:`.AtomGroup` class, which 
supports multiple coordinate sets, e.g. models from an NMR structure or 
snapshots from a molecular dynamics trajectory.
 
Instances of the class can be obtained by parsing a PDB file as follows:
    
>>> from prody import *
>>> ag = parsePDB('1aar')
>>> ag
<AtomGroup: 1aar (1218 atoms)>

In addition to :class:`.AtomGroup` class, following classes that act as 
pointers provide convenient access subset of data:

* :class:`.Selection` - Points to an arbitrary subset of atoms. See 
  :ref:`selections` and :ref:`selection-operations` for usage examples.

* :class:`.Segment` - Points to atoms that have the same segment name.

* :class:`.Chain` - Points to atoms in a segment that have the same chain 
  identifier.

* :class:`.Residue` - Points to atoms in a chain that have the same 
  residue number and insertion code.
                      
* :class:`.AtomMap` - Points to arbitrary subsets of atoms while 
  allowing for duplicates and missing atoms.  Indices of atoms are stored 
  in the order provided by the user.
    
* :class:`.Atom` - Points to a single atom

* :class:`.Bond` - Points to two connected atoms

Atom data fields
===============================================================================

:ref:`fields` defines an interface for handling data parsed from molecular
data files, in particular PDB files.  Aforementioned classes offer ``get``
and ``set`` functions for manipulating this data.  For example, the following 
prints residue names:

>>> print(ag.getResnames())
['MET' 'MET' 'MET' ..., 'HOH' 'HOH' 'HOH']

Atom flags
===============================================================================

:ref:`flags` module defines a way to mark atoms with certain properties, such
as atoms that are part of a **protein**.  Following example checks whether
all atoms of *ag* are protein atoms: 
    
>>> ag.isprotein
False

This indicates that there are some non-protein atoms, probably water atoms. 
We can easily make a count as follows:

>>> ag.numAtoms('protein')
1203
>>> ag.numAtoms('hetero')
15
>>> ag.numAtoms('water')
15


Atom selections
===============================================================================

:ref:`selections` offer a flexible and powerful way to access subsets of 
selections and is one of the most important features of ProDy.   The details 
of the selection grammar is described in :ref:`selections`.  Following examples
show how to make quick selections using the overloaded ``.`` operator:
    
>>> ag.chain_A # selects chain A
<Selection: 'chain A' from 1aar (608 atoms)>
>>> ag.calpha # selects alpha carbons
<Selection: 'calpha' from 1aar (152 atoms)>
>>> ag.resname_ALA # selects alanine residues
<Selection: 'resname ALA' from 1aar (20 atoms)>

It is also possible to combine selections with ``and`` and ``or`` operators:

>>> ag.chain_A_and_backbone
<Selection: 'chain A and backbone' from 1aar (304 atoms)>
>>> ag.acidic_or_basic
<Selection: 'acidic or basic' from 1aar (422 atoms)>

Using dot operator will behave like the logical ``and`` operator:
    
>>> ag.chain_A.backbone
<Selection: '(backbone) and (chain A)' from 1aar (304 atoms)>
  
For this to work, the first word following the dot operator must be a flag 
label or a field name, e.g. ``resname``, ``name``, ``apolar``, ``protein``, 
etc.  Underscores will be interpreted as white space, as obvious from the
previous examples.  The limitation of this is that parentheses, special 
characters cannot be used.

Functions
===============================================================================

Following functions can be used for permanent data storage:
    
  * :func:`.loadAtoms`
  * :func:`.saveAtoms`
  
Following function can be used to identify fragments in a group 
(:class:`.AtomGroup`) or subset (:class:`.Selection`) of atoms:
    
  * :func:`.findFragments`
  * :func:`.iterFragments`
  
Following function can be used check whether a word is reserved because
it is used internally by :mod:`.prody.atomic` classes:

  * :func:`.isReserved`
  * :func:`.getReservedWords`"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = ['Atomic', 'AtomGroup', 
           'HierView', 'Segment', 'Chain', 'Residue', 'Atom',
           'AtomPointer', 'AtomSubset',
           'Selection', 'AtomMap',
           'Bond', 'select', 'atomgroup', 'hierview', 'fields', 'flags']

from fields import ATOMIC_FIELDS

from atom import *
from bond import *
from flags import *
from chain import *
from subset import *
from atomic import *
from select import *
from atommap import *
from residue import *
from pointer import *
from segment import *
from hierview import *
from functions import *
from atomgroup import *
from selection import *

import flags
import atomic
import select
import atommap
import pointer
import hierview
import functions
import atomgroup
import selection

from atomgroup import SELECT
from chain import AAMAP, getSequence

__all__.extend(functions.__all__)
__all__.extend(select.__all__)
__all__.extend(flags.__all__)

from functions import isAtomic, isSubset

from select import checkSelstr, isKeyword, isSelectionMacro

atomic.SELECT = atomgroup.SELECT = selection.SELECT = SELECT = Select()
atomic.isSelectionMacro = isSelectionMacro
atomic.isKeyword = isKeyword
atomic.AtomMap = AtomMap
atomic.AtomGroup = AtomGroup
atomic.Selection = Selection

atomgroup.isReserved = isReserved
atomgroup.HierView = HierView

pointer.atommap = atommap
pointer.AtomMap = AtomMap
pointer.AtomGroup = AtomGroup
pointer.Selection = Selection

select.flags = flags
select.isReserved = isReserved
select.HierView = HierView

selection.HierView = HierView

import numpy as np

n_atoms = 10
ATOMGROUP = AtomGroup('Test')
ATOMGROUP.setCoords(np.random.random((n_atoms,3)))
ATOMGROUP.setNames(['CA']*n_atoms)
ATOMGROUP.setResnames(['GLY']*n_atoms)
ATOMGROUP.setResnums(np.arange(1,n_atoms+1))
ATOMGROUP.setChids(['A']*n_atoms)
ATOMGROUP.setAltlocs([' ']*n_atoms)
ATOMGROUP.setElements(['C']*n_atoms)
ATOMGROUP.setFlags('hetatm', [False]*n_atoms)
ATOMGROUP.setOccupancies([1]*n_atoms)
ATOMGROUP.setSecstrs(['H']*n_atoms)
ATOMGROUP.setSegnames(['PDB']*n_atoms)
ATOMGROUP.setAnisous(np.random.random((n_atoms,6)))
ATOMGROUP.setAnistds(np.random.random((n_atoms,6)))
ATOMGROUP.setIcodes([' ']*n_atoms)
ATOMGROUP.setTypes(['CH2']*n_atoms)
ATOMGROUP.setBetas([0]*n_atoms)
ATOMGROUP.setCharges([0]*n_atoms)
ATOMGROUP.setMasses([12]*n_atoms)
ATOMGROUP.setRadii([1.4]*n_atoms)

select.ATOMGROUP = ATOMGROUP
