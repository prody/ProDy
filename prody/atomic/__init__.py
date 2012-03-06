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
    
    
Atomic classes
===============================================================================

ProDy stores atomic data in instances of :class:`~.AtomGroup` class, which 
supports multiple coordinate sets, e.g. models from an NMR structure or 
snapshots from a molecular dynamics trajectory.
 
Instances of the class can be obtained by parsing a PDB file as follows:
    
>>> from prody import *
>>> ag = parsePDB('1aar')
>>> ag
<AtomGroup: 1aar (1218 atoms)>

All atomic data in :class:`~.AtomGroup` instances and comes
with other classes acting as pointers to provide convenient read/write access 
to such data.  These classes are:

* :class:`~.Atom` - Points to a single atom in an :class:`~.AtomGroup` 
  instance.                          

* :class:`~.Selection` - Points to an arbitrary subset of atoms. See 
  :ref:`selections` and :ref:`selection-operations` for usage examples.

* :class:`~.Segment` - Points to atoms that have the same segment name.

* :class:`~.Chain` - Points to atoms in a segment that have the same chain 
  identifier.

* :class:`~.Residue` - Points to atoms in a chain that have the same 
  residue number and insertion code.
                      
* :class:`~.AtomMap` - Points to arbitrary subsets of atoms while 
  allowing for duplicates and missing atoms.  Indices of atoms are stored 
  in the order provided by the user.
    
Atom selections
===============================================================================

Flexible and powerful atom selections is one of the most important features 
of ProDy.  The details of the selection grammar is described in 
:ref:`selections`. 

Using the flexibility of Python, atom selections are made much easier by
overriding the ``.`` operator, so the following are interpreted as selections:
    
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
  
For this to work, the first word following the dot operator must be a selection
keyword, e.g. ``resname``, ``name``, ``apolar``, ``protein``, etc. 
Underscores will be interpreted as white space, as obvious from the
previous examples.  The limitation of this is that parentheses, special 
characters cannot be used.     

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from fields import ATOMIC_ATTRIBUTES, ATOMIC_FIELDS

"""
Behavioral differences
-------------------------------------------------------------------------------

Atomic classes behave differently to indexing and to calls of certain built-in 
functions.  These differences are:

=========  ====================================================================
Class               Properties and differences
=========  ====================================================================
                       
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

"""

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = ['Atomic', 'AtomGroup', 
           'HierView', 'Segment', 'Chain', 'Residue', 'Atom',
           'AtomPointer', 'AtomSubset',
           'Selection', 'AtomMap',
           'Bond', 'loadAtoms', 'saveAtoms',
           'select', 'atomgroup', 'hierview']

from fields import ATOMIC_FIELDS

from atom import *
from bond import *
from chain import *
from subset import *
from atomic import *
from select import *
from atommap import *
from residue import *
from pointer import *
from segment import *
from hierview import *
from atomgroup import *
from selection import *

import atomic
import select
import pointer
import hierview
import atomgroup
import selection

from atomgroup import SELECT

from chain import AAA2A, getSequence

atomgroup.HierView = HierView
selection.HierView = HierView

#from atom import Atom
from functions import loadAtoms, saveAtoms

__all__ += select.__all__

from select import isMacro, isKeyword, isReserved
atomic.isMacro = isMacro
atomic.isKeyword = isKeyword
atomic.SELECT = atomgroup.SELECT = selection.SELECT = SELECT = Select()
atomgroup.isReserved = isReserved
pointer.AtomMap = AtomMap
pointer.AtomGroup = AtomGroup
subset.Selection = Selection

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
ATOMGROUP.setHeteros([False]*n_atoms)
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
