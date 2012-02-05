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

import atomic
from atomic import *
__all__ = atomic.__all__
