# -*- coding: utf-8 -*-
"""This module defines classes for handling atomic data. Read this page using
``help(atomic)``.

.. _atomic:

Atomic classes
^^^^^^^^^^^^^^

ProDy stores atomic data in instances of :class:`.AtomGroup` class, which
supports multiple coordinate sets, e.g. models from an NMR structure or
snapshots from a molecular dynamics trajectory.

Instances of the class can be obtained by parsing a PDB file as follows:

.. ipython:: python

   from prody import *
   ag = parsePDB('1aar')
   ag


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
^^^^^^^^^^^^^^^^

:ref:`fields` defines an interface for handling data parsed from molecular
data files, in particular PDB files.  Aforementioned classes offer ``get``
and ``set`` functions for manipulating this data.  For example, the following
prints residue names:

.. ipython:: python

   ag.getResnames()

Atom flags
^^^^^^^^^^

:ref:`flags` module defines a way to mark atoms with certain properties, such
as atoms that are part of a **protein**.  Following example checks whether
all atoms of *ag* are protein atoms:

.. ipython:: python

   ag.isprotein

This indicates that there are some non-protein atoms, probably water atoms.
We can easily make a count as follows:

.. ipython:: python

   ag.numAtoms('protein')
   ag.numAtoms('hetero')
   ag.numAtoms('water')


Atom selections
^^^^^^^^^^^^^^^

:ref:`selections` offer a flexible and powerful way to access subsets of
selections and is one of the most important features of ProDy.   The details
of the selection grammar is described in :ref:`selections`.  Following examples
show how to make quick selections using the overloaded ``.`` operator:

.. ipython:: python

   ag.chain_A  # selects chain A
   ag.calpha  # selects alpha carbons
   ag.resname_ALA  # selects alanine residues

It is also possible to combine selections with ``and`` and ``or`` operators:

.. ipython:: python

   ag.chain_A_and_backbone
   ag.acidic_or_basic

Using dot operator will behave like the logical ``and`` operator:

.. ipython:: python

   ag.chain_A.backbone

For this to work, the first word following the dot operator must be a flag
label or a field name, e.g. ``resname``, ``name``, ``apolar``, ``protein``,
etc.  Underscores will be interpreted as white space, as obvious from the
previous examples.  The limitation of this is that parentheses, special
characters cannot be used.

Functions
^^^^^^^^^

Following functions can be used for permanent data storage:

  * :func:`.loadAtoms`
  * :func:`.saveAtoms`

Following function can be used to identify fragments in a group
(:class:`.AtomGroup`) or subset (:class:`.Selection`) of atoms:

  * :func:`.findFragments`
  * :func:`.iterFragments`

Following function can be used to get an :class:`.AtomMap` that sorts atoms
based on a given property:

  * :func:`.sortAtoms`

Following function can be used check whether a word is reserved because
it is used internally by :mod:`.prody.atomic` classes:

  * :func:`.isReserved`
  * :func:`.listReservedWords`"""

import prody
LOGGER = prody.LOGGER
SETTINGS = prody.SETTINGS

__all__ = ['Atomic', 'AtomGroup',
           'HierView', 'Segment', 'Chain', 'Residue', 'Atom',
           'AtomPointer', 'AtomSubset',
           'Selection', 'AtomMap',
           'Bond', 'select', 'atomgroup', 'hierview', 'chain', 'fields', 'flags']

from .fields import ATOMIC_FIELDS

from .atom import *
from .bond import *
from .flags import *
from .chain import *
from .subset import *
from .atomic import *
from .select import *
from .atommap import *
from .residue import *
from .pointer import *
from .segment import *
from .hierview import *
from .functions import *
from .atomgroup import *
from .selection import *

from . import flags
from . import atomic
from . import select
from . import atommap
from . import pointer
from . import hierview
from . import functions
from . import atomgroup
from . import selection
from . import chain
from . import segment

from .chain import getSequence

__all__.extend(functions.__all__)
__all__.extend(select.__all__)
__all__.extend(flags.__all__)

from .select import checkSelstr, isSelectionMacro

atomic.SELECT = selection.SELECT = SELECT = Select()
atomic.isSelectionMacro = isSelectionMacro
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

chain.HierView = HierView

segment.HierView = HierView

import numpy as np

n_atoms = 10
ATOMGROUP = AtomGroup('Test')
ATOMGROUP.setCoords(np.random.random((n_atoms, 3)))
ATOMGROUP.setNames(['CA']*n_atoms)
ATOMGROUP.setResnames(['GLY']*n_atoms)
ATOMGROUP.setResnums(np.arange(1, n_atoms+1))
ATOMGROUP.setChids(['A']*n_atoms)
ATOMGROUP.setAltlocs([' ']*n_atoms)
ATOMGROUP.setElements(['C']*n_atoms)
ATOMGROUP.setFlags('hetatm', [False]*n_atoms)
ATOMGROUP.setOccupancies([1]*n_atoms)
ATOMGROUP.setSecstrs(['H']*n_atoms)
ATOMGROUP.setSegnames(['PDB']*n_atoms)
ATOMGROUP.setAnisous(np.random.random((n_atoms, 6)))
ATOMGROUP.setAnistds(np.random.random((n_atoms, 6)))
ATOMGROUP.setIcodes([' ']*n_atoms)
ATOMGROUP.setTypes(['CH2']*n_atoms)
ATOMGROUP.setBetas([0]*n_atoms)
ATOMGROUP.setCharges([0]*n_atoms)
ATOMGROUP.setMasses([12]*n_atoms)
ATOMGROUP.setRadii([1.4]*n_atoms)

select.ATOMGROUP = ATOMGROUP
