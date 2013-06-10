.. _hierview:

Hierarchical Views
===============================================================================

This part describes how to use hierarchical views.  We start by importing
everything from the ProDy package:

.. ipython:: python

   from prody import *

Hierarchical Views
-------------------------------------------------------------------------------

Then we parses a structure to get an :class:`.AtomGroup` instance which has a
plain view of atoms:

.. ipython:: python

   structure = parsePDB('3mkb')
   structure

A hierarchical view of the structure can be simply get by calling the
:meth:`.AtomGroup.getHierView` method:

.. ipython:: python

   hv = structure.getHierView()
   hv


Indexing
^^^^^^^^

Indexing :class:`.HierView` instances return :class:`.Chain`:

.. ipython:: python

   hv['A']
   hv['B']
   hv['Z'] # This will return None, which means chain Z does not exist

The length of the *hv* variable gives the number of chains in the structure:

.. ipython:: python

   len(hv)
   hv.numChains()

It is also possible to get a :class:`.Residue` by
directly indexing the :class:`.HierView` instance:

.. ipython:: python

   hv['A', 100]

Insertion codes can also be passed:

.. ipython:: python

   hv['A', 100, 'B']

But this does not return anything, since residue 100B does not exist.


Iterations
^^^^^^^^^^

One can iterate over :class:`.HierView` instances to get chains:

.. ipython:: python

   for chain in hv:
       chain

It is also possible to get a :func:`list` of chains simply as follows:

.. ipython:: python

   chains = list(hv)
   chains


Residues
^^^^^^^^

In addition, one can also iterate over all residues:

.. ipython:: python

   for i, residue in enumerate(hv.iterResidues()):
       if i == 4: break
       print(residue)


Chains
-------------------------------------------------------------------------------


.. ipython:: python

   chA = hv['A']
   chA

Length of the chain equals to the number of residues in it:

.. ipython:: python

   len(chA)
   chA.numResidues()


Indexing
^^^^^^^^

Indexing a :class:`.Chain` instance returns a :class:`.Residue` instance.

.. ipython:: python

   chA[1]

If a residue does not exist, ``None`` is returned:

.. ipython:: python

   chA[1000]
   chA[1, 'A'] # Residue 1 with insertion code A also does not exist

If residue with given integer number does not exist, ``None`` is returned.


Iterations
^^^^^^^^^^

Iterating over a chain yields residues:

.. ipython:: python

   for i, residue in enumerate(chA):
       if i == 4: break
       print(residue)

Note that water atoms, each constituting a residue, are also part of a chain
if they are labeled with that chain's identifier.

This enables getting a :func:`list` of residues simply as follows:

.. ipython:: python

   chA_residues = list(chA)
   chA_residues[:4]
   chA_residues[-4:]


Get data
^^^^^^^^

All methods defined for :class:`.AtomGroup` class are also defined for
:class:`.Chain` and :class:`.Residue` classes:

.. ipython:: python

   chA.getCoords()
   chA.getBetas()


Selections
^^^^^^^^^^

Finally, you can select atoms from a :class:`.Chain` instance:

.. ipython:: python

   chA_backbone = chA.select('backbone')
   chA_backbone
   chA_backbone.getSelstr()

As you see, the selection string passed by the user is augmented with
"chain" keyword and identifier automatically to provide internal
consistency:

.. ipython:: python

   structure.select(chA_backbone.getSelstr())

Residues
-------------------------------------------------------------------------------

.. ipython:: python

   chA_res1 = chA[1]
   chA_res1


Indexing
^^^^^^^^

:class:`.Residue` instances can be indexed to get individual atoms:

.. ipython:: python

   chA_res1['CA']
   chA_res1['CB']
   chA_res1['X'] # if atom does not exist, None is returned


Iterations
^^^^^^^^^^

Iterating over a residue instance yields :class:`Atom` instances:

.. ipython:: python

   for i, atom in enumerate(chA_res1):
       if i == 4: break
       print(atom)

This makes it easy to get a :func:`list` of atoms:

.. ipython:: python

   list(chA_res1)


Get data
^^^^^^^^

All methods defined for :class:`.AtomGroup` class are also defined for
:class:`.Residue` class:

.. ipython:: python

   chA_res1.getCoords()
   chA_res1.getBetas()


Selections
^^^^^^^^^^

Finally, you can select atoms from a :class:`.Residue` instance:

.. ipython:: python

   chA_res1_bb = chA_res1.select('backbone')
   chA_res1_bb
   chA_res1_bb.getSelstr()

Again, the selection string is augmented with the chain identifier and
residue number (:term:`resnum`).


Atoms
-------------------------------------------------------------------------------

The lowest level of the hierarchical view contains :class:`Atom` instances.

.. ipython:: python

   chA_res1_CA = chA_res1['CA']
   chA_res1_CA

*Get atomic data*

All methods defined for :class:`.AtomGroup` class are also defined for
:class:`.Atom` class with the difference that method names are singular
(except for coordinates):

.. ipython:: python

   chA_res1_CA.getCoords()
   chA_res1_CA.getBeta()


State Changes
-------------------------------------------------------------------------------

A :class:`.HierView` instance represents the state of an :class:`.AtomGroup`
instance at the time it is built.  When chain identifiers or residue numbers
change, the state that hierarchical view represents may not match the current
state of the atom group:

.. ipython:: python

   chA.setChid('X')
   chA
   hv['X'] # returns None, since hierarchical view is not updated
   hv.update() # this updates hierarchical view
   hv['X']

When this is the case, :meth:`.HierView.update` method can be used to update
hierarchical view.