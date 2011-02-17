.. currentmodule:: prody.select

.. _selection-examples:

*******************************************************************************
Atom Selection Examples
*******************************************************************************

Synopsis
===============================================================================

One of the most powerful features of ProDy its atom selection engine.
It enables users to select well defined atom subsets easily by passing
simple keywords or make sophisticated selections by using composite keyword
arguments. This example shows how to make selections.

User Input
-------------------------------------------------------------------------------

Structure of a protein in PDB format or simply as PDB id.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Using keywords
-------------------------------------------------------------------------------

>>> prot = parsePDB('1p38')
>>> prot
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

Most of the time single word keywords may be enough for you to get the set
of atoms that you want:

For selecting :term:`protein` atoms:

>>> protein = prot.select('protein')
>>> protein
<Selection: "protein" from 1p38 (2833 atoms; 1 coordinate sets, active set index: 0)>

The above shows, 2833 of 2962 atoms are protein atoms. 

For selecting Cα atoms you can use :term:`calpha`:

>>> calpha = prot.select('calpha')
>>> calpha 
<Selection: "calpha" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)>

The above shows that there are 351 amino acid residues.

Or, for selecting :term:`backbone` atoms:

>>> backbone = prot.select('backbone')
>>> backbone
<Selection: "backbone" from 1p38 (1404 atoms; 1 coordinate sets, active set index: 0)>

You can also combine these selection keywords:

>>> calpha_water = prot.select('calpha or water')
>>> calpha_water
<Selection: "calpha or water" from 1p38 (480 atoms; 1 coordinate sets, active set index: 0)>

For a list of keywords see :ref:`selections`. 

Select atoms by name
-------------------------------------------------------------------------------

We can select side-chain atoms as follows:

>>> side_chain_atoms = prot.select('protein and not name N CA C O')
>>> print side_chain_atoms
Selection "protein and not name N CA C O" from 1p38
>>> len(side_chain_atoms)
1429

Same selection could also be made using :term:`sidechain` keyword or 
:term:`backbone` keyword preceded by *not*:

>>> prot.select('sidechain')
<Selection: "sidechain" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

>>> prot.select('not backbone')
<Selection: "not backbone" from 1p38 (1558 atoms; 1 coordinate sets, active set index: 0)>

Oops, :term:`not` :term:`backbone` did not select the same number of atoms. 
Let's try to see why:

>>> print set(prot.select('not backbone').getResidueNames())
set(['CYS', 'ILE', 'VAL', 'GLN', 'LYS', 'HOH', 'PRO', 'THR', 'PHE', 'ASN', 'HIS', 'MET', 'ASP', 'LEU', 'ARG', 'TRP', 'ALA', 'GLU', 'TYR', 'SER'])

Note that we used built-in Python type :class:`set`.

As you can see atoms of **HOH** residues are also included in the selection.

Let's try:

>>> prot.select('not backbone and not water')
<Selection: "not backbone and not water" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

We also used :term:`water` term. This has now worked as :term:`sidechain` did.
This was to show that it is possible to select same set of atoms in a number 
of different ways. 

Select amino acids by type/name
-------------------------------------------------------------------------------

Let's say we want to select charged residues. We can use :term:`resname`
keyword followed by 3-letter residue names:  

>>> prot.select('resname ARG LYS HIS ASP GLU')
<Selection: "resname ARG LYS HIS ASP GLU" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

Or, we can use predefined keywords :term:`acidic` and :term:`basic`.

>>> charged = prot.select('acidic or basic')
>>> print charged
Selection "acidic or basic" from 1p38
>>> len(charged)
906
>>> set(charged.getResidueNames())
set(['HIS', 'ASP', 'LYS', 'GLU', 'ARG'])

Same selection could also be made using :term:`charged` keyword:

>>> prot.select('charged')
<Selection: "charged" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

More examples
-------------------------------------------------------------------------------

There is much more to what you can do with this flexible and fast atom 
selection engine, without the need for writing nested loops with comparisons 
or changing the source code. See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of :class:`~atomic.Selection`
    objects
  * :ref:`contacts` for selecting interacting atoms
