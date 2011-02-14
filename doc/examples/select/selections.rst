.. currentmodule:: prody.select

.. _selection-examples:

*******************************************************************************
Atom Selection Examples
*******************************************************************************

Synopsis
===============================================================================

Coming soon.

User Input
-------------------------------------------------------------------------------


ProDy Code
===============================================================================


**Select atoms by name**:

We can select side-chain atoms as follows:

>>> from prody import *
>>> prot = parsePDB('1p38')

>>> side_chain_atoms = prot.select('protein and not name N CA C O')
>>> print side_chain_atoms
Selection "protein and not name N CA C O" from 1p38
>>> len(side_chain_atoms)
1429

Same selection could also be made using :term:`sidechain` keyword or :term:`backbone` keyword
preceded by *not*:

>>> prot.select('sidechain')
<Selection: "sidechain" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

>>> prot.select('not backbone')
<Selection: "not backbone" from 1p38 (1558 atoms; 1 coordinate sets, active set index: 0)>

Oops, *not* :term:`backbone` did not select the same number of atoms. Let's try to
see why:

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

**Select amino acids by type/name**:

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

.. seealso::
   To see all of what you can do with atom selections, 
   go to :ref:`selections`. Also see :ref:`selection-operations` and 
   :ref:`contacts` for more advanced examples.

Full list of selection keywords 
are given in section :ref:`selections`. 
