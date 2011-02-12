.. currentmodule:: prody.select

.. _selection-operations:

*******************************************************************************
Operations on Selections 
*******************************************************************************



Before reading this section, familiarity with :ref:`selections` can be 
helpful.

Let's import all classes and functions from ProDy and parse coordinates from
the PDB structure 1p38:

>>> from prody import *
>>> prot = parsePDB('1p38')


Union
===============================================================================

Let's select β-carbon atoms for non-GLY amino acid residues, and 
α-carbons for GLYs in two steps:

>>> betas = prot.select('name CB and protein')
>>> print len(betas)
336
>>> gly_alphas = prot.select('name CA and resname GLY')
>>> print len(gly_alphas)
15

The above shows that the p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

>>> betas_gly_alphas = betas | gly_alphas
>>> print betas_gly_alphas
Selection "(name CB and pr...nd resname GLY)" from 1p38
>>> print len(betas_gly_alphas)
351

The selection string for the union of selections becomes:

>>> print betas_gly_alphas.getSelectionString()
(name CB and protein) or (name CA and resname GLY)

Note that it is also possible to yield the same selection using selection 
string ``(name CB and protein) or (name CA and resname GLY)``.


Intersection
===============================================================================

It is as easy to get the intersection of two selections. Let's find 
:term:`charged` and :term:`medium` size residues in a protein:

>>> charged = prot.select('charged')
>>> print charged
Selection "charged" from 1p38
>>> medium = prot.select('medium')
>>> print medium
Selection "medium" from 1p38

>>> medium_charged = medium & charged
>>> print medium_charged
Selection "(medium) and (charged)" from 1p38
>>> print medium_charged.getSelectionString()
(medium) and (charged)

Let's see which amino acids are considered :term:`charged` and :term:`medium`:

>>> print set(medium_charged.getResidueNames())
set(['ASP'])

What about amino acids that are :term:`medium` or :term:`charged`:

>>> print set((medium | charged).getResidueNames())
set(['CYS', 'ASP', 'VAL', 'LYS', 'PRO', 'THR', 'GLU', 'HIS', 'ARG', 'ASN'])


Inversion
===============================================================================

It is also possible to invert a selection:

>>> only_protein = prot.select('protein')
>>> print only_protein
Selection "protein" from 1p38
>>> only_non_protein = ~only_protein
>>> print only_non_protein
Selection "not (protein) " from 1p38

>>> water = prot.select('water')
>>> print water
Selection "water" from 1p38

The above shows that 1p38 does not contain any non-:term:`water` 
:term:`hetero` atoms.

Addition
===============================================================================

.. versionadded:: 0.2.1

Another operation defined on the :class:`Select` object is addition
(also on other :class:`~prody.atomic.AtomPointer` derived classes). 

This may be useful if you want to yield atoms in an :class:`AtomGroup` in a 
specific order.
Let's think of a simple case, where we want to output atoms in 1p38 in a 
specific order:

>>> protein = prot.select('protein')
>>> water = prot.select('water')
>>> water_protein = water + protein
>>> writePDB('1p38_water_protein.pdb', water_protein)
'1p38_water_protein.pdb'

In the resulting file, the :term:`water` atoms will precedes the 
:term:`protein` atoms.


Membership
===============================================================================

.. versionadded:: 0.5.3

Selections also allows membership test operations:

>>> backbone = prot.select('protein') 
>>> calpha = prot.select('calpha')

Is ``calpha`` a subset of ``backbone``?

>>> calpha in backbone
True

Or, is water in protein selection?

>>> water in protein
False

Other tests include:

>>> protein in prot
True
>>> backbone in prot
True
>>> prot in prot
True
>>> calpha in calpha
True

Equality
===============================================================================

.. versionadded:: 0.5.3

You can also check the equality of selections. Comparison will return
``True`` if both selections refer to the same atoms.

>>> calpha = prot.select('protein and name CA') 
>>> calpha2 = prot.select('calpha')
>>> calpha == calpha2
True

>>> all = prot.select('all')
>>> prot == all
True
