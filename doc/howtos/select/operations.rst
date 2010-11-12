.. currentmodule:: prody.select

.. _selops:

*******************************************************************************
Selection operations
*******************************************************************************

Before reading this section, familiarity with :ref:`selections` might be 
helpful.

Let's import all classes and functions from ProDy and parse coordinates from
PDB structure 1p38:

>>> from prody import *
>>> prot = parsePDB('1p38')


Union
===============================================================================

Let's select beta-carbon atoms of for non-GLY amino acid residues, and 
alpha-carbons for GLYs in two steps:

>>> betas = prot.select('name CB and protein')
"name CB and protein" from 1p38 (336 atoms; 1 coordinate sets, active set index: 0)
>>> somealphas = prot.select('name CA and resname GLY')
"name CA and resname GLY" from 1p38 (15 atoms; 1 coordinate sets, active set index: 0)

This shows that p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

>>> betas_somealphas = betas | somealphas
>>> print betas_somealphas
"(name CB and prote...resname GLY)" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)

The selection string for the union of selections become:

>>> print betas_somealphas.getSelectionString()
(name CB and protein) or (name CA and resname GLY)

Note that it is also possibe to reach the same selection using selection string
``(name CB and protein) or (name CA and resname GLY)``.


Intersection
===============================================================================

It is as easy to get the intersection of two selections. Let's find charged
and medium size residues in a protein:

>>> charged = prot.select('charged')
>>> print charged
"charged" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)
>>> medium = prot.select('medium')
>>> print medium
"medium" from 1p38 (751 atoms; 1 coordinate sets, active set index: 0)

>>> medium_charged = medium & charged
>>> print medium_charged
"(medium) and (charged)" from 1p38 (216 atoms; 1 coordinate sets, active set index: 0)
>>> print medium_charged.getSelectionString()
(medium) and (charged)

Let's see which amino acids are considered charged and medium:

>>> print set(medium_charged.getResidueNames())
set(['ASP'])

What about amino acids that are medium or charged:

>>> print set((medium | charged).getResidueNames())
set(['CYS', 'ASP', 'VAL', 'LYS', 'PRO', 'THR', 'GLU', 'HIS', 'ARG', 'ASN'])


Inversion
===============================================================================

It is also possible to invert a selection:

>>> only_protein = prot.select('protein')
>>> print only_protein
"protein" from 1p38 (2833 atoms; 1 coordinate sets, active set index: 0)
>>> only_non_protein = ~only_protein
>>> print only_non_protein
"not (protein) " from 1p38 (129 atoms; 1 coordinate sets, active set index: 0)

>>> water = prot.select('water')
>>> print water
"water" from 1p38 (129 atoms; 1 coordinate sets, active set index: 0)

This shows that, structure 1p38 does not contain any non-water hetero atoms.
