.. currentmodule:: prody.select

.. _selops:

*******************************************************************************
Operations on Selections 
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
>>> print len(betas)
336
>>> somealphas = prot.select('name CA and resname GLY')
>>> print len(somealphas)
15

This shows that p38 structure contains 15 GLY residues.

These two selections can be combined as follows:

>>> betas_somealphas = betas | somealphas
>>> print betas_somealphas
Selection "(name CB and pr...nd resname GLY)" from 1p38
>>> print len(betas_somealphas)
351

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
Selection "charged" from 1p38
>>> medium = prot.select('medium')
>>> print medium
Selection "medium" from 1p38

>>> medium_charged = medium & charged
>>> print medium_charged
Selection "(medium) and (charged)" from 1p38
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
Selection "protein" from 1p38
>>> only_non_protein = ~only_protein
>>> print only_non_protein
Selection "not (protein) " from 1p38

>>> water = prot.select('water')
>>> print water
Selection "water" from 1p38

This shows that, structure 1p38 does not contain any non-water hetero atoms.

Addition
===============================================================================

.. versionadded:: 0.2.1

Final operation defined on :class:`Selection` object (also on other 
:class:`AtomPointer` derived classes) is addition. 

This may be useful if you want to get atoms in an :class:`AtomGroup` in a 
specific order.
Let's think of a simple case, where we want to write atoms in 1p38 in a 
specific order.

>>> protein = prot.select('protein')
>>> water = prot.select('water')
>>> water_protein = water + protein
>>> writePDB('1p38_water_protein.pdb', water_protein)
'1p38_water_protein.pdb'

In the resulting file, water atoms will preceed protein atoms.
 
  

