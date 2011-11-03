.. currentmodule:: prody.atomic

.. _hierview:

*******************************************************************************
Hierarchical Views
*******************************************************************************

Synopsis
===============================================================================

This example shows how to get a hierarchical view (:class:`HierView`) of atoms 
in a Protein Data Bank structure. A :class:`HierView` instance contains
chains which contains residues which contains atoms.

Input
-------------------------------------------------------------------------------

A protein structure file in PDB format or a PDB identifier.

Output
-------------------------------------------------------------------------------

A data structure with a hierarchical view of the contents of a protein 
structure. 

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse a structure
-------------------------------------------------------------------------------

Parsing a structure returns an :class:`AtomGroup` instance which has a plain
view of atoms. 

>>> structure = parsePDB('3mkb')
>>> structure
<AtomGroup: 3mkb (4776 atoms; 1 coordinate sets, active set index: 0)>

One can access individual atoms by indexing the :class:`AtomGroup` instance:

>>> structure[10]
<Atom: CG from 3mkb (index 10; 1 coordinate sets, active set index: 0)>

Or, one can iterate over all atoms one by one:

>>> for atom in structure:
...     atom # doctest: +ELLIPSIS
<Atom: N from 3mkb (index 0; 1 coordinate sets, active set index: 0)>
...
<Atom: O from 3mkb (index 4775; 1 coordinate sets, active set index: 0)>

Hierarchical view
-------------------------------------------------------------------------------

A hierarchical view of the structure can be simply get by calling the
:meth:`~AtomGroup.getHierView` method:

>>> hv = structure.getHierView()
>>> hv
<HierView: AtomGroup 3mkb>

*Indexing*

Indexing a :class:`HierView` instance returns :class:`Chain` instances:

>>> hv['A']
<Chain: A from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>
>>> hv['B']
<Chain: B from 3mkb (1193 atoms; 1 coordinate sets, active set index: 0)>
>>> hv['Z'] # This will return None, which means chain Z does not exist

The length of the *hv* variable gives the number of chains in the structure:

>>> len(hv)
4
>>> hv.numChains()
4

It is also possible to get a :class:`Residue` by directly indexing the
:class:`HierView` instance:

>>> hv['A', 100]
<Residue: MET 100 from Chain A from 3mkb (8 atoms; 1 coordinate sets, active set index: 0)>

Insertion codes can also be passed:

>>> hv['A', 100, 'B']

But this does not return anything, since residue 100B does not exist.

*Iterations*

One can iterate over :class:`HierView` instances to get chains:

>>> for chain in hv:
...     chain # doctest: +ELLIPSIS
<Chain: A from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>
...
<Chain: D from 3mkb (1196 atoms; 1 coordinate sets, active set index: 0)>

It is also possible to get a :func:`list` of chains simply as follows:

>>> chains = list( hv )
>>> chains # doctest: +SKIP
[<Chain: A from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>, 
 <Chain: B from 3mkb (1193 atoms; 1 coordinate sets, active set index: 0)>, 
 <Chain: C from 3mkb (1189 atoms; 1 coordinate sets, active set index: 0)>, 
 <Chain: D from 3mkb (1196 atoms; 1 coordinate sets, active set index: 0)>]

*Iterate residues*

In addition, one can also iterate over all residues:

>>> for residue in hv.iterResidues():
...     residue # doctest: +ELLIPSIS
<Residue: ALA 1 from Chain A from 3mkb (5 atoms; 1 coordinate sets, active set index: 0)>
<Residue: PHE 2 from Chain A from 3mkb (11 atoms; 1 coordinate sets, active set index: 0)>
...
<Residue: HOH 475 from Chain D from 3mkb (1 atoms; 1 coordinate sets, active set index: 0)>
<Residue: HOH 492 from Chain D from 3mkb (1 atoms; 1 coordinate sets, active set index: 0)>

Chains
-------------------------------------------------------------------------------

>>> chA = hv['A']
>>> chA
<Chain: A from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>

Length of the chain equals to the number of residues in it:

>>> len(chA)
254
>>> chA.numResidues()
254

*Indexing*

Indexing a :class:`Chain` instance returns a :class:`Residue` instance.

>>> chA[1]
<Residue: ALA 1 from Chain A from 3mkb (5 atoms; 1 coordinate sets, active set index: 0)>

If a residue does not exist, ``None`` is returned:

>>> chA[1000]
>>> chA[1, 'A'] # Residue 1 with insertion code A also does not exist  



If residue with given integer number does not exist, ``None`` is returned. 

*Iterations*

Iterating over a chain yields residues:

>>> for residue in chA:
...     residue # doctest: +ELLIPSIS
<Residue: ALA 1 from Chain A from 3mkb (5 atoms; 1 coordinate sets, active set index: 0)>
<Residue: PHE 2 from Chain A from 3mkb (11 atoms; 1 coordinate sets, active set index: 0)>
...
<Residue: HOH 490 from Chain A from 3mkb (1 atoms; 1 coordinate sets, active set index: 0)>
<Residue: HOH 493 from Chain A from 3mkb (1 atoms; 1 coordinate sets, active set index: 0)>

Note that water atoms, each constituting a residue, are also part of a chain
if they are labeled with that chain's identifier.

This enables getting a :func:`list` of residues simply as follows:

>>> chA_residues = list(chA)
>>> chA_residues # doctest: +SKIP
[<Residue: ALA 1 from Chain A from 3mkb (5 atoms; 1 coordinate sets, active set index: 0)>,
 ...,
 <Residue: HOH 493 from Chain A from 3mkb (1 atoms; 1 coordinate sets, active set index: 0)>]

*Get atomic data*

All methods defined for :class:`AtomGroup` class are also defined for 
:class:`Chain` (and also :class:`Residue`) class:

>>> print( chA.getCoords() ) # doctest: +ELLIPSIS
[[ -2.139  17.026 -13.287]
 [ -1.769  15.572 -13.111]
 [ -0.296  15.257 -13.467]
 ...
 [ -5.843  17.181 -16.86 ]
 [-13.199  -9.21  -49.692]
 [ -0.459   0.378 -46.156]]
>>> print( chA.getBetas() )
[ 59.35  59.14  58.5  ...,  57.79  47.77  40.77]

*Select atoms*

Finally, you can select atoms from a :class:`Chain` instance:

>>> chA_backbone = chA.select('backbone')
>>> chA_backbone
<Selection: "(backbone) and (chain A)" from 3mkb (560 atoms; 1 coordinate sets, active set index: 0)>
>>> chA_backbone.getSelstr()
'(backbone) and (chain A)'

As you see, the selection string passed by the user is augmented with 
"chain" keyword and identifier automatically to provide internal
consistency:

>>> structure.select( chA_backbone.getSelstr() )
<Selection: "(backbone) and (chain A)" from 3mkb (560 atoms; 1 coordinate sets, active set index: 0)>
 

Residues
-------------------------------------------------------------------------------

>>> chA_res1 = chA[1]
>>> chA_res1
<Residue: ALA 1 from Chain A from 3mkb (5 atoms; 1 coordinate sets, active set index: 0)>

*Indexing*

:class:`Residue` instances can be indexed to get individual atoms:

>>> chA_res1['CA']
<Atom: CA from 3mkb (index 1; 1 coordinate sets, active set index: 0)>
>>> chA_res1['CB']
<Atom: CB from 3mkb (index 4; 1 coordinate sets, active set index: 0)>
>>> chA_res1['X'] # if atom does not exist, None is returned

*Iterations*

Iterating over a residue instance yields :class:`Atom` instances:

>>> for atom in chA_res1:
...     atom # doctest: +ELLIPSIS
<Atom: N from 3mkb (index 0; 1 coordinate sets, active set index: 0)>
...
<Atom: CB from 3mkb (index 4; 1 coordinate sets, active set index: 0)>

This makes it easy to get a :func:`list` of atoms:

>>> list( chA_res1 ) # doctest: +SKIP
[<Atom: N from 3mkb (index 0; 1 coordinate sets, active set index: 0)>,
 <Atom: CA from 3mkb (index 1; 1 coordinate sets, active set index: 0)>,
 <Atom: C from 3mkb (index 2; 1 coordinate sets, active set index: 0)>,
 <Atom: O from 3mkb (index 3; 1 coordinate sets, active set index: 0)>,
 <Atom: CB from 3mkb (index 4; 1 coordinate sets, active set index: 0)>]

*Get atomic data*

All methods defined for :class:`AtomGroup` class are also defined for 
:class:`Residue` class:

>>> print( chA_res1.getCoords() )
[[ -2.139  17.026 -13.287]
 [ -1.769  15.572 -13.111]
 [ -0.296  15.257 -13.467]
 [  0.199  14.155 -13.155]
 [ -2.752  14.639 -13.898]]
>>> print( chA_res1.getBetas() )
[ 59.35  59.14  58.5   59.13  59.02]

*Select atoms*

Finally, you can select atoms from a :class:`Residue` instance:

>>> chA_res1_bb = chA_res1.select('backbone')
>>> chA_res1_bb
<Selection: "(backbone) and ...A and resnum 1)" from 3mkb (4 atoms; 1 coordinate sets, active set index: 0)>
>>> chA_res1_bb.getSelstr()
'(backbone) and (chain A and resnum 1)'

Again, the selection string is augmented with the chain identifier and 
residue number ("resnum").

Atoms
-------------------------------------------------------------------------------

The lowest level of the hierarchical view contains :class:`Atom` instances.

>>> chA_res1_CA = chA_res1['CA']
>>> chA_res1_CA
<Atom: CA from 3mkb (index 1; 1 coordinate sets, active set index: 0)>


*Get atomic data*

All methods defined for :class:`AtomGroup` class are also defined for 
:class:`Atom` class with the difference that method names are singular 
(except for coordinates):

>>> print( chA_res1_CA.getCoords() )
[ -1.769  15.572 -13.111]
>>> print( chA_res1_CA.getBeta() )
59.14



Changes in state
-------------------------------------------------------------------------------

A :class:`HierView` instance represents the state of an :class:`AtomGroup` 
instance at the time it is built. When chain identifiers or residue 
numbers change, the state that hierarchical view represents may not
match the current state of the atom group:

>>> chA.setIdentifier('X')
>>> chA
<Chain: X from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>
>>> hv['X'] # returns None, since hierarchical view is not updated
>>> hv.update() # this updates hierarchical view
>>> hv['X']
<Chain: X from 3mkb (1198 atoms; 1 coordinate sets, active set index: 0)>


|questions|

|suggestions|


