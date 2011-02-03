.. currentmodule:: prody.compare

.. _deformation:

*******************************************************************************
Deformation and ANM analysis
*******************************************************************************

Synopsis
===============================================================================

This example shows how to calculate the deformation vector describing the 
change between two structures of a protein. 

User Input
-------------------------------------------------------------------------------

Two structures for the same protein.


ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *


Parse structures
-------------------------------------------------------------------------------


Let's parse two p38 MAP Kinase structures: 1p38 and 1zz2

>>> reference = parsePDB('1p38')
>>> mobile = parsePDB('1zz2')   # this is the one we want to superimpose

Match chains
-------------------------------------------------------------------------------

ProDy offers the function :func:`matchChains` to find matching chains
in two structures easily. We use it to find the chains for which we will 
calculate the deformation vector:

>>> matches = matchChains(reference, mobile)

:func:`matchChains` function returns a list. If there are no matching chains, 
list is empty, else the list contains a tuple for each pair of matching chains.

>>> print len(matches) 
1
>>> match = matches[0]

There is only one match in this case. First item is a subset of atoms from the 
first structure (*reference*). Second item is a subset of atoms from the 
second structure (*mobile*).

>>> ref_chain = match[0]
>>> mob_chain = match[1]

Matched atoms are returned in :class:`~prody.atomic.AtomMap` instances.
We can get information on matched subset of atoms by entering the variable 
name:

>>> ref_chain
<AtomMap: AtomGroup 1p38 -> AtomGroup 1zz2 (from 1p38; 337 atoms; 337 mapped; 0 unmapped; 1 coordinate sets, active set index: 0)>
>>> mob_chain
<AtomMap: AtomGroup 1zz2 -> AtomGroup 1p38 (from 1zz2; 337 atoms; 337 mapped; 0 unmapped; 1 coordinate sets, active set index: 0)>

Both :class:`~prody.atomic.AtomMap` instances refer to same number of atoms, 
and their name suggests how they were retrieved.

In addition, we can find out the sequence identity that the matched atoms 
(residues) share (third item in the tuple):

>>> match[2]
99.406528189910986

The fourth item in the tuple shows the coverage of the matching:

>>> match[3]
96

This is the percentage of matched residues with respect to the longer chain.
1p38 chain A contains 351 resiudes, 96% of it is 337 residues, which
is the number of atoms in the returned atom maps. 


Calculate RMSD and superpose
-------------------------------------------------------------------------------

We calculate the RMSD using :func:`~prody.measure.calcRMSD` function: 

>>> print calcRMSD(ref_chain, mob_chain) # doctest: +SKIP
72.93

Let's find the transformation that minimizes RMSD between these chains
using :func:`~prody.measure.calcTransformation` function:

>>> t = calcTransformation(mob_chain, ref_chain)

We apply this transformation to *mobile* structure (not to *mob_chain*, 
to preserve structures integrity).

>>> t.apply(mobile)
<AtomGroup: 1zz2 (2872 atoms; 1 coordinate sets, active set index: 0)>
>>> print calcRMSD(ref_chain, mob_chain) # doctest: +SKIP
1.86

Deformation vector
-------------------------------------------------------------------------------

Once matching chains are identified it is straightforward to calculate the
deformation vector using :func:`~prody.measure.calcDeformVector`

>>> defvec = calcDeformVector(ref_chain, mob_chain)
>>> print abs(defvec) # doctest: +SKIP
34.20


To show how RMSD and deformation vector are related, we can be calculate 
RMSD from the magnitude of the deformation vector:

>>> print (abs(defvec)**2 / len(ref_chain)) ** 0.5 # doctest: +SKIP
1.86

Array of numbers for this deformation can be obtained as follows

>>> arr = defvec.getArray() # arr is a NumPy array
>>> print arr # doctest: +SKIP

Following yields the normalized deformation vector

>>> defvecnormed = defvec.getNormed()
>>> print abs(defvecnormed) # doctest: +SKIP
1.0

Compare with ANM modes
-------------------------------------------------------------------------------

Let's get ANM model for the reference chain using 
:func:`~prody.dynamics.calcANM` (a shorthand function for ANM calculations):

>>> anm = calcANM(ref_chain)

Calculate overlap between slowest ANM mode and the deformation vector

>>> print anm[0] * defvecnormed # note that we used normalized deformation vector # doctest: +SKIP
-0.422463159673

We can do this for a set of ANM modes (slowest 6) as follows

>>> import numpy as np
>>> print np.array( anm[:6].getModes() ) * defvecnormed # doctest: +SKIP 
[-0.422463159673 -0.135595183602 0.489432673122 0.0312043675943
 -0.17064639491 -0.095208459294]

|questions|

|suggestions|
