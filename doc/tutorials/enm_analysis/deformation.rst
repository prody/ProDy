.. _deformation:


Deformation Analysis
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to calculate the deformation vector describing the
change between two structures of a protein.  Two structures of the same
protein in PDB format will be used.  A :class:`.Vector` instance that
contains the deformation vector describing the change in protein structure
will be calculated. This object will be compared to :class:`.ANM` modes.

Parse structures
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:


.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()


Let's parse two p38 MAP Kinase structures: 1p38 and 1zz2

.. ipython:: python

   reference = parsePDB('1p38')
   mobile = parsePDB('1zz2')   # this is the one we want to superimpose

Match chains
-------------------------------------------------------------------------------

ProDy offers the function :func:`.matchChains` to find matching chains
in two structures easily. We use it to find the chains for which we will
calculate the deformation vector:

.. ipython:: python

   matches = matchChains(reference, mobile)

:func:`.matchChains` function returns a list. If there are no matching chains,
list is empty, else the list contains a tuple for each pair of matching chains.

.. ipython:: python

   len(matches)
   match = matches[0]

There is only one match in this case. First item is a subset of atoms from the
first structure (*reference*). Second item is a subset of atoms from the
second structure (*mobile*).

.. ipython:: python

   ref_chain = match[0]
   mob_chain = match[1]

Matched atoms are returned in :class:`.AtomMap` instances.
We can get information on matched subset of atoms by entering the variable
name:

.. ipython:: python

   ref_chain
   mob_chain


Both :class:`.AtomMap` instances refer to same number of atoms,
and their name suggests how they were retrieved.

In addition, we can find out the sequence identity that the matched atoms
(residues) share (third item in the tuple):

.. ipython:: python

   match[2]

The fourth item in the tuple shows the coverage of the matching:

.. ipython:: python

   match[3]

This is the percentage of matched residues with respect to the longer chain.
1p38 chain A contains 351 resiudes, 96% of it is 337 residues, which
is the number of atoms in the returned atom maps.


RMSD and superpose
-------------------------------------------------------------------------------

We calculate the RMSD using :func:`.calcRMSD` function:

.. ipython:: python

   calcRMSD(ref_chain, mob_chain).round(2)


Let's find the transformation that minimizes RMSD between these chains
using :func:`.calcTransformation` function:

.. ipython:: python

   t = calcTransformation(mob_chain, ref_chain)

We apply this transformation to *mobile* structure (not to *mob_chain*,
to preserve structures integrity).

.. ipython:: python

   t.apply(mobile)
   calcRMSD(ref_chain, mob_chain).round(2)


Deformation vector
-------------------------------------------------------------------------------

Once matching chains are identified it is straightforward to calculate the
deformation vector using :func:`.calcDeformVector`

.. ipython:: python

   defvec = calcDeformVector(ref_chain, mob_chain)
   abs(defvec).round(3)


To show how RMSD and deformation vector are related, we can be calculate
RMSD from the magnitude of the deformation vector:

.. ipython:: python

   (abs(defvec)**2 / len(ref_chain)) ** 0.5


Array of numbers for this deformation can be obtained as follows

.. ipython:: python

   arr = defvec.getArray() # arr is a NumPy array
   arr.round(2)

Following yields the normalized deformation vector

.. ipython:: python

   defvecnormed = defvec.getNormed()
   abs(defvecnormed)

Compare with ANM modes
-------------------------------------------------------------------------------

Let's get ANM model for the reference chain using
:func:`.calcANM` (a shorthand function for ANM calculations):

.. ipython:: python

   anm = calcANM(ref_chain)[0]

Calculate overlap between slowest ANM mode and the deformation vector

.. ipython:: python

   (anm[0] * defvecnormed).round(2) # used normalized deformation vector

We can do this for a set of ANM modes (slowest 6) as follows

.. ipython:: python

   (array(list(anm[:6])) * defvecnormed).astype(float64).round(2)