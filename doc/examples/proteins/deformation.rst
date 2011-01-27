.. currentmodule:: prody.compare

.. _deformation:

*******************************************************************************
Deformation and ANM analysis
*******************************************************************************

Getting the deformation vector that describes the change in atomic positions
of a protein from one of it's structures to another is most of the time
tricky. ProDy implements :func:`matchChains` to make this an easy task.
This function finds matching chains and returns corresponding atoms form these
chains.


Let's parse two p38 MAP Kinase structures: 1p38 and 1zz2

>>> from prody import *
>>> reference = parsePDB('1p38')
>>> mobile = parsePDB('1zz2')   # this is the one we want to superimpose

Let's find matching chains in these structures

>>> matches = matchChains(reference, mobile)

:func:`matchChains` function returns a list. If there are no matching chains, 
list is empty else, list contains a tuple for each pair of matching chains.

>>> print len(matches) 
1

Let's get the first match (there is only one match in this case)

>>> match = matches[0]

First item is a selection of atoms from the first structure (reference)

>>> ref_chain = match[0]

Second item is a selection of atoms from the second structure (mobile)

>>> mob_chain = match[1]


RMSD and superimposition
===============================================================================

>>> print calcRMSD(ref_chain, mob_chain) # doctest: +SKIP
72.93

Let's find the transformation that minimizes RMSD between these chains

>>> t = calcTransformation(mob_chain, ref_chain)

We apply this transformation to mobile structure (not to mob_chain, 
to preserve structures integrity)

>>> t.apply(mobile)
<AtomGroup: 1zz2 (2872 atoms; 1 coordinate sets, active set index: 0)>
>>> print calcRMSD(ref_chain, mob_chain) # doctest: +SKIP
1.86
>>> print len(ref_chain)
337

Deformation vector
===============================================================================

Let's get the deformation vector

>>> defvec = calcDeformVector(ref_chain, mob_chain)
>>> print abs(defvec) # doctest: +SKIP
34.20

RMSD can be calculated from the magnitude of the deformation vector

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
===============================================================================

Let's get ANM model for the reference chain

>>> anm = calcANM(ref_chain)

Calculate overlap between slowest ANM mode and the deformation vector

>>> print anm[0] * defvecnormed # note that we used normalized deformation vector # doctest: +SKIP
-0.422463159673

We can do this for a set of ANM modes (slowest 6) as follows

>>> import numpy as np
>>> print np.array( anm[:6].getModes() ) * defvecnormed # doctest: +SKIP 
[-0.422463159673 -0.135595183602 0.489432673122 0.0312043675943
 -0.17064639491 -0.095208459294]

