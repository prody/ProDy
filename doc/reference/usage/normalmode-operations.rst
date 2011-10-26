.. currentmodule:: prody.dynamics

.. _normalmode-operations:

*******************************************************************************
Normal Mode Algebra 
*******************************************************************************

In this example we will compare modes from two ANMs for the same protein, 
but everything applies to comparison of ANMs and PCAs (as long as they contain 
same number of atoms).

Let's get started by getting ANM models for two related protein structures:

>>> from prody import *
>>> str_one = parsePDB('1p38')
>>> str_two = parsePDB('1r39')

**Find and align matching chains**

>>> matches = matchChains(str_one, str_two)
>>> match = matches[0]
>>> ch_one = match[0]
>>> ch_two = match[1]

>>> ch_two, t = superpose(ch_two, ch_one) # minimizes RMSD by moving ch_two onto ch_one
>>> # t is transformation, which is already applied to ch_two, so we don't use it
>>> rmsd = calcRMSD(ch_one, ch_two)
>>> print( '{0:.2f}'.format(rmsd) ) # We are just printing rmsd with some formatting
0.90


**Get ANM models for each chain**

>>> anm_one, ch_one = calcANM(ch_one)
>>> anm_two, ch_two = calcANM(ch_two)

>>> print( anm_one[0] )
Mode 1 from ANM 1p38

Let's rename these :class:`ANM` instances, so that they print short: 

>>> anm_one.setTitle('1p38_anm')
>>> anm_two.setTitle('1r39_anm')

This is how they print now:

>>> print( anm_one[0] )
Mode 1 from ANM 1p38_anm
>>> print( anm_two[0] )
Mode 1 from ANM 1r39_anm

Calculate overlap
===============================================================================

We need Numpy in this part:

>>> import numpy as np
>>> np.set_printoptions(precision=3) # This is to make sure same numbers are printed in your session

Multiplication of two :class:`Mode` instances returns dot product
of their eigenvectors. This dot product is the overlap or cosine correlation
between modes.

Let's calculate overlap for slowest modes:

>>> overlap = anm_one[0] * anm_two[0]
>>> print( '{0:.3f}'.format(overlap) ) 
-0.984

This show that the overlap between these two modes is 0.98, which is not 
surprising since ANM modes come from structures of the *same* protein.

To compare multiple modes, convert a list of modes to a :func:`numpy.array`:

>>> print( np.array(list(anm_one[:3])) * np.array(list(anm_two[:3])) )
[-0.98402119545 -0.98158348545 -0.991357811832]

This shows that slowest three modes are almost identical.

We could also generate a matrix of overlaps using :func:`numpy.outer`:

>>> outer = np.outer( np.array(list(anm_one[:3])),  np.array(list(anm_two[:3])) )
>>> print( outer.astype(np.float64).round(2) )
[[-0.98 -0.14 -0.  ]
 [ 0.15 -0.98  0.08]
 [ 0.01 -0.08 -0.99]]

This could also be printed in a pretty table format using :func:`printOverlapTable`:

>>> printOverlapTable(anm_one[:3], anm_two[:3])
Overlap Table
                      ANM 1r39_anm
                    #1     #2     #3
ANM 1p38_anm #1   -0.98  -0.14   0.00
ANM 1p38_anm #2   +0.15  -0.98  +0.08
ANM 1p38_anm #3   +0.01  -0.08  -0.99
<BLANKLINE>


**Scaling**

:class:`Mode` instances can be scaled, but after this operation they will
become :class:`Vector` instances:

>>> anm_one[0] * 10
<Vector: 10*(Mode 1 from ANM 1p38_anm)>

Linear combination
===============================================================================

It is also possible to linearly combine normal modes:

>>> anm_one[0] * 3 + anm_one[1] + anm_one[2] * 2
<Vector: 3*(Mode 1 from ANM 1p38_anm) + Mode 2 from ANM 1p38_anm + 2*(Mode 3 from ANM 1p38_anm)>

Or, we could use eigenvalues for linear combination:

>>> lincomb = anm_one[0] * anm_one[0].getEigenvalue() + anm_one[1] * anm_one[1].getEigenvalue()

It is the name of the :class:`Vector` instance that keeps track of operations.

>>> print( lincomb.getTitle() )  
0.148971269751*(Mode 1 from ANM 1p38_anm) + 0.24904210757*(Mode 2 from ANM 1p38_anm)

Approximate a deformation vector
===============================================================================

Let's get the deformation vector between *ch_one* and *ch_two*:

>>> defvec = calcDeformVector(ch_one, ch_two)
>>> defvec_magnitude = abs(defvec)
>>> print( '{0:.2f}'.format(defvec_magnitude) )
16.69

Let's see how deformation projects onto ANM modes:

>>> print( np.array(list(anm_one[:3])) * defvec )
[-5.60860594784 2.15393365959 -3.13701609199]

We can use these numbers to combine ANM modes:

>>> approximate_defvec = np.sum( (np.array(list(anm_one[:3])) * defvec) * np.array(list(anm_one[:3])) ) 
>>> print( approximate_defvec )
-5.60860594784*(Mode 1 from ANM 1p38_anm) + 2.15393365959*(Mode 2 from ANM 1p38_anm) + -3.13701609199*(Mode 3 from ANM 1p38_anm)

Let's deform 1r39 chain along this approximate deformation vector and see
how RMSD changes:

>>> ch_two.setCoordinates(ch_two.getCoordinates() - approximate_defvec.getArrayNx3())
>>> rmsd = calcRMSD(ch_one, ch_two)
>>> print( '{0:.2f}'.format(rmsd) )
0.82

RMSD decreases from 0.89 A to 0.82 A.
