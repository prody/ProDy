.. module:: prody.dynamics

.. _nmops:

*******************************************************************************
Normal Mode Operations
*******************************************************************************

In this example we will compare modes from two ANMs for the same protein, 
but everything applies to comparison of ANMs and PCAs (as long as they contain 
same number of atoms).

Let's get started by getting ANM models for two related protein structures:

>>> from prody import *
>>> str_one = parsePDB('1p38')
@> 1p38 downloaded (./1p38.pdb.gz)
@> 2962 atoms and 1 coordinate sets were parsed in 0.09s.
>>> str_two = parsePDB('1r39')
@> 1r39 downloaded (./1r39.pdb.gz)
@> 2912 atoms and 1 coordinate sets were parsed in 0.08s.

**Find and align matching chains**

>>> matches = findMatchingChains(str_one, str_two)
@> Checking 1p38: 1 chains are identified
@> Checking 1r39: 1 chains are identified
@> Comparing Chain A from 1p38 (len=351) and Chain A from 1r39 (len=345):
@> 	345 residues match with 99% sequence identity and 98% coverage based on residue numbers.
>>> match = matches[0]
>>> ch_one = match[0]
>>> ch_two = match[1]

>>> ch_two, t = superimpose(ch_two, ch_one) # minimizes RMSD by moving ch_two onto ch_one
>>> # t is transformation, which is already applied to ch_two, so we don't use it
>>> getRMSD(ch_one, ch_two)
0.89840163398680861


**Get ANM models for each chain**

>>> anm_one = getANM(ch_one)
@> Hessian was built in 1.59s.
@> 20 modes were calculated in 1.50s.
>>> anm_two = getANM(ch_two)
@> Hessian was built in 1.58s.
@> 20 modes were calculated in 1.50s.

>>> print anm_one[0]
NormalMode 0 from "(calpha and all) a...o 2797 2804)" from 1p38 (345 atoms; 1 coordinate sets, active set index: 0)

Let's rename these :class:`ANM` instances, so that they print short: 

>>> anm_one.setName('1p38_anm')
>>> anm_two.setName('1r39_anm')

This is how they print now:

>>> print anm_one[0]
NormalMode 0 from 1p38_anm
>>> print anm_two[0]
NormalMode 0 from 1r39_anm

Calculate overlap
===============================================================================

Multiplication of two :class:`NormalMode` instances returns dot product
of their eigenvectors. This dot product is the overlap or cosine correlation
between modes.

Let's calculate overlap for slowest modes:

>>> print anm_one[0] * anm_two[0]
-0.98402119545

This show that the overlap between these two modes is 0.98, which is not 
surprising since ANM modes come from structures of the *same* protein.

To compare multiple modes, convert a list of modes to a :func:`numpy.array`:

>>> print np.array(anm_one[:6]) * np.array(anm_two[:3])
[-0.98402119545 -0.98158348545 -0.991357811832]

This shows that slowest three modes are almost identical.

We could also generate a matrix of overlaps using :func:`numpy.outer`:

>>> print np.outer( np.array(anm_one[:3]),  np.array(anm_two[:3]) )
[[-0.98402119545 -0.144944616676 -0.00217115583244]
 [0.148366788279 -0.98158348545 0.0807736109529]
 [0.0104328721628 -0.084078114473 -0.991357811832]]

This could also be printed in a pretty table format using :func:`printOverlapTable`:

>>> printOverlapTable(anm_one[:3], anm_two[:3])
Overlap/Correlation Table
                    1r39_anm      
                #1     #2     #3  
1p38_anm #1   -0.98  -0.14  -0.00 
1p38_anm #2   +0.15  -0.98  +0.08 
1p38_anm #3   +0.01  -0.08  -0.99 

**Scaling**

:class:`NormalMode` instances can be scaled, but after this operation they will
become :class:`SimpleMode` instances:

>>> anm_one[0] * 10
<SimpleMode: 10*(NormalMode 0 from 1p38_anm)>

Linear combination
===============================================================================

It is also possible to linearly combine normal modes:

>>> anm_one[0] * 3 + anm_one[1] + anm_one[2] * 2
<SimpleMode: 3*(NormalMode 0 from 1p38_anm) + NormalMode 1 from 1p38_anm + 2*(NormalMode 2 from 1p38_anm)>

Or, we could use eigenvalues for linear combination:

>>> lincomb = anm_one[0] * anm_one[0].getEigenvalue() + anm_one[1] * anm_one[1].getEigenvalue()

It is the name of the :class:`SimpleMode` instance that keeps track of operations.

>>> print lincomb.getName()  
0.148971269751*(NormalMode 0 from 1p38_anm) + 0.24904210757*(NormalMode 1 from 1p38_anm)

Approximate a deformation vector
===============================================================================

Let's get the deformation vector between *ch_one* and *ch_two*:

>>> defvec = getDeformation(ch_one, ch_two)
>>> abs(defvec)
16.687069727870373

Let's see how deformation projects onto ANM modes:

>>> print np.array(anm_one[:3]) * defvec
[-5.60860594784 2.15393365959 -3.13701609199]

We can use these numbers to combine ANM modes:

>>> approximate_defvec = sum( (np.array(anm_one[:3]) * defvec) * np.array(anm_one[:3]) ) 
>>> print approximate_defvec
-5.60860594784*(NormalMode 0 from 1p38_anm) + 2.15393365959*(NormalMode 1 from 1p38_anm) + -3.13701609199*(NormalMode 2 from 1p38_anm)

Let's deform 1r39 chain along this approximate deformation vector and see
how RMSD changes:

>>> ch_two.setCoordinates(ch_two.getCoordinates() - approximate_defvec.get3dArray())
>>> getRMSD(ch_one, ch_two)
0.82096008703377343

RMSD decreases from 0.89 A to 0.82 A.
