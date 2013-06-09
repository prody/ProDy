.. _normalmode-operations:

Normal Mode Algebra
===============================================================================

This part shows how to use some handy features of :class:`.Mode` objects.

ANM Calculations
-------------------------------------------------------------------------------

We will compare modes from two ANMs for the same protein, but everything
applies to comparison of ANMs and PCAs (as long as they contain same number
of atoms).

Let's get started by getting ANM models for two related protein structures:

.. ipython:: python

   from prody import *
   str1 = parsePDB('1p38')
   str2 = parsePDB('1r39')

**Find and align matching chains**

.. ipython:: python

   matches = matchChains(str1, str2)
   match = matches[0]
   ch1 = match[0]
   ch2 = match[1]

Minimize RMSD by superposing ``ch2`` onto ``ch1``:

.. ipython:: python

   ch2, t = superpose(ch2, ch1)  # t is transformation, already applied to ch2
   calcRMSD(ch1, ch2)

**Get ANM models for each chain**

.. ipython:: python

   anm1, ch1 = calcANM(ch1)
   anm2, ch2 = calcANM(ch2)
   anm1[0]

Let's rename these :class:`.ANM` instances, so that they print short:

.. ipython:: python

   anm1.setTitle('1p38_anm')
   anm2.setTitle('1r39_anm')

This is how they print now:

.. ipython:: python

   anm1[0]
   anm2[0]


Calculate overlap
-------------------------------------------------------------------------------

We need Numpy in this part:

.. ipython:: python

   from numpy import *

Multiplication of two :class:`.Mode` instances returns dot product
of their eigenvectors. This dot product is the overlap or cosine correlation
between modes.

Let's calculate overlap for slowest modes:

.. ipython:: python

   overlap = anm1[0] * anm2[0]
   overlap

This show that the overlap between these two modes is 0.98, which is not
surprising since ANM modes come from structures of the *same* protein.

To compare multiple modes, convert a list of modes to a :func:`numpy.array`:

.. ipython:: python

   array(list(anm1[:3])) * array(list(anm2[:3]))

This shows that slowest three modes are almost identical.

We could also generate a matrix of overlaps using :func:`numpy.outer`:

.. ipython:: python

   outer_product = outer(array(list(anm1[:3])), array(list(anm2[:3])))
   outer_product

This could also be printed in a pretty table format using
:func:`.printOverlapTable`:

.. ipython:: python

   printOverlapTable(anm1[:3], anm2[:3])

**Scaling**

:class:`.Mode` instances can be scaled, but after this operation they will
become :class:`.Vector` instances:

.. ipython:: python

   anm1[0] * 10

Linear combination
-------------------------------------------------------------------------------

It is also possible to linearly combine normal modes:

.. ipython:: python

   anm1[0] * 3 + anm1[1] + anm1[2] * 2


Or, we could use eigenvalues for linear combination:

.. ipython:: python

   lincomb = anm1[0] * anm1[0].getEigval() + anm1[1] * anm1[1].getEigval()

It is the name of the :class:`.Vector` instance that keeps track of operations.

.. ipython:: python

   lincomb.getTitle()

Approximate a deformation vector
-------------------------------------------------------------------------------

Let's get the deformation vector between *ch1* and *ch2*:

.. ipython:: python

   defvec = calcDeformVector(ch1, ch2)
   abs(defvec)


Let's see how deformation projects onto ANM modes:

.. ipython:: python

   array(list(anm1[:3])) * defvec


We can use these numbers to combine ANM modes:

.. ipython:: python

   approximate_defvec = sum((array(list(anm1[:3])) * defvec) *
                            array(list(anm1[:3])))
   approximate_defvec

Let's deform 1r39 chain along this approximate deformation vector and see
how RMSD changes:

.. ipython:: python

   ch2.setCoords(ch2.getCoords() - approximate_defvec.getArrayNx3())
   calcRMSD(ch1, ch2)

RMSD decreases from 0.89 A to 0.82 A.