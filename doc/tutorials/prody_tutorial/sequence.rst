Sequence Analysis
===============================================================================

Evol component of ProDy package brought new fast, flexible, and efficient
features for handling multiple sequence alignments and analyzing sequence
evolution.  Here, we just give a brief introduction to these features.  For
more detailed examples, see :ref:`evol-tutorial`.


.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Access Pfam
-------------------------------------------------------------------------------

First, let's fetch an MSA file from `Pfam`_ database:

.. ipython:: python

   filename = fetchPfamMSA('pkinase', alignment='seed')
   filename


We downloaded the seed alignment for Protein Kinase (:pfam:`Pkinase`) family.

.. _Pfam: http://pfam.sanger.ac.uk/


Parse MSA
-------------------------------------------------------------------------------

As you might guess, we will parse this file using :func:`.parseMSA` function:

.. ipython:: python

   msa = parseMSA(filename)
   msa


Sequences
-------------------------------------------------------------------------------

You can access :class:`.Sequence` objects by indexing :class:`.MSA`:

.. ipython:: python

   seq = msa[0]
   seq
   print(seq)


You can also slice :class:`.MSA` objects and iterate over sequences:

.. ipython:: python

   for seq in msa[:4]:
       repr(seq)



Analysis
-------------------------------------------------------------------------------

Evol component includes several functions for calculating conservation and
coevolution properties of amino acids, which are shown in :ref:`evol-tutorial`.
Here, let's take a look at :func:`.calcMSAOccupancy` and
:func:`.showMSAOccupancy` functions:

.. ipython:: python

   occ = calcMSAOccupancy(msa, count=True)
   occ.min()

This shows that an amino acid is present only in one of the sequences in the
MSA.


.. ipython:: python

   @savefig prody_tutorial_sequence_occ.png width=6in
   showMSAOccupancy(msa, count=True);


You see that many residues are not present in all sequences. You will see
how to refine such MSA instances in :ref:`evol-tutorial`.