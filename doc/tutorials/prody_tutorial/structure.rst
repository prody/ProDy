Structure Analysis
===============================================================================

Measure geometric properties
-------------------------------------------------------------------------------

Functions for analyzing molecular structure can be found in :mod:`.measure`
module. For example, you can calculate phi (φ) and psi (ψ) for the 10th
residue, or the radius of gyration of the protein as follows:

.. ipython:: python
   :suppress:

   from prody import *; from pylab import *; ion(); p38 = parsePDB('1p38')

.. ipython:: python

   print(calcPhi(p38[10,]).round(2)) # rounded numbers for printing purposes
   print(calcPsi(p38[10,]).round(2))
   print(calcGyradius(p38).round(2))


Compare and align structures
-------------------------------------------------------------------------------

You can also compare different structures using some of the methods in
:mod:`.proteins` module.  Let's parse another p38 MAP kinase structure

.. ipython:: python

   bound = parsePDB('1zz2')

You can find similar chains in structure 1p38 and 1zz2 using
:func:`.matchChains` function:


.. ipython:: python

   apo_chA, bnd_chA, seqid, overlap = matchChains(p38, bound)[0]
   repr(apo_chA)
   repr(bnd_chA)
   print(int(seqid)) # precent sequence identity between two chains
   print(int(overlap)) # percent overlap between two chains

Matching Cα atoms are selected and returned as :class:`.AtomMap` instances.
We can use them to calculate RMSD and superpose structures.

.. ipython:: python

   print(calcRMSD(bnd_chA, apo_chA))
   bnd_chA, transformation = superpose(bnd_chA, apo_chA)
   print(calcRMSD(bnd_chA, apo_chA))


.. ipython:: python

   showProtein(p38);
   @savefig prody_tutorial_structure_compare.png width=4in
   showProtein(bound);


Writing PDB files
-------------------------------------------------------------------------------

PDB files can be written using the :func:`.writePDB` function.
The function accepts objects containing or referring to atomic data.

Writing selected atoms:

.. ipython:: python

   writePDB('1p38_calphas.pdb', p38.select('calpha'))


Writing a chain:

.. ipython:: python

   chain_A = p38['A']
   writePDB('1p38_chain_A.pdb', chain_A)


As you may have noticed, this function returns the file name after it is
successfully written.  This is a general behavior for ProDy output functions.
For more PDB writing examples see :ref:`writepdb`.