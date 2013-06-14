.. _biomolt:

Building Biomolecules
===============================================================================

Some PDB files contain coordinates for a monomer of a functional/biological
multimer (biomolecule).  ProDy offers functions to build structures of
biomolecules using the header data from the PDB file.  We will use PDB file
that contains the coordinates for a monomer of a biological
multimeric protein and the transformations in the header section to
generate the multimer coordinates.  Output will be an :class:`~.AtomGroup`
instance that contains the multimer coordinates.


We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Build a Multimer
-------------------------------------------------------------------------------

Let's build the dimeric form of :pdb:`3enl` of :wiki:`enolase`:

.. ipython:: python

   monomer, header = parsePDB('3enl', header=True)
   monomer

Note that we passed ``header=True`` argument to parse header data in addition
to coordinates.

.. ipython:: python

   showProtein(monomer);
   @savefig structure_analysis_biomolt_monomer.png width=4in
   legend();


Let's get the dimer coordinates using :func:`.buildBiomolecules` function:

.. ipython:: python

   dimer = buildBiomolecules(header, monomer)
   dimer


This function takes biomolecular tarnsformations from the *header* dictionary
(item with key ``'biomoltrans'``) and applies them to the
*monomer*.

.. ipython:: python

   showProtein(dimer);
   @savefig structure_analysis_biomolt_dimer.png width=4in
   legend();


The *dimer* object now has two chains:

.. ipython:: python

   list(dimer.iterChains())


Build a Tetramer
-------------------------------------------------------------------------------


Let's build the tetrameric form of :pdb:`1k4c` of
:wiki:`KcsA_potassium_channel`:

.. ipython:: python

   monomer, header = parsePDB('1k4c', header=True)
   monomer

.. ipython:: python

   showProtein(monomer);
   @savefig structure_analysis_biomolt_monomer.png width=4in
   legend();

Note that we do not want to replicate potassium ions, so we will exclude them:

.. ipython:: python

   potassium = monomer.name_K
   potassium
   without_K = ~ potassium
   without_K

.. ipython:: python

   tetramer = buildBiomolecules(header, without_K)
   tetramer

Now, let's append potassium ions to the tetramer:

.. ipython:: python

   potassium.setChids('K')
   kcsa = tetramer + potassium.copy()
   kcsa.setTitle('KcsA')


Here is a view of the tetramer:

.. ipython:: python

   showProtein(kcsa);
   @savefig structure_analysis_biomolt_dimer.png width=4in
   legend();

Let's get a list of all the chains:

.. ipython:: python

   list(kcsa.iterChains())
