.. _compare-chains:

Structure Comparison
===============================================================================

This section shows how to find identical or similar protein chains in two
structures files and align them.

:mod:`.proteins` module contains functions for matching and mapping
chains. Results can be used for RMSD fitting and PCA analysis.

Output will be :class:`.AtomMap` instances that can be used as input
to ProDy classes and functions.


Match chains
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Matching chains is useful when comparing two structures.  We will find
matching chains in two different HIV :wiki:`Reverse Transcriptase` structures.


First we define a function that prints information on paired (matched) chains:

.. ipython:: python

   def printMatch(match):
       print('Chain 1     : {}'.format(match[0]))
       print('Chain 2     : {}'.format(match[1]))
       print('Length      : {}'.format(len(match[0])))
       print('Seq identity: {}'.format(match[2]))
       print('Seq overlap : {}'.format(match[3]))
       print('RMSD        : {}\n'.format(calcRMSD(match[0], match[1])))

Now let's parse bound RT structure :pdb:`1vrt` and unbound structure
:pdb:`1dlo`:

.. ipython:: python

   bound = parsePDB('1vrt')
   unbound = parsePDB('1dlo')

Let's verify that these structures are not aligned:

.. ipython:: python

   showProtein(unbound, bound);
   @savefig structure_analysis_compare_before.png width=4in
   legend();


We find matching chains as follows:

.. ipython:: python

   matches = matchChains(bound, unbound)
   for match in matches:
       printMatch(match)

This resulted in two matches. Chains A and B of two structures are paired.
These chains contain only Cα atoms:

.. ipython:: python

   match[0][0].iscalpha
   match[0][1].iscalpha


For a structural alignment based on both chains, we merge these matches as
follows:

.. ipython:: python

   bound_ca = matches[0][0] + matches[1][0]
   bound_ca
   unbound_ca = matches[0][1] + matches[1][1]
   unbound_ca

Let's calculate RMSD:

.. ipython:: python

   calcRMSD(bound_ca, unbound_ca)

We find the transformation that minimizes RMSD between these two
selections and apply it to unbound structure:

.. ipython:: python

   calcTransformation(unbound_ca, bound_ca).apply(unbound);
   calcRMSD(bound_ca, unbound_ca)

Let's see the aligned structures now:

.. ipython:: python

   showProtein(unbound, bound);
   @savefig structure_analysis_compare_after.png width=4in
   legend();

By default, :func:`.matchChains` function matches Cα atoms.
*subset* argument allows for matching larger numbers of atoms.
We can match backbone atoms as follows:

.. ipython:: python

   matches = matchChains(bound, unbound, subset='bb')
   for match in matches:
       printMatch(match)


Or, we can match all atoms as follows:

.. ipython:: python

   matches = matchChains(bound, unbound, subset='all')
   for match in matches:
       printMatch(match)


Map onto a chain
-------------------------------------------------------------------------------

Mapping is different from matching. When chains are matched, all matching
atoms are returned as :class:`.AtomMap` instances. When atoms
are mapped onto a *chain*, missing atoms are replaced by dummy atoms. The
length of the mapping is equal to the length of *chain*. Mapping is used
particularly useful in assembling coordinate data in analysis of heterogeneous
datasets (see :ref:`pca`).

Let's map bound structure onto unbound chain A (subunit p66):


.. ipython:: python

   def printMapping(mapping):
       print('Mapped chain     : {}'.format(mapping[0]))
       print('Target chain     : {}'.format(mapping[1]))
       print('Mapping length   : {}'.format(len(mapping[0])))
       print('# of mapped atoms: {}'.format(mapping[0].numMapped()))
       print('# of dummy atoms : {}'.format(mapping[0].numDummies()))
       print('Sequence identity: {}'.format(mapping[2]))
       print('Sequence overlap : {}\n'.format(mapping[3]))

   unbound_hv = unbound.getHierView()
   unbound_A = unbound_hv['A']
   mappings = mapOntoChain(bound, unbound_A)
   for mapping in mappings:
       printMapping(mapping)

:func:`.mapOntoChain` mapped only Cα atoms. *subset* argument allows for
matching larger numbers of atoms. We can map backbone atoms as follows:

.. ipython:: python

   mappings = mapOntoChain(bound, unbound_A, subset='bb')
   for mapping in mappings:
       printMapping(mapping)

Or, we can map all atoms as follows:

.. ipython:: python

   mappings = mapOntoChain(bound, unbound_A, subset='all')
   for mapping in mappings:
       printMapping(mapping)