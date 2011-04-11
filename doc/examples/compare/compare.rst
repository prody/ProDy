.. currentmodule:: prody.compare

.. _compare-chains:

*******************************************************************************
Compare chains
*******************************************************************************

Synopsis
=============================================================================

:mod:`~prody.compare` module contains functions for matching and mapping 
chains. Results can be used for RMSD fitting and PCA analysis. 

Input
-------------------------------------------------------------------------------

Two structures of the same protein in PDB file format.

Output
-------------------------------------------------------------------------------

Output is :class:`~prody.atomic.AtomMap` instances that can be used as input
to other classes and functions.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Match chains
-------------------------------------------------------------------------------

Matching chains is useful when comparing two chains.
Let's find matching chains in two different HIV-RT structures:
    
>>> bound = parsePDB('1vrt')
>>> unbound = parsePDB('1dlo')
>>> matches = matchChains(bound, unbound)
>>> for match in matches:  
...     print 'Chain 1          :', match[0]      
...     print 'Chain 2          :', match[1]
...     print 'Length           :', len(match[0])
...     print 'Sequence identity:', match[2]
...     print 'Sequence overlap :', match[3]
...     print 'RMSD             :', calcRMSD(match[0], match[1])
...     print ''
Chain 1          : AtomMap Chain B from 1vrt -> Chain B from 1dlo
Chain 2          : AtomMap Chain B from 1dlo -> Chain B from 1vrt
Length           : 400
Sequence identity: 99.2518703242
Sequence overlap : 96
RMSD             : 110.45149192
<BLANKLINE>
Chain 1          : AtomMap Chain A from 1vrt -> Chain A from 1dlo
Chain 2          : AtomMap Chain A from 1dlo -> Chain A from 1vrt
Length           : 524
Sequence identity: 99.0458015267
Sequence overlap : 94
RMSD             : 142.084163869
<BLANKLINE>


We can do a structural alignment using both chains as follows:

First we add the :class:`~prody.atomic.AtomGroup` instances for both 
structures: 

>>> bound_ca = matches[0][0] + matches[1][0]
>>> bound_ca
<AtomMap: (AtomMap Chain B from 1vrt -> Chain B from 1dlo) + (AtomMap Chain A from 1vrt -> Chain A from 1dlo) (from 1vrt; 924 atoms; 924 mapped; 0 unmapped; 1 coordinate sets, active set index: 0)>
>>> unbound_ca = matches[0][1] + matches[1][1]
>>> unbound_ca
<AtomMap: (AtomMap Chain B from 1dlo -> Chain B from 1vrt) + (AtomMap Chain A from 1dlo -> Chain A from 1vrt) (from 1dlo; 924 atoms; 924 mapped; 0 unmapped; 1 coordinate sets, active set index: 0)>
>>> calcRMSD(bound_ca, unbound_ca)
129.34348658001386

Then we find the transformation that minimizes RMSD between these two
selections and apply it to unbound structure:

>>> calcTransformation(unbound_ca, bound_ca).apply(unbound)
<AtomGroup: 1dlo (7691 atoms; 1 coordinate sets, active set index: 0)>
>>> calcRMSD(bound_ca, unbound_ca)
6.0020747465625366

By default, :func:`matchChains` function matches Cα atoms. 
*subset* argument allows for matching larger numbers of atoms. 
We can match backbone atoms as follows:

>>> matches = matchChains(bound, unbound, subset='bb')
>>> for match in matches:  
...     print 'Chain 1          :', match[0]      
...     print 'Chain 2          :', match[1]
...     print 'Length           :', len(match[0])
...     print 'Sequence identity:', match[2]
...     print 'Sequence overlap :', match[3]
...     print 'RMSD             :', calcRMSD(match[0], match[1])
...     print ''
Chain 1          : AtomMap Chain B from 1vrt -> Chain B from 1dlo
Chain 2          : AtomMap Chain B from 1dlo -> Chain B from 1vrt
Length           : 1600
Sequence identity: 99.2518703242
Sequence overlap : 96
RMSD             : 1.71102621571
<BLANKLINE>
Chain 1          : AtomMap Chain A from 1vrt -> Chain A from 1dlo
Chain 2          : AtomMap Chain A from 1dlo -> Chain A from 1vrt
Length           : 2096
Sequence identity: 99.0458015267
Sequence overlap : 94
RMSD             : 7.78386812028
<BLANKLINE>


Or, we can match all atoms as follows:

>>> matches = matchChains(bound, unbound, subset='all')
>>> for match in matches:  
...     print 'Chain 1          :', match[0]      
...     print 'Chain 2          :', match[1]
...     print 'Length           :', len(match[0])
...     print 'Sequence identity:', match[2]
...     print 'Sequence overlap :', match[3]
...     print 'RMSD             :', calcRMSD(match[0], match[1])
...     print ''
Chain 1          : AtomMap Chain B from 1vrt -> Chain B from 1dlo
Chain 2          : AtomMap Chain B from 1dlo -> Chain B from 1vrt
Length           : 3225
Sequence identity: 99.2518703242
Sequence overlap : 96
RMSD             : 2.20947196284
<BLANKLINE>
Chain 1          : AtomMap Chain A from 1vrt -> Chain A from 1dlo
Chain 2          : AtomMap Chain A from 1dlo -> Chain A from 1vrt
Length           : 4159
Sequence identity: 99.0458015267
Sequence overlap : 94
RMSD             : 7.83814068858
<BLANKLINE>


Map onto a chain
-------------------------------------------------------------------------------

Mapping is different from matching. When chains are matched, all matching
atoms are returned as :class:`~prody.atomic.AtomMap` instances. When atoms
are mapped onto a *chain*, missing atoms are replaced by dummy atoms. The 
length of the mapping is equal to the length of *chain*. Mapping is used
particularly useful in assembling coordinate data in analysis of heterogeneous
datasets (see :ref:`pca`).

Let's map bound structure onto unbound chain A (subunit p66):
    
>>> unbound_hv = unbound.getHierView()
>>> unbound_A = unbound_hv['A'] 
>>> mappings = mapOntoChain(bound, unbound_A)
>>> for mapping in mappings:  
...     print 'Mapped chain       :', mapping[0]      
...     print 'Target chain       :', mapping[1]
...     print 'Mapping length     :', len(mapping[0])
...     print '# of mapped atoms  :', mapping[0].getNumOfMapped()
...     print '# of unmapped atoms:', mapping[0].getNumOfUnmapped()
...     print 'Sequence identity  :', mapping[2]
...     print 'Sequence overlap   :', mapping[3]
Mapped chain       : AtomMap Chain B from 1vrt -> Chain A from 1dlo
Target chain       : AtomMap Chain A from 1dlo -> Chain B from 1vrt
Mapping length     : 556
# of mapped atoms  : 524
# of unmapped atoms: 32
Sequence identity  : 99
Sequence overlap   : 94

:func:`mapOntoChain` mapped only Cα atoms. *subset* argument allows for
matching larger numbers of atoms. We can map backbone atoms as follows:

>>> mappings = mapOntoChain(bound, unbound_A, subset='bb')
>>> for mapping in mappings:  
...     print 'Mapped chain       :', mapping[0]      
...     print 'Target chain       :', mapping[1]
...     print 'Mapping length     :', len(mapping[0])
...     print '# of mapped atoms  :', mapping[0].getNumOfMapped()
...     print '# of unmapped atoms:', mapping[0].getNumOfUnmapped()
...     print 'Sequence identity  :', mapping[2]
...     print 'Sequence overlap   :', mapping[3]
Mapped chain       : AtomMap Chain B from 1vrt -> Chain A from 1dlo
Target chain       : AtomMap Chain A from 1dlo -> Chain B from 1vrt
Mapping length     : 2224
# of mapped atoms  : 2096
# of unmapped atoms: 128
Sequence identity  : 99
Sequence overlap   : 94

Or, we can map all atoms as follows:

>>> mappings = mapOntoChain(bound, unbound_A, subset='all') 
>>> for mapping in mappings:  
...     print 'Mapped chain       :', mapping[0]      
...     print 'Target chain       :', mapping[1]
...     print 'Mapping length     :', len(mapping[0])
...     print '# of mapped atoms  :', mapping[0].getNumOfMapped()
...     print '# of unmapped atoms:', mapping[0].getNumOfUnmapped()
...     print 'Sequence identity  :', mapping[2]
...     print 'Sequence overlap   :', mapping[3]
Mapped chain       : AtomMap Chain B from 1vrt -> Chain A from 1dlo
Target chain       : AtomMap Chain A from 1dlo -> Chain B from 1vrt
Mapping length     : 4370
# of mapped atoms  : 4159
# of unmapped atoms: 211
Sequence identity  : 99
Sequence overlap   : 94

|questions|

|suggestions|
