.. currentmodule:: prody.select

.. _selection-examples:

*******************************************************************************
Atom Selection Examples
*******************************************************************************

Synopsis
===============================================================================

One of the most powerful features of ProDy its atom selection engine.
It enables users to select well defined atom subsets easily by passing
simple keywords or make sophisticated selections by using composite keyword
arguments. This example shows how to make selections.

User Input
-------------------------------------------------------------------------------

Structure of a protein in PDB format or simply as PDB id.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Using keywords
-------------------------------------------------------------------------------

>>> prot = parsePDB('1p38')
>>> prot
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

Most of the time single word keywords may be enough for you to get the set
of atoms that you want:

For selecting protein atoms:

>>> protein = prot.select('protein')
>>> protein
<Selection: "protein" from 1p38 (2833 atoms; 1 coordinate sets, active set index: 0)>

The above shows, 2833 of 2962 atoms are protein atoms. 

For selecting Cα atoms you can use calpha:

>>> calpha = prot.select('calpha')
>>> calpha 
<Selection: "calpha" from 1p38 (351 atoms; 1 coordinate sets, active set index: 0)>

The above shows that there are 351 amino acid residues.

Or, for selecting backbone atoms:

>>> backbone = prot.select('backbone')
>>> backbone
<Selection: "backbone" from 1p38 (1404 atoms; 1 coordinate sets, active set index: 0)>

You can also combine these selection keywords:

>>> calpha_water = prot.select('calpha or water')
>>> calpha_water
<Selection: "calpha or water" from 1p38 (480 atoms; 1 coordinate sets, active set index: 0)>

For a list of keywords see :ref:`selections`. 

Select atoms by name
-------------------------------------------------------------------------------

We can select side-chain atoms as follows:

>>> side_chain_atoms = prot.select('protein and not name N CA C O')
>>> print( side_chain_atoms )
Selection "protein and not name N CA C O" from 1p38
>>> len(side_chain_atoms)
1429

Same selection could also be made using sidechain keyword or 
backbone keyword preceded by *not*:

>>> prot.select('sidechain')
<Selection: "sidechain" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

>>> prot.select('not backbone')
<Selection: "not backbone" from 1p38 (1558 atoms; 1 coordinate sets, active set index: 0)>

Oops, not backbone did not select the same number of atoms. 
Let's try to see why:

>>> print( set(prot.select('not backbone').getResidueNames()) )
set(['CYS', 'ILE', 'VAL', 'GLN', 'LYS', 'HOH', 'PRO', 'THR', 'PHE', 'ASN', 'HIS', 'MET', 'ASP', 'LEU', 'ARG', 'TRP', 'ALA', 'GLU', 'TYR', 'SER'])

Note that we used built-in Python type :class:`set`.

As you can see atoms of **HOH** residues are also included in the selection.

Let's try:

>>> prot.select('not backbone and not water')
<Selection: "not backbone and not water" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

We also used water keyword. This has now worked as sidechain did.
This was to show that it is possible to select same set of atoms in a number 
of different ways. 

Select amino acids by type/name
-------------------------------------------------------------------------------

Let's say we want to select charged residues. We can use resname
keyword followed by 3-letter residue names:  

>>> prot.select('resname ARG LYS HIS ASP GLU')
<Selection: "resname ARG LYS HIS ASP GLU" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

Or, we can use predefined keywords acidic and basic.

>>> charged = prot.select('acidic or basic')
>>> print( charged )
Selection "acidic or basic" from 1p38
>>> len(charged)
906
>>> set(charged.getResidueNames())
set(['HIS', 'ASP', 'LYS', 'GLU', 'ARG'])

Same selection could also be made using charged keyword:

>>> prot.select('charged')
<Selection: "charged" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

More examples
-------------------------------------------------------------------------------

There is much more to what you can do with this flexible and fast atom 
selection engine, without the need for writing nested loops with comparisons 
or changing the source code. See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of :class:`~atomic.Selection`
    objects
  * :ref:`contacts` for selecting interacting atoms


.. doctest::
    :hide:
    
    >>> # Testing the selection parser
    >>> # All tricky selection strings should be contained here
    >>> from prody import *
    >>> import numpy as np
    >>> c = parsePDB('1zz2')
    >>> p = c.copy('protein')
    >>> len(p)
    2716
    >>> len(set(p.getResidueNames()))
    20
    >>> i = c.copy('resname B11')
    >>> len(i)
    33
    >>> len(c.select('(x < 5)'))
    319
    >>> len(c.select('(x <  sqrt(-(- sq(-(-5)))))'))
    319
    >>> len(c.select('(x < 5) and protein and within 10 of water'))
    301
    >>> len(c.select('backbone and within 5 of not protein'))
    649
    >>> len(c.select('backbone and within 5 of not index < 2716'))
    649
    >>> c.select('backbone and same residue as not index < 2716')
    >>> len(c.select('backbone and same residue as within 5 of not protein'))
    1052
    >>> len(p.select('within 4 of inhibitor', inhibitor=i))
    50
    >>> len(c.select('protein and within 4 of resname B11', inhibitor=i))
    50
    >>> len(c.select('exwithin 4 of resname B11', inhibitor=i))
    55
    >>> len(p.select('calpha and (same residue as within 4 of inhibitor)', inhibitor=i))
    20
    >>> len(p.select('backbone and within 5 of somepoint', somepoint=np.array((25, 73, 13))))
    18
    >>> len(p.select('backbone and within 5 of index 1172'))
    21
    >>> len(p.select('backbone and not within 5 of index 1172'))
    1326
    >>> p.select('backbone and not within 5 of not index 1172')
    >>> len(p.select('backbone and within 5 of not index 1172'))
    1347
    >>> len(p.select('backbone and sqrt((x - 25)**2 + (y - 74)**2 + (z - 13)**2) <= 5'))
    12
    >>> len(p.select('sqrt((x - 25)**2 + (y - 74)**2 + (z - 13)**2) <= 5'))
    26
    >>> a = c[1173]
    >>> point = a.getCoordinates()
    >>> len(p.select('within 5 of index 1173'))
    29
    >>> len(c.select('within 5 of index 1173'))
    29
    >>> len(p.select('within 5 of point', point=point))
    29
    >>> len(c.select('within 5 of point', point=point))
    29
    >>> len(c.select('sqrt((x - {0[0]:.3f})**2 + (y - {0[1]:.3f})**2 + (z - {0[2]:.3f})**2) <= 5'.format(point)))
    29
    >>> contacts = Contacts(p)
    >>> print( len(contacts.select(5, i)) == len(p.select('within 5 of inhibitor', inhibitor=i)) )
    True
    >>> print( len(contacts.select(5, np.array((25, 73, 13)))) == len(p.select('within 5 of point', point=np.array((25, 73, 13)))) )
    True
    >>> len(p.select('ca'))
    337
    >>> len(p.select('name "(CA)"'))
    337
    >>> len(p.select('resname "[A-K].*"'))
    1272
    >>> len(p.select('not resname "[A-K].*"'))
    1444
    >>> len(p.select('name "(C|O).*"'))
    2240
    >>> len(p.select('not name "(C|O).*"'))
    476
    >>> len(p.select('x > -0.5 * -5'))
    2540
