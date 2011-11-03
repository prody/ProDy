.. currentmodule:: prody.atomic

.. _atomgroup:

*******************************************************************************
Constructing an :class:`AtomGroup`
*******************************************************************************

Synopsis
===============================================================================

This example shows how to construct an :class:`AtomGroup` instance from scratch. 
It is particularly useful for those who intend to design a molecular data file 
parser that returns parsed data in an :class:`AtomGroup` instance.

Input
-------------------------------------------------------------------------------

Atomic data, such as coordinates, atom names, residue names, etc.


Output
-------------------------------------------------------------------------------

Output is an :class:`~prody.atomic.AtomGroup` instance that stores atomic data
and can be used as input to functions and classes for dynamics analysis.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package and the NumPy package:

>>> from prody import *
>>> import numpy as np

Instantiate an AtomGroup
-------------------------------------------------------------------------------

>>> wtr1 = AtomGroup('Water')
>>> wtr1
<AtomGroup: Water (0 atoms; 0 coordinate sets)>


Set coordinates
-------------------------------------------------------------------------------

The best way to start constructing an :class:`AtomGroup` is by setting the
coordinates first. Number of atoms will be automatically set according to
the size of the coordinate data array:

>>> coords = np.array( [ [1, 0, 0], [0, 0, 0], [0, 0, 1] ] )
>>> print( coords )
[[1 0 0]
 [0 0 0]
 [0 0 1]]
>>> wtr1.setCoords( coords )
>>> wtr1
<AtomGroup: Water (3 atoms; 1 coordinate sets, active set index: 0)>

Set attributes
-------------------------------------------------------------------------------

Attributes must be passed in a list or an array whose size is the same
as the number of atoms.

>>> wtr1.setNames( ['H', 'O', 'H'] )
>>> wtr1.setResnums( [1, 1, 1] )
>>> wtr1.setResnames( ['WAT', 'WAT', 'WAT'] )

Accessing data will return a copy of the data:

>>> print( wtr1.getNames() )
['H' 'O' 'H']

Individual atoms
-------------------------------------------------------------------------------

Individual atoms are represented by instance of :class:`Atom`.

**Iteration**

Atoms in an :class:`AtomGroup` can be iterated over

>>> for a in wtr1: a
... 
<Atom: H from Water (index 0; 1 coordinate sets, active set index: 0)>
<Atom: O from Water (index 1; 1 coordinate sets, active set index: 0)>
<Atom: H from Water (index 2; 1 coordinate sets, active set index: 0)>

**Indexing**

Atoms in an atom group can be accessed via indexing:

>>> a = wtr1[0]
>>> a
<Atom: H from Water (index 0; 1 coordinate sets, active set index: 0)>
>>> print( a.getCoords() )
[ 1.  0.  0.]

Coordinate sets
-------------------------------------------------------------------------------

Let's add another coordinate set to the atom group:

>>> wtr1.addCoordset( np.array( [ [0, 1, 0], [0, 0, 0], [0, 0, 1] ] ) )
>>> wtr1
<AtomGroup: Water (3 atoms; 2 coordinate sets, active set index: 0)>

Note that number of coordinate sets is now 2, but active coordinate set index
is still 0. Active coordinate set incex can be changed for :class:`AtomGroup`

>>> a.setACSI(1)
>>> a
<Atom: H from Water (index 0; 2 coordinate sets, active set index: 1)>

Changing active coordinate set for an atom group, does not affect the active 
coordinate set of the atom group:

>>> wtr1
<AtomGroup: Water (3 atoms; 2 coordinate sets, active set index: 0)>

Coordinates for the atom group will be returned from the active coordinate set

>>> print( a.getCoords() )
[ 0.  1.  0.]

**Iterations**

Coordinate sets can also be iterated over for :class:`Atom` and 
:class:`AtomGroup` instances:

>>> for xyz in a.iterCoordsets(): print( xyz )
... 
[ 1.  0.  0.]
[ 0.  1.  0.]

Clone atom groups
-------------------------------------------------------------------------------

Now let's make another copy of this water:

>>> wtr2 = wtr1.copy()
>>> wtr2
<AtomGroup: Copy of Water (3 atoms; 2 coordinate sets, active set index: 0)>

**Translate clone**

Let's translate the coordinates of wtr2 so that it does not overlap with wtr1

>>> wtr2.setCoords( wtr2.getCoords() + 2 )
>>> print( wtr2.getCoords() )
[[ 3.  2.  2.]
 [ 2.  2.  2.]
 [ 2.  2.  3.]]

Above operation only translated the coordinate set at index 0

>>> wtr2.setACSI(1)
>>> print( wtr2.getCoords() )
[[ 0.  1.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  1.]]
>>> wtr2.setCoords( wtr2.getCoords() + 2 ) # translate the second coordinate set as well

**Change clone attributes**

Before we merge wtr1 and wtr2, let's change resid's of wtr2:

>>> wtr2.setResnums( [2, 2, 2] )
>>> print( wtr2.getResnums() )
[2 2 2]

We can do this in an alternate way too:

>>> wtr2.select('all').setResnums(2)
>>> print( wtr2.getResnums() )
[2 2 2]

Note that the following won't work:

>>> wtr2.setResnums(2)
Traceback (most recent call last):
  File "/usr/lib/python2.6/doctest.py", line 1248, in __run
    compileflags, 1) in test.globs
  File "<doctest __main__[29]>", line 1, in <module>
    wtr2.resids = 2
  File "....", line 424, in set_resnums
    if len(resids) != self._n_atoms:
TypeError: object of type 'int' has no len()

Merge atom groups
-------------------------------------------------------------------------------

Let's merge two water atom groups:

>>> wtrs = wtr1 + wtr2
>>> wtrs
<AtomGroup: Water + Copy of Water (6 atoms; 2 coordinate sets, active set index: 0)>
>>> print( wtrs.getCoords() )
[[ 1.  0.  0.]
 [ 0.  0.  0.]
 [ 0.  0.  1.]
 [ 3.  2.  2.]
 [ 2.  2.  2.]
 [ 2.  2.  3.]]
>>> print( wtrs.getNames() )
['H' 'O' 'H' 'H' 'O' 'H']
>>> print( wtrs.getResnums() )
[1 1 1 2 2 2]

.. note::
   This hints to why :class:`AtomGroup` instead of Molecule is used. The entire 
   content of a PDB file is not a molecule in a strict sense, even when so it 
   is not complete always. If it was, we don't store bond information anyhow. 
   We merely store coordinate, etc. data on some atoms in an :class:`AtomGroup`
   instance.

Hierarchical view
-------------------------------------------------------------------------------

Hierarchical views of atom groups are represented by :class:`HierView`.

Residues (and also chains) in an atom group can also be iterated over

>>> for res in wtrs.getHierView().iterResidues(): res
<Residue: WAT 1 from Chain  from Water + Copy of Water (3 atoms; 2 coordinate sets, active set index: 0)>
<Residue: WAT 2 from Chain  from Water + Copy of Water (3 atoms; 2 coordinate sets, active set index: 0)>

Finally, it's is possible to change the name of *wtrs* from 
"Water + Copy of Water" to something shorter:

>>> wtrs.setTitle('2Waters')
>>> wtrs
<AtomGroup: 2Waters (6 atoms; 2 coordinate sets, active set index: 0)>

|questions|

|suggestions|
