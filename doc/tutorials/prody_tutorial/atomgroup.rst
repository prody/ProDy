.. _atomgroup:

Atom Groups
===============================================================================

Below example shows how to build an :class:`.AtomGroup` from scratch.  We start
by importing everything from the ProDy package and the NumPy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Building an Atom Group
-------------------------------------------------------------------------------

The best way to start constructing an :class:`.AtomGroup` is by setting the
coordinates first. Number of atoms will be automatically set according to
the size of the coordinate data array:

.. ipython:: python

   wtr1 = AtomGroup('Water')
   wtr1
   coords = array([[1, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=float)
   coords
   wtr1.setCoords(coords)
   wtr1

Attributes
^^^^^^^^^^

Attributes must be passed in a list or an array whose size is the same
as the number of atoms.

.. ipython:: python

   wtr1.setNames(['H', 'O', 'H'])
   wtr1.setResnums([1, 1, 1])
   wtr1.setResnames(['WAT', 'WAT', 'WAT'])

Accessing data will return a copy of the data:

.. ipython:: python

   wtr1.getNames()

Atoms
^^^^^

Atoms are represented by instance of :class:`.Atom`.

Iteration
"""""""""

Atoms in an :class:`.AtomGroup` can be iterated over

.. ipython:: python

   for a in wtr1: a


Indexing
""""""""

Atoms in an atom group can be accessed via indexing:

.. ipython:: python

   a = wtr1[0]
   a
   a.getCoords()


Coordinate sets
^^^^^^^^^^^^^^^

Let's add another coordinate set to the atom group:

.. ipython:: python

   wtr1.addCoordset(array([[0, 1, 0], [0, 0, 0], [0, 0, 1]], dtype=float))
   wtr1


Note that number of coordinate sets is now 2, but active coordinate set index
is still 0. Active coordinate set incex can be changed for :class:`.AtomGroup`

.. ipython:: python

   a.setACSIndex(1)
   a

Changing active coordinate set for an atom group, does not affect the active
coordinate set of the atom group:

.. ipython:: python

   wtr1

Coordinates for the atom group will be returned from the active coordinate set

.. ipython:: python

   a.getCoords()

Iterations
""""""""""

Coordinate sets can also be iterated over for :class:`.Atom` and
:class:`.AtomGroup` instances:

.. ipython:: python

   for xyz in a.iterCoordsets(): xyz

Copying and Merging
^^^^^^^^^^^^^^^^^^^

Now let's make another copy of this water:

.. ipython:: python

   wtr2 = wtr1.copy()
   wtr2

Translate copy
""""""""""""""

Let's translate the coordinates of wtr2 so that it does not overlap with wtr1

.. ipython:: python

   wtr2.setCoords(wtr2.getCoords() + 2)
   wtr2.getCoords()

Above operation only translated the coordinate set at index 0

.. ipython:: python

   wtr2.setACSIndex(1)
   wtr2.getCoords()
   wtr2.setCoords(wtr2.getCoords() + 2)  # translate the 2nd coordset as well

Change attributes
"""""""""""""""""

Before we merge wtr1 and wtr2, let's change resid's of wtr2:

.. ipython:: python

   wtr2.setResnums( [2, 2, 2] )
   wtr2.getResnums()

We can do this in an alternate way too:

.. ipython:: python

   wtr2.select('all').setResnums(2)
   wtr2.getResnums()


Note that the following won't work:

.. ipython:: python

   wtr2.setResnums(2)

Merge two copies
""""""""""""""""

Let's merge two water atom groups:

.. ipython:: python

   wtrs = wtr1 + wtr2
   wtrs
   wtrs.getCoords()
   wtrs.getNames()
   wtrs.getResnums()

Hierarchical views
^^^^^^^^^^^^^^^^^^

Hierarchical views of atom groups are represented by :class:`.HierView`.

Residues (and also chains) in an atom group can also be iterated over

.. ipython:: python

   for res in wtrs.getHierView().iterResidues(): res

Renaming an atom group
^^^^^^^^^^^^^^^^^^^^^^

Finally, it's is possible to change the name of *wtrs* from
"Water + Water" to something shorter:

.. ipython:: python

   wtrs.setTitle('2Waters')
   wtrs


.. _attributes:

Storing data in AtomGroup
-------------------------------------------------------------------------------

Now let's get an atom group from a PDB file:

.. ipython:: python

   structure = parsePDB('1p38')

In addition to what's in a PDB file, you can store arbitrary atomic attributes
in :class:`.AtomGroup` objects.

Set a new attribute
^^^^^^^^^^^^^^^^^^^

For the purposes of this example, we will manufacture atomic data by
dividing the residue number of each atom by 10:

.. ipython:: python

   myresnum = structure.getResnums() / 10.0

We will add this to the atom group using :meth:`.AtomGroup.setData`
method by passing a name for the attribute and the data:

.. ipython:: python

   structure.setData('myresnum', myresnum)

We can check if a custom atomic attribute exists using
:meth:`.AtomGroup.isDataLabel` method:

.. ipython:: python

   structure.isDataLabel('myresnum')


Access subset of data
^^^^^^^^^^^^^^^^^^^^^

Custom attributes can be accessed from selections:

.. ipython:: python

   calpha = structure.calpha
   calpha.getData('myresnum')


Make selections
^^^^^^^^^^^^^^^

Custom atomic attributes can be used in selections:

.. ipython:: python

   mysel = structure.select('0 < myresnum and myresnum < 10')
   mysel

This gives the same result as the following selection:

.. ipython:: python

   structure.select('0 < resnum and resnum < 100') == mysel


Save attributes
^^^^^^^^^^^^^^^

It is not possible to save custom attributes in PDB files, but
:func:`.saveAtoms` function can be used them to save in disk for later use:

.. ipython:: python

   saveAtoms(structure)

Let's load it using :func:`.loadAtoms` function:

.. ipython:: python

   structure = loadAtoms('1p38.ag.npz')
   structure.getData('myresnum')


Delete an attribute
^^^^^^^^^^^^^^^^^^^

Finally, when done with an attribute, it can be deleted using
:meth:`.AtomGroup.delData` method:

.. ipython:: python

   structure.delData('myresnum')
