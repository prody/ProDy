.. _attributes:

Storing data in AtomGroup
===============================================================================

Synopsis
-------------------------------------------------------------------------------

This example shows how to store arbitrary atomic attributes in an
:class:`.AtomGroup` instance. Input is protein structure in PDB format and
other atomic data.

Parse structure
-------------------------------------------------------------------------------


.. ipython:: python

   from prody import *
   from pylab import *
   ion()

Let's parse a structure of p38 MAP kinase:

.. ipython:: python

   ag = parsePDB('1p38')
   ag

Set a new attribute
-------------------------------------------------------------------------------

For the purposes of this example, we will manufacture atomic data by
dividing the residue number of each atom by 10:

.. ipython:: python

   myresnum = ag.getResnums() / 10.0

We will add this to the atom group using :meth:`.AtomGroup.setData`
method by passing a name for the attribute and the data:

.. ipython:: python

   ag.setData('myresnum', myresnum)

We can check if a custom atomic attribute exists using
:meth:`.AtomGroup.isDataLabel` method:

.. ipython:: python

   ag.isDataLabel('myresnum')


Access subset of data
-------------------------------------------------------------------------------

Custom attributes can be accessed from selections:

.. ipython:: python

   calpha = ag.calpha
   calpha.getData('myresnum')


Make selections
-------------------------------------------------------------------------------

Custom atomic attributes can be used in selections:

.. ipython:: python

   mysel = ag.select('0 < myresnum and myresnum < 10')
   mysel

This gives the same result as the following selection:

.. ipython:: python

   ag.select('0 < resnum and resnum < 100') == mysel


Save attributes
-------------------------------------------------------------------------------

It is not possible to save custom attributes in PDB files, but
:func:`~.saveAtoms` function can be used them to save in disk for later use:

.. ipython:: python

   saveAtoms(ag)

Let's load it using :func:`~.loadAtoms` function:

.. ipython:: python

   ag = loadAtoms('1p38.ag.npz')
   ag.getData('myresnum')


Delete an attribute
-------------------------------------------------------------------------------

Finally, when done with an attribute, it can be deleted using
:meth:`.AtomGroup.delData` method:

.. ipython:: python

   ag.delData('myresnum')
