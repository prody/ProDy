.. _attributes:

*******************************************************************************
Storing data in AtomGroup instances
*******************************************************************************

Synopsis
===============================================================================

This example shows how to store arbitrary atomic attributes in an 
:class:`~.AtomGroup` instance. 


Input
-------------------------------------------------------------------------------

A protein structure in PDB format and other atomic data.


Output
-------------------------------------------------------------------------------

Output is an :class:`~.AtomGroup` instance that stores atomic data
and can be used as input to functions and classes for dynamics analysis.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package and the NumPy package:

>>> from prody import *

Read protein structure
-------------------------------------------------------------------------------

Let's parse a structure of p38 MAP kinase:

>>> ag = parsePDB('1p38')
>>> ag
<AtomGroup: 1p38 (2962 atoms)>

Set a new attribute
-------------------------------------------------------------------------------

For the purposes of this example, we will manufacture atomic data by
dividing the residue number of each atom by 10:

>>> myresnum = ag.getResnums() / 10.0

We will add this to the atom group using :meth:`.AtomGroup.setData`
method by passing a name for the attribute and the data:

>>> ag.setData('myresnum', myresnum)

We can check if a custom atomic attribute exists using 
:meth:`.AtomGroup.isData` method:

>>> ag.isData('myresnum')
True


Access data from selections
-------------------------------------------------------------------------------

Custom attributes can be accessed from selections:

>>> calpha = ag.calpha
>>> print( calpha.getData('myresnum') ) # doctest: +ELLIPSIS
[  0.4   0.5   0.6   0.7   0.8   0.9   1.    1.1   1.2   1.3   1.4   1.5
   1.6   1.7   1.8   1.9   2.    2.1   2.2   2.3   2.4   2.5   2.6   2.7
  ...
  34.   34.1  34.2  34.3  34.4  34.5  34.6  34.7  34.8  34.9  35.   35.1
  35.2  35.3  35.4]



Use attribute in atom selections
-------------------------------------------------------------------------------

Custom atomic attributes can be used in selections:

>>> ag.select('0 < myresnum and myresnum < 10')
<Selection: "0 < myresnum and myresnum < 10" from 1p38 (788 atoms)>

This gives the same result as the following selection:

>>> ag.select('0 < resnum and resnum < 100') == ag.select('0 < myresnum and myresnum < 10') 
True


Save attributes
-------------------------------------------------------------------------------

It is not possible to save custom attributes in PDB files, but 
:func:`~.saveAtoms` function can be used them to save in disk for later use:

>>> saveAtoms(ag)
'1p38.ag.npz'

Let's load it using :func:`~.loadAtoms` function:

>>> ag = loadAtoms('1p38.ag.npz')
>>> ag.getData('myresnum')
array([  0.4,   0.4,   0.4, ...,  77.1,  77.3,  77.6])


Delete an attribute
-------------------------------------------------------------------------------

Finally, when done with an attribute, it can be deleted using 
:meth:`.AtomGroup.delData` method:

>>> ag.delData('myresnum')
array([  0.4,   0.4,   0.4, ...,  77.1,  77.3,  77.6])

|questions|

|suggestions|
