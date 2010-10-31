.. module:: prody.proteins

.. _parsepdb:

*******************************************************************************
Parse PDB coordinate/header data
*******************************************************************************

PDB files can be parsed using :func:`parsePDB` function. This function
returns an :class:`prody.proteins.atomgroup.AtomGroup` instance that contains
coordinate data, and optionally a :class:`dict` instance that contains part of 
the data from header section.


.. literalinclude:: parsepdb.py
