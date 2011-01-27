.. currentmodule:: prody.proteins

.. _parsepdb:

*******************************************************************************
Parse PDB data
*******************************************************************************

PDB files can be parsed using :func:`parsePDB` function.

|more| See also :ref:`pdbparser-performance`.

>>> from prody import *
 
Coordinate data
===============================================================================

You can parse PDB files by passing a filename (gzipped files are handled):

>>> fetchPDB('1p38') # doctest: +SKIP
'1p38.pdb.gz'
>>> atoms = parsePDB('1p38.pdb.gz')
>>> atoms
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

PDBParser will tell you what was parsed and how long it took. The following 
will be printed to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

This excludes the time spent on reading the file.
The time to evaluate coordinate lines and build an atomgroup is measured.


Such info messages can be turned of using the function :func:`~prody.ProDySetVerbosity`.

>>> ProDySetVerbosity('warning')

PDB files can be parsed by passing only an identifier. A PDB file will be 
downloaded if necessary (see also :ref:`fetchpdb`).

>>> atoms = parsePDB('1mkp')
>>> atoms
<AtomGroup: 1mkp (1183 atoms; 1 coordinate sets, active set index: 0)>

Parser method returns an :class:`~prody.atomic.AtomGroup` instance.

:more: See our report on :ref:`pdbparser-performance`.


Header data
===============================================================================

If you also need header data from the PDB file, type in as follows:

>>> atoms, header = parsePDB('1mkp', header=True)

Header data is returned in a :class:`dict`. Printing its keys will show what
was parsed.

>>> header['experiment']
'X-RAY DIFFRACTION'
>>> header['resolution']
'2.35 ANGSTROMS'
>>> print header.keys() # doctest: +SKIP
['reference', 'classification', 'compounds', 'resolution', 'title', 
'source', 'experiment', 'authors', 'identifier', 'deposition_date']

