.. _pdbfiles:

PDB files
===============================================================================

This examples demonstrates how to use the flexible PDB fetcher,
:func:`.fetchPDB`. Valid inputs are PDB identifier, e.g :pdb:`2k39`, or a list
of PDB identifiers, e.g. ``["2k39", "1mkp", "1etc"]``.
Compressed PDB files (:file:`pdb.gz`) will be saved to the current working
directory or a target folder.


.. _fetchpdb:

Fetch PDB files
-------------------------------------------------------------------------------

Single file
^^^^^^^^^^^

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *


The function will return a filename if the download is successful.

.. ipython:: python

   filename = fetchPDB('1p38')
   filename

Multiple files
^^^^^^^^^^^^^^

This function also accepts a list of PDB identifiers:

.. ipython::

   In [1]: mkdir temp
   In [2]: filenames = fetchPDB(['1p38', '1r39', '@!~#'], folder='temp')
   In [3]: filenames

For failed downloads, ``None`` will be returned (or the list will contain
``None`` item).


Also note that in this case we passed a folder name. Files are saved in
this folder, after it is created if it did not exist.


ProDy will give you a report of download results and return a list of
filenames. The report will be printed on the screen, which in this case would
be::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> @!~# is not a valid identifier.
  @> 1r39 downloaded (./1r39.pdb.gz)
  @> PDB download completed (1 found, 1 downloaded, 1 failed).


.. _parsepdb:


Parse PDB files
-------------------------------------------------------------------------------

ProDy offers a fast and flexible PDB parser, :func:`.parsePDB`.
Parser can be used to read well defined subsets of atoms, specific chains or
models (in NMR structures) to boost the performance. This example shows how to
use the flexible parsing options.

Three types of input are accepted from user:

  * PDB file path, e.g. ``"../1MKP.pdb"``
  * compressed (gzipped) PDB file path, e.g. ``"1p38.pdb.gz"``
  * PDB identifier, e.g. :pdb:`2k39`


Output is an :class:`.AtomGroup` instance that stores atomic data
and can be used as input to functions and classes for dynamics analysis.

Parse a file
^^^^^^^^^^^^

You can parse PDB files by passing a filename (gzipped files are handled).
We do so after downloading a PDB file (see :ref:`fetchpdb` for more
information):

.. ipython:: python

   fetchPDB('1p38') # doctest: +SKIP
   atoms = parsePDB('1p38.pdb.gz')
   repr(atoms)

Parser returns an :class:`.AtomGroup` instance.

Also note that the time it took to parse the file is printed on
the screen. This includes the time that it takes to evaluate
coordinate lines and build an :class:`.AtomGroup` instance and
excludes the time spent on reading the file from disk.

Use an identifier
^^^^^^^^^^^^^^^^^

PDB files can be parsed by passing simply an identifier. arser will look for a
PDB file that matches the given identifier in the current working directory.
If a matching file is not found, ProDy will downloaded it from PDB FTP server
automatically and saved it in the current working directory.

.. ipython:: python

   atoms = parsePDB('1mkp')
   repr(atoms)


Subsets of atoms
^^^^^^^^^^^^^^^^

Parser can be used to parse backbone or Cα atoms:

.. ipython:: python

   backbone = parsePDB('1mkp', subset='bb')
   repr(backbone)
   calpha = parsePDB('1mkp', subset='ca')
   repr(calpha)


Specific chains
^^^^^^^^^^^^^^^

Parser can be used to parse a specific chain from a PDB file:

.. ipython:: python

   chA = parsePDB('3mkb', chain='A')
   repr(chA)
   chC = parsePDB('3mkb', chain='C')
   repr(chC)

Multiple chains can also be parsed in the same way:

.. ipython:: python

   chAC = parsePDB('3mkb', chain='AC')
   repr(chAC)


Specific models
^^^^^^^^^^^^^^^

Parser can be used to parse a specific model from a file:

.. ipython:: python

   model1 = parsePDB('2k39', model=10)
   repr(model1)

Alternate locations
^^^^^^^^^^^^^^^^^^^

When a PDB file contains alternate locations for some of the atoms, by default
alternate locations with indicator ``A`` are parsed.

.. ipython:: python

   altlocA = parsePDB('1ejg')
   repr(altlocA)

Specific alternate locations can be parsed as follows:

.. ipython:: python

   altlocB = parsePDB('1ejg', altloc='B')
   repr(altlocB)

Note that in this case number of atoms are different between the two atom
groups. This is because the residue types of atoms with alternate locations
are different.

Also, all alternate locations can be parsed as follows:

.. ipython:: python

   all_altlocs = parsePDB('1ejg', altloc=True)
   repr(all_altlocs)

Note that this time parser returned three coordinate sets. One for each
alternate location indicator found in this file (A, B, C). When parsing
multiple alternate locations, parser will expect for the same residue type
for each atom with an alternate location. If residue names differ, a warning
message will be printed.

Composite arguments
^^^^^^^^^^^^^^^^^^^

Parser can be used to parse coordinates from a specific model for a subset of
atoms of a specific chain:

.. ipython:: python

   composite = parsePDB('2k39', model=10, chain='A', subset='ca')
   repr(composite)

Header data
^^^^^^^^^^^

PDB parser can be used to extract header data in a :class:`dict` from PDB
files as follows:

.. ipython:: python

   atoms, header = parsePDB('1mkp', header=True)
   list(header)
   header['experiment']
   header['resolution']

It is also possible to parse only header data by passing `model=0` as an
argument:

.. ipython:: python

   header = parsePDB('1mkp', header=True, model=0)

or using :func:`.parsePDBHeader` function:

.. ipython:: python

   header = parsePDBHeader('1mkp')


.. _writepdb:

Write PDB file
-------------------------------------------------------------------------------

PDB files can be written using :func:`.writePDB` function. This
example shows how to write PDB files for :class:`.AtomGroup`
instances and subsets of atoms.

Write all atoms
^^^^^^^^^^^^^^^

All atoms in an :class:`.AtomGroup` can be written in PDB format
as follows:

.. ipython:: python

   writePDB('MKP3.pdb', atoms)

Upon successful writing of PDB file, filename is returned.

Write a subset
^^^^^^^^^^^^^^

It is also possible to write subsets of atoms in PDB format:

.. ipython:: python

   alpha_carbons = atoms.select('calpha')
   writePDB('1mkp_ca.pdb', alpha_carbons)
   backbone = atoms.select('backbone')
   writePDB('1mkp_bb.pdb', backbone)
