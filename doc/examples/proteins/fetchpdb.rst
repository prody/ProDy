.. currentmodule:: prody.proteins

.. _fetchpdb:

*******************************************************************************
Fetch PDB files
*******************************************************************************

Synopsis
===============================================================================

This examples demonstrates how to use the flexible PDB fetcher, 
:func:`fetchPDB`. 

Input
-------------------------------------------------------------------------------

Valid inputs are:

  * PDB identifier, e.g ``"2k39"``
  * list of PDB identifiers, e.g. ``["2k39", "1mkp", "1etc"]`` 
  
Output
-------------------------------------------------------------------------------

Compressed PDB files (:file:`pdb.gz`) saved in the current working directory 
or a target folder.
  
ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Single file
-------------------------------------------------------------------------------

The function will return a filename if the download is successful.
 
>>> filename = fetchPDB('1p38')
>>> print( filename ) # doctest: +SKIP
1p38.pdb.gz

Multiple files
-------------------------------------------------------------------------------

This function also accepts a list of PDB identifiers:

>>> filenames = fetchPDB(['1p38', '1r39', '@!~#'], folder='temp')
>>> print( filenames ) # doctest: +SKIP
['1p38.pdb.gz', '1r39.pdb.gz', None]

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

|questions|

|suggestions|

