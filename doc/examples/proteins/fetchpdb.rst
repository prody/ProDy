.. currentmodule:: prody.proteins

.. _fetchpdb:

*******************************************************************************
Fetch PDB files
*******************************************************************************

Quick access to PDB structures is essential especially when working in
an interactive (ProDy) session. This example shows how to fetch PDB structures 
with their PDB identifiers using :func:`fetchPDB` function.  

>>> from prody import *

The function will return a filename if the download is succesfull.
 
>>> filename = fetchPDB('1p38')
>>> print filename
1p38.pdb.gz


This function also accepts a list of PDB identifiers.

>>> filenames = fetchPDB(['1p38', '1r39', '@!~#'], folder='.')
>>> print filenames
['1p38.pdb.gz', '1r39.pdb.gz', None]

For failed downloads, ``None`` will be returned (or the list will contain ``None`` item).

ProDy will give you a report of download results and return you a list of 
filenames. The report will be printed to the screen. in this case::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> @!~# is not a valid identifier.
  @> 1r39 downloaded (./1r39.pdb.gz)
  @> PDB download completed (1 found, 1 downloaded, 1 failed).

Such info messages can be turned of using the function :func:`~prody.ProDySetVerbosity`.

>>> ProDySetVerbosity('warning')

We can remove downloaded files as follows:

>>> import os
>>> os.remove('1r39.pdb.gz') 
>>> os.remove('1p38.pdb.gz')
