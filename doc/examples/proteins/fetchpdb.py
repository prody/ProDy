#!/usr/bin/env python
"""fetchPDB - fetch PDB files with their PDB identifiers

Quick access to PDB structures is essential especially when working in
an interactive (ProDy) session. This example shows how to fetch PDB structures 
with their PDB identifiers using :func:`fetchPDB` function.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# The function will return a filename if the download is succesfull.
 
filename = fetchPDB('1p38')
# @> 1p38 downloaded (./1p38.pdb.gz)
print filename
# ./1p38.pdb.gz


# This function also accepts a list of PDB identifiers.

filenames = fetchPDB(['1p38', '1r39', '@!~#'], folder='.')
# @> 1p38 (./1p38.pdb.gz) is found in the target directory.
# @> @!~# is not a valid identifier.
# @> 1r39 downloaded (./1r39.pdb.gz)
# @> PDB download completed (1 found, 1 downloaded, 1 failed).

# It will give you a report of download results and return you a list of 
# filenames.

print filenames
# ['./1p38.pdb.gz', './1r39.pdb.gz', None]

# For failed downloads, ``None`` will be returned.
