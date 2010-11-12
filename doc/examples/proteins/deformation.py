#!/usr/bin/env python
"""getDeformVector - 


"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# Let's parse two p38 MAP Kinase structures: 1p38 and 1zz2

reference = parsePDB('1p38')
mobile = parsePDB('1zz2')   # this is the one we want to superimpose

# Let's find matching chains in these structures
matches = matchChains(reference, mobile)

# matchChains() function returns a list
# if there are no matching chains, list is empty
# else, list contains a tuple for each pair of matching chains

# Let's get the first match (there is only one match in this case)
match = matches[0]

# First item is a selection of atoms from the first structure (reference)
ref_chain = match[0]
# Second item is a selection of atoms from the second structure (mobile)
mob_chain = match[1]


# RMSD and superimposition
#==============================================================================

print getRMSD(ref_chain, mob_chain)
# this prings 72.93 (angstroms)

# Let's find the transformation that minimizes RMSD between these chains
t = getTransformation(mob_chain, ref_chain)

# We apply this transformation to mobile structure (not to mob_chain, to preserve structures integrity)
t.apply(mobile)

print getRMSD(ref_chain, mob_chain)
# this prints 1.86 (angstroms)

print len(ref_chain)
# this prints 337, all being alpha carbons

# Deformation vector
#==============================================================================

# Let's get the deformation vector
defvec = getDeformVector(ref_chain, mob_chain)

# abs(defvec) returns the magnitude of the deformation
print abs(defvec)
# this prints 34.20 (angstroms)

# RMSD can be calculated from the magnitude of the deformation vector
print (abs(defvec)**2 / len(ref_chain)) ** 0.5
# this prints 1.86 (angstroms)

# Array of numbers for this deformation can be obtained as follows
arr = defvec.getArray() # arr is a NumPy array
print arr

# Following yields the normalized deformation vector
defvecnormed = defvec.getNormed()
print abs(defvecnormed)


# Compare with ANM modes
#==============================================================================
import numpy as np
# Let's get ANM model for the reference chain
anm = getANM(ref_chain)

# Calculate overlap between slowest ANM mode and the deformation vector
print anm[0] * defvecnormed # note that we used normalized deformation vector

# We can do this for a set of ANM modes (slowest 6) as follows
print np.array( anm[:6].getModes() ) * defvecnormed
