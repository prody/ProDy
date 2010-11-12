#!/usr/bin/env python
"""Constructing AtomGroups

This example shows how to construct an AtomGroup instance from scratch.
It is particularly useful for those who intends to design a molecular data 
file parser.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *
import numpy as np


# Instantiate an AtomGroup
#==============================================================================
 
# At instantiation, a name must be given to the AtomGroup. It may be an empty
# sting, i.e. ""
wtr1 = AtomGroup('Water')

# Printing an AtomGroup instance will show you some usefull information.
print wtr1

# Set/get atom attributes
#==============================================================================

# Number of atoms in an AtomGroup is inferred when first time attributes
# of atoms are set
 
# Setting coordinates
wtr1.setCoordinates( np.array( [ [1, 0, 0], [0, 0, 0], [0, 0, 1] ] ) )
wtr1

# In this case we set coordinates using a 3x3 array. Number of atoms
# will automatically be set to 3.

# Setting other attributes
wtr1.setAtomNames( ['H', 'O', 'H'] )
wtr1.setResidueNumbers( [1, 1, 1] )
wtr1.setResidueNames( ['WAT', 'WAT', 'WAT'] )

# Data must always be passed in a list or an array. The length of the container
# must be equal to the number of atoms.

# Accessing data will return a *copy* of the data
print wtr1.getAtomNames()

# Individual atoms
#==============================================================================

# Individual atoms are represented by instance of :class:`prody.proteins.atom.Atom`.

# **Iteration**

# Atoms in an atom group can be iterated over
for a in wtr1: print a

# **Indexing**

# Atoms in an atom group can be accessed via indexing
a = wtr1[0]
# This returns an AtomGroup instance, which provides access to arbitrary 
# groups of atoms in an atom group
print a
print a.getCoordinates()


# Coordinate sets
#==============================================================================

# Let's add a coordinate set to the atom group
wtr1.addCoordset( np.array( [ [0, 1, 0], [0, 0, 0], [0, 0, 1] ] ) )
# Note that number of coordinate sets is now 2, but active coordinate set index
# is still 0
print wtr1
# Active coordinate set incex can be changed for AtomGroups
a.setActiveCoordsetIndex(1)
print a

# Changing active coordinate set for an atom group, does not affect the 
# active coordinate set of the atom group
print wtr1
# Coordinates for the atom group will be returned from the active coordinate set
print a.getCoordinates()

# **Iterations**

# coordinate sets can also be iterated over for atoms and atom groups
for xyz in a.iterCoordsets(): xyz

# Clone atom groups
#==============================================================================

# Now let's make another copy of this water
wtr2 = wtr1.copy()
print wtr2

# **Translate clone**

# let's translate the coordinates of wtr2 so that it does not overlap with wtr1
wtr2.setCoordinates( wtr2.getCoordinates() + 2 )
print wtr2.getCoordinates()
# above operation only translated the coordinate set at index 0
wtr2.setActiveCoordsetIndex(1)
print wtr2.getCoordinates()
wtr2.setCoordinates( wtr2.getCoordinates() + 2 ) # translate the second coordinate set as well

# **Change clone attributes**

# before we merge wtr1 and wtr2, let's change resid's of wtr2
wtr2.setResidueNumbers( [2, 2, 2] )
print wtr2.getResidueNumbers()
# we can do this in an alternate way too
wtr2.select('all').setResidueNumbers(2)
print wtr2.getResidueNumbers()
# note that the following won't work
# wtr2.setResidueNumbers(2)


# Merge atom groups
#==============================================================================

# let's merge two water atom groups
wtrs = wtr1 + wtr2
print wtrs
print wtrs.getCoordinates()
print wtrs.getAtomNames()
print wtrs.getResidueNumbers()

"""
   This hints to why :class:`AtomGroup` instead of Molecule is used. The entire 
   content of a PDB file is not a molecule in strict sense, even when so it is 
   not complete always. If it was, we don't store bond information anyhow. 
   We merely store coordinate, etc. data on some atoms in an AtomGroup.
"""

# Hierarchical view
#==============================================================================

# Hierarchical views of atom groups are represented by 
# :class:`prody.proteins.hierview.HierView`.

# Residues (and also chains) in an atom group can also be iterated over
for res in wtrs.getHierView().iterResidues(): print res

# It's is also possible to change the name of AtomGroup "Water + Copy of Water"
wtrs.setName('2Waters')
print wtrs

