#!/usr/bin/env python
"""PCA of NMR models

Ubiquitin models are analyzed for dominant modes of structural variation.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# We parse only CA atoms in this case 
ubi = parsePDB('2k39', subset='calpha')

# We use residues 1 to 70. 71 and above are disordered
ubi = ubi.copy('resnum < 71')

# Start an ensemble instance
ensemble = Ensemble('Ubiquitin NMR ensemble')


# Set the reference coordinates (this sets coordinates from model 1 as the reference)
ensemble.setCoordinates( ubi.getCoordinates() )

# This gets coordinate sets from all models and adds them to the ensemble
ensemble.addCoordset( ubi.getCoordsets() ) 
    
# Perform an iterative superimposition
ensemble.iterimpose()

# Let's do some plotting
import matplotlib.pyplot as pl
pl.hist(ensemble.getRMSDs())

# Perform PCA
pca = PCA('ubi PCA')
pca.buildCovariance(ensemble)
pca.calcModes()

# Write principal modes into an NMD file for NMWiz
writeNMD('ubi_pca.nmd', pca[:3], ubi)

# Let's print fraction of variance for top raking 4 PCs (or essential modes)
# These numbers are listed in the table S3
for mode in pca[:4]:
    print mode.getFractOfVariance()

# We set the active coordinate set to 79, which is the one that is closest 
# to the mean structure (note that indices start from 1 in Python)
ubi.setActiveCoordsetIndex(78)

# Get ANM instance for the active coordset (model 79)
anm = getANM(ubi)
anm.setName('ubi ANM')

# Calculate overlaps between ANM and PCA modes
printOverlapTable(pca[:4], anm[:4])
# These numbers are included in Table 1

# Let's do some cleaning
import os
os.remove('ubi_pca.nmd')
