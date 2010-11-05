#!/usr/bin/env python
"""EDA of MD trajectories

An MD trajectory is analyzed. 

Requirements:

- PSF, PDB, and DCD files
- MDAnalysis Python package

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# Import MDAnalysis
import MDAnalysis

# Instantiate a Universe for the simulated system
universe = MDAnalysis.Universe('protein.psf', 'protein.dcd')

# Select atoms of interest
calpha = universe.selectAtoms('name CA')

# Get coordinates of CA atoms
ca_coords = universe.dcd.timeseries(universe_ca, format='fac')

# Instantiate an emseble
ensemble = Ensemble('MD snapshots')
# Add all coordinate sets to ensemble
ensemble.addCoordset(ca_coords)
# Set reference coordinates 
ensemble.setCoordinates(ca_coords[0])
# Perform iterative sueprimposition
ensemble.iterimpose()

# Instantiate EDA and perform calculations
eda = EDA('EDA')
eda.buildCovariance(ensemble)
eda.calcModes()



# Write essential modes into an NMD file for NMWiz
# for this we will need to parse the protein structure as well
prot = parsePDB('protein.pdb')
prot_ca = prot.select('calpha')
writeNMD('md_eda.nmd', eda[:3], prot_ca)

# Let's print fraction of variance for top raking 4 essential modes
for mode in eda[:4]:
    print mode.getFractOfVariance()
