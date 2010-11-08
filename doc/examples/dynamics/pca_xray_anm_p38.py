#!/usr/bin/env python
"""PCA of X-ray structures

p38 structures are analyzed and principal modes are compared with ANM modes.

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# We make a list of PDB identifiers. 
# This is NOT the same list that was used in our paper.
# Some structures superceding others lead to changes in PDB identifiers.
# The list was updated to include corresponding PDB identifiers.
pdbids = ['1A9U', '1BL6', '1BL7', '1BMK', '1DI9', '1IAN', '1KV1', '1KV2', '1LEW', '1LEZ', 
          '1M7Q', '1OUK', '1OUY', '1OVE', '1OZ1', '1P38', '1R39', '1R3C', '1W7H', '1W82', 
          '1W83', '1W84', '1WBN', '1WBO', '1WBS', '1WBT', '1WBV', '1WBW', '1WFC', '1YQJ', 
          '1YW2', '1YWR', '1ZYJ', '1ZZ2', '1ZZL', '2BAJ', '2BAK', '2BAL', '2BAQ', '2EWA', 
          '2FSL', '2FSM', '2FSO', '2FST', '2GFS', '2GHL', '2GHM', '2GTM', '2GTN', '2I0H', 
          '2NPQ', '2OKR', '2OZA', '3HVC', '3MH0', '3MH3', '3MH2', '2PUU', '3MGY', '3MH1', 
          '2QD9', '2RG5', '2RG6', '2ZAZ', '2ZB0', '2ZB1', '3BV2', '3BV3', '3BX5', '3C5U', 
          '3L8X', '3CTQ', '3D7Z', '3D83', '2ONL'] 

# Fetch PDB files for these identifiers into folder p38
# fetchPDB function returns a list of filenames
pdbfiles = fetchPDB(pdbids, folder='p38')

# ProDy will print some output to the screen
# Sometimes it is useful to have them in a file
# So, we will start a logfile 
ProDyStartLogfile('p38_pca_anm_calculations') 

# Parse the reference PDB file
ref_structure = parsePDB('p38/1p38.pdb.gz')

# Note that we want to regenerate data  published in the PNAS paper
# hence, we use the same set of atoms (this must be given in the supplementary table)
# copy function makes a copy of selected atoms
ref_structure = ref_structure.copy('resnum 5 to 31 36 to 114 122 to 169 185 to 351 and calpha')

# Eename the reference structure (this way screen logs print shorter)
ref_structure.setName('p38 reference')

# Start an ensemble instance
ensemble = Ensemble('p38 X-ray ensemble')

# Select chain A of the reference structure
ref_chain = ref_structure.getHierView().getChain('A')


# Let's start a list to keep PDB filenames for which mapping to reference failed 
failures = []

# For each PDB file we find matching chain and add it to the ensemble
for pdbfile in pdbfiles:
    
    # Parse next PDB file
    # Note that we are parsing only alpha carbons (it's faster)
    current_structure = parsePDB(pdbfile, subset='calpha')
    
    # For 1yw2, residue numbers are shifted by 1000
    # We reset the numbers to ensure mapping works
    if current_structure.getName() == '1yw2':
        protein = current_structure.select('protein')
        protein.setResidueNumbers( protein.getResidueNumbers() - 1000)
    
    # Get mapping to the reference chain
    current_mapping = mapAtomsToChain(current_structure, ref_chain)
    
    # If mapping fails, add to the failures list
    if not current_mapping:
        failures.append(pdbfile)
        continue
    
    # Get the atom mapping
    current_atommap = current_mapping[0][0]
    
    # Rename the atom mapping
    current_atommap.setName( current_structure.getName())

    # Add the atommap (mapped coordinates) to the ensemble
    # Note that some structures do not completely map
    # so we pass weights (1 for mapped atoms, 0 for unmapped atoms)
    ensemble.addCoordset(current_atommap, weights=current_atommap.getMappedFlags())    

# Let's see if any mapping failed
print failures

# Set the reference coordinates
ensemble.setCoordinates(ref_structure) 
    
# Perform an iterative superimposition
ensemble.iterimpose()

# Let's do some plotting
import matplotlib.pyplot as pl
pl.hist(ensemble.getRMSDs())
# Note that this is slightly different from Fig. S1B. In the paper, different
# sets of atoms were used in RMSD calculation

# Perform PCA
pca = PCA('p38 xray')
pca.buildCovariance(ensemble)
pca.calcModes()

# Write principal modes into an NMD file for NMWiz
writeNMD('p38_principal_modes.nmd', pca[:3], ref_structure)

# Let's print fraction of variance for top raking 6 PCs (or principal modes)
# These numbers are listed in Table 1
for mode in pca[:3]:
    print mode.getFractOfVariance()


# Get ANM instance for reference structure
anm = getANM(ref_structure)
anm.setName('p38 ANM')

# Calculate overlaps between ANM and PCA modes
printOverlapTable(pca[:3], anm[:3])
# These numbers are included in Table 1

# Let's do some cleaning
import os
os.remove('p38_principal_modes.nmd')
#for fn in pdbfiles:
#    os.remove(fn)
#os.rmdir('p38')


# Most of the screen output will be stored in p38_pca_anm_calculations.log file
# It may be useful if you want to see how each individual structure mapped to 
# the reference chain 

# If you want to keep working in the same Python session and close the logfile, try:
ProDyCloseLogfile('p38_pca_anm_calculations')
