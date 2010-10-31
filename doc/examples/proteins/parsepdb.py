
from prody import *
 
# Parsing Coordinates
#==============================================================================

# You can parse PDB files by passing a filename (gzipped files are handled)
atoms = parsePDB('1p38.pdb.gz')
# @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.
# PDBParser will tell you what was parsed and how long it took
# This excludes the time spent on reading the file
# The time to evaluate coordinate lines and build an atomgroup is measured

# Parser method returns an AtomGroup instance.
# printing AtomGroup shows some information, such as number of atoms
print atoms


# PDB files can be parser by passing only an identifier.
atoms = parsePDB('1p38')
# @> 1p38 (./1p38.pdb.gz) is found in the target directory.
# @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

atoms = parsePDB('1mkp')
# @> 1mkp downloaded (./1mkp.pdb.gz)
# @> PDBParser: 1183 atoms and 1 coordinate sets were parsed in 0.03s.

# A PDB file will be downloaded if necessary.

# Parsing Header Data
#==============================================================================

# If you also need header data from the PDB file, type in as follows:

atoms, header = parsePDB('1mkp', header=True)
# @> 1mkp downloaded (./1mkp.pdb.gz)
# @> PDBParser: 1183 atoms and 1 coordinate sets were parsed in 0.03s.

# Header data is returned in a dictionary. Printing its keys will show what
# was parsed.

header['experiment']
# 'X-RAY DIFFRACTION'
header['resolution']
# '2.35 ANGSTROMS'
print header.keys()
# ['reference', 'classification', 'compounds', 'resolution', 'title', 
# 'source', 'experiment', 'authors', 'identifier', 'deposition_date']
