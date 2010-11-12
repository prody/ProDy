.. currentmodule:: prody.select

.. _contacts:

*******************************************************************************
Intermolecular Contacts
*******************************************************************************

ProDy selection engine has a powerful feature that enables identifying 
intermolecular contacts very easily. To show this we will identify protein atoms
interacting with an inhibitor.  

Let's start with loading ProDy and parse a PDB file that contains a bound ligand.

>>> from prody import *
>>> pdb = parsePDB('1zz2')

``1zz2`` contains an inhibitor bound p38 MAP kinase structure. Residue name of 
inhibitor is ``B11``. We will make copies of inhibitor and protein. These
copies will imitate atom groups that are parsed from separate files.

>>> inhibitor = pdb.copy('resname B11')
>>> inhibitor
<AtomGroup: Copy of 1zz2 selection "resname B11" (33 atoms; 1 coordinate sets, active set index: 0)>
>>> protein = pdb.copy('protein')
>>> protein
<AtomGroup: Copy of 1zz2 selection "protein" (2716 atoms; 1 coordinate sets, active set index: 0)>

We see that inhibitor molecule contains 33 atoms.

Now, let's say that we want protein atoms that are within 4 angstrom of the 
inhibitor.

>>> contacts = protein.select('within 4 of inhibitor', inhibitor=inhibitor)
>>> contacts
<Selection: "within 4 of inhibitor" from Copy of 1zz2 selection "protein" (50 atoms; 1 coordinate sets, active set index: 0)>

We found that 50 protein atoms are contacting with the inhibitor.
Now, let's try something more sophisticated. We select alpha carbons of
residues that have at least one atom interacting with the inhibitor:

>>> contacts_ca = protein.select('calpha and (same residue as within 4 of inhibitor)', inhibitor=inhibitor)
>>> contacts_ca
<Selection: "calpha and (sam...4 of inhibitor)" from Copy of 1zz2 selection "protein" (20 atoms; 1 coordinate sets, active set index: 0)>

This shows that, 20 residues have atoms interacting with the inhibitor.

Similarly, one can give arbitrary coordinate arrays for identifying atoms
in a spherical region. Let's find backbone atoms within 5 angstroms of point 
(10, 10, 10):

>>> import numpy as np # We will need to pass a Numpy array
>>> sel = protein.select('backbone within 5 of somepoint', somepoint=np.array((10, 10, 10)))

