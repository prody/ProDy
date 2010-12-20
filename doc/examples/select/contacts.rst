.. currentmodule:: prody.select

.. _contacts:

*******************************************************************************
Intermolecular Contacts
*******************************************************************************

Composite Contact Selections
===============================================================================

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
<Selection: "index 227 to 22...56 to 1356 1358" from Copy of 1zz2 selection "protein" (50 atoms; 1 coordinate sets, active set index: 0)>

We found that 50 protein atoms are contacting with the inhibitor.
Now, let's try something more sophisticated. We select alpha carbons of
residues that have at least one atom interacting with the inhibitor:

>>> contacts_ca = protein.select('calpha and (same residue as within 4 of inhibitor)', inhibitor=inhibitor)
>>> contacts_ca
<Selection: "index 225 to 22...51 to 1351 1359" from Copy of 1zz2 selection "protein" (20 atoms; 1 coordinate sets, active set index: 0)>

This shows that, 20 residues have atoms interacting with the inhibitor.

Similarly, one can give arbitrary coordinate arrays for identifying atoms
in a spherical region. Let's find backbone atoms within 5 angstroms of point 
(25, 73, 13):

>>> import numpy as np # We will need to pass a Numpy array
>>> sel = protein.select('backbone and within 5 of somepoint', somepoint=np.array((25, 73, 13)))


Faster Contact Selections
===============================================================================

For repeated and faster contact identification :class:`Contacts` class is
recommended.

>>> # We pass the protein as argument
>>> protein_contacts = Contacts(protein)
>>> # The following corresponds to "within 5 of inhibitor"
>>> protein_contacts.select(5, inhibitor)
<Selection: "index 226 227 2... 1358 1359 1362" from Copy of 1zz2 selection "protein" (93 atoms; 1 coordinate sets, active set index: 0)>

This method is 20 times faster than the one in the previous part, but it is
limited to selecting only contacting atoms (other selection arguments cannot be 
passed). Again, it should be noted that :class:`Contacts` does not update the 
KDTree that it uses, so it should be used if protein coordinates does not change 
between selections. 
