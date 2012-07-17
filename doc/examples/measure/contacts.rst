.. _contacts:

*******************************************************************************
Intermolecular Contacts
*******************************************************************************

Synopsis
===============================================================================

This examples shows how to identify intermolecular contacts, e.g. protein
atoms interacting with a bound inhibitor.  A structure of a protein-ligand 
complex in PDB format will be used.  Output will be :class:`~.Selection` 
instances that points to atoms matching the contact criteria given by the user. 
:class:`~.Selection` instances can be used as input to other
functions for further analysis.

Simple contact selections
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

ProDy selection engine has a powerful feature that enables identifying 
intermolecular contacts very easily. We will see this by identifying protein 
atoms interacting with an inhibitor.

We start with parsing a PDB file that contains a protein and a bound ligand.

>>> pdb = parsePDB('1zz2')

``1zz2`` contains an inhibitor bound p38 MAP kinase structure. Residue name of 
inhibitor is ``B11``. Protein atoms interacting with the inhibitor can simply 
be identified as follows:

>>> contacts = pdb.select('protein and within 4 of resname B11')
>>> contacts
<Selection: 'protein and wit... of resname B11' from 1zz2 (50 atoms)>

``'protein and within 4 of resname B11'`` is interpreted as select protein
atoms that are within 4 A of residue whose name is B11. This selects
protein atoms that within 4 A of the inhibitor. 

Contacts between different atom groups
===============================================================================

In some cases, the protein and the ligand may be in separate files. 
We will imitate this case by making copies of protein and ligand.

>>> inhibitor = pdb.select('resname B11').copy()
>>> inhibitor
<AtomGroup: 1zz2 Selection 'resname B11' (33 atoms)>
>>> protein = pdb.select('protein').copy()
>>> protein
<AtomGroup: 1zz2 Selection 'protein' (2716 atoms)>

We see that inhibitor molecule contains 33 atoms.

Now we have two different atom groups, and we want protein atoms that are 
within 4 Å of the inhibitor.

>>> contacts = protein.select('within 4 of inhibitor', inhibitor=inhibitor)
>>> contacts
<Selection: 'index 227 230 2... 1354 1356 1358' from 1zz2 Selection 'protein' (50 atoms)>

We found that 50 protein atoms are contacting with the inhibitor.
In this case, we passed the atom group *inhibitor* as a keyword argument 
to the selection function. Note that the keyword must match that is used 
in the selection string. 


Composite contact selections
===============================================================================

Now, let's try something more sophisticated. We select Cα atoms of
residues that have at least one atom interacting with the inhibitor:

>>> contacts_ca = protein.select('calpha and (same residue as within 4 of inhibitor)', inhibitor=inhibitor)
>>> contacts_ca
<Selection: 'index 225 232 2... 1328 1351 1359' from 1zz2 Selection 'protein' (20 atoms)>

In this case, ``'calpha and (same residue as within 4 of inhibitor)'`` is 
interpreted as select Cα atoms of residues that have at least
one atom within 4 A of any inhibitor atom.

This shows that, 20 residues have atoms interacting with the inhibitor.

Spherical atom selections
===============================================================================

Similarly, one can give arbitrary coordinate arrays as keyword arguments to 
identify atoms in a spherical region. Let's find backbone atoms within 5 
Å of point (25, 73, 13):

>>> import numpy as np # We will need to pass a Numpy array
>>> sel = protein.select('backbone and within 5 of somepoint', somepoint=np.array((25, 73, 13)))


Fast contact selections
===============================================================================

For repeated and faster contact identification :class:`~.Contacts` class is
recommended.

>>> # We pass the protein as argument
>>> protein_contacts = Contacts(protein)
>>> # The following corresponds to "within 5 of inhibitor"
>>> protein_contacts.select(4, inhibitor)
<Selection: 'index 227 230 2... 1354 1356 1358' from 1zz2 Selection 'protein' (50 atoms)>

This method is 20 times faster than the one in the previous part, but it is
limited to selecting only contacting atoms (other selection arguments cannot be 
passed). Again, it should be noted that :class:`~.Contacts` does not update the 
KDTree that it uses, so it should be used if protein coordinates does not change 
between selections. 
