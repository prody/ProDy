.. _pca-xray-calculations:

*******************************************************************************
PCA of X-ray structures: Calculations
*******************************************************************************

Synopsis
===============================================================================

This is the first part of a lengthy ProDy example. The aim is to repeat the 
calculations for p38 MAP kinase (MAPK) that were published in [AB09]_. 
In this part, we perform the calculations using the same p38 MAPK dataset.

Input
-------------------------------------------------------------------------------

A set of structures of the same protein. Specifying a list of PDB identifiers 
is sufficient.

Output
-------------------------------------------------------------------------------

A :class:`~.PCA` instance that stores the covariance matrix and principal modes
describing the dominant changes in the dataset. The :class:`~.PCA` instance
and principal modes (:class:`~.Mode`) can be used as input for the functions in 
:mod:`~.dynamics` module.


ProDy Code
===============================================================================
  
We start by importing everything from the ProDy package:

>>> from prody import *

Gather dataset
-------------------------------------------------------------------------------

We use a list of PDB identifiers for the structures that we want to 
include in our analysis.

>>> pdbids = ['1A9U', '1BL6', '1BL7', '1BMK', '1DI9', '1IAN', '1KV1', '1KV2', '1LEW', '1LEZ', 
...           '1M7Q', '1OUK', '1OUY', '1OVE', '1OZ1', '1P38', '1R39', '1R3C', '1W7H', '1W82', 
...           '1W83', '1W84', '1WBN', '1WBO', '1WBS', '1WBT', '1WBV', '1WBW', '1WFC', '1YQJ', 
...           '1YW2', '1YWR', '1ZYJ', '1ZZ2', '1ZZL', '2BAJ', '2BAK', '2BAL', '2BAQ', '2EWA', 
...           '2FSL', '2FSM', '2FSO', '2FST', '2GFS', '2GHL', '2GHM', '2GTM', '2GTN', '2I0H', 
...           '2NPQ', '2OKR', '2OZA', '3HVC', '3MH0', '3MH3', '3MH2', '2PUU', '3MGY', '3MH1', 
...           '2QD9', '2RG5', '2RG6', '2ZAZ', '2ZB0', '2ZB1', '3BV2', '3BV3', '3BX5', '3C5U', 
...           '3L8X', '3CTQ', '3D7Z', '3D83', '2ONL']


Note that we used a list of identifiers that are different from what was 
listed in the supporting material of [AB09]_. 
Since the paper was published, the Protein Data Bank has refined some
of the structures  and changed their identifiers. 
These changes are reflected to the above list.
  
Also note that, it is possible to update this list to include all of the p38
structures currently available in the PDB using the 
:func:`~.blastPDB` function as follows: 
 
>>> p38_sequence ='''GLVPRGSHMSQERPTFYRQELNKTIWEVPERYQNLSPVGSGAYGSVCAAFDTKTGHRVAVKKLSRPFQS
... IIHAKRTYRELRLLKHMKHENVIGLLDVFTPARSLEEFNDVYLVTHLMGADLNNIVKCQKLTDDHVQFLIYQILRGLKYI
... HSADIIHRDLKPSNLAVNEDCELKILDFGLARHTDDEMTGYVATRWYRAPEIMLNWMHYNQTVDIWSVGCIMAELLTGRT
... LFPGTDHIDQLKLILRLVGTPGAELLKKISSESARNYIQSLAQMPKMNFANVFIGANPLAVDLLEKMLVLDSDKRITAAQ
... ALAHAYFAQYHDPDDEPVADPYDQSFESRDLLIDEWKSLTYDEVISFVPPPLDQEEMES''' 
>>> # blast_record = blastPDB('''p38_sequence''')
>>> # pdbids = blast_record.getHits() # uncomment this and previous line to update PDB list

We use the same set of structures to reproduce the results.
After we listed the PDB identifiers, we obtain them using 
:func:`~.fetchPDB` function as follows:
 
>>> pdbfiles = fetchPDB(pdbids, folder='pdbfiles', compressed=False)
  
``pdbfiles`` variable contains the list of filenames of obtained PDB structures.

Set reference chain
-------------------------------------------------------------------------------

The next step is setting one of the p38 structures as the reference
structure. We use 1p38 chain A. Note that we won't use
all of the resolved residues in this structure. We select only those residues
which are resolved in at least 90% of the dataset. 

>>> ref_structure = parsePDB('pdbfiles/1p38.pdb')
>>> ref_structure = ref_structure.copy('resnum 5 to 31 36 to 114 122 to 169 185 to 351 and calpha')
>>> ref_structure.setTitle('p38 reference')

Select chain A from the reference structure

>>> ref_chain = ref_structure.getHierView().getChain('A')
>>> ref_chain
<Chain: A from p38 reference (321 residues, 321 atoms)>

We use the :func:`~.parsePDB` function to parse a PDB file.
This returns a :class:`~.AtomGroup` instance. We make a copy
of α-carbon atoms of select residues for analysis.   

|more| See :ref:`selections` for making selections.

Prepare ensemble
-------------------------------------------------------------------------------

X-ray structural ensembles are heterogenous, i.e. different structures
have different sets of unresolved residues. Hence, it is not straightforward
to analyzed them as it would be for NMR models (see :ref:`pca-nmr`). 

ProDy has special functions and classes for facilitating efficient analysis
of the PDB X-ray data. In this example we use :func:`~.mapOntoChain` 
function which returns an :class:`~.AtomMap` instance.

|more| See :ref:`atommaps` for more details.   

Start a logfile to save screen output: 

>>> startLogfile('p38_pca') 

Instantiate an :class:`~.PDBEnsemble` object:
  
>>> ensemble = PDBEnsemble('p38 X-ray')
  
Set the reference coordinates:

>>> ensemble.setCoords(ref_chain) 
      
For each PDB file, we find the matching chain and add it to the ensemble:

>>> for pdbfile in pdbfiles:
...     # Parse next PDB file. (only alpha carbons, since it's faster)
...     structure = parsePDB(pdbfile, subset='calpha')
...     # Get mapping to the reference chain
...     mappings = mapOntoChain(structure, ref_chain)
...     atommap = mappings[0][0]
...     # Add the atommap (mapped coordinates) to the ensemble
...     # Note that some structures do not completely map (missing residues)
...     # so we pass weights (1 for mapped atoms, 0 for unmapped atoms)
...     ensemble.addCoordset(atommap, weights=atommap.getMappedFlags())    

>>> ensemble
<PDBEnsemble: p38 X-ray (75 conformations; 321 atoms)>
>>> len(ensemble) == len(pdbfiles)
True

Perform an iterative superimposition:

>>> ensemble.iterpose()

Close the logfile (file content shows how chains were paired/mapped):

>>> closeLogfile('p38_pca')

Save coordinates
-------------------------------------------------------------------------------

We use :class:`~.PDBEnsemble` to store coordinates of the X-ray 
structures. The :class:`~.PDBEnsemble` instances do not store any 
other atomic data. If we want to write aligned coordinates into a file, we 
need to pass the coordinates to an :class:`~.AtomGroup` instance.
Then we use :func:`~.writePDB` function to save coordinates:

>>> xray_coords = ref_structure.copy()
>>> xray_coords.delCoordset(0) # Delete existing coordinate set
>>> xray_coords.addCoordset( ensemble.getCoordsets() )
>>> writePDB('p38_xray_coors.pdb', xray_coords)
'p38_xray_coors.pdb'


PCA calculations
-------------------------------------------------------------------------------

Once the coordinate data is prepared, it is straightforward to perform the 
:class:`~.PCA` calculations:

>>> pca = PCA('p38 xray')           # Instantiate a PCA instance
>>> pca.buildCovariance(ensemble)   # Build covariance for the ensemble
>>> pca.calcModes()                 # Calculate modes (20 of the by default)

**Approximate method**

In the following we are using singular value decomposition for faster 
and more memory efficient calculation of principal modes:

>>> pca_svd = PCA('p38 svd')
>>> pca_svd.performSVD(ensemble)

The resulting eigenvalues and eigenvectors may show small differences due to
missing atoms in the datasets:

>>> '%.2f' % abs(pca_svd.getEigenvalues()[:20] - pca.getEigenvalues()).max()
'0.40'
>>> '%.3f' % abs(calcOverlap(pca, pca_svd).diagonal()[:20]).min()
'0.998'

Note that building and diagonalizing the covariance matrix is the preferred
method for heterogeneous ensembles. For NMR models or MD trajectories SVD 
method may be preferred over covariance method.

ANM calculations
-------------------------------------------------------------------------------

To perform :class:`~.ANM` calculations:

>>> anm = ANM('1p38')             # Instantiate a ANM instance
>>> anm.buildHessian(ref_chain)   # Build Hessian for the reference chain  
>>> anm.calcModes()               # Calculate slowest non-trivial 20 modes 

Save your work
-------------------------------------------------------------------------------

Calculated data can be saved in a ProDy internal format
to use in a later session or to share it with others.

If you are in an interactive Python session, and wish to continue without
leaving your session, you do not need to save the data. Saving data is useful
if you want to use it in another session or at a later time, or if you want
to share it with others.

>>> saveModel(pca)
'p38_xray.pca.npz'
>>> saveModel(anm)
'1p38.anm.npz'
>>> saveEnsemble(ensemble)
'p38_X-ray.ens.npz'
>>> writePDB('p38_ref_chain.pdb', ref_chain)
'p38_ref_chain.pdb'

We use the :func:`~.saveModel` and :func:`~.saveEnsemble` functions to save 
calculated data. In :ref:`pca-xray-analysis`, we will use the 
:func:`~.loadModel` and :func:`~.loadEnsemble` functions to load the data.

See Also
===============================================================================

This example is continued in :ref:`pca-xray-analysis` 

|questions|

|suggestions|
