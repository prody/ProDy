.. _tutorial:

.. module:: prody

*******************************************************************************
Tutorial
*******************************************************************************

ProDy can be used in two ways:

  * by writing scripts for automated tasks
  * interactively in a Python interpreter
  
First part of the tutorial shows how to use one of the scripts that is 
distributed with ProDy. Following parts assume user is in an interactive
Python shell, but should be useful for those who may want to write their own
scripts for automatizing certain analysis tasks. 

Upon completion of the tutorial, user is referred to :ref:`examples` section 
that gives more specific and comprehensive examples. Also, :ref:`howtos` section 
that elaborates on important features of ProDy is recommended for further 
reading.

Using ProDy scripts
===============================================================================

ProDy scripts come with source distribution which can be downloaded from 
http://pypi.python.org/pypi/ProDy/ (tar.gz file). Latest versions of these 
scripts can also be downloaded from https://github.com/abakan/ProDy/tree/master/scripts/.

:ref:`scripts` section shows the help text that comes with these scripts.

Using these scripts is straight forward (on linux). Let's assume you are
in the same folder as ANM script (also using a linux machine). 
Running the following command will perform calculations for p38 MAP kinase 
structure, and will write eigenvalues/vectors in plain text and NMD formats::

  ./anm 1p38
  
In this example, default parameters are used (``cutoff=15.`` and ``gamma=1.``).
Also, all alpha carbons of the protein structure 1p38 are used.

In the following case, cutoff distance is changed to 14 angstroms, 
alpha carbons of residues with numbers smaller than 340 are used, 
and output files are prefixed with "p38anm"::

  ./anm -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

Output file :file:`p38_anm.nmd` file can be visualized using |nmwiz|. 

Using ProDy interactively
===============================================================================

In the rest of the tutorial, it is assumed that you are typing commands in an 
interaction Python session. If you are using Python IDE, try using 
`IPython <http://ipython.scipy.org/>`_ that comes fancy coloring and handy features.


Import all functions and classes from ProDy as follows:

>>> from prody import *

Parse a PDB file
-------------------------------------------------------------------------------

Let's start with reading the contents of a PDB file. ProDy offers a fast PDB 
file parser (:func:`prody.proteins.parsePDB`). For this function to work, 
PDB file does not need to exist on your machine, so long as you have internet 
connection. ProDy will fetch the PDB file if a valid identifier is provided.

>>> prot = parsePDB('1p38')
@> 1p38 downloaded (./1p38.pdb.gz)
@> 2962 atoms and 1 coordinate sets were parsed in 0.08s.


1p38 is an unliganded structure of p38 MAP kinase. :file:`1p38.pdb.gz` has been 
downloaded, and coordinates were parsed into an :class:`prody.proteins.AtomGroup`
instance.

.. note::
   ProDy prints some log information to the console. The level of verbosity can
   be adjusted using :func:`ProDySetVerbosity` function.

To get some information on the AtomGroup instance, you can type variable name
and press *Enter*.

>>> prot
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>


There are more ways to access ProteinDataBank (http://www.pdb.org/) content.
See :ref:`parsepdb`, :ref:`fetchpdb`, and :ref:`blastpdb` examples.

List of functions for accessing protein data can also be found in :ref:`prodb`

Select subset of atoms
-------------------------------------------------------------------------------

ProDy :class:`prody.proteins.AtomGroup` instances offer powerful atom
selection capabilities that are comparable to that of 
`VMD <http://www.ks.uiuc.edu/Research/vmd/>`_. Full list of selection keywords 
are given in section :ref:`selections`. Here, only a few examples are shown:

**Select protein atoms**:

>>> only_protein = prot.select('protein')
>>> only_protein
<Selection: "protein" from 1p38 (2833 atoms; 1 coordinate sets, active set index: 0)>

This shows that 2833 of 2962 atoms are protein atoms.

**Select atoms by name**:

>>> side_chain_atoms = prot.select('protein and not name N CA C O')
>>> print side_chain_atoms
<Selection: "protein and not name N CA C O" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>
>>> len(side_chain_atoms)
1429

Same selection could also be made using ``sidechain`` or ``not backbone`` keywords:

>>> prot.select('sidechain')
<Selection: "sidechain" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

>>> prot.select('not backbone')
<Selection: "not backbone" from 1p38 (1558 atoms; 1 coordinate sets, active set index: 0)>

Oops, ``not backbone`` did not select the same number of atoms. Let's try to
see why:

>>> print set(prot.select('not backbone').getResidueNames())
set(['CYS', 'ILE', 'VAL', 'GLN', 'LYS', 'HOH', 'PRO', 'THR', 'PHE', 'ASN', 'HIS', 'MET', 'ASP', 'LEU', 'ARG', 'TRP', 'ALA', 'GLU', 'TYR', 'SER'])

Note that we used built-in Python type :class:`set`.

As you can see atoms of **HOH** residues are also included in the selection.

Let's try:

>>> prot.select('not backbone and not water')
<Selection: "not backbone and not water" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>

This has now worked as "sidecain" did.

**Select amino acids by type/name**:

>>> charged = prot.select('acidic or basic')
>>> print charged
<Selection: "acidic or basic" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>
>>> len(charged)
906
>>> set(charged.getResidueNames())
set(['ARG', 'ASP', 'GLU', 'HIS', 'LYS'])

Same selection could also be made using ``charged`` keyword:

>>> prot.select('charged')
<Selection: "charged" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

Or, we could use residue names expilicitly:

>>> prot.select('resname ARG LYS HIS ASP GLU')
<Selection: "resname ARG LYS HIS ASP GLU" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

.. seealso::
   For more information, tips and tricks see :ref:`selections` and :ref:`selops`.

Hierarchical view of atoms in a PDB 
-------------------------------------------------------------------------------

:class:`prody.proteins.AtomGroup` instances has a flat view of atoms in PDB
files, but it is possible to get a hierarchical view (:class:`prody.proteins.HierView`) 
of them:

>>> hv = prot.getHierView()

Now, one can iterate over chains and residues:

>>> for chain in hv:
>>>     print chain
>>> ...
Chain A from 1p38

>>> for res in hv.iterResidues():
>>>     print res
>>> ...
GLU 4 from Chain A from 1p38
ARG 5 from Chain A from 1p38
PRO 6 from Chain A from 1p38
THR 7 from Chain A from 1p38
...

Write a PDB file
-------------------------------------------------------------------------------

PDB files can be written using :func:`prody.proteins.writePDB` function.
This function accepts objects containing or referring to atomic data.

Writing a selection:

>>> calphas = prot.select('calpha')
>>> writePDB('1p38_calphas.pdb', calphas)
'1p38_calphas.pdb'

Write a chain:

>>> chain_A = hv.getChain('A')
>>> writePDB('1p38_chain_A.pdb', chain_A)
'1p38_chain_A.pdb'

As you might have noticed, this function returns the file name after it is
successfully written.


Perform ANM calculations
-------------------------------------------------------------------------------

Let's perform an ANM analysis for chain A alpha carbon atoms. ANM instances
are instantiated using a name:

>>> anm = ANM('p38 ANM anlaysis')

Hessian matrix can be built for any set of atoms. In this case, we will 
use selection that contains alpha carbon atoms. 

>>> anm.buildHessian(calphas)
@> Hessian was built in 1.62s.

Normal modes are calculated by calling :meth:`prody.dynamics.ANM.calcModes`. 
This will calculate 20 modes by default.

>>> anm.calcModes()
@> 20 modes were calculated in 1.52s.

This omits modes with zero eigenvalues. For moreinformation on ANM methods 
see :func:`prody.dynamics.ANM`.

Individual modes can be accessed by indexing ANM instance:

>>> slowest_mode = anm[0]
>>> print slowest_mode
Mode 1 from p38 ANM analysis

Note that indices in Python start from 0. 0th mode is the 1st non-zero mode,
in this case.

The following function (:func:`prody.dynamics.writeNMD`) writes ANM results 
in NMD format. NMD files can be viewed using |vmd| plugin |nmwiz|. 

>>> writeNMD('p38anm.nmd', anm[:6], calphas) 
