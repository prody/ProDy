.. _tutorial:

.. currentmodule:: prody

*******************************************************************************
Tutorial
*******************************************************************************

ProDy can be used:

* interactively in a Python shell,
* as a command line program via :ref:`scripts`,
* or as a toolkit for developing new software.

For those who are new to Python, `A Primer on Python for Life Science Researchers 
<http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.0030199>`_
or the more comprehensive `The Python Tutorial 
<http://docs.python.org/tutorial/>`_ is recommended.

ProDy documentation is organized in three main sections:

* :ref:`tutorial`, which starts shortly, contains brief usage examples of 
  some ProDy features and is the best place to start with learning ProDy.
* :ref:`examples` section contains comprehensive usage examples, which may be 
  applied to other cases after small modifications.
* :ref:`reference` section describes all ProDy classes and functions. The 
  tutorial and example pages link to relevant parts of this section. 
  
  In interactive Python sessions, reference documentation can be accessed 
  using built-in Python :func:`help` function: 

  >>> help(parsePDB) # doctest: +SKIP

  This function can be used to print the description of functions, classes, 
  and class methods on the screen. You might need to type ``q`` to quit from 
  help. Also, if you are using interactive Python shell (IPython), you can 
  get help as follows::
  
    prot ?


Finally, this page and others in the :ref:`examples` and :ref:`reference` 
sections contain ProDy code snippets. These snippets can be retrieved using 
the "Show Code Snippets" button on the right hand side panel. ProDy code 
will be displayed in a popup window. The code in the new page can be directly 
copied into a file. Click on the text, press :kbd:`Ctrl+A` and then 
:kbd:`Ctrl+C` to have the text in your clipboard. 

.. image:: /_static/codesnippets.png
   :align: center
   :alt: Getting ProDy code snippets.

First part of the tutorial shows how to use one of the :ref:`scripts` that is 
distributed with ProDy. Following parts assume user is in an interactive
Python shell. 

.. contents:: The ProDy Tutorial
   :local:
   :backlinks: none


Command line usage
===============================================================================

ProDy scripts come with source distribution (:file:`tar.gz` file in 
:ref:`getprody`). Latest versions of these scripts can also be obtained from 
https://github.com/abakan/ProDy/tree/master/scripts/.

Using these scripts is straight forward (on Linux). Let's assume you are
in the same folder as the :file:`anm.py` script. Running the following 
command will perform ANM calculations for p38 MAP kinase structure, and will 
write eigenvalues/vectors in plain text and :term:`NMD` formats::

  ./python anm.py 1p38
  
In this example, default parameters (``cutoff=15.`` and ``gamma=1.``)
and all alpha carbons of the protein structure 1p38 are used.

In the following, cutoff distance is changed to 14 angstroms, 
alpha carbons of residues with numbers smaller than 340 are used, 
and output files are prefixed with "p38_anm"::

  ./python anm.py -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

Output file :file:`p38_anm.nmd` file can be visualized using |nmwiz|. 

:ref:`scripts` section shows more detailed information.


Interactive usage
===============================================================================

One of our aims is making it ProDy suitable for interactive usage
by designing flexible functions and classes and giving them easy to remember 
names in a consistent manner. 

For best interactive usage experience, 
we strongly recommend `IPython <http://ipython.scipy.org/>`_ or a similar 
interactive shell over the standard Python shell. It offers a nice coloring 
scheme, handy features, such as tab completion, and convenient
integration of Numpy and Matplotlib (see http://matplotlib.sourceforge.net).

In the rest of this tutorial, we assume that commands are typed in a 
Python shell. We start with importing all functions and classes from 
ProDy into the current namespace as follows:

>>> from prody import *

Parse a PDB file
===============================================================================

Let's start with reading the contents of a PDB file. ProDy offers a fast PDB 
file parser function :func:`~proteins.parsePDB` (:ref:`pdbparser-performance`). 
It is sufficient to pass a PDB identifier to get all atomic data from a file.
The parser will seek the PDB file in the current working directory, and if
it is not found it will download it from the PDB FTP server automatically.

>>> prot = parsePDB('1p38')

1p38 is an unbound structure of p38 MAP kinase. :file:`1p38.pdb.gz` has been 
downloaded, and coordinates were parsed into an :class:`~prody.atomic.AtomGroup`
instance.

To get some information on the :class:`~prody.atomic.AtomGroup` instance, you 
can type variable name and hit :kbd:`enter` key.

>>> prot
<AtomGroup: 1p38 (2962 atoms; 1 coordinate sets, active set index: 0)>

All atomic data can be accessed using ``get`` methods:

>>> print prot.getResidueNames()
['GLU' 'GLU' 'GLU' ..., 'HOH' 'HOH' 'HOH']
>>> print prot.getCoordinates() # doctest: +SKIP
[[ 28.492   3.212  23.465]
 [ 27.552   4.354  23.629]
 [ 26.545   4.432  22.489]
 ..., 
 [ 18.872   8.33   36.716]
 [-22.062  21.632  42.029]
 [  1.323  30.027  65.103]]

More examples
-------------------------------------------------------------------------------

:func:`~proteins.parsePDB` function is very flexible and can be extremely
efficient depending on what you want to extract from a PDB file. You can
parse specific chains, alternate locations, subsets of atoms, or models
from a file. A detailed usage example can be found in :ref:`parsepdb`.

There are more ways to access Protein Data Bank (http://www.pdb.org/) 
content via ProDy. See the :ref:`fetchpdb` for downloading and :ref:`blastpdb` 
for blast searching PDB files.

Descriptions of all functions for accessing and handling protein data can be 
found in :mod:`~prody.proteins` module reference documentation.


Select atoms
===============================================================================

ProDy :class:`~atomic.AtomGroup` objects have a plain view of
atoms, but offer a powerful atom selector. You can get well defined subsets of 
atoms by passing simple keywords as arguments or make sophisticated selections
using composite arguments. Selection keywords and grammar is very much similar
to that of `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_. Let's see some 
examples here:

Select protein atoms:

>>> protein = prot.select('protein')
>>> protein
<Selection: "protein" from 1p38 (2833 atoms; 1 coordinate sets, active set index: 0)>

Using :term:`protein` keyword we selected 2833 atoms out of 2962 atoms. 
:meth:`~prody.atomic.Atomic.select` method returned a :class:`~atomic.Selection` 
instance.

Select atoms by :term:`name`:

>>> backbone = prot.select('protein and not name N CA C O')
>>> backbone
<Selection: "protein and not name N CA C O" from 1p38 (1429 atoms; 1 coordinate sets, active set index: 0)>
>>> len(backbone)
1429

We could also use :term:`backbone` to make the same selection. 

Select amino acids by type/name (:term:`resname`):

>>> prot.select('resname ARG LYS HIS ASP GLU')
<Selection: "resname ARG LYS HIS ASP GLU" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>

Or, we can use predefined keywords :term:`acidic` and :term:`basic`.

>>> charged = prot.select('acidic or basic')
>>> charged
<Selection: "acidic or basic" from 1p38 (906 atoms; 1 coordinate sets, active set index: 0)>
>>> set(charged.getResidueNames())
set(['HIS', 'ASP', 'LYS', 'GLU', 'ARG'])

Here is a more sophisticated composite selection:

>>> print protein.getCoordinates().mean(0) # doctest: +SKIP
[  1.005  17.533  40.052]
>>> prot.select('protein and name CA CB and same residue as ((x-1)**2 + (y-17.5)**2 + (z-40.0)**2)**0.5 < 10')
<Selection: "protein and nam...)**2)**0.5 < 10" from 1p38 (66 atoms; 1 coordinate sets, active set index: 0)>

In this case, we first calculated the geometric center of the protein atoms.
Then, we selected carbon alpha and beta atoms of residues that have at least 
one atom within 10 angstrom away from the geometric center.

More examples
-------------------------------------------------------------------------------

You can do much more with this flexible and fast atom selection engine,
without the need for writing nested loops with comparisons or changing the 
source code. See the following for more information:

  * :ref:`selection-examples` for a detailed usage example
  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of :class:`~atomic.Selection`
    objects
  * :ref:`contacts` for selecting interacting atoms

Hierarchical views
===============================================================================

It is always useful to have a hierarchical view of protein structure in hand.
You do not get this by default when you parse a PDB file, but it is not 
missing from ProDy. Calling :meth:`~atomic.AtomGroup.getHierView` method of an 
:class:`~atomic.AtomGroup` instance will return you a 
hierarchical view of the atoms (:class:`~atomic.HierView`):

>>> hv = prot.getHierView()

Now, one can iterate over chains and residues:

>>> for chain in hv: print chain
Chain A
>>> for res in hv.iterResidues(): print res
GLU 4
ARG 5
PRO 6
...
HOH 771
HOH 773
HOH 776

A more detailed usage example can be found in :ref:`hierview`.

Write a PDB file
===============================================================================

PDB files can be written using :func:`~prody.proteins.writePDB` function.
This function accepts objects containing or referring to atomic data.

Write selected atoms:

>>> writePDB('1p38_calphas.pdb', prot.select('calpha'))
'1p38_calphas.pdb'

Write a chain:

>>> chain_A = hv.getChain('A')
>>> writePDB('1p38_chain_A.pdb', chain_A)
'1p38_chain_A.pdb'

As you might have noticed, this function returns the file name after it is
successfully written.

PCA calculations
===============================================================================

We show how to perform principal component analysis (:class:`~dynamics.PCA`) 
of a set of NMR models for ubiquitin (PDB ID: 2k39).

Parse and align the coordinate data:

>>> ubi = parsePDB('2k39', subset='calpha')
>>> ubi = ubi.copy('resnum < 71')
>>> alignCoordsets(ubi)

Perform the PCA:

>>> pca = PCA('Ubiquitin')
>>> pca.buildCovariance(ubi)
>>> pca.calcModes()

Print the fraction of variance for top raking 4 PCs:

>>> for mode in pca[:4]:
...     print mode.getFractOfVariance() # doctest: +SKIP
0.299016803492
0.0959780950608
0.0647918823066
0.058247703612

>>> writeNMD('ubi_pca.nmd', pca[:3], ubi)
'ubi_pca.nmd'

This was a short example for an easy case. :ref:`pca` section
contains more comprehensive examples for heterogeneous datasets. 

ANM calculations
===============================================================================

Anisotropic network model (:class:`~prody.dynamics.ANM`) analysis can be 
performed in two ways:

The longer and more controlled way:

>>> calphas = ubi.select('calpha') # select alpha carbons
>>> anm = ANM('ubi') # instantiate ANM object
>>> anm.buildHessian(calphas) # build Hessian matrix for selected atoms 
>>> anm.calcModes() # calculate normal modes
  
This gives more control to the user. For example, instead of building the
Hessian matrix using ProDy, user can provide a pre-calculated matrix 
using :meth:`~dynamics.ANM.setHessian` and diagonalize it using ProDy for further
analysis. 

Alternative shorter way, which may be suitable for interactive sessions:

>>> anm = calcANM(ubi, selstr='calpha')

Individual modes can be accessed by indexing ANM instance:

>>> slowest_mode = anm[0]
>>> print slowest_mode
Mode 1 from ANM Copy of 2k39 selection "resnum < 71"
>>> slowest_mode.getEigenvalue() # doctest: +SKIP
1.7142408905432185

Note that indices in Python start from zero (0). 
0th mode is the 1st non-zero mode, in this case.

The following function (:func:`~dynamics.writeNMD`) writes ANM results 
in NMD format. NMD files can be viewed using |vmd| plugin NMWiz. 

>>> writeNMD('p38_anm.nmd', anm[:6], calphas) 
'p38_anm.nmd'

For more information on elastic network model calculations see
:ref:`enm` section.

Comparative analysis
===============================================================================

ProDy comes with many built-in functions to facilitate a comparative analysis
of experimental and theoretical data. For example, using 
:func:`printOverlapTable` function you can see the agreement between 
experimental (PCA) modes and theoretical (ANM) modes calculated above:

>>> printOverlapTable(pca[:4], anm[:4])
Overlap Table
                  ANM Copy of 2k39 selection "resnum < 71"
                     #1     #2     #3     #4
PCA Ubiquitin #1   +0.00  +0.14  +0.07  -0.20
PCA Ubiquitin #2   +0.22  -0.35  +0.17  +0.51
PCA Ubiquitin #3   -0.13  -0.67  -0.15  -0.24
PCA Ubiquitin #4   -0.29  +0.13  -0.11  +0.02
<BLANKLINE>

For a more comprehensive example usage of analysis functions see 
:ref:`pca-xray`.

ProDy verbosity
===============================================================================

Finally, you might have noted that ProDy prints some information to the console
after parsing a file or doing some calculations. For example, PDB parser will 
print what was parsed and how long it took to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

The level of verbosity can be adjusted using :func:`changeVerbosity` function.
If you do not want any logs printed on the screen, you can enter 
``changeVerbosity(None)``.

Questions or suggestions
===============================================================================

If you have any questions or suggestions please contact us in one of the
following ways:

|questions|

|suggestions|


