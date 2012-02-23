.. _tutorial:

.. currentmodule:: prody

*******************************************************************************
Tutorial
*******************************************************************************

Foreword
===============================================================================

ProDy is designed for structure-based analysis of protein dynamics. It is 
suitable for a variety of uses from everyday parsing and analysis of protein 
structures to developing or prototyping software for structure-based analysis.  

Dynamics from experiments
-------------------------------------------------------------------------------

ProDy allows for quantitative analysis of *large* heterogeneous experimental 
structural datasets.  It takes, for example, under 6 seconds to parse, align, 
and analyze 75 p38 MAP kinase structures on a personal laptop computer 
(see :ref:`pca-xray-calculations`).  Such experimental datasets may contain 
sequence homologs, orthologs, mutant, or ligand bound forms of a protein.  
Dominant patterns in structural variability are extracted by principal 
component analysis (PCA) of the ensemble, and can be compared with 
theoretically predicted conformational dynamics using ProDy.

Dynamics from theory
-------------------------------------------------------------------------------

On the theoretical side, protein dynamics is predicted using normal mode 
analysis (NMA) based on elastic network models (ENMs) or is extracted from 
molecular dynamics (MD) trajectories using essential dynamics analysis (EDA).  
Numerous helper functions enable comparative analysis of thus obtained 
experimental and theoretical data, and visualize the principal changes 
in conformations that are accessible in different functional states.

Input for ProDy
-------------------------------------------------------------------------------

The input is for calculations is the set of atomic coordinates of the query 
protein in PDB file format, or simply the PDB id or single letter amino acid 
sequence of the protein.  Fast and flexible ProDy parsers are used to retrieve 
data from the `wwPDB <http://www.wwpdb.org/>`_ FTP servers and extract the
coordinate data and other relevant information. 

Additionally, PCA/NMA/EDA results calculated by other software can be parsed
and analyzed using rich analysis and plotting functions of ProDy.  User just 
needs to provide normal mode data in plain text files.  See the tutorial for 
an example.

Visualization
-------------------------------------------------------------------------------

Analysis and plotting capabilities of ProDy are complemented VMD plugin 
:ref:`NMWiz`.  NMWiz can be used to visualize and animate normal mode data 
calculated using ProDy or other software.  See :ref:`NMWiz` documentation for 
details. 

How to Use
===============================================================================

How to use ProDy
-------------------------------------------------------------------------------

ProDy can be used:

* interactively in a Python shell,
* as a command line program via :ref:`commands`,
* from within VMD via :ref:`nmwiz`,
* or as a toolkit for developing new software.

Python for beginners
-------------------------------------------------------------------------------

Users who are new to programming or Python are referred to the following 
resources for an introduction to programming in Python:

* `The Python Tutorial <http://docs.python.org/tutorial/>`_
* `Python Scientific Lecture Notes <http://scipy-lectures.github.com/>`_
* `A Primer on Python for Life Science Researchers
  <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.0030199>`_


Using documentation
-------------------------------------------------------------------------------

ProDy documentation is organized in three main sections:

* Tutorial, contains brief usage examples of some ProDy features and is the 
  best place to start to learn ProDy.
* :ref:`examples` section contains comprehensive usage examples, which may be 
  applied to other cases after small modifications.
* :ref:`reference` section describes all ProDy classes and functions, with some
  usage examples. Features are divided into following main modules:

  * :mod:`~.atomic` - efficient handling of atomic data
  * :mod:`~.dynamics` - analysis and modeling of protein dynamics
  * :mod:`~.ensemble` - analysis of ensembles and trajectories
  * :mod:`~.measure` - analysis of geometrical properties 
  * :mod:`~.proteins` - structure analysis, IO, and data retrieval
  * :mod:`~.select` - atom selections
  * :mod:`~.trajectory` - trajectory IO
  
In interactive sessions, the reference documentation can be accessed 
using the built-in Python function :func:`help`:: 

  help(atomic)
  help(select)
  help(parsePDB)

This function prints the description of functions, classes, and class methods 
to the screen. Note that you might need to type ``q`` to exit from 
help. If you are using the interactive Python shell (IPython), you can also 
receive help by typing::
  
  select ?

Copying examples
-------------------------------------------------------------------------------

Documentation contains ProDy code snippets.  These snippets can be reformatted 
using the :guilabel:`Show Code Snippets` button on the right hand side panel. 
The code will be displayed in the same window, and can be copied directly into 
a file. Click on the text, press :kbd:`Ctrl+A` and then :kbd:`Ctrl+C` to have 
the text in your clipboard. To return to the documentation click the 
:guilabel:`Show documentation` button at the top.

.. image:: /_static/codesnippets.png
   :align: center
   :alt: Getting ProDy code snippets.

.. contents:: The ProDy Tutorial
   :local:
   :backlinks: none


Interactive usage
-------------------------------------------------------------------------------

One of our aims is making ProDy suitable for interactive usage by designing 
flexible functions and classes and giving them easy to remember names in a 
consistent manner. 

For best interactive usage experience, we strongly recommend that you use 
`IPython <http://ipython.scipy.org/>`_ or a similar interactive shell instead 
of the standard Python shell.  The IPython shell, for example, provides 
user-friendly features, such as dynamic introspection and help, and also 
optionally convenient integration of Numpy and `Matplotlib 
<http://matplotlib.sourceforge.net>`_.

In the rest of this tutorial, it is assumed that the user is typing commands 
in a Python shell. To begin the Tutorial, import all the functions and classes 
from ProDy into the current namespace as follows:

>>> from prody import *

There are other ways to import ProDy contents. You may use 
``import prody as pd`` and prefix all functions calls with ``pd.``, 
if you prefer not to overcrowd your namespace.
Alternatively, if you want to use functions in a specific module, 
:mod:`~.proteins` let's say, you can use ``from prody.proteins import *``. 
You should, however, avoid using ``from prody.proteins.pdbfile import *``, 
because location of methods in submodules may change without notice.

.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np
   prot = parsePDB('1p38')

   plt.close('all')

Protein structure
===============================================================================

Protein structure files in :file:`.pdb` format are the standard input for 
ProDy.  PDB files are parsed using :func:`~.parsePDB` (see 
:ref:`pdbparser-performance` and :ref:`pdbparser-performance-2` for 
benchmarks).  It is sufficient to pass a PDB identifier to read the file, and 
the parser will download it automatically if needed.

>>> prot = parsePDB('1p38')

In the above line, :file:`1p38.pdb` is downloaded and coordinates and atomic 
data are parsed from the file. 
  
Managing resources
-------------------------------------------------------------------------------

You can tell ProDy where to get PDB files from or where to store downloaded 
files:

  * one of the `wwPDB <http://www.wwpdb.org/>`_ FTP servers in US, Europe or
    Japan can be picked for downloads using :func:`~.setWWPDBFTPServer`
  * a local PDB mirror can be be set for faster access to files using
    :func:`~.setPDBMirrorPath` 
  * a local folder can be set for storing downloaded files for later access
    using :func:`~.setPDBLocalFolder` 

Note that when these functions are used, ProDy will save your settings as 
:file:`.prodyrc` in your home folder.


Quick visualization
-------------------------------------------------------------------------------

:file:`1p38.pdb` contains an unbound structure of the p38 MAP kinase.
If you have `Matplotlib`_ installed, you can take 
a quick look at what you parsed using :func:`~.showProtein` function:  


>>> import matplotlib.pyplot as plt
>>> plt.figure(figsize=(5,4)) # doctest: +SKIP
>>> showProtein(prot) # doctest: +SKIP
>>> plt.legend(prop={'size': 10}) # doctest: +SKIP

.. plot::
   :context:
   
   import matplotlib.pyplot as plt
   plt.figure(figsize=(5,4))
   showProtein(prot)
   plt.legend(prop={'size': 10})

   
.. plot::
   :nofigs: 
   :context: 
   
   plt.close('all')

:func:`~.parsePDB` returns data in an :class:`~.AtomGroup` 
instance. 


Atomic data
-------------------------------------------------------------------------------

To get information on an :class:`~.AtomGroup` instance, 
type in the variable name and hit :kbd:`enter` key:

>>> prot
<AtomGroup: 1p38 (2962 atoms)>

The above shows that atom group object contains 2962 atoms. 
All atomic data from this object can be retrieved using ``get`` methods. 
For example:

>>> print( prot.getResnames() )
['GLU' 'GLU' 'GLU' ..., 'HOH' 'HOH' 'HOH']
>>> print( prot.getCoords() ) # doctest: +ELLIPSIS
[[ 28.492   3.212  23.465]
 [ 27.552   4.354  23.629]
 ...
 [-22.062  21.632  42.029]
 [  1.323  30.027  65.103]]
 
The list of methods for getting and setting atomic data is provided in
:class:`~.AtomGroup` reference documentation. 

**Indexing**:

An individual :class:`~.Atom` can be accessed by indexing atom group 
instances:

>>> atom = prot[0]
>>> atom
<Atom: N from 1p38 (index 0)>

Not that all ``get/set`` functions defined for :class:`~.AtomGroup` 
instances are also defined for :class:`~.Atom` instances, using singular
form of the function name.  

>>> atom.getResname()
'GLU'

**Slicing**:

It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

>>> prot[::2]
<Selection: "index 0:2962:2" from 1p38 (1481 atoms)>

Hierarchical view
-------------------------------------------------------------------------------

You can also access specific chains or residues in an atom group.  Indexing
by a single letter identifier will return a :class:`~.Chain` instance: 

>>> prot['A']
<Chain: A from 1p38 (480 residues, 2962 atoms)>

Indexing atom group with a chain identifier and a residue number will return
:class:`~.Residue` instance:

>>> prot['A', 100]
<Residue: ASN 100 from Chain A from 1p38 (8 atoms)>

See :ref:`atomic` for details of indexing atom groups and :ref:`hierview`
for more on hierarchical views.


Writing PDB files
-------------------------------------------------------------------------------

PDB files can be written using the :func:`~.writePDB` function.
The function accepts objects containing or referring to atomic data.

Writing selected atoms:

>>> writePDB('1p38_calphas.pdb', prot.select('calpha'))
'1p38_calphas.pdb'

Writing a chain:

>>> chain_A = prot['A']
>>> writePDB('1p38_chain_A.pdb', chain_A)
'1p38_chain_A.pdb'

As you may have noticed, this function returns the file name after it is
successfully written.  This is a general behavior for ProDy output functions.
For more PDB writing examples see :ref:`writepdb`.

More examples
-------------------------------------------------------------------------------

:func:`~.parsePDB` function is very flexible and can be extremely
efficient depending on what you want to extract from a PDB file.  It can be 
used to parse specific chains, models, alternate locations, or well-defined 
subsets of atoms from a file.  A detailed usage example can be found in 
:ref:`parsepdb`.  

ProDy can parse other file types, including :file:`.psf` and :file:`.pqr` files.
All of the functions for accessing and handling protein structural data are 
described in :mod:`~.proteins` module reference documentation.
Also, :ref:`fetchpdb` and :ref:`blastpdb` examples show other ways to 
access the Protein Data Bank (|pdb|) content.

For more details on atomic objects see :ref:`atomic`.  
:class:`~.AtomGroup` instances can be build from scratch or 
parsers for other file types (e.g. mol2) can be developed. The example in 
:ref:`atomgroup` can be helpful to this aim.


Atom selections
===============================================================================

:class:`~.AtomGroup` instances have a plain view of atoms for efficiency, 
but they are coupled with a powerful atom selection engine.  You can get well 
defined atom subsets by passing simple keywords or make rather sophisticated 
selections using composite statements.  Selection keywords and grammar is very 
much similar to those found in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.  
Some examples are shown here:

Keyword selections
-------------------------------------------------------------------------------

>>> protein = prot.select('protein')
>>> protein
<Selection: "protein" from 1p38 (2833 atoms)>

Using the "protein" keyword we selected 2833 atoms out of 2962 atoms. 
:meth:`~.Atomic.select` method returned a :class:`~.Selection` 
instance.  Note that all ``get`` and ``set`` methods defined for
the :class:`~.AtomGroup` class are also defined for 
:class:`~.Selection` class. For example:

>>> print( protein.getResnames() )
['GLU' 'GLU' 'GLU' ..., 'ASP' 'ASP' 'ASP']

Select by name/type
-------------------------------------------------------------------------------

We select backbone atoms by passing atom names following "name" keyword:

>>> backbone = prot.select('protein and name N CA C O')
>>> backbone
<Selection: "protein and name N CA C O" from 1p38 (1404 atoms)>
>>> len(backbone)
1404

We can also use "backbone" to make the same selection. 

We select acidic and basic residues by using residue names with 
"resname" keyword:

>>> prot.select('resname ARG LYS HIS ASP GLU')
<Selection: "resname ARG LYS HIS ASP GLU" from 1p38 (906 atoms)>

Alternatively, we can use predefined keywords "acidic" and "basic".

>>> charged = prot.select('acidic or basic')
>>> charged
<Selection: "acidic or basic" from 1p38 (906 atoms)>
>>> set(charged.getResnames())
set(['HIS', 'ASP', 'LYS', 'GLU', 'ARG'])

Composite selections
-------------------------------------------------------------------------------

Let's try a more sophisticated selection.  We first calculate the geometric 
center of the protein atoms. Then, we select the Cα and Cβ atoms of residues 
that have at least one atom within 10 Å away from the geometric center.

>>> print( protein.getCoords().mean(0).round(3) ) # doctest: +ELLIPSIS
[  1.005  17.533  40.052]
>>> prot.select('protein and name CA CB and same residue as ((x-1)**2 + (y-17.5)**2 + (z-40.0)**2)**0.5 < 10')
<Selection: "protein and nam...)**2)**0.5 < 10" from 1p38 (66 atoms)>

Selection operations
-------------------------------------------------------------------------------

:class:`~.Selection` instances can be

>>> ca = prot.select('name CA') 
>>> cb = prot.select('name CB')
>>> ca | cb
<Selection: "(name CA) or (name CB)" from 1p38 (687 atoms)>
>>> ca & cb

Selections simplified
-------------------------------------------------------------------------------

In interactive sessions, typing in ``.select('backbone')`` or even 
``.select('bb')`` may be time consuming.  An alternative to this is using
dot operator:

>>> prot.protein
<Selection: "protein" from 1p38 (2833 atoms)>

You can use dot operator multiple times:

>>> prot.protein.backbone
<Selection: "(backbone) and (protein)" from 1p38 (1404 atoms)>

This may go on and on:

>>> prot.protein.backbone.resname_ALA.calpha
<Selection: "(calpha) and ((...and (protein)))" from 1p38 (26 atoms)>


More examples
-------------------------------------------------------------------------------

There is much more to what you can do with this flexible and fast atom 
selection engine, without the need for writing nested loops with comparisons 
or changing the source code.  See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of 
    :class:`~.Selection` objects
  * :ref:`contacts` for selecting interacting atoms


ProDy verbosity
-------------------------------------------------------------------------------

Finally, you might have noted that ProDy prints some information to the console
after parsing a file or doing some calculations. For example, PDB parser will 
print what was parsed and how long it took to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

This behavior is useful in interactive sessions, but may be problematic for
automated tasks as the messages are printed to stderr. The level of verbosity 
can be adjusted using :func:`setVerbosity` function, and 
``setVerbosity(None)`` will stop all information messages.


Dynamics analysis
===============================================================================

PCA calculations
-------------------------------------------------------------------------------

We show how to perform principal component analysis (:class:`~.PCA`) 
of a set of NMR models for ubiquitin (PDB ID: 2k39).

Parse and align the coordinate data:

>>> ubi = parsePDB('2k39', subset='calpha')
>>> ubi_selection = ubi.select('resnum < 71')
>>> ubi_ensemble = Ensemble(ubi_selection)
>>> ubi_ensemble.iterpose()

Perform the PCA:

>>> pca = PCA('Ubiquitin')
>>> pca.buildCovariance(ubi_ensemble)
>>> pca.calcModes()

Print the fraction of variance for top raking 4 PCs:

>>> for mode in pca[:4]:
...     print( mode.getFractOfVariance().round(2) ) # doctest: +ELLIPSIS
0.13
0.09
0.08
0.07

PCA data can be saved on disck using :func:`~.saveModel`
function:

>>> saveModel(pca)
'Ubiquitin.pca.npz'

This functions writes data in binary format, so is an efficient way of 
storing data permanently.  In a later session, this data can be loaded using 
:func:`~.loadModel` function.

NMD files 
-------------------------------------------------------------------------------

The :func:`~.writeNMD` function writes PCA results 
in NMD format. NMD files can be viewed using the :ref:`nmwiz` VMD plugin.

>>> writeNMD('ubi_pca.nmd', pca[:3], ubi_selection)
'ubi_pca.nmd'

Write files 
-------------------------------------------------------------------------------

Additionally, results can be written in plain text files for analysis with
other programs using the :func:`~.writeArray` function:

>>> writeArray('ubi_pca_modes.txt', pca.getArray(), format='%8.3f')
'ubi_pca_modes.txt'


ANM calculations
-------------------------------------------------------------------------------

Anisotropic network model (:class:`~.ANM`) analysis can be 
performed in two ways:

The shorter way, which may be suitable for interactive sessions:

>>> anm, atoms = calcANM(ubi_selection, selstr='calpha')

The longer and more controlled way:

>>> anm = ANM('ubi') # instantiate ANM object
>>> anm.buildHessian(ubi_selection) # build Hessian matrix for selected atoms 
>>> anm.calcModes() # calculate normal modes
>>> saveModel(anm)
'ubi.anm.npz'


:ref:`anm` provides a more detailed discussion of ANM calculations. 
The above longer way gives more control to the user. For example, instead of 
building the Hessian matrix using uniform force constant and cutoff distance, 
customized force constant functions (see :ref:`gamma`) or a pre-calculated matrix 
(see :meth:`~.ANM.setHessian`) may be used. 

Individual :class:`~.Mode` instances can be accessed by 
indexing the :class:`~.ANM` instance:

>>> slowest_mode = anm[0]
>>> print( slowest_mode )
Mode 1 from ANM ubi
>>> print( slowest_mode.getEigenvalue().round(3) )
1.714

Note that indices in Python start from zero (0). 
0th mode is the 1st non-zero mode in this case.

The :func:`~.writeNMD` function writes ANM results 
in NMD format. NMD files can be viewed using the :ref:`nmwiz` VMD plugin. 

>>> writeNMD('p38_anm.nmd', anm[:6], ubi_selection) 
'p38_anm.nmd'

For more information on elastic network model calculations see
:ref:`enm` section.

Comparative analysis
-------------------------------------------------------------------------------

ProDy comes with many built-in functions to facilitate a comparative analysis
of experimental and theoretical data. For example, using 
:func:`~.printOverlapTable` function you can see the agreement between 
experimental (PCA) modes and theoretical (ANM) modes calculated above:

>>> printOverlapTable(pca[:4], anm[:4])
Overlap Table
                            ANM ubi
                     #1     #2     #3     #4
PCA Ubiquitin #1   -0.21  +0.30  -0.17  -0.47
PCA Ubiquitin #2   +0.01  +0.72  +0.08  +0.05
PCA Ubiquitin #3   +0.31  +0.11  +0.18  +0.19
PCA Ubiquitin #4   +0.11  -0.02  -0.17  -0.39
<BLANKLINE>

Output above shows that PCA mode 2 and ANM mode 2 for ubiquitin show the 
highest overlap (cosine-correlation). 

.. plot::
   :context:
   :nofigs:
   
   pca = loadModel('Ubiquitin.pca.npz')
   anm = loadModel('ubi.anm.npz')

We can also make a plot of this table using :func:`~.showOverlapTable`
function:

.. plot::
   :include-source:
   :context:
   
   plt.figure( figsize=(5,4) )
   showOverlapTable(pca[:4], anm[:4])
   
.. plot::
   :nofigs:
   :context:
   
   plt.close('all')

This was a short example for a simple case. :ref:`pca` section contains more 
comprehensive examples for heterogeneous datasets. :ref:`pca-xray-analysis` 
shows more analysis function usage examples and :ref:`dynamics` module 
documentation lists all of the analysis functions. 

External data 
-------------------------------------------------------------------------------

Normal mode data from other NMA, EDA, or PCA programs can be parsed using
:func:`~.parseModes` function for ProDy analysis. 

In this case, we will parse ANM modes for p38 MAP Kinase calculated using 
`ANM server <http://ignmtest.ccbb.pitt.edu/cgi-bin/anm/anm1.cgi>`_  as the 
external software.  We use :download:`oanm.eigvals <doctest/oanm_eigvals.txt>` 
and :download:`oanm.slwevs <doctest/oanm_slwevs.txt>` files from the ANM 
server. 

You can either download these files to your current working directory from here
or obtain them for another protein from the ANM server.

>>> nma = parseModes(normalmodes='oanm_slwevs.txt', 
...                  eigenvalues='oanm_eigvals.txt', 
...                  nm_usecols=range(1,21), 
...                  ev_usecols=[1], ev_usevalues=range(6,26))
>>> nma
<NMA: oanm_slwevs (20 modes, 351 atoms)>
>>> nma.setTitle('1p38 ANM')
>>> slowmode = nma[0]
>>> print( slowmode.getEigenvalue().round(2) )
0.18

.. plot::
   :context:
   :nofigs:
   
   nma = parseModes(normalmodes='oanm_slwevs.txt', 
                    eigenvalues='oanm_eigvals.txt', 
                    nm_usecols=range(1,21), ev_usecols=[1], 
                    ev_usevalues=range(6,26))
   nma.setTitle('1p38 ANM')
   slowmode = nma[0]

Plotting data 
-------------------------------------------------------------------------------

If you have `Matplotlib`_, you can use ProDy functions whose name start with
``show`` to plot data:

.. plot::
   :include-source:
   :context:
   
   plt.figure( figsize=(5,4) )
   showSqFlucts( slowmode )
   
.. plot::
   :nofigs:
   :context:
   
   plt.close('all')
   
      
:ref:`pca-xray-plotting` shows more plotting examples and 
:ref:`dynamics` module documentation lists all of the plotting functions. 


|questions|

|suggestions|
