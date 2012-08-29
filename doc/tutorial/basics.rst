.. _prody-basics:

*******************************************************************************
ProDy Basics
*******************************************************************************

In the rest of this tutorial, it is assumed that the user is typing commands 
in a Python shell.  ProDy will automatically download and save some files,
so you may want start a Python shell from inside of a directory that you may 
make for this tutorial::

  $ mkdir pdtut
  $ cd pdtut
  $ ipython --pylab
  
if you are using or::
  
  $ ipython -pylab
  
if you are using an older version of `IPython <ipython.org>`_.

On Windows, after you make the directory, you can make a shift-right click in 
it in Windows Explorer and then select :menuselection:`Open command window here`
option. 


Import from ProDy
===============================================================================

To begin the Tutorial, import all the functions and classes from ProDy into the
current namespace as follows:

>>> from prody import *

There are other ways to import ProDy contents.  You may use 
``import prody as pd`` and prefix all functions calls with ``pd.``, if you 
prefer not to overcrowd your namespace.  Alternatively, if you want to use 
functions in a specific module, :mod:`.proteins` let's say, you can use 
``from prody.proteins import *``.  You should, however, avoid using 
``from prody.proteins.pdbfile import *``, because location of methods in 
submodules may change without notice.

.. plot::
   :nofigs: 
   :context: 
    
   from prody import *
   import matplotlib.pyplot as plt
   import numpy as np
   structure = parsePDB('1p38')

   plt.close('all')
   
Protein structure
===============================================================================

Protein structure files in :file:`.pdb` format are the standard input for 
ProDy.  PDB files are parsed using :func:`.parsePDB` (see 
:ref:`pdbparser-performance` and :ref:`pdbparser-performance-2` for 
benchmarks).  It is sufficient to pass a PDB identifier to read the file, and 
the parser will download it automatically if needed.

>>> structure = parsePDB('1p38')

In the above line, :file:`1p38.pdb` is downloaded and coordinates and atomic 
data are parsed from the file. 
  
Managing resources
-------------------------------------------------------------------------------

You can tell ProDy where to get PDB files from or where to store downloaded 
files:

  * one of the `wwPDB <http://www.wwpdb.org/>`_ FTP servers in US, Europe or
    Japan can be picked for downloads using :func:`.setWWPDBFTPServer`
  * a local PDB mirror can be be set for faster access to files using
    :func:`.setPDBMirrorPath` 
  * a local folder can be set for storing downloaded files for later access
    using :func:`.setPDBLocalFolder` 

Note that when these functions are used, ProDy will save your settings as 
:file:`.prodyrc` in your home folder.


Quick visualization
-------------------------------------------------------------------------------

:file:`1p38.pdb` contains an unbound structure of the p38 MAP kinase.
If you have `Matplotlib <http://matplotlib.sourceforge.net>`_ installed, you 
can take a quick look at what you parsed using :func:`.showProtein` function:  


>>> import matplotlib.pyplot as plt
>>> plt.figure(figsize=(5,4)) # doctest: +SKIP
>>> showProtein(structure) # doctest: +SKIP
>>> plt.legend(prop={'size': 10}) # doctest: +SKIP

.. plot::
   :context:
   
   import matplotlib.pyplot as plt
   plt.figure(figsize=(5,4))
   showProtein(structure)
   plt.legend(prop={'size': 10})

   
.. plot::
   :nofigs: 
   :context: 
   
   plt.close('all')

More examples
-------------------------------------------------------------------------------

:func:`.parsePDB` function is very flexible and can be extremely
efficient depending on what you want to extract from a PDB file.  It can be 
used to parse specific chains, models, alternate locations, or well-defined 
subsets of atoms from a file.  A detailed usage example can be found in 
:ref:`parsepdb`.  

ProDy can parse other file types, including :file:`.psf` and :file:`.pqr` files.
All of the functions for accessing and handling protein structural data are 
described in :mod:`.proteins` module reference documentation.
Also, :ref:`fetchpdb` and :ref:`blastpdb` examples show other ways to 
access the Protein Data Bank (|pdb|) content.

For more details on atomic objects see :ref:`atomic`.  
:class:`.AtomGroup` instances can be build from scratch or 
parsers for other file types (e.g. mol2) can be developed. The example in 
:ref:`atomgroup` can be helpful to this aim.

Atomic Data
===============================================================================

:func:`.parsePDB` returns data in an :class:`.AtomGroup` instance.  
To get information on an :class:`.AtomGroup` instance, type in the 
variable name and hit :kbd:`enter` key:

>>> structure
<AtomGroup: 1p38 (2962 atoms)>

The above shows that atom group object contains 2962 atoms. 
All atomic data from this object can be retrieved using ``get`` methods. 
For example:

>>> print( structure.getResnames() )
['GLU' 'GLU' 'GLU' ..., 'HOH' 'HOH' 'HOH']
>>> print( structure.getCoords() ) # doctest: +ELLIPSIS
[[ 28.492   3.212  23.465]
 [ 27.552   4.354  23.629]
 ...
 [-22.062  21.632  42.029]
 [  1.323  30.027  65.103]]
 
The list of methods for getting and setting atomic data is provided in
:class:`.AtomGroup` reference documentation. 

**Indexing**:

An individual :class:`.Atom` can be accessed by indexing atom group 
instances:

>>> atom = structure[0]
>>> atom
<Atom: N from 1p38 (index 0)>

Note that all ``get/set`` functions defined for :class:`.AtomGroup` 
instances are also defined for :class:`.Atom` instances, using singular
form of the function name.  

>>> atom.getResname()
'GLU'

**Slicing**:

It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

>>> structure[::2]
<Selection: 'index 0:2962:2' from 1p38 (1481 atoms)>

Hierarchical view
-------------------------------------------------------------------------------

You can also access specific chains or residues in an atom group.  Indexing
by a single letter identifier will return a :class:`.Chain` instance: 

>>> structure['A']
<Chain: A from 1p38 (480 residues, 2962 atoms)>

Indexing atom group with a chain identifier and a residue number will return
:class:`.Residue` instance:

>>> structure['A', 100]
<Residue: ASN 100 from Chain A from 1p38 (8 atoms)>

See :ref:`atomic` for details of indexing atom groups and :ref:`hierview`
for more on hierarchical views.


Writing PDB files
-------------------------------------------------------------------------------

PDB files can be written using the :func:`.writePDB` function.
The function accepts objects containing or referring to atomic data.

Writing selected atoms:

>>> writePDB('1p38_calphas.pdb', structure.select('calpha'))
'1p38_calphas.pdb'

Writing a chain:

>>> chain_A = structure['A']
>>> writePDB('1p38_chain_A.pdb', chain_A)
'1p38_chain_A.pdb'

As you may have noticed, this function returns the file name after it is
successfully written.  This is a general behavior for ProDy output functions.
For more PDB writing examples see :ref:`writepdb`.


Atom Selections
===============================================================================

:class:`.AtomGroup` instances have a plain view of atoms for efficiency, 
but they are coupled with a powerful atom selection engine.  You can get well 
defined atom subsets by passing simple keywords or make rather sophisticated 
selections using composite statements.  Selection keywords and grammar is very 
much similar to those found in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.  
Some examples are shown here:

Keyword selections
-------------------------------------------------------------------------------

>>> protein = structure.select('protein')
>>> protein
<Selection: 'protein' from 1p38 (2833 atoms)>

Using the "protein" keyword we selected 2833 atoms out of 2962 atoms. 
:meth:`~.Atomic.select` method returned a :class:`.Selection` instance.  
Note that all ``get`` and ``set`` methods defined for the :class:`.AtomGroup` 
class are also defined for :class:`.Selection` class. For example:

>>> print( protein.getResnames() )
['GLU' 'GLU' 'GLU' ..., 'ASP' 'ASP' 'ASP']

Select by name/type
-------------------------------------------------------------------------------

We select backbone atoms by passing atom names following "name" keyword:

>>> backbone = structure.select('protein and name N CA C O')
>>> backbone
<Selection: 'protein and name N CA C O' from 1p38 (1404 atoms)>
>>> len(backbone)
1404

We can also use "backbone" to make the same selection. 

We select acidic and basic residues by using residue names with 
"resname" keyword:

>>> structure.select('resname ARG LYS HIS ASP GLU')
<Selection: 'resname ARG LYS HIS ASP GLU' from 1p38 (906 atoms)>

Alternatively, we can use predefined keywords "acidic" and "basic".

>>> charged = structure.select('acidic or basic')
>>> charged
<Selection: 'acidic or basic' from 1p38 (906 atoms)>
>>> set(charged.getResnames())
set(['HIS', 'ASP', 'LYS', 'GLU', 'ARG'])

Composite selections
-------------------------------------------------------------------------------

Let's try a more sophisticated selection.  We first calculate the geometric 
center of the protein atoms using :func:`.calcCenter` function.  Then, we 
select the Cα and Cβ atoms of residues that have at least one atom within 
10 Å away from the geometric center.

>>> center = calcCenter(protein).round(3)
>>> print( center )
[  1.005  17.533  40.052]
>>> structure.select('protein and name CA CB and same residue as ((x-1)**2 + (y-17.5)**2 + (z-40.0)**2)**0.5 < 10')
<Selection: 'protein and nam...)**2)**0.5 < 10' from 1p38 (66 atoms)>

Alternatively, this selection could be done as follows:

>>> structure.select('protein and name CA CB and same residue as within 10 of center', center=center)
<Selection: 'index 576 579 5... 1687 1707 1710' from 1p38 (66 atoms)>

Selection operations
-------------------------------------------------------------------------------

:class:`.Selection` instances can used with bitwise operators:

>>> ca = structure.select('name CA') 
>>> cb = structure.select('name CB')
>>> ca | cb
<Selection: '(name CA) or (name CB)' from 1p38 (687 atoms)>
>>> ca & cb

Selections simplified
-------------------------------------------------------------------------------

In interactive sessions, an alternative to typing in ``.select('protein')`` 
or ``.select('backbone')`` is using dot operator:

>>> structure.protein
<Selection: 'protein' from 1p38 (2833 atoms)>

You can use dot operator multiple times:

>>> structure.protein.backbone
<Selection: '(backbone) and (protein)' from 1p38 (1404 atoms)>

This may go on and on:

>>> structure.protein.backbone.resname_ALA.calpha
<Selection: '(calpha) and ((...and (protein)))' from 1p38 (26 atoms)>


More examples
-------------------------------------------------------------------------------

There is much more to what you can do with this flexible and fast atom 
selection engine, without the need for writing nested loops with comparisons 
or changing the source code.  See the following pages:

  * :ref:`selections` for description of all selection keywords
  * :ref:`selection-operations` for handy features of :class:`.Selection`
  * :ref:`contacts` for selecting interacting atoms
  

Analyze Structures
===============================================================================

Measure geometric properties
-------------------------------------------------------------------------------

ProDy offers several functions for analyzing molecular structure in 
:mod:`~prody.measure` module. For example, you can calculate phi (φ) and psi 
(ψ) for the 10th residue, or the radius of gyration of the protein as follows:

>>> print calcPhi(structure[10,]).round(2)
-115.54
>>> print calcPsi(structure[10,]).round(2)
147.49
>>> print calcGyradius(structure).round(2)
22.06

Most functions start with ``calc`` prefix.  You can type in ``calc`` and press
tab to see what is available in an IPython session.

Compare chains
-------------------------------------------------------------------------------

You can also compare different structures using some of the methods in
:mod:`.proteins` module.  Let's parse another p38 MAP kinase structure

>>> bound = parsePDB('1zz2')

You can find similar chains in structure 1p38 and 1zz2 using 
:func:`.matchChains` function:

>>> apo_chA, bnd_chA, seqid, overlap = matchChains(structure, bound)[0]
>>> apo_chA
<AtomMap: Chain A from 1p38 -> Chain A from 1zz2 from 1p38 (337 atoms)>
>>> bnd_chA
<AtomMap: Chain A from 1zz2 -> Chain A from 1p38 from 1zz2 (337 atoms)>
>>> print int(seqid) # precent sequence identity between two chains
99
>>> print int(overlap) # percent overlap between two chains
96

Matching Cα atoms are selected and returned as :class:`.AtomMap` instances.
We can use them to calculate RMSD and superpose structures.

>>> print calcRMSD(bnd_chA, apo_chA)
72.9302308695
>>> bnd_chA, transformation = superpose(bnd_chA, apo_chA)
>>> print calcRMSD(bnd_chA, apo_chA)
1.86280149087

>>> plt.figure(figsize=(5,4)) # doctest: +SKIP
>>> showProtein(structure) # doctest: +SKIP
>>> showProtein(bound) # doctest: +SKIP
>>> plt.legend(prop={'size': 10}) # doctest: +SKIP

.. plot::
   :context:
   
   import matplotlib.pyplot as plt
   bound = parsePDB('1zz2')
   matchAlign(structure, bound)
   plt.figure(figsize=(5,4))
   showProtein(structure)
   showProtein(bound)
   plt.legend(prop={'size': 10})

   
.. plot::
   :nofigs: 
   :context: 
   
   plt.close('all')


ProDy Verbosity
===============================================================================

Finally, you might have noted that ProDy prints some information to the console
after parsing a file or doing some calculations. For example, PDB parser will 
print what was parsed and how long it took to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

This behavior is useful in interactive sessions, but may be problematic for
automated tasks as the messages are printed to stderr.  The level of verbosity 
can be controlled using :func:`.confProDy` function, and calling it as 
``confProDy(verbosity='none')`` will stop all information messages permanently.

|questions|

|suggestions|
