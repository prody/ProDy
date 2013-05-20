.. _prody-basics:

ProDy Basics
===============================================================================


We start with importing everything from ProDy package:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()


Functions and classes are named such that they should not create a conflict
with any other package.  In this part we will familiarize with different
categories of functions and methods.


File Parsers
-------------------------------------------------------------------------------

Let's start with parsing a protein structure and then keep working on that
in this part.  File parser function names are prefixed with ``parse``.
You can get a list of parser functions by pressing :kbd:`TAB` after typing
in ``parse``:

.. ipython::

   @verbatim
   In [2]: parse<TAB>
   parseArray         parseHeatmap       parseNMD           parsePDBStream
   parseSTRIDE        parseDCD           parseMSA           parsePDB
   parsePQR           parseSparseMatrix  parseDSSP          parseModes
   parsePDBHeader     parsePSF


When using :func:`parsePDB`, an identifier will be sufficient, if your
machine is connected to the internet.  If corresponding file is not found
in the current working directory, it will be downloaded from PDB servers.

Let's parse structure :pdb:`1p38` of p38 MAP kinase (MAPK):

.. ipython:: python

   p38 = parsePDB('1p38') # returns an AtomGroup object
   p38 # typing in variable name will give some information


Similar to listing parser function names, we can use tab completion to
introspect ``p38`` object:

.. ipython::

   @verbatim
   In [2]: p38.num<TAB>
   p38.numAtoms      p38.numChains     p38.numFragments  p38.numSegments
   p38.numBonds      p38.numCoordsets  p38.numResidues


This action printed a list of methods with `num` prefix. Let's use some of
them to get information on the structure:


.. ipython::

   @doctest
   In [2]: p38.numAtoms()
   Out[2]: 2962

   @doctest
   In [2]: p38.numCoordsets() # returns number of models
   Out[2]: 1

   @doctest
   In [2]: p38.numResidues()  # water molecules also count as residues
   Out[2]: 480


Analysis Functions
-------------------------------------------------------------------------------

Similarly to parsers, analysis functions have a ``calc`` prefix:

.. ipython::

   @verbatim
   In [2]: calc<TAB>
   calcADPAxes          calcCrossProjection  calcMSF              calcRMSF
   calcADPs             calcCumulOverlap     calcOccupancies      calcRankorder
   calcANM              calcDeformVector     calcOmega            calcShannonEntropy
   calcAngle            calcDihedral         calcOverlap          calcSqFlucts
   calcCenter           calcDistance         calcPerturbResponse  calcSubspaceOverlap
   calcCollectivity     calcFractVariance    calcPhi              calcTempFactors
   calcCovOverlap       calcGNM              calcProjection       calcTransformation
   calcCovariance       calcGyradius         calcPsi
   calcCrossCorr        calcMSAOccupancy     calcRMSD


Let's read documentation of :func:`calcGyradius` function and use it to
calculate the radius of gyration of p38 MAPK structure:

.. ipython::

   In [1]: calcGyradius ?

   In [1]: calcGyradius(p38)


Plotting Functions
-------------------------------------------------------------------------------

Likewise, plotting function names have ``plot`` prefix and here is a list
of them:

.. ipython::

   @verbatim
   In [2]: show<TAB>
   showContactMap       showEllipsoid        showNormedSqFlucts   showScaledSqFlucts
   showCrossCorr        showFractVars        showOccupancies      showShannonEntropy
   showCrossProjection  showHeatmap          showOverlap          showSqFlucts
   showCumulFractVars   showMSAOccupancy     showOverlapTable
   showCumulOverlap     showMode             showProjection
   showDiffMatrix       showMutinfoMatrix    showProtein

We can use :func:`showProtein` function to make a quick plot of p38 structure:

.. ipython:: python

   @savefig prody_tutorial_basics_protein.png width=4in
   showProtein(p38);

This of course does not compare to any visualization software that you
might be familiar with, but it comes handy to see what you are dealing with.



Protein Structures
-------------------------------------------------------------------------------

Protein structures (:file:`.pdb` files) will be the standard input for most
*ProDy* calculations, so it is good to familiarize with ways to access and
manage PDB files.

First of all, *ProDy* downloads compressed PDB files when needed.  If you
prefer saving decompressed files, you can use :func:`.fetchPDB` function as
follows:

.. ipython:: python

  fetchPDB('1p38', compressed=False)

Note that ProDy functions that fetch files or output files return filename
upon successful completion of the task.  You can use this behavior to
minimize the code you write as follows:


.. ipython:: python

  parsePDB(fetchPDB('1p38', compressed=False)) # same as p38 parsed above

Secondly, ProDy can manage local mirror of PDB server or a local PDB folders,
as well as using a server close to your physical location for downloads:

  * One of the `wwPDB`_ FTP servers in US, Europe or Japan can be picked for
    downloads using :func:`.wwPDBServer`.

  * A local PDB mirror can be set for faster access to files using
    :func:`.pathPDBMirror`.


  * A local folder can be set for storing downloaded files for future access
    using :func:`.pathPDBFolder`.

If you are in the Americas now, you can choose the PDB server in the US
as follows:

.. ipython:: python

   wwPDBServer('us')

If you would like to have a central folder, such as :file:`~Downloads/pdb`,
for storing downloaded PDB files (you will need to make it), do as follows:

.. ipython:: python

   mkdir /home/abakan/Downloads/pdb;
   pathPDBFolder('/home/abakan/Downloads/pdb')

Note that when these functions are used, ProDy will save your settings
in :file:`.prodyrc` file stored in your home folder.

.. _wwPDB: http://www.wwpdb.org/

..
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

Atom Groups
-------------------------------------------------------------------------------

:func:`.parsePDB` returns structure data in an :class:`.AtomGroup` instance,
such as the ``p38`` variable we parsed above:

.. ipython:: python

   p38


The above shows that atom group object contains 2962 atoms.  Data
from this object can be retrieved using ``get`` methods.  For example:


.. ipython:: python

   p38.getResnames()
   p38.getCoords()


The get a list of all methods use tab completion, i.e. ``p38.<TAB>``.  We
will learn more about atom groups in the following chapters.

Indexing
^^^^^^^^

An individual :class:`.Atom` can be accessed by indexing atom group
instances:

.. ipython:: python

   atom = p38[0]
   atom

Note that all ``get/set`` functions defined for :class:`.AtomGroup`
instances are also defined for :class:`.Atom` instances, using singular
form of the function name.

.. ipython:: python

   atom.getResname()

An individual :class:`.Atom` can be accessed by indexing atom group
instances:

.. ipython:: python

   atom = p38[0]
   atom

Note that all ``get/set`` functions defined for :class:`.AtomGroup`
instances are also defined for :class:`.Atom` instances, using singular
form of the function name.

.. ipython:: python

   atom.getResname()

Slicing
^^^^^^^

It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]



It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]



It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]



It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]



It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]



It is also possible to get a slice of an atom group, for example we can get
every other atom as follows:

.. ipython:: python

   p38[::2]


Hierarchical view
^^^^^^^^^^^^^^^^^

You can also access specific chains or residues in an atom group.  Indexing
by a single letter identifier will return a :class:`.Chain` instance:

.. ipython:: python

   p38['A']

Indexing atom group with a chain identifier and a residue number will return
:class:`.Residue` instance:
..

.. ipython:: python

   p38['A', 100]


See :ref:`atomic` for details of indexing atom groups and :ref:`hierview`
for more on hierarchical views.



ProDy Verbosity
-------------------------------------------------------------------------------

Finally, you might have noted that ProDy prints some information to the console
after parsing a file or doing some calculations. For example, PDB parser will
print what was parsed and how long it took to the screen::

  @> 1p38 (./1p38.pdb.gz) is found in the target directory.
  @> PDBParser: 2962 atoms and 1 coordinate sets were parsed in 0.08s.

This behavior is useful in interactive sessions, but may be problematic for
automated tasks as the messages are printed to stderr.  The level of verbosity
can be controlled using :func:`.confProDy` function, and calling it as
``confProDy(verbosity='none')`` will stop all information messages permanently.

