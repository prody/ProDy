.. _scripts:

*******************************************************************************
Scripts
*******************************************************************************

Command line usage
===============================================================================

ProDy scripts come with both the source and the binary distributions.
On Linux, when installing ProDy from source, the scripts are placed into a 
default folder that is included in the environment variable :envvar:`PATH`, 
e.g. :file:`/usr/local/bin/`. 
On Windows, installer places the scripts into the :file:`Scripts` folder 
under the path to the corresponding Python distribution, 
e.g. :file:`C:\Python27\Scripts` if you used Python 2.7. 
You may need to add this path to the environment variable :envvar:`PATH` 
yourself. 
  

ProDy scripts cane used as command line programs. Running the following 
command will perform ANM calculations for the p38 MAP kinase structure, and 
will write eigenvalues/vectors in plain text and :term:`NMD` formats::

  $ anm.py 1p38
  
In the above example, the default parameters (``cutoff=15.`` and ``gamma=1.``)
and all of the Cα atoms of the protein structure 1p38 are used.

In the example below, the *cutoff* distance is changed to 14 Å, 
and the Cα atoms of residues with numbers smaller than 340 are used, 
the output files are prefixed with :file:`p38_anm`::

  $ anm.py -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

The output file :file:`p38_anm.nmd` can be visualized using NMWiz (|nmwiz|). 

.. _scripts-anm:


:file:`anm.py`
===============================================================================

Perform ANM calculations and output the results in plain text, NMD, and 
graphical formats.

Usage:
-------------------------------------------------------------------------------

Running :command:`anm.py -h` displays:

.. literalinclude:: anm.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`anm.py --examples` displays:

.. literalinclude:: anm_eg.txt


.. _scripts-gnm:

:file:`gnm.py`
===============================================================================

Perform GNM calculations and output the results in plain text and graphical 
formats. 
 
 
Usage:
-------------------------------------------------------------------------------

Running :command:`gnm.py -h` displays:

.. literalinclude:: gnm.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`gnm.py --examples` displays:

.. literalinclude:: gnm_eg.txt


.. _scripts-pca:

:file:`pca.py`
===============================================================================

Perform PCA calculations and output the results in plain text, NMD formats,
and graphical formats.

Usage:
-------------------------------------------------------------------------------

Running :command:`pca.py -h` displays:

.. literalinclude:: pca.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`pca.py --examples` displays:

.. literalinclude:: pca_eg.txt


.. _scripts-alignmodels:

:file:`alignmodels.py`
===============================================================================

Align models in a PDB file.

Usage:
-------------------------------------------------------------------------------

Running :command:`alignmodels.py -h` displays:

.. literalinclude:: alignmodels.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`alignmodels.py --examples` displays:

.. literalinclude:: alignmodels_eg.txt


.. _scripts-biomolecule:

:file:`biomolecule.py`
===============================================================================
 
Generate biomolecule structure using the transformation from the header 
section of the PDB file.
 
Usage:
-------------------------------------------------------------------------------
 
Running :command:`biomolecule.py -h` displays:

.. literalinclude:: biomolecule.txt

Examples:
-------------------------------------------------------------------------------
 
Running :command:`biomolecule.py --examples` displays:

.. literalinclude:: biomolecule_eg.txt

.. _scripts-blastpdb:

:file:`blastpdb.py`
===============================================================================

Search Protein Data Bank for structures matching a user given sequence.


Usage:
-------------------------------------------------------------------------------

Running :command:`blastpdb.py -h` displays:

.. literalinclude:: blastpdb.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`blastpdb.py --examples` displays:

.. literalinclude:: blastpdb_eg.txt

.. _scripts-fetchpdb:

:file:`fetchpdb.py`
===============================================================================

Download PDB for given identifiers.
 
Usage:
-------------------------------------------------------------------------------

Running :command:`fetchpdb.py -h` displays:

.. literalinclude:: fetchpdb.txt

Examples:
-------------------------------------------------------------------------------

Running :command:`fetchpdb.py --examples` displays:

.. literalinclude:: fetchpdb_eg.txt

.. _scripts-pdbselect:

:file:`pdbselect.py`
===============================================================================

Extract a selection of atoms from a PDB file.
 
Usage:
-------------------------------------------------------------------------------
 
Running :command:`pdbselect.py -h` displays:

.. literalinclude:: pdbselect.txt

Examples:
-------------------------------------------------------------------------------
 
Running :command:`pdbselect.py --examples` displays:

.. literalinclude:: pdbselect_eg.txt
