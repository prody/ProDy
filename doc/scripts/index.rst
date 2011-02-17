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

This script is used to perform ANM calculations and output the 
eigenvalues/eigenvectors in plain text and NMD formats, i.e. running 
:command:`anm.py 1p38` will perform calculations for PDB structure 1p38. 
The resulting NMD files can be visualized using NMWiz.

Running :command:`anm.py -h` displays:

.. literalinclude:: anm.txt


.. _scripts-gnm:

:file:`gnm.py`
===============================================================================

This script is used to perform GNM calculations and output the 
eigenvalues/eigenvectors in plain text format, i.e. running 
:command:`gnm.py 1p38` will perform calculations for PDB structure 1p38. 
 
Running :command:`gnm.py -h` displays:

.. literalinclude:: gnm.txt


:file:`pca.py`
===============================================================================

This script is used to perform PCA calculations and output the 
eigenvalues/eigenvectors in plain text and NMD formats, i.e. running 
``pca.py 2k39`` will perform calculations for NMR models in structure 
2k39. The resulting NMD files can be visualized using NMWiz.

Running :command:`pca.py -h` displays:

.. literalinclude:: pca.txt


.. _scripts-alignmodels:

:file:`alignmodels.py`
===============================================================================

This script can be used to align models in a PDB file.

Running :command:`alignmodels.py -h` displays:

.. literalinclude:: alignmodels.txt


.. _scripts-biomolecule:

:file:`biomolecule.py`
===============================================================================
 
This script can be used to generate biomolecule structure using
the transformation in header section of the PDB file.
 
Running :command:`biomolecule.py -h` displays:

.. literalinclude:: biomolecule.txt


.. _scripts-blastpdb:

:file:`blastpdb.py`
===============================================================================

This script can be used to download PDB files matching a user given sequence.

Running :command:`blastpdb.py -h` displays:

.. literalinclude:: blastpdb.txt


.. _scripts-fetchpdb:

:file:`fetchpdb.py`
===============================================================================

This script can be used to download PDB for given identifiers.
 
Running :command:`fetchpdb.py -h` displays:

.. literalinclude:: fetchpdb.txt


.. _scripts-pdbselect:

:file:`pdbselect.py`
===============================================================================

This script is used to extract a selection of atoms from a PDB file, i.e. 
running :command:`pdbselect.py 1p38 "protein and name CA"`
will write Cα atoms in :file:`1p38_selected.pdb` file.
 
Running :command:`pdbselect.py -h` displays:

.. literalinclude:: pdbselect.txt
