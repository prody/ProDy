.. _commands:

*******************************************************************************
ProDy Commands
*******************************************************************************

Command line usage
===============================================================================

ProDy command line programs, or scripts, come with both the source and the 
binary distributions.  On Linux, when installing ProDy from source, the scripts 
are placed into a default folder that is included in the environment variable 
:envvar:`PATH`, e.g. :file:`/usr/local/bin/`.  On Windows, installer places the
scripts into the :file:`Scripts` folder under the path to the corresponding 
Python distribution, e.g. :file:`C:\\Python27\\Scripts` if you used Python 2.7. 
You may need to add this path to the environment variable :envvar:`PATH` 
yourself. 

A list of ProDy commands can be obtained by running :command:`prody`::
  
  $ prody 
  
This will display available commands and short descriptions:

.. literalinclude:: prody.txt

To get more information on a specific command, type in command name, e.g.
:command:`prody anm`.

Example
-------------------------------------------------------------------------------

Running the following command will perform ANM calculations for the p38 MAP 
kinase structure, and will write eigenvalues/vectors in plain text and 
:term:`NMD` formats::

  $ prody anm 1p38
  
In the above example, the default parameters (``cutoff=15.`` and ``gamma=1.``)
and all of the Cα atoms of the protein structure 1p38 are used.

In the example below, the *cutoff* distance is changed to 14 Å, 
and the Cα atoms of residues with numbers smaller than 340 are used, 
the output files are prefixed with :file:`p38_anm`::

  $ prody anm -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

The output file :file:`p38_anm.nmd` can be visualized using NMWiz (|nmwiz|). 

.. _command-anm:


ProDy commands
===============================================================================

anm
-------------------------------------------------------------------------------

Perform ANM calculations and output the results in plain text, NMD, and 
graphical formats.

Running :command:`prody anm -h` displays:

.. literalinclude:: prody_anm.txt

Running :command:`prody anm --examples` displays:

.. literalinclude:: prody_anm_eg.txt


.. _command-gnm:

gnm
-------------------------------------------------------------------------------

Perform GNM calculations and output the results in plain text NMD, and 
graphical formats.
 
 
Running :command:`prody gnm -h` displays:

.. literalinclude:: prody_gnm.txt

Running :command:`prody gnm --examples` displays:

.. literalinclude:: prody_gnm_eg.txt


.. _command-pca:

pca
-------------------------------------------------------------------------------

Perform PCA calculations and output the results in plain text, NMD, and 
graphical formats.


Running :command:`prody pca -h` displays:

.. literalinclude:: prody_pca.txt

Download :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`.

Running :command:`prody pca --examples` displays:

.. literalinclude:: prody_pca_eg.txt


.. _command-eda:

eda
-------------------------------------------------------------------------------

Perform EDA calculations and output the results in plain text, NMD, and 
graphical formats.


Running :command:`prody eda -h` displays:

.. literalinclude:: prody_eda.txt

Download :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`.

Running :command:`prody eda --examples` displays:

.. literalinclude:: prody_eda_eg.txt


.. _command-align:

align
-------------------------------------------------------------------------------

Align models in a PDB file.

Running :command:`prody align -h` displays:

.. literalinclude:: prody_align.txt

Running :command:`prody align --examples` displays:

.. literalinclude:: prody_align_eg.txt


.. _command-biomol:

biomol
-------------------------------------------------------------------------------
 
Generate biomolecule structure using the transformation from the header 
section of the PDB file.
 
Running :command:`prody biomol -h` displays:

.. literalinclude:: prody_biomol.txt

Running :command:`prody biomol --examples` displays:

.. literalinclude:: prody_biomol_eg.txt

.. _command-blast:

blast
-------------------------------------------------------------------------------

Search Protein Data Bank for structures matching a user given sequence.

Running :command:`prody blast -h` displays:

.. literalinclude:: prody_blast.txt

Running :command:`prody blast --examples` displays:

.. literalinclude:: prody_blast_eg.txt

.. _command-catdcd:

catdcd
-------------------------------------------------------------------------------

Concatenate, slice, and/or reduce DCD files.

Running :command:`prody catdcd -h` displays:

.. literalinclude:: prody_catdcd.txt

Download :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`.

Running :command:`prody catdcd --examples` displays:

.. literalinclude:: prody_catdcd_eg.txt

.. _command-fetch:

fetch
-------------------------------------------------------------------------------

Download PDB for given identifiers.
 
Running :command:`prody fetch -h` displays:

.. literalinclude:: prody_fetch.txt

Running :command:`prody fetch --examples` displays:

.. literalinclude:: prody_fetch_eg.txt

.. _command-select:

select
-------------------------------------------------------------------------------

Extract a selection of atoms from a PDB file.

Running :command:`prody select -h` displays:

.. literalinclude:: prody_select.txt

Running :command:`prody select --examples` displays:

.. literalinclude:: prody_select_eg.txt
