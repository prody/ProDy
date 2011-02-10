.. _scripts:

*******************************************************************************
Scripts
*******************************************************************************

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
and all Cα atoms of the protein structure 1p38 are used.

In the following, cutoff distance is changed to 14 Å, 
Cα atoms of residues with numbers smaller than 340 are used, 
and output files are prefixed with "p38_anm"::

  ./python anm.py -c 14 -s "calpha resnum < 340" -p p38_anm 1p38

Output file :file:`p38_anm.nmd` file can be visualized using |nmwiz|. 

.. _scripts-anm:

:file:`anm.py`
===============================================================================

This script can be used to perform ANM calculations and output 
eigenvalues/eigenvectors in plain text and NMD formats, i.e. running ``python gnm.py 1p38``
will perform calculations for PDB structure 1p38. Resulting NMD files can be 
visualized using |nmwiz|.

Running ``python anm.py -h`` (or :command:`python anm -h`) prints::

  Usage: anm.py [options] <pdb>

  ProDy - Anisotropic Network Model

  Options:
    -h, --help            show this help message and exit
    -n INT, --nmodes=INT  number of non-zero eigenvalues/vectors to calculate
                          (20)
    -c FLOAT, --cutoff=FLOAT
                          cutoff distance (15.0)
    -g FLOAT, --gamma=FLOAT
                          spring constant (1.0)
    -p STRING, --prefix=STRING
                          prefix for output files (anm)
    -s STRING, --select=STRING
                          selection string (protein and name CA)
    --silent              omit verbose information (False)
    -e, --examples        show usage examples

.. _scripts-gnm:

:file:`gnm.py`
===============================================================================

This script can be used to perform GNM calculations and output 
eigenvalues/eigenvectors in plain text format, i.e. running ``python gnm.py 1p38``
will perform calculations for PDB structure 1p38. 
 
Running ``python gnm.py -h`` prints::

  Usage: gnm.py [options] <pdb>

  ProDy - Gaussian Network Model

  Options:
    -h, --help            show this help message and exit
    -n INT, --nmodes=INT  number of non-zero eigenvalues/vectors to calculate
                          (20)
    -c FLOAT, --cutoff=FLOAT
                          cutoff distance (10.0)
    -g FLOAT, --gamma=FLOAT
                          spring constant (1.0)
    -p STRING, --prefix=STRING
                          prefix for output files (gnm)
    -s STRING, --select=STRING
                          selection string (protein and name CA)
    --silent              omit verbose information (False)
    -e, --examples        show usage examples

.. _scripts-pdbselect:

:file:`pdbselect.py`
===============================================================================

This script can be used to extract a selection of atoms from a PDB file, i.e. 
running ``python pdbselect.py 1p38 selected.pdb "protein and name CA"``
will write Cα atoms in :file:`selected.pdb` file.
 
Running ``python pdbselect.py -h`` prints::

  ProDy - PDBSelect
  usage: pdbselect.py <input> <output> <selection>
