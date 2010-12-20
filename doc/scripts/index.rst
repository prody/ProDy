.. _scripts:

*******************************************************************************
Scripts
*******************************************************************************

There are three stand-alone scripts distributed with ProDy in the 
:file:`scripts` folder. Alternatively, they can be obtained from 
https://github.com/abakan/ProDy/tree/master/scripts/.  
These scripts are written by Lidio Meireles.

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
will write alpha carbon atoms in :file:`selected.pdb` file.
 
Running ``python pdbselect.py -h`` prints::

  ProDy - PDBSelect
  usage: pdbselect.py <input> <output> <selection>
