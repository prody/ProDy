.. _scripts:

*******************************************************************************
Scripts
*******************************************************************************

There are three stand-alone scripts distributed with ProDy in the 
:file:`scripts` folder. These scripts are written by Lidio Meireles.

:file:`anm.py`
===============================================================================

This script can be used to perform ANM calculations and output 
eigenvalues/eigenvectors in plain text and NMD formats. NMD files can be 
visualized using |nmwiz|.
 
Running :command:`python anm.py -h` prints::

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
                          prefix for output files (anm_)
    -s STRING, --select=STRING
                          selection string (protein and name CA)
    --silent              omit verbose information (False)
    -e, --examples        show usage examples

:file:`gnm.py`
===============================================================================

This script can be used to perform GNM calculations and output 
eigenvalues/eigenvectors in plain text format.
 
Running :command:`python gnm.py -h` prints::

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
                          prefix for output files (gnm_)
    -s STRING, --select=STRING
                          selection string (protein and name CA)
    --silent              omit verbose information (False)
    -e, --examples        show usage examples


:file:`pdbselect.py`
===============================================================================

This script can be used to extract a selection of atoms from a PDB file.
 
Running :command:`python pdbselect.py -h` prints::

  ProDy - PDBSelect
  usage: pdbselect.py <input> <output> <selection>
