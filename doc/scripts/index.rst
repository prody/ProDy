.. _scripts:

*******************************************************************************
Scripts
*******************************************************************************

There are three stand-alone scripts distributed with ProDy. These scripts
are contributed by Lidio Meireles.

:file:`anm.py`
===============================================================================

This script can be used to perform ANM calculations and output 
eigenvalues/eigenvectors in plain text and NMD formats. NMD files can be 
visualized using NMWiz.
 
Running :command:`anm.py -h` prints::

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

:file:`pdbselect.py`
===============================================================================
