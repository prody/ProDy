.. _prody-biomol:

*******************************************************************************
prody biomol
*******************************************************************************

Usage
===============================================================================

Running :command:`prody biomol -h` displays::

  usage: prody biomol [-h] [--quiet] [--examples] [-p STR] [-b INT] pdb
  
  positional arguments:
    pdb                   PDB identifier or filename
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -p STR, --prefix STR  prefix for output files (default: pdb_biomol_)
    -b INT, --biomol INT  index of the biomolecule, by default all are generated

Examples
===============================================================================

Running :command:`prody biomol --examples` displays::

  Generate biomolecule coordinates:
  
  $ prody biomol 2bfu
