.. _prody-fetch:

*******************************************************************************
prody fetch
*******************************************************************************

Download PDB files for given identifiers.

Usage
===============================================================================

Running :command:`prody fetch -h` displays::

  usage: prody fetch [-h] [--quiet] [--examples] [-d PATH] [-f FILE] [-z]
                     pdb [pdb ...]
  
  positional arguments:
    pdb                   PDB identifier(s)
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -d PATH, --dir PATH   target directory for saving PDB files
    -f FILE, --file FILE  file that contains PDB identifiers
    -z, --gzip            write compressed PDB file

Examples
===============================================================================

Running :command:`prody fetch --examples` displays::

  Download PDB file(s) by specifying identifiers:
  
  $ prody fetch 1mkp 1p38
