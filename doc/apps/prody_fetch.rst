.. _prody-fetch:

*******************************************************************************
prody fetch
*******************************************************************************

Usage
===============================================================================

Running :command:`prody fetch -h` displays::

  usage: prody fetch [-h] [--quiet] [--examples] [-d PATH] [-z] pdb [pdb ...]
  
  positional arguments:
    pdb                  PDB identifier(s) or a file that contains them
  
  optional arguments:
    -h, --help           show this help message and exit
    --quiet              suppress info messages to stderr
    --examples           show usage examples and exit
    -d PATH, --dir PATH  target directory for saving PDB file(s)
    -z, --gzip           write compressed PDB file(S)

Examples
===============================================================================

Running :command:`prody fetch --examples` displays::

  Download PDB file(s) by specifying identifiers:
  
    $ prody fetch 1mkp 1p38
