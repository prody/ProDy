.. _prody-select:

*******************************************************************************
prody select
*******************************************************************************

Extract a selection of atoms from a PDB file.

Usage
===============================================================================

Running :command:`prody select -h` displays::

  usage: prody select [-h] [--quiet] [--examples] [-o STR] pdb selstr
  
  positional arguments:
    pdb                   PDB identifier or filename
    selstr                atom selection string
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -o STR, --output STR  output filanem (default: 'pdb_selected.pdb')

Examples
===============================================================================

Running :command:`prody select --examples` displays::

  This command selects specified atoms and writes them in a PDB file.
  
  Fetch PDB 2bfu and write backbone atoms in a file:
  
    $ prody select 2bfu "backbone"
