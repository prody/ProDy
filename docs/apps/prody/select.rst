.. _prody-select:

prody select
====================

Usage
--------------------

Running :command:`prody select -h` displays::

  usage: prody select [-h] [--quiet] [--examples] [-o STR] [-p STR] [-x STR]
                      select pdb [pdb ...]
  
  positional arguments:
    select                atom selection string
    pdb                   PDB identifier(s) or filename(s)
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  output options:
    -o STR, --output STR  output PDB filename (default: pdb_selected.pdb)
    -p STR, --prefix STR  output filename prefix (default: PDB filename)
    -x STR, --suffix STR  output filename suffix (default: _selected)

Examples
--------------------

Running :command:`prody select --examples` displays::

  This command selects specified atoms and writes them in a PDB file.
  
  Fetch PDB files 1p38 and 1r39 and write backbone atoms in a file:
  
    $ prody select backbone 1p38 1r39
  
