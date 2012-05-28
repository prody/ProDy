.. _prody-contacts:

*******************************************************************************
prody contacts
*******************************************************************************

Usage
===============================================================================

Running :command:`prody contacts -h` displays::

  usage: prody contacts [-h] [--quiet] [--examples] [-s SELSTR] [-r FLOAT]
                        [-t STR] [-p STR] [-x STR]
                        target ligand [ligand ...]
  
  positional arguments:
    target                target PDB identifier or filename
    ligand                ligand PDB identifier(s) or filename(s)
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -s SELSTR, --select SELSTR
                          selection string for target
    -r FLOAT, --radius FLOAT
                          contact radius (default: 4.0)
    -t STR, --extend STR  output same residue residue, chain, or segment as the
                          contacts
    -p STR, --prefix STR  output filename prefix (default: target filename)
    -x STR, --suffix STR  output filename suffix (default: _contacts)

Examples
===============================================================================

Running :command:`prody contacts --examples` displays::

  Identify contacts of a target structure with one or more ligands.
  
  Fetch PDB structure 1zz2, save PDB files for individual ligands, and
  identify
  contacting residues of the target protein:
  
      $ prody select -o lig_B11 1zz2 "resname B11"
      $ prody select -o lig_BOG 1zz2 "resname BOG"
      $ prody contacts -r 4.0 -t residue -s protein 1zz2 lig_B11.pdb
  lig_BOG.pdb
  
