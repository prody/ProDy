.. _prody-catdcd:

*******************************************************************************
prody catdcd
*******************************************************************************

Concatenate, slice, and/or select DCD files.

Usage
===============================================================================

Running :command:`prody catdcd -h` displays::

  usage: prody catdcd [-h] [--quiet] [--examples] [-s SELSTR] [-o FILENAME] [-n]
                      [--psf PSF | --pdb PDB] [--first INT] [--last INT]
                      [--stride INT] [--align SELSTR]
                      dcd [dcd ...]
  
  positional arguments:
    dcd                   DCD filename(s) (all must have same number of atoms)
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -s SELSTR, --select SELSTR
                          atom selection (default: "all")
    -o FILENAME, --output FILENAME
                          output filename (default: "trajectory.dcd")
    -n, --num             print the number of frames in each file and exit
    --psf PSF             PSF filename (must have same number of atoms as DCDs)
    --pdb PDB             PDB filename (must have same number of atoms as DCDs)
    --first INT           the first frame to be written to the output file
                          (default: 0, first frame)
    --last INT            the last frame to be written to the output file
                          (default: -1, last frame)
    --stride INT          number of frames to skip when writing (default: 1,
                          skip none)
    --align SELSTR        atom selection for aligning frames, one of PSF or PDB
                          files must be provided

Examples
===============================================================================

Running :command:`prody catdcd --examples` displays::

  Concatenate two DCD files and output all atmos:
  
    $ prody catdcd mdm2.dcd mdm2sim2.dcd
  
  Concatenate two DCD files and output backbone atoms:
  
    $ prody catdcd mdm2.dcd mdm2sim2.dcd --pdb mdm2.pdb -s bb
