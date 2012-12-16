.. _prody-catdcd:

*******************************************************************************
prody catdcd
*******************************************************************************

Usage
===============================================================================

Running :command:`prody catdcd -h` displays::

  usage: prody catdcd [-h] [--quiet] [--examples] [-s SEL] [-o FILE] [-n]
                      [--psf PSF] [--pdb PDB] [--first INT] [--last INT]
                      [--stride INT] [--align SEL]
                      dcd [dcd ...]
  
  positional arguments:
    dcd                   DCD filename(s) (all must have same number of atoms)
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -s SEL, --select SEL  atom selection (default: all)
    -o FILE, --output FILE
                          output filename (default: trajectory.dcd)
    -n, --num             print the number of frames in each file and exit
    --psf PSF             PSF filename (must have same number of atoms as DCDs)
    --pdb PDB             PDB filename (must have same number of atoms as DCDs)
    --first INT           index of the first output frame, default: 0
    --last INT            index of the last output frame, default: -1
    --stride INT          number of steps between output frames, default: 1
    --align SEL           atom selection for aligning frames, a PSF or PDB file
                          must be provided, if a PDB is provided frames will be
                          superposed onto PDB coordinates

Examples
===============================================================================

Running :command:`prody catdcd --examples` displays::

  Concatenate two DCD files and output all atmos:
  
    $ prody catdcd mdm2.dcd mdm2sim2.dcd
  
  Concatenate two DCD files and output backbone atoms:
  
    $ prody catdcd mdm2.dcd mdm2sim2.dcd --pdb mdm2.pdb -s bb
  
