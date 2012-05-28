.. _prody-pca:

*******************************************************************************
prody pca
*******************************************************************************

Perform PCA/EDA calculations and output the results in plain text, NMD, and 
graphical formats.

Download example :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`.

Usage
===============================================================================

Running :command:`prody pca -h` displays::

  usage: prody pca [-h] [--quiet] [--examples] [-n INT] [-s SELSTR] [-a]
                   [-o PATH] [-e] [-r] [-q] [-v] [-z] [-t STR] [-j] [-p STR]
                   [-f STR] [-d STR] [-x STR] [-A] [-R] [-Q] [-J STR] [-F STR]
                   [-D INT] [-W FLOAT] [-H FLOAT] [--psf PSF | --pdb PDB]
                   [--aligned]
                   coords
  
  positional arguments:
    coords                PDB or DCD filename
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    --psf PSF             PSF filename
    --pdb PDB             PDB filename
    --aligned             trajectory is already aligned
  
  parameters:
    -n INT, --number-of-modes INT
                          number of non-zero eigenvalues/vectors to calculate
                          (default: 10)
    -s SELSTR, --select SELSTR
                          atom selection (default: "protein and name CA or
                          nucleic and name P C4' C2")
  
  output:
    -a, --all-output      write all outputs
    -o PATH, --output-dir PATH
                          output directory (default: ".")
    -e, --eigenvs         write eigenvectors/values
    -r, --cross-correlations
                          write cross-correlations
    -q, --square-fluctuations
                          write square-fluctuations
    -v, --covariance      write covariance matrix
    -z, --npz             write compressed ProDy data file
    -t STR, --extend STR  output NMD files for the model extended to "backbone"
                          ("bb") or "all" atoms of the residue, model must have
                          one node per residue
    -j, --projection      write projections onto PCs
  
  output options:
    -p STR, --file-prefix STR
                          prefix for output files (default: pdb_pca)
    -f STR, --number-format STR
                          delimiter (default: "%12g")
    -d STR, --delimiter STR
                          delimiter (default: " ")
    -x STR, --extension STR
                          file extension (default: .txt)
  
  figures:
    -A, --all-figures     save all figures
    -R, --cross-correlations-figure
                          save cross-correlations
    -Q, --square-fluctuations-figure
                          save square-fluctuations
    -J STR, --projection-figure STR
                          save projections onto specified subspaces, e.g. "1,2"
                          for projections onto PCs 1 and 2; "1,2 1,3" for
                          projections onto PCs 1,2 and 1, 3; "1 1,2,3" for
                          projections onto PCs 1 and 1, 2, 3
  
  figure options:
    -F STR, --figure-format STR
                          figure format, one of eps, pdf, png, ps, raw, rgba,
                          svg, svgz (default: pdf)
    -D INT, --resolution INT
                          figure resolution (dpi) (default: 300)
    -W FLOAT, --width FLOAT
                          figure width (inch) (default: 8.0)
    -H FLOAT, --height FLOAT
                          figure height (inch) (default: 6.0)

Examples
===============================================================================

Running :command:`prody pca --examples` displays::

  This command performs PCA (or EDA) calculations for given multi-model
  PDB structure or DCD format trajectory file and outputs results in NMD
  format.  If a PDB identifier is given, structure file will be
  downloaded from the PDB FTP server.  DCD files may be accompanied with
  PDB or PSF files to enable atoms selections.
  
  Fetch pdb 2k39, perform PCA calculations, and output NMD file:
  
    $ prody pca 2k39
  
  Fetch pdb 2k39 and perform calculations for backbone of residues up to
  71, and save all output and figure files:
  
    $ prody pca 2k39 --select "backbone and resnum < 71" -a -A
  
  Perform EDA of MDM2 trajectory:
  
    $ prody eda mdm2.dcd
  
  Perform EDA for backbone atoms:
  
    $ prody eda mdm2.dcd --pdb mdm2.pdb --select backbone
