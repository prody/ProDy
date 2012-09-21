.. _prody-pca:

*******************************************************************************
prody pca
*******************************************************************************

Usage
===============================================================================

Running :command:`prody pca -h` displays::

  usage: prody pca [-h] [--quiet] [--examples] [-n INT] [-s SEL] [-a] [-o PATH]
                   [-e] [-r] [-q] [-v] [-z] [-t STR] [-j] [-p STR] [-f STR]
                   [-d STR] [-x STR] [-A] [-R] [-Q] [-J STR] [-F STR] [-D INT]
                   [-W FLOAT] [-H FLOAT] [--psf PSF | --pdb PDB] [--aligned]
                   dcd
  
  positional arguments:
    dcd                   file in DCD or PDB format
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    --psf PSF             PSF filename
    --pdb PDB             PDB filename
    --aligned             trajectory is already aligned
  
  parameters:
    -n INT, --number-of-modes INT
                          number of non-zero eigenvectors (modes) to calculate
                          (default: 10)
    -s SEL, --select SEL  atom selection (default: "protein and name CA or
                          nucleic and name P C4' C2")
  
  output:
    -a, --all-output      write all outputs
    -o PATH, --output-dir PATH
                          output directory (default: ".")
    -e, --eigenvs         write eigenvalues/vectors
    -r, --cross-correlations
                          write cross-correlations
    -q, --square-fluctuations
                          write square-fluctuations
    -v, --covariance      write covariance matrix
    -z, --npz             write compressed ProDy data file
    -t STR, --extend STR  write NMD file for the model extended to "backbone"
                          ("bb") or "all" atoms of the residue, model must have
                          one node per residue
    -j, --projection      write projections onto PCs
  
  output options:
    -p STR, --file-prefix STR
                          output file prefix (default: "pdb_pca")
    -f STR, --number-format STR
                          number output format (default: "%12g")
    -d STR, --delimiter STR
                          number delimiter (default: " ")
    -x STR, --extension STR
                          numeric file extension (default:".txt")
  
  figures:
    -A, --all-figures     save all figures
    -R, --cross-correlations-figure
                          save cross-correlations figure
    -Q, --square-fluctuations-figure
                          save square-fluctuations figure
    -J STR, --projection-figure STR
                          save projections onto specified subspaces, e.g. "1,2"
                          for projections onto PCs 1 and 2; "1,2 1,3" for
                          projections onto PCs 1,2 and 1, 3; "1 1,2,3" for
                          projections onto PCs 1 and 1, 2, 3
  
  figure options:
    -F STR, --figure-format STR
                          pdf (default: pdf)
    -D INT, --dpi INT     figure resolution (dpi) (default: 300)
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
