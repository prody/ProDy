.. _evol-refine:

*******************************************************************************
evol refine
*******************************************************************************

Usage
===============================================================================

Running :command:`evol refine -h` displays::

  usage: evol refine [-h] [--quiet] [--examples] [-l STR] [-s FLOAT] [-c FLOAT]
                     [-r FLOAT] [-k] [-o STR] [-f STR] [-z]
                     msa
  
  positional arguments:
    msa                   MSA filename to be refined
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  refinement options:
    -l STR, --label STR   sequence label, UniProt ID code, or PDB and chain
                          identifier
    -s FLOAT, --seqid FLOAT
                          identity threshold for selecting unique sequences
    -c FLOAT, --colocc FLOAT
                          column (residue position) occupancy
    -r FLOAT, --rowocc FLOAT
                          row (sequence) occupancy
    -k, --keep            keep columns corresponding to residues not resolved in
                          PDB structure, applies label argument is a PDB
                          identifier
  
  output options:
    -o STR, --outname STR
                          output filename, default is msa filename with _refined
                          suffix
    -f STR, --format STR  output MSA file format, default is same as input
    -z, --compressed      gzip refined MSA output

Examples
===============================================================================

Running :command:`evol refine --examples` displays::

  Sequence coevolution analysis involves several steps that including
  retrieving data and refining it for calculations.  These steps are
  illustrated below for RnaseA protein family.
  
  Search Pfam database:
  
    $  evol search 2w5i
  
  Download Pfam MSA file:
  
    $  evol fetch RnaseA
  
  Refine MSA file:
  
    $ evol refine RnaseA_full.slx -l RNAS1_BOVIN --seqid 0.98 --rowocc 0.8
  
  Checking occupancy:
  
    $ evol occupancy RnaseA_full.slx -l RNAS1_BOVIN -o col -S
  
  Conservation analysis:
  
    $ evol conserv RnaseA_full_refined.slx
  
  Coevolution analysis:
  
    $ evol coevol RnaseA_full_refined.slx -S -c apc
  
  Rank order analysis:
  
    $ evol rankorder RnaseA_full_refined_mutinfo_corr_apc.txt -p 2w5i_1-121.pdb --seq-sep 3
  
