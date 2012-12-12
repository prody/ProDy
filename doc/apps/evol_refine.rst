.. _evol-refine:

*******************************************************************************
evol refine
*******************************************************************************

Usage
===============================================================================

Running :command:`evol refine -h` displays::

  usage: evol refine [-h] [--quiet] [--examples] [-l STR] [-s FLOAT] [-c FLOAT]
                     [-r FLOAT] [-o STR] [-f STR] [-z]
                     msa
  
  positional arguments:
    msa                   MSA filename to be refined
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  refinement options:
    -l STR, --label STR   sequence label, UniProt ID code or PDB and chain
                          identifier
    -s FLOAT, --seqid FLOAT
                          identity threshold for selecting unique sequences
    -c FLOAT, --colocc FLOAT
                          column (residue position) occupancy
    -r FLOAT, --rowocc FLOAT
                          row (sequence) occupancy
  
  output options:
    -o STR, --outname STR
                          output filename, default is msa filename with _refined
                          suffix
    -f STR, --format STR  output MSA file format, default is same as input
    -z, --compressed      gzip refined MSA output

Examples
===============================================================================

Running :command:`evol refine --examples` displays::

  This application refines MSA by removing gapped columns (residue
  positions) and rows (sequences).  Following example will save entropy
  data and plot using default options:
  
      $ evol refine piwi.slx -l GTHB2_ONCKE
