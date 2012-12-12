.. _evol-merge:

*******************************************************************************
evol merge
*******************************************************************************

Usage
===============================================================================

Running :command:`evol merge -h` displays::

  usage: evol merge [-h] [--quiet] [--examples] [-o STR] [-f STR] [-z]
                    msa [msa ...]
  
  positional arguments:
    msa                   MSA filenames to be merged
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  output options:
    -o STR, --outname STR
                          output filename, default is first input filename with
                          _merged suffix
    -f STR, --format STR  output MSA file format, default is same as first input
                          MSA
    -z, --compressed      gzip merged MSA output

Examples
===============================================================================

Running :command:`evol merge --examples` displays::

  This application merges a number of MSAs into one large MSA. The
  merging of
  seqences in done based on common labels appearing across all the input
  MSAs The
  following example show how to merge two MSAs:
  
      $ evol merge piwi.slx -l GTHB2_ONCKE
