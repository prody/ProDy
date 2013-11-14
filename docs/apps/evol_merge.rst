.. _evol-merge:

evol merge
====================

Usage
--------------------

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
--------------------

Running :command:`evol merge --examples` displays::

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
  
