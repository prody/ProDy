.. _evol-conserv:

evol conserv
====================

Usage
--------------------

Running :command:`evol conserv -h` displays::

  usage: evol conserv [-h] [--quiet] [--examples] [-n] [-g] [-p STR] [-f STR]
                      [-S] [-H FLOAT] [-W FLOAT] [-F STR] [-D INT]
                      msa
  
  positional arguments:
    msa                   refined MSA file
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  calculation options:
    -n, --no-ambiguity    treat amino acids characters B, Z, J, and X as non-
                          ambiguous
    -g, --gaps            do not omit gap characters
  
  output options:
    -p STR, --prefix STR  output filename prefix, default is msa filename with
                          _conserv suffix
    -f STR, --number-format STR
                          number output format (default: %12g)
  
  figure options:
    -S, --save-plot       save conservation plot
    -H FLOAT, --height FLOAT
                          figure height (inch) (default: 6)
    -W FLOAT, --width FLOAT
                          figure width (inch) (default: 8)
    -F STR, --figure-format STR
                          figure file format, one of raw, png, ps, svgz, eps,
                          pdf, rgba, svg (default: pdf)
    -D INT, --dpi INT     figure resolution (dpi) (default: 300)

Examples
--------------------

Running :command:`evol conserv --examples` displays::

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
  
