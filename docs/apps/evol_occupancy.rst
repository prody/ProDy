.. _evol-occupancy:

evol occupancy
====================

Usage
--------------------

Running :command:`evol occupancy -h` displays::

  usage: evol occupancy [-h] [--quiet] [--examples] [-o STR] [-p STR] [-l STR]
                        [-f STR] [-S] [-X STR] [-Y STR] [-T STR] [-D INT]
                        [-W FLOAT] [-F STR] [-H FLOAT]
                        msa
  
  positional arguments:
    msa                   MSA file
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  calculation options:
    -o STR, --occ-axis STR
                          calculate row or column occupancy or both., one of
                          row, col, both (default: row)
  
  output options:
    -p STR, --prefix STR  output filename prefix, default is msa filename with
                          _occupancy suffix
    -l STR, --label STR   index for column based on msa label
    -f STR, --number-format STR
                          number output format (default: %12g)
  
  figure options:
    -S, --save-plot       save occupancy plot/s
    -X STR, --xlabel STR  specify xlabel
    -Y STR, --ylabel STR  specify ylabel
    -T STR, --title STR   figure title
    -D INT, --dpi INT     figure resolution (dpi) (default: 300)
    -W FLOAT, --width FLOAT
                          figure width (inch) (default: 8)
    -F STR, --figure-format STR
                          figure file format, one of svgz, svg, pdf, raw, rgba,
                          ps, eps, png (default: pdf)
    -H FLOAT, --height FLOAT
                          figure height (inch) (default: 6)

Examples
--------------------

Running :command:`evol occupancy --examples` displays::

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
  
