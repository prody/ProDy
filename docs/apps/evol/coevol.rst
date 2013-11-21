.. _evol-coevol:

evol coevol
====================

Usage
--------------------

Running :command:`evol coevol -h` displays::

  usage: evol coevol [-h] [--quiet] [--examples] [-n] [-c STR] [-m STR] [-t]
                     [-p STR] [-f STR] [-S] [-L FLOAT] [-U FLOAT] [-X STR]
                     [-T STR] [-D INT] [-H FLOAT] [-W FLOAT] [-F STR]
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
    -c STR, --correction STR
                          also save corrected mutual information matrix data and
                          plot, one of apc, asc
    -m STR, --normalization STR
                          also save normalized mutual information matrix data
                          and plot, one of sument, minent, maxent, mincon,
                          maxcon, joint
  
  output options:
    -t, --heatmap         save heatmap files for all mutual information matrices
    -p STR, --prefix STR  output filename prefix, default is msa filename with
                          _coevol suffix
    -f STR, --number-format STR
                          number output format (default: %12g)
  
  figure options:
    -S, --save-plot       save coevolution plot
    -L FLOAT, --cmin FLOAT
                          apply lower limits for figure plot
    -U FLOAT, --cmax FLOAT
                          apply upper limits for figure plot
    -X STR, --xlabel STR  specify xlabel, by default will be applied on ylabel
    -T STR, --title STR   figure title
    -D INT, --dpi INT     figure resolution (dpi) (default: 300)
    -H FLOAT, --height FLOAT
                          figure height (inch) (default: 6)
    -W FLOAT, --width FLOAT
                          figure width (inch) (default: 8)
    -F STR, --figure-format STR
                          figure file format, one of svgz, rgba, png, pdf, eps,
                          svg, ps, raw (default: pdf)

Examples
--------------------

Running :command:`evol coevol --examples` displays::

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
  
