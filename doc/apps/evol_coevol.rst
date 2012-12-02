.. _evol-coevol:

*******************************************************************************
evol coevol
*******************************************************************************

Usage
===============================================================================

Running :command:`evol coevol -h` displays::

  usage: evol coevol [-h] [--quiet] [--examples] [-n] [-g] [-c STR] [-m STR]
                     [-p STR] [-f STR] [-S] [-L FLOAT] [-U FLOAT] [-F STR]
                     [-D INT] [-W FLOAT] [-H FLOAT]
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
    -c STR, --correction STR
                          also save corrected mutual information matrix data and
                          plot, one of apc, asc
    -m STR, --normalization STR
                          also save normalized mutual information matrix data
                          and plot, one of sument, minent, maxent, mincon,
                          maxcon, joint
  
  output options:
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
    -F STR, --figure-format STR
                          figure file format, one of svgz, ps, svg, eps, raw,
                          rgba, pdf, png (default: pdf)
    -D INT, --dpi INT     figure resolution (dpi) (default: 300)
    -W FLOAT, --width FLOAT
                          figure width (inch) (default: 8)
    -H FLOAT, --height FLOAT
                          figure height (inch) (default: 6)

Examples
===============================================================================

Running :command:`evol coevol --examples` displays::

  This application calculates mutual information between MSA postions
  for a refined multiple sequence alignment.  Following example will
  save coevolution data and plot using default options:
  
      $ evol coevol piwi_refined.slx -S
  
  Following example will save coevolution data and plot for all
  correction and normalizations:
  
      $ evol coevol piwi_refined.slx -S -c apc -c asc -m sument -m minent -m maxent -m mincon -m maxcon -m joint
