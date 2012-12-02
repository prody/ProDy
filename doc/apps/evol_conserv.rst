.. _evol-conserv:

*******************************************************************************
evol conserv
*******************************************************************************

Usage
===============================================================================

Running :command:`evol conserv -h` displays::

  usage: evol conserv [-h] [--quiet] [--examples] [-n] [-g] [-p STR] [-f STR]
                      [-S] [-F STR] [-D INT] [-W FLOAT] [-H FLOAT]
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

Running :command:`evol conserv --examples` displays::

  This application calculates conservation using Shannon entropy for a \
  refined multiple sequence alignment.  Following example will save
  entropy data and plot using default options:
  
      $ evol conserv piwi_refined.slx -S
