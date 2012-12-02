.. _evol-filter:

*******************************************************************************
evol filter
*******************************************************************************

Usage
===============================================================================

Running :command:`evol filter -h` displays::

  usage: evol filter [-h] [--quiet] [--examples] (-s | -e | -c) [-F] [-o STR]
                     [-f STR] [-z]
                     msa word [word ...]
  
  positional arguments:
    msa                   MSA filename to be filtered
    word                  word to be compared to sequence label
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  filtering method (required):
    -s, --startswith      sequence label starts with given words
    -e, --endswith        sequence label ends with given words
    -c, --contains        sequence label contains with given words
  
  filter option:
    -F, --full-label      compare full label with word(s)
  
  output options:
    -o STR, --outname STR
                          output filename, default is msa filename with _refined
                          suffix
    -f STR, --format STR  output MSA file format, default is same as input
    -z, --compressed      gzip refined MSA output

Examples
===============================================================================

Running :command:`evol filter --examples` displays::

  This application filters sequences in an MSA based on label data.
  Following example will filter human sequences:
  
      $ evol filter piwi_seed.slx HUMAN -e
