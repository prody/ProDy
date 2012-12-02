.. _evol-search:

*******************************************************************************
evol search
*******************************************************************************

Usage
===============================================================================

Running :command:`evol search -h` displays::

  usage: evol search [-h] [--quiet] [--examples] [-b] [-s] [-g] [-e FLOAT]
                     [-t INT] [-o STR] [-d STR]
                     query
  
  positional arguments:
    query                 protein UniProt ID or sequence, a PDB identifier, or a
                          sequence file, where sequence have no gaps and 12 or
                          more characters
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  sequence search options:
    -b, --searchBs        search Pfam-B families
    -s, --skipAs          do not search Pfam-A families
    -g, --ga              use gathering threshold
    -e FLOAT, --evalue FLOAT
                          e-value cutoff, must be less than 10.0
    -t INT, --timeout INT
                          timeout in seconds for blocking connection attempt
                          (default: 30)
  
  output options:
    -o STR, --outname STR
                          name for output file, default is standard output
    -d STR, --delimiter STR
                          delimiter for output data columns (default: )

Examples
===============================================================================

Running :command:`evol search --examples` displays::

  This application searches Pfam database with a UniProt ID or protein
  sequence or filename containing the sequence.  Some specific search
  options are included for sequence search.  Minimum length of query
  sequence should be 12 and should not contain gaps.  If outname is
  specified it will output the results obtained in a file or the output
  will be directed to standard output.
  
  Search Pfam with PDB and chain identifier and output results to
  screen:
  
      $ evol search 1mkpA
  
  Search Pfam with UniProt ID and write output into a file:
  
      $ evol search P08581 --outname families.txt
  
  Search Pfam with a sequence with search options:
  
      $ evol search PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCALCSLHSIGKIGGAQNRSYSKLLCGLLAERLRISPDRVYINYYDMNAANVGWNNSTFA --evalue 2 --searchBs
  
