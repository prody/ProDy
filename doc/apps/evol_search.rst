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
  
