.. _evol-fetch:

evol fetch
====================

Usage
--------------------

Running :command:`evol fetch -h` displays::

  usage: evol fetch [-h] [--quiet] [--examples] [-a STR] [-f STR] [-o STR]
                    [-i STR] [-g STR] [-t INT] [-d PATH] [-p STR] [-z]
                    acc
  
  positional arguments:
    acc                   Pfam accession or ID
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  download options:
    -a STR, --alignment STR
                          alignment type, one of full, seed, ncbi, metagenomics
                          (default: full)
    -f STR, --format STR  Pfam supported MSA format, one of selex, fasta,
                          stockholm (default: selex)
    -o STR, --order STR   ordering of sequences, one of tree, alphabetical
                          (default: tree)
    -i STR, --inserts STR
                          letter case for inserts, one of upper, lower (default:
                          upper)
    -g STR, --gaps STR    gap character, one of dashes, dots, mixed (default:
                          dashes)
    -t INT, --timeout INT
                          timeout for blocking connection attempts (default: 60)
  
  output options:
    -d PATH, --outdir PATH
                          output directory (default: .)
    -p STR, --outname STR
                          output filename, default is accession and alignment
                          type
    -z, --compressed      gzip downloaded MSA file

Examples
--------------------

Running :command:`evol fetch --examples` displays::

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
  
