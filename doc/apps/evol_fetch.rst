.. _evol-fetch:

*******************************************************************************
evol fetch
*******************************************************************************

Usage
===============================================================================

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
                          timeout for blocking connection attempts (default: 5)
  
  output options:
    -d PATH, --outdir PATH
                          output directory (default: .)
    -p STR, --outname STR
                          output filename, default is accession and alignment
                          type
    -z, --compressed      gzip downloaded MSA file

Examples
===============================================================================

Running :command:`evol fetch --examples` displays::

  Given a Pfam ID or accession, this program fetches the MSA of that
  family. Supported alignment options are full, seed, ncbi or
  metagenomics and alignment formats are selex, stockholm or fasta. The
  output MSA is downloaded and saved in the specified or default '.'
  directory.
  
  Fetch PFAM ID Cys_knot:
  
      $ evol fetch Cys_knot
  
  Fetch PFAM accession with specific parameters:
  
      $ evol fetch PF00007 --compressed --format fasta --outname mymsa
