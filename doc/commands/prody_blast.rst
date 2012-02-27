.. _prody-blast:

*******************************************************************************
prody blast
*******************************************************************************

Blast Protein Data Bank for structures matching a user given sequence.

Usage
===============================================================================

Running :command:`prody blast -h` displays::

  usage: prody blast [-h] [--quiet] [--examples] [-i FLOAT] [-o FLOAT] [-d PATH]
                     [-z]
                     seq
  
  positional arguments:
    seq                   sequence or file in fasta format
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -i FLOAT, --identity FLOAT
                          percent sequence identity (default: 90.0)
    -o FLOAT, --overlap FLOAT
                          percent sequence overlap (default: 90.0)
    -d PATH, --dir PATH   download uncompressed PDB files to given path
    -z, --gzip            write compressed PDB file

Examples
===============================================================================

Running :command:`prody blast --examples` displays::

  Blast search PDB for the first sequence in a fasta file:
  
    $ prody blast seq.fasta -i 70
  
  Blast search PDB for the sequence argument:
  
    $ prody blast MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGR
  TLSDYNIQKESTLHLVLRLRGG
  
  Blast search PDB for avidin structures, download files, and align all
  files onto the 2avi structure:
  
    $ prody blast -d . ARKCSLTGKWTNDLGSNMTIGAVNSRGEFTGTYITAVTATSNEIKESPL
  HGTQNTINKRTQPTFGFTVNWKFSESTTVFT
    $ prody align 2avi.pdb *pdb
