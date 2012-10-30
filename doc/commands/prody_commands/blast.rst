.. _prody-blast:

*******************************************************************************
prody blast
*******************************************************************************

Usage
===============================================================================

Running :command:`prody blast -h` displays::

  usage: prody blast [-h] [--quiet] [--examples] [-i FLOAT] [-o FLOAT] [-d PATH]
                     [-z] [-f STR] [-e FLOAT] [-l INT] [-s INT] [-t INT]
                     sequence
  
  positional arguments:
    sequence              sequence or file in fasta format
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
    -i FLOAT, --identity FLOAT
                          percent sequence identity (default: 90.0)
    -o FLOAT, --overlap FLOAT
                          percent sequence overlap (default: 90.0)
    -d PATH, --output-dir PATH
                          download uncompressed PDB files to given directory
    -z, --gzip            write compressed PDB file
  
  Blast Parameters:
    -f STR, --filename STR
                          a filename to save the results in XML format
    -e FLOAT, --expect FLOAT
                          blast search parameter
    -l INT, --hit-list-size INT
                          blast search parameter
    -s INT, --sleep-time INT
                          how long to wait to reconnect for results (sleep time
                          is doubled when results are not ready)
    -t INT, --timeout INT
                          when to give up waiting for results

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
