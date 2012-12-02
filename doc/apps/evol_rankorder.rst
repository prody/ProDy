.. _evol-rankorder:

*******************************************************************************
evol rankorder
*******************************************************************************

Usage
===============================================================================

Running :command:`evol rankorder -h` displays::

  usage: evol rankorder [-h] [--quiet] [--examples] [-z] [-d STR] [-p STR]
                        [-m STR] [-l STR] [-n INT] [-q INT] [-t FLOAT] [-u]
                        [-o STR]
                        mutinfo
  
  positional arguments:
    mutinfo               mutual information matrix
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  input options:
    -z, --zscore          apply zscore for identifying top ranked coevolving
                          pairs
    -d STR, --delimiter STR
                          delimiter used in mutual information matrix file
    -p STR, --pdb STR     PDB file that contains same number of residues as the
                          mutual information matrix, output residue numbers will
                          be based on PDB file
    -m STR, --msa STR     MSA file used for building the mutual info matrix,
                          output residue numbers will be based on the most
                          complete sequence in MSA if a PDB file or sequence
                          label is not specified
    -l STR, --label STR   label in MSA file for output residue numbers
  
  output options:
    -n INT, --num-pairs INT
                          number of top ranking residue pairs to list (default:
                          100)
    -q INT, --seq-sep INT
                          report coevolution for residue pairs that are
                          sequentially separated by input value (default: 3)
    -t FLOAT, --min-dist FLOAT
                          report coevolution for residue pairs whose CA atoms
                          are spatially separated by at least the input value,
                          used when a PDB file is given and --use-dist is true
                          (default: 10.0)
    -u, --use-dist        use structural separation to report coevolving pairs
    -o STR, --outname STR
                          output filename, default is mutinfo_rankorder.txt

Examples
===============================================================================

Running :command:`evol rankorder --examples` displays::

  This application identifies that top ranking pairs of residues that
  coevolve based on their mutual information.  By default coevolution is
  reported for pairs that are at least 3 residues apart in sequence. A
  z-score normalization can be applied to the mutinfo matrix to identify
  coevolving pairs.  The following examples show how to use with default
  as well as additional options:
  
      $ evol rankorder piwi_refined_mutinfo.txt -z
  
      $ evol rankorder piwi_refined_mutinfo.txt --msa piwi_refined.slx --label AGO6_ARATH
