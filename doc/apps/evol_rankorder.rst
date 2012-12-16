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
  
