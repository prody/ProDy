.. _evolapps:

EVOL Application
===============================================================================

Synopsis
-------------------------------------------------------------------------------

The part shows how to:

  * Use evol applications to obtain conservation and coevolution plots without 
    the python API 

 
Evol Applications 
-------------------------------------------------------------------------------

Evol Applications have similar fuctionality as the python API. We can ``search``
Pfam, ``fetch`` from Pfam and also ``refine MSA``, ``merge`` two or more MSA 
and calculate ``conservation`` and ``coevolution`` properties and also 
``rankorder`` results from mutual information to get top-ranking pairs. 

All ``evol`` functions and their options can be obtained using the -h option. 
We should be in /prody/scripts directory to run the following commands::

    python evol -h
    python evol search -h
    python evol search 2W5IB
    python evol fetch PF00074

Using the above we can search and fetch MSA. Next we can refine the MSA:: 

    python evol refine -h
    python evol refine PF00074_full.slx -l RNAS1_BOVIN -s 0.98 -r 0.8

Next we can calculate conservation using shannon entropy and coevolution using
mutual information with correction and also save the plots.:: 

    python evol conserv PF00074_full_refined.slx -S
    python evol coevol PF00074_full_refined.slx -S -F png -c apc -cmin 0.0

We can rank order the residues with highest covariance and apply filters like
reporting only those pairs that are at a separation of at least 5 residues 
sequentially or are 15 Ang apart in structure. The residues may be numbered 
based on PDB:: 

    python evol rankorder -h
    python evol rankorder PF00074_full_refined_mutinfo_corr_apc.txt -q 5 -p 
    2W5IBI_1-121.pdb
    python evol rankorder PF00074_full_refined_mutinfo_corr_apc.txt -u -t 15 -p 2W5IB_1-121.pdb

   

