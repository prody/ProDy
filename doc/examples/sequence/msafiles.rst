.. _msafiles:

*******************************************************************************
Pfam Database and MSA Files
*******************************************************************************

Synopsis
===============================================================================

The following examples shows how to do the following:

  * search Pfam database to identify pfam accession no and other details
  * fetch the MSA of the pfam using acession no 
  * parse the MSA, filter, slice MSA, write the MSA


Search Pfam 
===============================================================================

This example demonstrates how to search Pfam database with a given query, 
:func:`.searchPfam`.  Valid inputs are UniProt ID, e.g. ``"PIWI_ARCFU"``, or 
PDB file, e.g. ``"3luc"`` or ``"3lucA"`` with chain identifier. Input can also 
be a protein sequence or a file containing the sequence, but sequence should 
not contain gaps and should be at least 12 characters long.

Matching Pfam accession (one or more) as keys will map to a dictionary that 
contains locations (alignment start, end, evalue etc), pfam family type, 
accession and id.
 
  
UniProt ID search
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

>>> from prody import *

The function will return a dictionary if successful.
 
>>> matches = searchPfam('PIWI_ARCFU')
>>> for family, in matches: # doctest: +SKIP
...     print family, matches[family]  
  PF02171 {'locations': [{'end': '406', 'bitscore': '215.70', 'hmm_start': '2', 
 'ali_end': '405', 'ali_start': '111', 'start': '110', 'evalue': '6.8e-61', 
 'hmm_end': '304'}], 'type': 'Pfam-A', 'accession': 'PF02171', 'id': 'Piwi'}


Sequence search
-------------------------------------------------------------------------------

This function also accepts a protein sequence:

>>> sequence = ('PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGS'
... 'SEPCALCSLHSIGKIGGAQNRSYSKLLCGLLAERLRISPDRVYINYYDMNAANVGWNNSTFA')
>>> matches = searchPfam(sequence) # doctest: +SKIP


For sequence searches, we can pass additional parameters to :func:`searchPfam`
like *search_b* which will search pfam B and *skip_a* that will not search 
pfamA database. Additional parameters include *ga* that uses gathering 
threshold instead of evalue,  *evalue* cutoff can also be specified and 
*timeout* that can be set higher especially when searching larger 
sequences, default is ``timeout=30`` seconds.

>>> matches = searchPfam(sequence, search_b=True, evalue=2.0) # doctest: +SKIP


Retrieve MSA files
===============================================================================

This example demonstrates how to search Pfam database with a given query using  
:func:`.fetchPfamMSA`. Valid inputs are Pfam ID, e.g. ``"Piwi"``, or Pfam
accession, e.g. ``"PF02171"`` obtained from :func:`.searchPfam`.  Alignment 
type can be ``'full'`` (default), ``"seed"``, ``"ncbi"`` or ``"metagenomics"``.

>>> fetchPfamMSA('piwi', alignment='seed')
'piwi_seed.sth'

A compressed file can be downloaded by setting ``compressed=True``. 
The ``format`` of the MSA can be of ``"selex"``, (default), ``"stockholm"`` or
``"fasta"``.  This will return the path of the downloaded MSA file. 
The ``output`` name can be specified, for by default it will have 
``"accession/ID_alignment.format"``.

Note that in this case we passed a folder name, the downloaded file is saved 
in this folder, after it is created if it did not exist. Also bigger timeouts
are necessary for larger families. Some other parameters like ``gap``, 
``order`` or ``inserts`` can be set, as shown in the following example. 

>>> fetchPfamMSA('PF02171', compressed=True, gaps='mixed', inserts='lower', 
... order='alphabetical', format='fasta')
'PF02171_full.fasta.gz'
    
>>> msafile = 'piwi_seed.sth'


Parsing MSA files
===============================================================================

This shows how to use the :class:`.MSAFile` or :func:`.parseMSA` to read the 
MSA file. :func:`.parseMSA` returns a :class:`.MSA` object. 

Reading using :class:`.MSAFile` yields an MSAFile object. Iterating over the 
object will yield sequence id, sequence, residue start and end indices:

>>> msafobj = MSAFile(msafile)
>>> for seq in msafobj: # doctest: +ELLIPSIS 
...     print(seq)
('YQ53_CAEEL', 'DILVGIAR.EKKP...NLAKRGRNNYK', 650, 977)
('Q21691_CAEEL', 'TIVFGIIA.EKRP...NLAKRGHNNYK', 673, 1001)
('AGO6_ARATH', 'FILCILPERKTSD...LAAAQVAQFTK', 541, 851)
...
('O02095_CAEEL', 'QLLFFVVK..SRY...RYSQRGAMVLA', 574, 878)
('Q19645_CAEEL', 'PFVLFISD..DVP...ELAKRGTGLYK', 674, 996)
('O62275_CAEEL', 'TFVFIITD.DSIT...EYAKRGRNLWN', 594, 924)

Reading using :func:`.parseMSA` yields an :class:`.MSA` object.  We can parse 
compressed files, but reading uncompressed files are much faster as shown.

>>> msa = parseMSA('PF02171_full.fasta.gz')

>>> fetchPfamMSA('PF02171', format='fasta')
'PF02171_full.fasta'
>>> msa = parseMSA('PF02171_full.fasta')


Filtering and Slicing
===============================================================================

This shows how to use the :class:`.MSAFile` object or :class:`.MSA` object to 
refine MSA using filters and slices. 

*Filtering sequences*
    
Any function that takes label and sequence arguments and returns a boolean 
value can be used for filtering the sequences.  A sequence will be yielded 
if the function returns **True**.  In the following example, sequences from
organism *ARATH* are filtered:
    
>>> msafobj = MSAFile(msafile, filter=lambda lbl, seq: 'ARATH' in lbl)
>>> for seq in msafobj: # doctest: +ELLIPSIS 
...     print(seq)
('AGO6_ARATH', 'FIL...FTK', 541, 851)
('AGO4_ARATH', 'FIL...FMK', 577, 885)
('AGO10_ARATH', 'LLL...YLE', 625, 946)

*Slicing sequences*
    
A list of integers can be used to slice sequences as follows.  This enables 
selective parsing of the MSA file. 
    
>>> msafobj = MSAFile(msafile, slice=list(range(10)) + list(range(394,404)))
>>> for seq in msafobj: # doctest: +ELLIPSIS 
...     print(seq)
('YQ53_CAEEL', 'DILVGIAR.ELAKRGRNNYK', 650, 977)
('Q21691_CAEEL', 'TIVFGIIA.ELAKRGHNNYK', 673, 1001)
('AGO6_ARATH', 'FILCILPERKAAAQVAQFTK', 541, 851)
(...)
('O02095_CAEEL', 'QLLFFVVK..YSQRGAMVLA', 574, 878)
('Q19645_CAEEL', 'PFVLFISD..LAKRGTGLYK', 674, 996)
('O62275_CAEEL', 'TFVFIITD.DYAKRGRNLWN', 594, 924)

Slicing can also be done using :class:`.MSA`. The :class:`.MSA`. object offers 
other functionalities like querying, indexing, slicing row and columns and 
refinement. 

*Querying*
    
You can query whether a sequence in contained in the instance using the 
UniProt identifier of the sequence as follows:


>>> msa = parseMSA(msafile)
>>> 'YQ53_CAEEL' in msa
True
 
    
Indexing MSA objects
===============================================================================
    
Retrieve a sequence at a given index:
    
>>> msa[0] # doctest: +ELLIPSIS
<Sequence: YQ53_CAEEL (piwi_seed[0]; length 404; 328 residues and 76 gaps)>
    
Retrieve a sequence by UniProt ID:
    
>>> msa['YQ53_CAEEL'] # doctest: +ELLIPSIS
<Sequence: YQ53_CAEEL (piwi_seed[0]; length 404; 328 residues and 76 gaps)>
    
Slice an MSA instance:
    
>>> msa[:2]
<MSA: piwi_seed' (2 sequences, 404 residues)>
    
Slice using a list of UniProt IDs:
    
>>> msa[:2] == msa[['YQ53_CAEEL', 'Q21691_CAEEL']]
True
    
Retrieve a character or a slice of a sequence:

>>> msa[0,0]
<Sequence: YQ53_CAEEL (length 1; 1 residues and 0 gaps)>
>>> msa[0,0:10]
<Sequence: YQ53_CAEEL (length 10; 9 residues and 1 gaps)>
    
Slice MSA rows and columns:
    
>>> msa[:10,20:40]
<MSA: piwi_seed' (10 sequences, 20 residues)>

    
Writing MSA files
===============================================================================

:func:`.writeMSA` can be used to write MSA. It takes filename as input 
which should contain appropriate extension that can be ``".slx"`` or 
``".sth"`` or  `".fasta"`` or format should be specified as ``"SELEX"``, 
``"FASTA"`` or ``"Stockholm"``. Input MSA should be :class:`.MSAFile` or 
:class:`.MSA` object. Filename can contain ``".gz"`` extension, in which case 
a compressed file will be written. 
Returns the name of the MSA file that is written. 

>>> writeMSA('sliced_MSA.gz', msa, format='SELEX')
'sliced_MSA.gz'
>>> filename = writeMSA('sliced_MSA.fasta', msafobj)


See Also
===============================================================================

See :mod:`~.prody.sequence` module for all sequence analysis functions. 

|questions|

|suggestions|

.. sectionauthor:: Anindita Dutta
