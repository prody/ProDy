.. _msafiles:

*******************************************************************************
Accessing Pfam data and handling MSA files
*******************************************************************************

Synopsis
===============================================================================

The following examples shows how to do the following:

  * search Pfam database to identify pfam accession no and other details
  * fetch the MSA of the pfam using acession no 
  * parse the MSA, filter, slice MSA, write the MSA

search Pfam database
===============================================================================

This example demonstrates how to search Pfam database with a given query, 
:func:`.searchPfam`.  Valid inputs are UniProt ID, e.g. ``"P69332"``, or 
PDB file, e.g. ``"1mkp"`` or ``"1mkpA"`` with chain identifier. Input can also 
be a protein sequence or a file containing the sequence, but sequence should 
not contain gaps and should be at least 12 characters long.

Matching Pfam accession (one or more) as keys will map to a dictionary that 
contains locations (alignment start, end, evalue etc), pfam family type, 
accession and id.
  
UniProt ID 
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

The function will return a dictionary if successful.
 
>>> matches = searchPfam('P69332')
>>> for family, match in matches: # doctest: +SKIP
...     print family, match  
 PF00001 {'locations': [{'end': '291', 'bitscore': '161.50', 'hmm_start': '1', 
'ali_end': '291', 'ali_start': '50', 'start': '50', 'evalue': '2.1e-44', 
'hmm_end': '257'}], 'type': 'Pfam-A', 'accession': 'PF00001', 'id': '7tm_1'}

Sequence and Parameters
===============================================================================

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

fetch MSA from Pfam
===============================================================================

This example demonstrates how to search Pfam database with a given query using  
:func:`.fetchPfamMSA`. Valid inputs are Pfam ID, e.g. ``"7tm_1"``, or Pfam
accession, e.g. ``"PF00001"`` obtained from :func:`.searchPfam`.  Alignment 
type can be ``'full'`` (default), ``"seed"``, ``"ncbi"`` or ``"metagenomics"``.

A compressed file can be downloaded by setting ``compressed=True``. 
The ``format`` of the MSA can be of ``"selex"``, (default), ``"stockholm"`` or
``"fasta"``.  This will return the path of the downloaded MSA file. 
The ``output`` name can be specified, for by default it will have 
``"accession/ID_alignment.format"``.

>>> fetchPfamMSA('PF00001', compressed=True, format='fasta', timeout=60)
'PF00001_full.fasta.gz'

Note that in this case we passed a folder name, the downloaded file is saved 
in this folder, after it is created if it did not exist. Also bigger timeouts
are necessary for larger families. Some other parameters like ``gap``, 
``order`` or ``inserts`` can be set, as shown in the following example. 

>>> fetchPfamMSA('piwi', 'seed', gaps='mixed', inserts='lower', 
... order='alphabetical')
'piwi_seed.slx'
    
>>> msafile = 'piwi_seed.slx'

Parse, Modify and Write MSAs
===============================================================================

These examples show how to use :class:`.MSAFile` object and 
:class:`.MSA` object to parse, refine the MSA and write the MSA. 

Parse MSAs
===============================================================================

This shows how to use the :class:`.MSAFile` or :func:`.parseMSA` to read the 
MSA file. 

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

>>> msa = parseMSA('PF00001_full.fasta.gz')

>>> msa = parseMSA('PF00001_full.fasta')

Filter or Slice MSA
===============================================================================

This shows how to use the :class:`.MSAFile` object or :class:`.MSA` object to 
refine MSA using filters and slices. 

*Filtering sequences*
    
Any function that takes label and sequence arguments and returns a boolean 
value can be used for filtering the sequences.  A sequence will be yielded 
if the function returns **True**.  In the following example, sequences from
organism *ARATH* are filtered:
    
>>> msa = MSAFile(msafile, filter=lambda lbl, seq: 'ARATH' in lbl)
>>> for seq in msa: # doctest: +ELLIPSIS 
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
    
*Indexing and slicing*
    
Retrieve a sequence at a given index:
    
>>> msa[0] # doctest: +ELLIPSIS
('YQ53_CAEEL', 'DIL...YK', 650, 977)
    
Retrieve a sequence by UniProt ID:
    
>>> msa['YQ53_CAEEL'] # doctest: +ELLIPSIS
('YQ53_CAEEL', 'DIL...YK', 650, 977)
    
Slice an MSA instance:
    
>>> msa[:2]
<MSA: piwi_seed' (2 sequences, 404 residues)>
    
Slice using a list of UniProt IDs:
    
>>> msa[:2] == msa[['YQ53_CAEEL', 'Q21691_CAEEL']]
True
    
Retrieve a character or a slice of a sequence:

>>> msa[0,0]
'D'
>>> msa[0,0:10]
'DILVGIAR.E'
    
Slice MSA rows and columns:
    
>>> msa[:10,20:40]
<MSA: piwi_seed' (10 sequences, 20 residues)>
    
write MSA
===============================================================================

:func:`.writeMSA` can be used to write MSA. It takes filename as input 
which should contain appropriate extension that can be ``".slx"`` or 
``".sth"`` or  `".fasta"`` or format should be specified as ``"SELEX"``, 
``"FASTA"`` or ``"Stockholm"``. Input MSA should be :class:`.MSAFile` or 
:class:`.MSA` object. Filename can contain ``".gz"`` extension, in which case 
a compressed file will be written. 
Returns the name of the MSA file that is written. 

>>> writeMSA('sliced_MSA.gz', msa, format='SELEX')
    'test.gz'
>>> filename = writeMSA('sliced_MSA.fasta', msafobj)


See Also
===============================================================================

See :mod:`.sequence` module for all sequence analysis functions. 

|questions|

|suggestions|
