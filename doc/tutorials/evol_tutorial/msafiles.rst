.. _msafiles:

Pfam Database and MSA Files
===============================================================================

Synopsis
-------------------------------------------------------------------------------

The part shows how to:

  * search Pfam database to identify family accession numbers and information
  * fetch the MSA of the Pfam using accession no
  * parse the MSA, filter, slice MSA, write the MSA


.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on

Search Pfam
-------------------------------------------------------------------------------

This example demonstrates how to search Pfam database with a given query,
:func:`.searchPfam`.  Valid inputs are UniProt ID, e.g. :uniprot:`PIWI_ARCFU`,
or PDB identifier, e.g. :pdb:`3luc` or ``"3lucA"`` with chain identifier.
Input can also be a protein sequence or a file containing the sequence,
but sequence should not contain gaps and should be at least 12 characters long.

Matching Pfam accession (one or more) as keys will map to a dictionary that
contains locations (alignment start, end, evalue etc), pfam family type,
accession and id.



We query Pfam using the :func:`.searchPfam`. with a UniProt ID.

.. ipython:: python

   matches = searchPfam('PIWI_ARCFU')
   matches


This function also accepts a protein sequence:

.. ipython:: python

   sequence = ('PMFIVNTNVPRASVPDGFLSELTQQLAQATGKPPQYIAVHVVPDQLMAFGGSSEPCA'
   'LCSLHSIGKIGGAQNRSYSKLLCGLLAERLRISPDRVYINYYDMNAANVGWNNSTFA')
   matches = searchPfam(sequence)
   matches


For sequence searches, we can pass additional parameters to :func:`.searchPfam`
like *search_b* which will search pfam B and *skip_a* that will not search
pfamA database. Additional parameters include *ga* that uses gathering
threshold instead of e-value, *evalue* cutoff can also be specified and
*timeout* that can be set higher especially when searching larger
sequences, default is ``timeout=60`` seconds.

.. ipython:: python

   matches = searchPfam(sequence, search_b=True, evalue=2.0)


Retrieve MSA files
-------------------------------------------------------------------------------

This example demonstrates how to search Pfam database with a given query using
:func:`.fetchPfamMSA`. Valid inputs are Pfam ID, e.g. :pfam:`Piwi`, or Pfam
accession, e.g. :pfam:`PF02171` obtained from :func:`.searchPfam`.  Alignment
type can be ``"full'`` (default), ``"seed"``, ``"ncbi"`` or ``"metagenomics"``.

.. ipython:: python

   fetchPfamMSA('piwi', alignment='seed')
   msafile = 'piwi_seed.sth'

A compressed file can be downloaded by setting ``compressed=True``.
The ``format`` of the MSA can be of ``"selex"`` (default), ``"stockholm"`` or
``"fasta"``.  This will return the path of the downloaded MSA file.
The ``output`` name can be specified, for by default it will have
``"accession/ID_alignment.format"``.

Note that in this case we passed a folder name, the downloaded file is saved
in this folder, after it is created if it did not exist. Also bigger timeouts
are necessary for larger families. Some other parameters like ``gap``,
``order`` or ``inserts`` can be set, as shown in the following example.

.. ipython:: python

   fetchPfamMSA('PF02171', compressed=True, gaps='mixed', inserts='lower',
   order='alphabetical', format='fasta', timeout=40)



Parsing MSA files
-------------------------------------------------------------------------------

This shows how to use the :class:`.MSAFile` or :func:`.parseMSA` to read the
MSA file. :func:`.parseMSA` returns a :class:`.MSA` object.

Reading using :class:`.MSAFile` yields an MSAFile object. Iterating over the
object will yield an object of :class:`.Sequence` from which labels, sequence
can be obtained.

.. ipython:: python

   msafobj = MSAFile(msafile)
   msafobj
   msa_seq_list = list(msafobj)
   msa_seq_list[0]

Reading using :func:`.parseMSA` yields an :class:`.MSA` object.  We can parse
compressed files, but reading uncompressed files are much faster.

.. ipython:: python

   msa = parseMSA('PF02171_full.fasta.gz')
   msa
   msa = parseMSA(fetchPfamMSA('PF02171', format='fasta'))
   msa


Filtering and Slicing
-------------------------------------------------------------------------------

This shows how to use the :class:`.MSAFile` object or :class:`.MSA` object to
refine MSA using filters and slices.

Filtering
^^^^^^^^^

Any function that takes label and sequence arguments and returns a boolean
value can be used for filtering the sequences.  A sequence will be yielded
if the function returns **True**.  In the following example, sequences from
organism *ARATH* are filtered:

.. ipython:: python

   msafobj = MSAFile(msafile, filter=lambda lbl, seq: 'ARATH' in lbl)
   for seq in msafobj:
       print(seq.getLabel())

Slicing
^^^^^^^

A list of integers can be used to slice sequences as follows.  This enables
selective parsing of the MSA file.

.. ipython:: python

   msafobj = MSAFile(msafile, slice=list(range(10)) + list(range(374,384)))
   list(msafobj)[0]


Slicing can also be done using :class:`.MSA`. The :class:`.MSA` object offers
other functionalities like querying, indexing, slicing row and columns and
refinement.

Querying
^^^^^^^^

You can query whether a sequence in contained in the instance using the
UniProt identifier of the sequence as follows:

.. ipython:: python

   msa = parseMSA(msafile)
   'YQ53_CAEEL' in msa


Indexing MSA objects
-------------------------------------------------------------------------------

Retrieving a sequence at a given index, or by id will give an object of
:class:`.Sequence`:

.. ipython:: python

   msa = parseMSA(msafile)
   msa[0]

Retrieve a sequence by UniProt ID:

.. ipython:: python

   msa['YQ53_CAEEL']

Slice an MSA instance to give a new :class:`.MSA`. object :

.. ipython:: python

   new_msa = msa[:2]
   new_msa

Slice using a list of UniProt IDs:

.. ipython:: python

   msa[:2] == msa[['TAG76_CAEEL', 'O16720_CAEEL']]

Retrieve a character or a slice of a sequence:

.. ipython:: python

   msa[0,0]
   msa[0,0:10]

Slice MSA rows and columns:

.. ipython:: python

   msa[:10,20:40]

Writing MSA files
-------------------------------------------------------------------------------

:func:`.writeMSA` can be used to write MSA. It takes filename as input
which should contain appropriate extension that can be ``".slx"`` or
``".sth"`` or  `".fasta"`` or format should be specified as ``"SELEX"``,
``"FASTA"`` or ``"Stockholm"``. Input MSA should be :class:`.MSAFile` or
:class:`.MSA` object. Filename can contain ``".gz"`` extension, in which case
a compressed file will be written.
Returns the name of the MSA file that is written.

.. ipython:: python

   writeMSA('sliced_MSA.gz', msa, format='SELEX')
   filename = writeMSA('sliced_MSA.fasta', msafobj)
   filename