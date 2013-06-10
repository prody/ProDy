.. _msafiles:

MSA Files
===============================================================================

This part shows how to parse, refine, filter, slice, and write MSA files.


.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on

Let's get Pfam MSA file for protein family that contains :uniprot:`PIWI_ARCFU`:

.. ipython::
   :verbatim:

   In [1]: searchPfam('PIWI_ARCFU').keys()
   Out[1]: ['PF02171']
   In [2]: fetchPfamMSA('PF02171', alignment='seed')
   Out[2]: 'PF02171_seed.sth'

Parsing MSA files
-------------------------------------------------------------------------------

This shows how to use the :class:`.MSAFile` or :func:`.parseMSA` to read the
MSA file.

Reading using :class:`.MSAFile` yields an MSAFile object. Iterating over the
object will yield an object of :class:`.Sequence` from which labels, sequence
can be obtained.

.. ipython:: python

   msafile = 'PF02171_seed.sth'
   msafobj = MSAFile(msafile)
   msafobj
   msa_seq_list = list(msafobj)
   msa_seq_list[0]

:func:`.parseMSA` returns an :class:`.MSA` object.  We can parse
compressed files, but reading uncompressed files are much faster.

.. ipython::
   :verbatim:

   In [1]: fetchPfamMSA('PF02171', compressed=True)
   Out[1]: 'PF02171_full.sth.gz'

   In [2]: parseMSA('PF02171_full.sth.gz')
   Out[2]: <MSA: PF02171_full (2067 sequences, 1392 residues)>

   In [3]: fetchPfamMSA('PF02171', format='fasta')
   Out[3]: 'PF02171_full.fasta.gz'

   In [3]: parseMSA('PF02171_full.fasta.gz')
   Out[3]: <MSA: PF02171_full (2067 sequences, 1392 residues)>


Iterating over a file will yield sequence id, sequence, residue start and
end indices:

.. ipython:: python

   msa = MSAFile(`PF02171_seed.sth`)
   for seq in msa:
       seq

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
       seq.getLabel()

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


MSA objects
-------------------------------------------------------------------------------

Indexing
^^^^^^^^

Retrieving a sequence at a given index, or by id will give an object of
:class:`.Sequence`:

.. ipython:: python

   msa = parseMSA(msafile)
   seq = msa[0]
   seq
   str(seq)

Retrieve a sequence by UniProt ID:

.. ipython:: python

   msa['YQ53_CAEEL']


Querying
^^^^^^^^

You can query whether a sequence in contained in the instance using the
UniProt identifier of the sequence as follows:

.. ipython:: python

   'YQ53_CAEEL' in msa

Slicing
^^^^^^^


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


Merging MSAs
-------------------------------------------------------------------------------

:func:`.mergeMSA` can be used to merge two or more MSAs. Based on their labels
only those sequences that appear in both MSAs are retained, and concatenated
horizontally to give a joint or merged MSA. This can be useful while evaluating
covariance patterns for proteins with multiple domains or protein-protein
interactions. The example shows merging for the multi-domain receptor
:pdb:`3KG2` containing pfam domains :pfam:`PF01094` and :pfam:`PF00497`.

.. ipython::
   :verbatim:

   In [1]: fetchPfamMSA('PF01094', format='fasta', timeout=120)
   Out[1]: 'PF01094_full.fasta'
   In [2]: fetchPfamMSA('PF00497', format='fasta', timeout=120)
   Out[2]: 'PF00497_full.fasta'

Let's parse and merge the two files:

.. ipython:: python

   msa1 = parseMSA('PF01094_full.fasta')
   msa1
   msa2 = parseMSA('PF00497_full.fasta')
   msa2
   merged = mergeMSA(msa1, msa2)
   merged

Merged MSA contains 984 sequences.

Writing MSAs
-------------------------------------------------------------------------------

:func:`.writeMSA` can be used to write MSA. It takes filename as input
which should contain appropriate extension that can be ``".slx"`` or
``".sth"`` or  ``".fasta"`` or format should be specified as ``"SELEX"``,
``"Stockholm"`` or ``"FASTA"``. Input MSA should be :class:`.MSAFile` or
:class:`.MSA` object. Filename can contain ``".gz"`` extension, in which case
a compressed file will be written.


.. ipython:: python

   writeMSA('sliced_MSA.gz', msa, format='SELEX')
   writeMSA('sliced_MSA.fasta', msafobj)

:func:`.writeMSA` returns the name of the MSA file that is written.
