.. _pfamaccess:

Pfam Access
===============================================================================

Synopsis
-------------------------------------------------------------------------------

The tutorial shows how to:

  * search Pfam database to identify family accession numbers and information
  * fetch the MSA of the Pfam using accession no

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
sequences, default is ``timeout=30`` seconds.

.. ipython:: python

   matches = searchPfam(sequence, search_b=True, evalue=2.0)


Retrieve MSA files
-------------------------------------------------------------------------------

This example demonstrates how to search Pfam database with a given query using
:func:`.fetchPfamMSA`. Valid inputs are Pfam ID, e.g. :pfam:`Piwi`, or Pfam
accession, e.g. :pfam:`PF02171` obtained from :func:`.searchPfam`.  Alignment
type can be ``"full'`` (default), ``"seed"``, ``"ncbi"`` or ``"metagenomics"``
or ``"rp15"`` or ``"rp35"`` or ``"rp55"`` or ``"rp75"``.

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
