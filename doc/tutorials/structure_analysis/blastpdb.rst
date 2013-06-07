.. _blastpdb:


Blast Search PDB
===============================================================================

This examples demonstrates how to use Protein Data Bank blast search function,
:func:`.blastPDB`.

:func:`.blastPDB` is a utility function which can be used to check if
structures matching a sequence exists in PDB or to identify a set of related
structures for :ref:`pca`.

We will used amino acid sequence of a protein, e.g.
``ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLF...YDIVKMKKSNISPNFNFMGQLLDFERTL``

The :func:`.blastPDB` function accepts sequence as a Python :func:`str`.

Output will be :class:`.PDBBlastRecord` instance that stores PDB hits and
returns to the user those sharing sequence identity above a user specified
value.

Blast search
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *

Let's search for structures similar to that of MKP-3, using its sequence:

.. ipython:: python

   blast_record = blastPDB('''
   ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPISDHWSQNLSQFFPEA
   ISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDAYDIVKMKKSNISPNFNFMGQLLDFERTL''')

:func:`.blastPDB` function returns a :class:`.PDBBlastRecord`. Let's retrieve
hits from this record:

Best match
-------------------------------------------------------------------------------

To get the best match, :meth:`.PDBBlastRecord.getBest` method can be used:

.. ipython:: python

   best = blast_record.getBest()
   best['pdb_id']
   best['percent_identity']


PDB hits
-------------------------------------------------------------------------------

.. ipython:: python

   hits = blast_record.getHits()
   list(hits)

This results in only MKP-3 itself, since percent_identity argument was set
to 90 by default:

.. ipython:: python

   hits = blast_record.getHits(percent_identity=50)
   list(hits)
   hits = blast_record.getHits(percent_identity=40)
   list(hits)


This resulted in 7 hits, including structures of MKP-2, MKP-4, and MKP-5
More information on a hit can be obtained as follows:

.. ipython:: python

   hits['1zzw']['percent_identity']
   hits['1zzw']['align-len']
   hits['1zzw']['identity']

Download hits
-------------------------------------------------------------------------------

PDB hits can be downloaded using :func:`.fetchPDB` function::

  filenames = fetchPDB(hits.keys())
  filenames