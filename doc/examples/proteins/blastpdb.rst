.. currentmodule:: prody.proteins

.. _blastpdb:

*******************************************************************************
Blast search PDB
*******************************************************************************

Synopsis
===============================================================================

This examples demonstrates how to use Protein Data Bank blast search function, 
:func:`blastPDB`. 

:func:`blastPDB` is a utility function which can be used to check if 
structures matching a sequence exists in PDB or to access a set of related 
structures for ensemble analysis (i.e. ref:`pca`). 

User Input
-------------------------------------------------------------------------------

Protein amino acid sequence is the only import from the user, e.g. 
``ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPIS...
YDIVKMKKSNISPNFNFMGQLLDFERTL``

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Perform search
-------------------------------------------------------------------------------

Let's search for structures similar to that of MKP-3, using its sequence

>>> from prody import *
>>> blast_record = blastPDB('''
... ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPIS
... DHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDA
... YDIVKMKKSNISPNFNFMGQLLDFERTL''')

:func:`blastPDB` function returns a :class:`PDBlastRecord`. Let's retrieve 
hits from this record:

Best match
-------------------------------------------------------------------------------

To get the best match, :meth:`PDBlastRecord.getBest` method can be used:

>>> best = blast_record.getBest()
>>> print best # doctest: +SKIP
{'percent_identity': 100.0, 
'sbjct_end': 144, 
'pdb_title': 'Chain A, Crystal Structure Of Pyst1 (Mkp3)', 
'query': u'ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPISDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDAYDIVKMKKSNISPNFNFMGQLLDFERTL', 
'match': u'ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPISDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDAYDIVKMKKSNISPNFNFMGQLLDFERTL', 
'sbjct': u'ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPISDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDAYDIVKMKKSNISPNFNFMGQLLDFERTL', 
'identities': 144, 
'query_start': 1, 
'query_end': 144, 
'chain_id': 'A', 
'pdb_id': '1mkp', 
'sbjct_start': 1, 
'score': 762.0, 
'gaps': 0, 
'expect': 7.6441600000000005e-82, 
'percent_coverage': 100.0, 
'positives': 144, 
'align_length': 144, 
'bits': 298.13}

PDB hits
-------------------------------------------------------------------------------

>>> hits = blast_record.getHits()
>>> print hits.keys()
['1mkp']

This results in only MKP-3 itself, since percent_identity argument was set 
to 90 by default

>>> hits = blast_record.getHits(percent_identity=50)
>>> print hits.keys()
['2hxp', '3lj8', '1mkp']

>>> hits = blast_record.getHits(percent_identity=40)
>>> print hits.keys()
['3lj8', '1mkp', '1zzw', '2g6z', '2hxp', '3ezz', '1m3g', '2oud']

This resulted in 7 hits, including structures of MKP-2, MKP-4, and MKP-5
More information on a hit can be obtained as follows:

>>> print hits['1zzw']['percent_identity']
47.2222222222
>>> print hits['1zzw']['align_length']
138
>>> print hits['1zzw']['identities']
68

Download hits
-------------------------------------------------------------------------------

PDB hits can be downloaded using :func:`fetchPDB` 
function.

>>> filenames = fetchPDB(hits.keys())
>>> print filenames # doctest: +SKIP
['./1mkp.pdb.gz', './1zzw.pdb.gz', './2g6z.pdb.gz', './2hxp.pdb.gz', 
'./3ezz.pdb.gz', './1m3g.pdb.gz', './2oud.pdb.gz']



|questions|

|suggestions|

