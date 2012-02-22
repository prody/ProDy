.. _blastpdb:

*******************************************************************************
Blast search PDB
*******************************************************************************

Synopsis
===============================================================================

This examples demonstrates how to use Protein Data Bank blast search function, 
:func:`~.blastPDB`. 

:func:`~.blastPDB` is a utility function which can be used to check if 
structures matching a sequence exists in PDB or to identify a set of related 
structures for ensemble analysis (i.e. :ref:`pca`). 

Input
-------------------------------------------------------------------------------

Amino acid sequence of a protein, e.g. 
``ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLF...YDIVKMKKSNISPNFNFMGQLLDFERTL``

The :func:`~.blastPDB` function accepts sequence as a Python :func:`str`. 

Output
-------------------------------------------------------------------------------
 
Output is an :class:`~.PDBBlastRecord` instance that stores PDB hits and returns
to the user those sharing sequence identity above a user specified value. 

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Perform search
-------------------------------------------------------------------------------

Let's search for structures similar to that of MKP-3, using its sequence:

>>> from prody import *
>>> blast_record = blastPDB('''
... ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPIS
... DHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDA
... YDIVKMKKSNISPNFNFMGQLLDFERTL''')

:func:`~.blastPDB` function returns a :class:`~.PDBBlastRecord`. Let's retrieve 
hits from this record:

Best match
-------------------------------------------------------------------------------

To get the best match, :meth:`.PDBBlastRecord.getBest` method can be used::

  best = blast_record.getBest()
  print best['pdb_id']
  # 1mkp
  print best['percent_identity']
  # 100.0
 
PDB hits
-------------------------------------------------------------------------------

::

  hits = blast_record.getHits()
  print( hits.keys() )
  # ['1mkp']

This results in only MKP-3 itself, since percent_identity argument was set 
to 90 by default::

  hits = blast_record.getHits(percent_identity=50)
  print( hits.keys() )
  # ['1m3g', '2hxp', '3lj8', '3ezz', '1mkp']

  hits = blast_record.getHits(percent_identity=40)
  print( hits.keys() )
  # ['3lj8', '1mkp', '1zzw', '2g6z', '2hxp', '3ezz', '1m3g', '2oud']

This resulted in 7 hits, including structures of MKP-2, MKP-4, and MKP-5
More information on a hit can be obtained as follows::

  print( hits['1zzw']['percent_identity'] )
  49.2753623188
  print( hits['1zzw']['align-len'] )
  138
  print( hits['1zzw']['identity'] )
  68

Download hits
-------------------------------------------------------------------------------

PDB hits can be downloaded using :func:`~.fetchPDB` function::

  filenames = fetchPDB(hits.keys())
  print( filenames ) # doctest: +SKIP
  # ['./1mkp.pdb.gz', './1zzw.pdb.gz', './2g6z.pdb.gz', './2hxp.pdb.gz', 
  #  './3ezz.pdb.gz', './1m3g.pdb.gz', './2oud.pdb.gz']



|questions|

|suggestions|

