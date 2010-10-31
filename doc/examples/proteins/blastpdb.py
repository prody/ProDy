#!/usr/bin/env python
"""blastPDB - search PDB structures matching a amino acid sequence

This example shows how to search PDB structures matching an amino acid
sequence using :func:`blastPDB` function. :func:`blastPDB` is a utility 
function which can be used to check if structures matching a sequence exists 
in PDB or to access a set of related structures for ensemble analysis (i.e. 
essential dynamics analysis). 

:func:`blastPDB` returns a :class:`PDBlastRecord`, whose method can return the
user the best match (:meth:`PDBlastRecord.getBest`) or hits that sharing
share a sequence identity better than a user given value 
(:meth:`PDBlastRecord.getHits`).  

"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'

from prody import *

# Let's search for structures similar to that of MKP-3, using its sequence
blast_record = blastPDB('''
ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPIS
DHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMNDA
YDIVKMKKSNISPNFNFMGQLLDFERTL''')

# blastPDB function returns a BlastRecord
# Let's parse hits from the record
hits = blast_record.getHits()
print hits.keys()
# ['1mkp']

# If we wanted to get the best hit, we would use getBest() method
best = blast_record.getBest()
print best


# This results in only MKP-3 itself, since percent_identity argument was set 
# to 90 by default
hits = blast_record.getHits(percent_identity=50)
print hits.keys()
# ['2hxp', '1mkp']

hits = blast_record.getHits(percent_identity=40)
print hits.keys()
# ['1mkp', '1zzw', '2g6z', '2hxp', '3ezz', '1m3g', '2oud']
# This resulted in 7 hits, including structures of MKP-2, MKP-4, and MKP-5
# More information on a hit can be obtained as follows
print hits['1zzw']['percent_identity']
# 47.2222222222
print hits['1zzw']['align_length']
# 138
print hits['1zzw']['identities']
# 68

# PDB hits can be downloaded using :func:`prody.proteins.proteindata.fetchPDB` 
# function.

filenames = fetchPDB(hits.keys())
print filenames
# ['./1mkp.pdb.gz', './1zzw.pdb.gz', './2g6z.pdb.gz', './2hxp.pdb.gz', 
# './3ezz.pdb.gz', './1m3g.pdb.gz', './2oud.pdb.gz']

# Let's delete these files
import os
for fn in filenames:
    os.remove(fn)

