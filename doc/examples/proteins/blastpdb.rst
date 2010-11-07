.. currentmodule:: prody.proteins

.. _blastpdb:

*******************************************************************************
Blast search PDB content
*******************************************************************************

This example shows how to search for PDB structures matching an amino acid
sequence using :func:`blastPDB` function. :func:`blastPDB` is a utility 
function which can be used to check if structures matching a sequence exists 
in PDB or to access a set of related structures for ensemble analysis (i.e. 
essential dynamics analysis). 

:func:`blastPDB` returns a :class:`PDBlastRecord`, whose method can return the
user the best match (:meth:`PDBlastRecord.getBest`) or hits that sharing
share a sequence identity better than a user given value 
(:meth:`PDBlastRecord.getHits`).  

.. literalinclude:: blastpdb.py
