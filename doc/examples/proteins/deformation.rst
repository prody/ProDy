.. currentmodule:: prody.compare

.. _deformation:

*******************************************************************************
Deformation and ANM analysis
*******************************************************************************

Getting the deformation vector that describes the change in atomic positions
of a protein from one of it's structures to another is most of the time
tricky. ProDy implements :func:`findMatchingChains` to make this an easy task.
This function finds matching chains and returns corresponding atoms form these
chains.



.. literalinclude:: deformation.py
