.. currentmodule:: prody.atomic

.. _atommaps:

How AtomMap's work
===============================================================================
    
:class:`AtomMap` class adds great flexibility to manipulating atomic data.

First let's see how an instance of :class:`Selection` 
(or :class:`Chain`\/:class:`Residue`) works. Below table shows 
indices for a selection of atoms in an :class:`AtomGroup` and values returned when 
:meth:`~Selection.getNames`, :meth:`~Selection.getResnames` and
:meth:`~Selection.getResnums` methods are called.

.. csv-table:: **Atom Subset** 
   :header: "Indices", ":meth:`getNames`", ":meth:`getResnames`", ":meth:`getResnums`"

   0, N, PHE, 1
   1, CA, PHE, 1
   2, C, PHE, 1
   3, O, PHE, 1
   4, CB, PHE, 1
   5, CG, PHE, 1
   6, CD1, PHE, 1
   7, CD2, PHE, 1
   8, CE1, PHE, 1
   9, CE2, PHE, 1
   10, CZ, PHE, 1

:class:`Selection` instances keep indices ordered and do not allow duplicate values, hence
their use is limited. In an :class:`AtomMap`, indices do not need to be sorted,
duplicate indices may exist, even "DUMMY" atoms are allowed.

Let's say we instantiate the following AtomMap::
    
    >>> amap = AtomMap(atomgroup, indices=[0, 1, 3, 8, 8, 9, 10], mapping=[5, 6, 7, 0, 1, 2, 3])



The size of the AtomMap based on this mapping is 8, since the larger mapping is 7.

Calling the same functions for this AtomMap instance would result in the following:

.. csv-table:: **Atom Map**
   :header: "Mapping", "Indices", ":meth:`getNames`", ":meth:`getResnames`", ":meth:`getResnums`", ":meth:`getMappedFlags`", ":meth:`getUnmappedFlags`"

   0, 8, CE1, PHE, 1, 1, 0
   1, 8, CE1, PHE, 1, 1, 0
   2, 9, CE2, PHE, 1, 1, 0
   3, 10, CZ, PHE, 1, 1, 0
   4, , , , 0, 0, 1
   5, 0, N, PHE, 1, 1, 0
   6, 1, CA, PHE, 1, 1, 0
   7, 3, O, PHE, 1, 1, 0
   
For unmapped atoms, numeric attributes are set to 0, others to empty string,
i.e. ``""``.

AtomMaps are used by functions that compare protein chains.

.. seealso::
   :ref:`pca-xray` and :ref:`pca-dimer` examples make use of AtomMaps.
   

