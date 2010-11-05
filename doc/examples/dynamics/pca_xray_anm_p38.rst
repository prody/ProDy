.. module:: prody.dynamics

.. _egp38:

*******************************************************************************
PCA/ANM of p38 X-ray Structures
*******************************************************************************

This example repeats the calculations for p38 MAP kinase that was 
published in [AB09]_.

Following classes and functions are used in the example:

Classes:
  * :class:`ANM`
  * :class:`PCA`
  * :class:`Ensemble`
Functions:
  * :func:`prody.proteins.fetchPDB`
  * :func:`prody.proteins.mapAtomsToChain`
  
This example also use :class:`prody.proteins.AtomMap` instances. Further reading
:ref:`atommaps` may help to understand how an X-ray ensemble is constructed.
  
.. literalinclude:: pca_xray_anm_p38.py
