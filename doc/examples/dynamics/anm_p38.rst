.. module:: prody.dynamics

.. _eganm:

*******************************************************************************
Anisotropic Network Model analysis
*******************************************************************************

This example shows how to perform ANM calculations. ANM results are saved
in NMD format, which can be visualized using |nmwiz| |vmd| plugin. 

Following ProDy classes and functions are used in the example:

Classes:
  * :class:`ANM`
Functions:
  * :func:`prody.proteins.parsePDB`
  * :func:`getVMDpath`
  * :func:`viewNMDinVMD`
  * :func:`writeNMD`

.. literalinclude:: anm_p38.py
