.. _external-matrix:

*******************************************************************************
Using an External Matrices
*******************************************************************************

Synopsis
===============================================================================

This example shows how to use matrices from external software in ANM and GNM
analysis.

Input
-------------------------------------------------------------------------------

Hessian/Kirchhoff matrix stored in a text file.

Output
-------------------------------------------------------------------------------

:class:`~.ANM` or :class:`~.GNM` instances that can be used for analysis of protein
dynamics.

ProDy Code
===============================================================================

We start by importing everything from the ProDy package:

>>> from prody import *

Parse Hessian matrix
-------------------------------------------------------------------------------

The input file that contains the Hessian matrix has the following format 
(:download:`oanm_hes.txt </doctest/oanm_hes.txt>`)::

       1       1    9.958948135375977e+00
       1       2   -3.788214445114136e+00
       1       3    6.236155629158020e-01
       1       4   -7.820609807968140e-01
       1       5    1.050322428345680e-01
       1       6   -3.992616236209869e-01
       1       7   -7.818332314491272e-01
       1       8   -1.989762037992477e-01
       1       9   -3.619094789028168e-01
       1      10   -5.224789977073669e-01
       ...

:func:`~.parseSparseMatrix` can be used for parsing the above file:

>>> hessian = parseSparseMatrix('oanm_hes.txt', symmetric=True)
>>> hessian.shape
(1053, 1053)

ANM calculations
-------------------------------------------------------------------------------

Rest of the calculations can be performed as follows:

>>> anm = ANM('Using external Hessian')
>>> anm.setHessian(hessian)
>>> anm.calcModes()
>>> anm
<ANM: Using external Hessian (20 modes, 351 nodes)>

For more information, see :ref:`anm`.

Parse Kirchhoff matrix
-------------------------------------------------------------------------------

The input file that contains the Kirchhoff matrix has the following format
(:download:`ognm_kirchhoff.txt </doctest/ognm_kirchhoff.txt>`)::

        3316
       1       1       5.00
       1       2      -1.00
       1       3      -1.00
       1       4      -1.00
       1      91      -1.00
       1     343      -1.00
       2       2      10.00
       2       3      -1.00
       2       4      -1.00
       ...

>>> kirchhoff = parseSparseMatrix('ognm_kirchhoff.txt', symmetric=True, 
...                               skiprows=1)
>>> kirchhoff.shape
(351, 351)

GNM calculations
-------------------------------------------------------------------------------

Rest of the GNM calculations can be performed as follows:

>>> gnm = GNM('Using external Kirchhoff')
>>> gnm.setKirchhoff(kirchhoff)
>>> gnm.calcModes()
>>> gnm
<GNM: Using external Kirchhoff (20 modes, 351 nodes)>

For more information, see :ref:`gnm`.

|questions|

|suggestions|
