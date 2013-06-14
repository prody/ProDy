.. _external-matrix:

Using an External Matrix
===============================================================================

This example shows how to use matrices from external software in :class:`.ANM`
or :class:`.GNM` analysis of protein dynamics.

Parse Hessian
-------------------------------------------------------------------------------

We start by importing everything from the ProDy package:

.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()  # turn interactive mode on

The input file that contains the Hessian matrix has the following format
(:download:`oanm_hes.txt <enm_analysis_files/oanm_hes.txt>`)::

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

:func:`.parseSparseMatrix` can be used for parsing the above file:

.. ipython:: python

   hessian = parseSparseMatrix('enm_analysis_files/oanm_hes.txt',
                               symmetric=True)
   hessian.shape


ANM calculations
-------------------------------------------------------------------------------

Rest of the calculations can be performed as follows:

.. ipython:: python

   anm = ANM('Using external Hessian')
   anm.setHessian(hessian)
   anm.calcModes()
   anm

For more information, see :ref:`anm`.

Parse Kirchhoff
-------------------------------------------------------------------------------

The input file that contains the Kirchhoff matrix has the following format
(:download:`enm_analysis_files/ognm_kirchhoff.txt`)::

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

.. ipython:: python

   kirchhoff = parseSparseMatrix('enm_analysis_files/ognm_kirchhoff.txt',
                                 symmetric=True, skiprows=1)
   kirchhoff.shape


GNM calculations
-------------------------------------------------------------------------------

Rest of the GNM calculations can be performed as follows:

.. ipython:: python

   gnm = GNM('Using external Kirchhoff')
   gnm.setKirchhoff(kirchhoff)
   gnm.calcModes()
   gnm


For more information, see :ref:`gnm`.
