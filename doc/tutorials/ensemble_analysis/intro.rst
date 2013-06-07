Introduction
===============================================================================

This tutorial shows how to analyze ensembles of experimental structures.


.. [AB09] Bakan A, Bahar I. The intrinsic dynamics of enzymes
   plays a  dominant role in determining the structural
   changes induced upon inhibitor binding. *Proc Natl Acad Sci U S A.*
   **2009** 106(34):14349-54.


Required Programs
-------------------------------------------------------------------------------

Latest version of `ProDy`_ and `Matplotlib`_ are required.

.. _ProDy: http://csb.pitt.edu/ProDy/getprody.html
.. _Matplotlib: http://matplotlib.org/

Recommended Programs
-------------------------------------------------------------------------------

`IPython`_ and `Scipy`_ are recommended for this tutorial.


.. _IPython: http://ipython.org/
.. _Scipy: http://scipy.org/


Getting Started
-------------------------------------------------------------------------------

To follow this tutorial, you will need the following files which can be
downloaded from :ref:`tutorials`.

.. files.txt will be automatically generated

.. literalinclude:: files.txt


We assume that you will follow this tutorial by typing commands in an
interactive Python session. We recommend that you use IPython. First,
we will make necessary imports from ProDy and Matplotlib packages:

.. ipython:: python

   from prody import *
   from pylab import *
   ion()
