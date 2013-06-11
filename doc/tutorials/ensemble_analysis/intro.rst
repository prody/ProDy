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


We recommend that you will follow this tutorial by typing commands in an
IPython session, e.g.::

  $ ipython

or with pylab environment::

  $ ipython --pylab


First, we will make necessary imports from ProDy and Matplotlib
packages.

.. ipython:: python

   from prody import *
   from pylab import *
   ion()

We have included these imports in every part of the tutorial, so that
code copied from the online pages is complete. You do not need to repeat
imports in the same Python session.