Introduction
===============================================================================

This tutorial describes how to use elastic network models, in particular
:ref:`gnm` and :ref:`anm`, for studying protein dynamics.


Required Programs
-------------------------------------------------------------------------------

Latest version of `ProDy`_ and `Matplotlib`_ required.

.. _ProDy: http://csb.pitt.edu/ProDy/getprody.html
.. _Matplotlib: http://matplotlib.org/

Recommended Programs
-------------------------------------------------------------------------------

`IPython`_ is highly recommended for interactive usage.

.. _IPython: http://ipython.org/


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
