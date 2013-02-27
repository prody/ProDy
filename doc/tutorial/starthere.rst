.. _start-here:

.. currentmodule:: prody

*******************************************************************************
Start Here
*******************************************************************************

How to Use ProDy
===============================================================================

ProDy can be used:

  * interactively in a Python shell,
  * as a command line program via :ref:`prody-apps`,
  * from within VMD via :ref:`nmwiz`,
  * or as a toolkit for developing new software.

Python for Beginners
===============================================================================

Users who are new to programming or Python are referred to the following 
resources for an introduction to programming in Python:

* `The Python Tutorial <http://docs.python.org/tutorial/>`_
* `Python Scientific Lecture Notes <http://scipy-lectures.github.com/>`_
* `A Primer on Python for Life Science Researchers
  <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.0030199>`_


Using Documentation
===============================================================================

ProDy documentation is organized in three main sections:

* Tutorial, contains brief usage examples of some ProDy features and is the 
  best place to start to learn ProDy.
* :ref:`examples` section contains comprehensive usage examples, which may be 
  applied to other cases after small modifications.
* :ref:`reference` section describes all ProDy classes and functions, with some
  usage examples. Features are divided into following main modules:

  * :mod:`.atomic` - efficient handling of atomic data
  * :mod:`.dynamics` - analysis and modeling of protein dynamics
  * :mod:`.ensemble` - analysis of arbitrary conformational ensembles
  * :mod:`.measure` - measure properties, transform coordinates 
  * :mod:`.proteins` - structure analysis, IO, and data retrieval
  * :mod:`.trajectory` - trajectory IO
  * :mod:`.utilities` - helper functions and classes
  
In interactive sessions, the reference documentation can be accessed 
using the built-in Python function :func:`help`:: 

  help(atomic)
  help(select)
  help(parsePDB)

This function prints the description of functions, classes, and class methods 
to the screen. Note that you might need to type ``q`` to exit from 
help. If you are using the interactive Python shell (IPython), you can also 
receive help by typing::
  
  select ?

Copying Examples
===============================================================================

Documentation contains ProDy code snippets.  These snippets can be reformatted 
using the :guilabel:`Show Code Snippets` button on the right hand side panel. 
The code will be displayed in the same window, and can be copied directly into 
a file. Click on the text, press :kbd:`Ctrl+A` and then :kbd:`Ctrl+C` to have 
the text in your clipboard. To return to the documentation click the 
:guilabel:`Show documentation` button at the top.

.. image:: /_static/codesnippets.png
   :align: center
   :alt: Getting ProDy code snippets.


Interactive Usage
===============================================================================

One of our aims is making ProDy suitable for interactive usage by designing 
flexible functions and classes and giving them easy to remember names in a 
consistent manner. 

For best interactive usage experience, we strongly recommend that you use 
`IPython <http://ipython.scipy.org/>`_ or a similar interactive shell instead 
of the standard Python shell.  The IPython shell, for example, provides 
user-friendly features, such as dynamic introspection and help, and also 
optionally convenient integration of Numpy and `Matplotlib 
<http://matplotlib.sourceforge.net>`_.
