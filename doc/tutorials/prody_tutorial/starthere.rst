How to Start
===============================================================================

Using ProDy
-------------------------------------------------------------------------------

ProDy can be used in a number of ways:

  #. interactively in a Python shell,
  #. as a command line program via :ref:`prody-apps`,
  #. from within VMD via :ref:`nmwiz`,
  #. or as a toolkit for developing new software.


Python for beginners
^^^^^^^^^^^^^^^^^^^^

Familiarity with Python programming language will help when using *ProDy*.  If
you are new to Python, or to programming, you may start with one of the
following tutorials:

  * `The Python Tutorial <http://docs.python.org/tutorial/>`_
  * `Python Scientific Lecture Notes <http://scipy-lectures.github.com/>`_
  * `A Primer on Python for Life Science Researchers
    <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.0030199>`_



Interactive Usage
-------------------------------------------------------------------------------

In the rest of this tutorial, we assume that you will be typing commands in a
Python shell.  *ProDy* will automatically download PDB files and save them to
current working directory, so you may want start Python from inside of a
directory that you make for this tutorial::

  $ mkdir prody_tutorial
  $ cd prody_tutorial


Start Python shell
^^^^^^^^^^^^^^^^^^

For best interactive usage experience, we strongly recommend that you use
`IPython`_ instead of the standard Python shell. IPython shell provides many
user-friendly features, such as dynamic introspection and help, and also
convenient integration of `Numpy`_ and `Matplotlib`_.


If you have installed IPython, type in::

  $ ipython

If you also installed Matplotlib, use::

  $ ipython --pylab

``--pylab`` option will import Matplotlib and Numpy automatically, and is
equivalent to the following:

.. ipython:: python

   from pylab import *
   ion()  # turn interactive mode on


If you don't have IPython yet, use::

  $ python


On Windows, after you make the directory, make a :kbd:`Shift+right click` in it
in Windows Explorer and then select :menuselection:`Open command window here`
option.  Then start ``C:\Python27\python.exe``.  Alternatively, you may
run :program:`IDLE (Python GUI)` or :program:`Python (command line)` from the
start menu.


Import from ProDy
^^^^^^^^^^^^^^^^^

We import all *ProDy* functions and classes into the current namespace as
follows:

.. ipython:: python

   from prody import *

There are other ways to import *ProDy* contents.  You may use ``import prody as
pd`` and prefix all functions calls with ``pd.``, if you prefer not to
overcrowd the target namespace.  Alternatively, if you want to use contents of
a specific module, such as :mod:`proteins`, you can use
``from prody.proteins import *``.  You should, however, avoid using
``from prody.proteins.pdbfile import *``, because location of methods
in submodules may change without notice.

Using Documentation
-------------------------------------------------------------------------------

ProDy documentation is quite comprehensive and you can access it in a number of
different ways.  In interactive sessions, API reference can be accessed using
the built-in Python function :func:`help`::

   help(select)   # help on select module
   help(fetchPDB) # help on parsePDB function

This function prints documentation on screen, and you will need to type ``q``
to exit from help view.  If you are using the interactive Python shell
(IPython), you can also get help using ``?``:

.. ipython::

   In [1]: fetchPDB ?

Searching documentation
^^^^^^^^^^^^^^^^^^^^^^^

You can search entire documentation, including manual and tutorial pages,
by typing in a keyword, function, or class name.  Try searching for
*selections* to get to :ref:`selections`, for example.


.. image:: /_static/toolbox.png
   :align: center
   :alt: Searching ProDy documentation

Copying code snippets
^^^^^^^^^^^^^^^^^^^^^

When reading online documentation, you can use :guilabel:`Show code`
button on the right hand side panel to display only code snippets.
From this view, you can copy code directly into a file, i.e. click
:guilabel:`Select` and then :kbd:`Ctrl+C` to have the text in your clipboard.
To return to the documentation click the :guilabel:`Close` button.

.. image:: /_static/showcode.png
   :align: center
   :alt: Showing code examples



.. _IPython: http://ipython.scipy.org/
.. _Numpy: http://www.numpy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net

