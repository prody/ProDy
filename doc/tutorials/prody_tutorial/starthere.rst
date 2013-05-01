How to Start
===============================================================================

Using ProDy
-------------------------------------------------------------------------------

ProDy can be used in a number of different ways:

  #. interactively in a Python shell,
  #. as a command line program via :ref:`prody-apps`,
  #. from within VMD via :ref:`nmwiz`,
  #. or as a toolkit for developing new software.


Python for beginners
^^^^^^^^^^^^^^^^^^^^

Familiarity with Python programming language will help when using ProDy.  If
you are new to Python, or to programming, you may also start with one of the
following tutorials:

  * `The Python Tutorial <http://docs.python.org/tutorial/>`_
  * `Python Scientific Lecture Notes <http://scipy-lectures.github.com/>`_
  * `A Primer on Python for Life Science Researchers
    <http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.0030199>`_


Interactive Usage
-------------------------------------------------------------------------------

One of our aims is making ProDy suitable for interactive usage by designing 
flexible functions and classes and giving them easy to remember names in a 
consistent manner. 

For best interactive usage experience, we strongly recommend that you use
`IPython`_ instead of the standard Python shell. The IPython shell provides
user-friendly features, such as dynamic introspection and help, and also
optionally convenient integration of `Numpy`_ and `Matplotlib`_.


.. _IPython: http://ipython.scipy.org/
.. _Numpy: http://www.numpy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net


Start Python shell
^^^^^^^^^^^^^^^^^^

If you have installed IPython, type in::

  $ ipython

If you also installed Matplotlib, use::

  $ ipython --pylab

``--pylab`` option will import Matplotlib and Numpy automatically.

If you don't have IPython yet, use::

  $ python


On Windows, you will need to run :program:`IDLE (Python GUI)` or 
:program:`Python (command line)`.

Using Documentation
-------------------------------------------------------------------------------

ProDy documentation is quite comprehensive and you can access it in a number of
different ways.  In interactive sessions, API reference can be accessed using
the built-in Python function :func:`help`::

  help(atomic) # help on atomic module
  help(select) # help on select module
  help(parsePDB) # help on parsePDB function

This function prints documentation on screen, and you will need to type ``q``
to exit from help view.  If you are using the interactive Python shell
(IPython), you can also receive help by typing::
  
  select ?

Copying code snippets
^^^^^^^^^^^^^^^^^^^^^

When reading online documentation, you can use :guilabel:`Show Code Snippets`
button on the right hand side panel to display only code snippets.  From this
view, you can copy code directly into a file, i.e. click on the text, press
:kbd:`Ctrl+A` and then :kbd:`Ctrl+C` to have the text in your clipboard. To 
return to the documentation click the :guilabel:`Show documentation` button
at the top.

.. image:: /_static/codesnippets.png
   :align: center
   :alt: Getting ProDy code snippets.
   :scale: 80 %
