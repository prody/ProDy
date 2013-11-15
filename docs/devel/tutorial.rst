.. _write-tutorial:

.. currentmodule:: prody

Writing Tutorials
=================

.. contents::
   :local:


This is a short guide for writing ProDy tutorials that are published as part
of online documentation pages, and also as individual downloadable PDF files.

Tutorial Setup
--------------

First go to :file:`doc` folder in ProDy package and generate necessary files
for your tutorial using :program:`start-tutorial.sh` script::

  $ cd doc
  $ ./start-tutorial.sh
  Enter tutorial title: ENM Analysis using ProDy
  Enter a short title: ENM Analysis
  Enter author name: Ahmet Bakan

  Tutorial folders and files are prepared, see tutorials/enm_analysis

This will generate following folder and files::

  $ cd tutorials/enm_analysis/
  $ ls -lgo
  -rw-r--r-- 1  328 Apr 30 16:48 conf.py
  -rw-r--r-- 1  395 Apr 30 16:48 index.rst
  -rw-r--r-- 1  882 Apr 30 16:48 intro.rst
  -rw-r--r-- 1 1466 Apr 30 16:48 Makefile
  lrwxrwxrwx 1   13 Apr 30 16:48 _static -> ../../_static


Note that short title will be used as filename and part of the URL of the
online documentation pages.

If tutorial logo/image that you want to use is different from ProDy logo,
update the following line in :file:`conf.py`::

  tutorial_logo = u'enm.png'     # default is ProDy logo
  tutorial_prody_version = u''        # default is latest ProDy version

Also, note ProDy version if the tutorial is developed for a specific release.


Style and Organization
----------------------

ProDy documentation and tutorials are written using `reStructuredText`_,
an easy-to-read/write file format.  See `reStructuredText Primer`_ for a
quick introduction.

reStructuredText is stored in plain-text files with :file:`.rst` extension,
and converted to HTML and PDF pages using `Sphinx`_.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _reStructuredText Primer: http://sphinx-doc.org/rest.html
.. _Sphinx: http://sphinx-doc.org/

:file:`index.rst` and :file:`intro.rst` files are automatically generated.
:file:`index.rst` file should include title and table of contents of the
tutorial.  Table of contents is just a list of :file:`.rst` files that are
part of the tutorial.  They be listed in the order that they should appear
in the final PDF file::

  .. _enm-analysis:

  .. use "enm-analysis" to refer to this file, i.e. :ref:`enm-analysis`

  *******************************************************************************
  ENM Analysis using ProDy
  *******************************************************************************


  .. add .rst files to `toctree` in the order that you want them

  .. toctree::
     :glob:
     :maxdepth: 2

     intro


Add more :file:`.rst` files as needed.  See other tutorials in
:file:`doc/tutorials` folder as examples.


Input/Output Files
------------------

All files needed to follow the tutorial should be stored in
:file:`tutorial_name_files` folder.  There is usually no need to provide
PDB files, as ProDy automatically downloads them when needed.  Optionally,
output files can also be provided.

.. note::
 	 Small input and output files that contain textual information may
 	 be included in the :program:`git` repository, but please avoid including
 	 large files in particular those that contain binary data.


Including Code
--------------

Python code in tutorials should be included using `IPython Sphinx directive`_.
In the beginning of each :file:`.rst` file, you should make necessary imports
as follows::

  .. ipython:: python

     from prody import *
     from matplotlib.pylab import *
     ion()

This will convert to the following:


.. ipython:: python

   from prody import *
   from matplotlib.pylab import *
   ion()


Then you can add the code for the tutorial::


  .. ipython:: python

     pdb = parsePDB('1p38')

.. ipython:: python

   pdb = parsePDB('1p38')



.. _IPython Sphinx Directive: http://ipython.org/ipython-doc/dev/development/ipython_directive.html

Including Figures
-----------------

IPython directive should also be used for including figures::



  .. ipython:: python

     @savefig tutorial_name_figure_name.png width=4in
     plot(range(10))

     @savefig tutorial_name_figure_two.png width=4in
     plot(range(100)); # used ; to suppress output

``@savefig`` decorator was used to save the figure.

.. note::

   Figure names needs to be unique within the tutorial and should be prefixed
   with the tutorial name.


Note that in the second :func:`~matplotlib.pyplot.plot` call, we used a
semicolon to suppress the output of the function.


If you want to make modifications to the figure, save it after the last
modification::


  .. ipython:: python

     plot(range(10));
     grid();
     xlabel('X-axis')
     @savefig tutorial_name_figure_three.png width=4in
     ylabel('Y-axis')


Testing Code
------------

If there is any particular code output that you want to test, you can use
``@doctest`` decorator as follows::


  .. ipython::

     @doctest
     In [1]: 2 + 2
     Out[1]: 4


.. ipython::

   @doctest
   In [1]: 2 + 2
   Out[1]: 4

Failing to produce the correct output will prevent building the documentation.


Publishing Tutorial
-------------------

To see how your :file:`.rst` files convert to HTML format, use the following
command::

  $ make html

You will find HTML files in :file:`_build/html` folder.

Once your tutorial is complete and looks good in HTML (no code execution
problems), following commands can be used to generate a PDF file and
tutorial file achieves::

  $ make pdf
  $ make files

ProDy online documentation will contain these files as well as tutorial pages
in HTML format.
