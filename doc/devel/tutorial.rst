.. _write-tutorial:

.. currentmodule:: prody

*******************************************************************************
Writing Tutorials
*******************************************************************************

This is a short guide for writing ProDy tutorials that are published as part 
of online documentation pages, and also as individual downloadable PDF files.

Tutorial Setup
===============================================================================

First go to :file:`doc` folder in ProDy package and generate necessary files 
for your tutorial.  You will need a descriptive name for the tutorial which
will be used as filename and part of the URL of online documentation pages. 
For example, a tutorial on sampling protein conformers could have the name
:file:`sampling_conformers`.  We run the following command::

	$ cd doc
	$ ./start-tutorial.sh
  Enter tutorial title: ENM Analysis using ProDy
  Enter a short title: ENM Analysis
  Enter author name: Ahmet Bakan

  Tutorial folders and files are prepared, see tutorials/enm_analysis

This will generate following folder and files::

  $ cd tutorials/sampling_conformers/
  $ ls -lgo
  -rw-r--r-- 1  328 Apr 30 16:48 conf.py
  -rw-r--r-- 1  395 Apr 30 16:48 index.rst
  -rw-r--r-- 1  882 Apr 30 16:48 intro.rst
  -rw-r--r-- 1 1466 Apr 30 16:48 Makefile
  lrwxrwxrwx 1   13 Apr 30 16:48 _static -> ../../_static
	
	
If tutorial logo/image that you want to use is different from ProDy logo, 
update the following line in :file:`conf.py`::

	tutorial_logo = u'enm.png'     # default is ProDy logo
  tutorial_prody_version = u''        # default is latest ProDy version

Also, note ProDy version if the tutorial is developed for a specific release.
  

Style and Organization
===============================================================================

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
===============================================================================

All files needed to follow the tutorial should be stored in 
:file:`tutorial_name_files` folder.  There is usually no need to provide 
PDB files, as ProDy automatically downloads them when needed.  Optionally, 
output files can also be provided.

.. note::
 	 Small input and output files that contain textual information may
 	 be included in the :program:`git` repository, but please avoid including
 	 large files in particular those that contain binary data.


Testing Code
===============================================================================

Test code in your tutorial pages, i.e. :file:`.rst` files by running the 
following command::

  $ make test
  
All code snippets should produce expected results.


Publishing Tutorial
===============================================================================

Once your tutorial is ready, following commands can be used to generate
a PDF file and input/output file achieves::

  $ make pdf
  $ make files

ProDy online documentation will contain these files as well as tutorial pages 
in HTML format.
