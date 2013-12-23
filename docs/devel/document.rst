.. _document:

.. currentmodule:: prody

Documenting ProDy
=================

.. contents::
   :local:


ProDy documentation is written using reStructuredText_ markup and prepared
using Sphinx_.  You may install Sphinx using :program:`easy_install`, i.e.
``easy_install -U Sphinx``, or using package manager on your Linux machine.

Building Manual
---------------

ProDy Manual in HTML and PDF formats can be build as follows::

  $ cd docs
  $ make html
  $ make pdf

If all documentation strings and pages are properly formatted according to
reStructuredText_ markup, documentation pages should compile without any
warnings. Note that to build PDF files, you need to install :program:`latex`
and :program:`pdflatex` programs.


**Read the Docs**

A copy of ProDy manual is hosted on `Read the Docs <https://readthedocs.org/>`_
and can be viewed at http://prody.readthedocs.org/. Read the Docs is configured
to build manual pages for the ``devel`` branch (latest) and the recent stable
versions. The user name for Read the Docs is ``prody``.


Building Website
----------------

ProDy-website source is hosted at https://github.com/prody/ProDy-website
This project contains tutorial files and the home pages for ProDy and other
related software.

**Latest version**


To build website on ProDy server, start with pulling changes::

  $ cd ProDy-website
  $ git pull

Running the following command will build HTML pages for the latest stable
release of ProDy::

  $ make html

HTML pages for manual and all tutorials are build as a single project,
which allows for referencing from manual to tutorials.

PDF files for the manual and tutorials, and also download files are build
as follows::

  $ make pdf

PDF and TGZ/ZIP files are copied to appropriate places after they are built.

**Development version**

Finally, HTML and PDF pages for the development version can be built as
follows::

  $ make devel

Again, this will copy HTML and PDF files to appropriate places, and a link
to these files will be provided from he homepage.


.. _Sphinx: http://sphinx.pocoo.org/
.. _reStructuredText: http://docutils.sf.net/rst.html
.. _Google Analytics: http://www.google.com/analytics/