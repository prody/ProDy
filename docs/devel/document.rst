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
versions.


Building Website
----------------

ProDy-website source is hosted at https://github.com/prody/ProDy-website
This project contains tutorial files and the home pages for other software
in the


.. _Sphinx: http://sphinx.pocoo.org/
.. _reStructuredText: http://docutils.sf.net/rst.html
.. _Google Analytics: http://www.google.com/analytics/