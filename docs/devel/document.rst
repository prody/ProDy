.. _document:

.. currentmodule:: prody

Documenting ProDy
=================

.. contents::
   :local:


ProDy documentation is written using reStructuredText_ markup and prepared
using Sphinx_.  You may install Sphinx using easy_install, i.e.
``easy_install -U Sphinx``, or using package manager on your Linux machine.

Building Documentation
----------------------

ProDy HTML documentation can be build as follows::

  $ cd doc
  $ make clean
  $ make html

If all documentation strings and pages are properly formatted according to
reStructuredText_ markup, documentation pages should compile without any
warnings.


**New release**


When making a new release, also build PDF documentation using the following
command::

  $ make latexpdf

To be able to run this command, you will need :program:`latex` installed on
your machine.

**Development version**


For documenting work in progress, increment version number in
:file:`doc/conf.py` and append ``-dev``, i.e. ``version = '1.3.1-dev'``.
HTML documentation for development version should go in :file:`latest` folder
in where the documentation is hosted.


.. _Sphinx: http://sphinx.pocoo.org/
.. _reStructuredText: http://docutils.sf.net/rst.html
.. _Google Analytics: http://www.google.com/analytics/