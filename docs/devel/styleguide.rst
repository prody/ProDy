.. _styleguide:

Style Guide for ProDy
=====================

Introduction
------------

:pep:`8`, the *Style Guide for Python Code*, is adopted in the development of
ProDy package.  Contributions to ProDy shall follow :pep:`8` and the
specifications and additions provided in this addendum.


Code Layout
-----------

**Indentation**

Use 4 spaces per indentation level in source code (:file:`.py`) and never use
tabs as a substitute.

In documentation files (:file:`.rst`), use 2 spaces per indentation level.


**Maximum line length**

Limit all lines to a maximum of 79 characters in both source code and
documentation files.  Exceptions may be made when tabulating data in
documentation files and strings.  The length of lines in a paragraph
may be much less than 79 characters if the line ends align better with
the first line, as in this paragraph.


**Encodings**

In cases where an encoding for a :file:`.py` file needs to be specified,
such as when characters like α, β, or Å are used in docstrings, use UTF-8
encoding, i.e. start the file with the following line::

  # -*- coding: utf-8 -*-


**Imports**

In addition to :pep:`8#imports` recommendations regarding imports, the
following should be applied:

  * relative intra-ProDy imports are discouraged, use
    ``from prody.atomic import AtomGroup`` not ``from atomic import AtomGroup``
  * always import from second top level module, use
    ``from prody.atomic import AtomGroup`` and not
    ``from prody.atomic.atomgroup import AtomGroup``,
    because file names may change or files that grow too big may be split
    into smaller modules, etc.

Here is a series of properly formatted imports following a module documentation
string::

  """This module defines a function to calculate something interesting."""

  import os.path
  from collections import defaultdict
  from time import time

  import numpy as np

  from prody.atomic import AtomGroup
  from prody.measure import calcRMSD
  from prody.tools import openFile
  from prody import LOGGER, SETTINGS

  __all__ = ['calcSomething']


Whitespace in Documentation
---------------------------

In addition to recommendations regarding whitespace use in Python code
(:pep:`8#whitespace-in-expressions-and-statements`), two whitespace
characters should follow a period in documentation files and strings
to help reading documentation in terminal windows and text editors.


Naming Conventions
------------------

ProDy naming conventions aim at making the library suitable for interactive
sessions, i.e. easy to remember and type.


**Class names**

Naming style for classes is ``CapitalizedWords`` (or ``CapWords``, or
``CamelCase``).  Abbreviations and/or truncated names should be used to
keep class names short.  Some class name examples are:

  * :class:`.ANM` for Anisotropic Network Model
  * :class:`.HierView` for Hierarchical View


**Exception names**

Prefer using a suitable standard-library exception over defining a new
one. If you absolutely need to define one, use the class naming convention.
Use the suffix "Error" for exception names, when exception is an error:

  * :exc:`.SelectionError`, the only exception defined in ProDy package


**Method and function names**

Naming style for methods and functions is ``mixedCase``, that differs from
``CapWords`` by initial lowercase character.  Starting with a lowercase
(no shift key) and using no underscore characters decreases the number of
key strokes by half in many cases in interactive sessions.

Method and function names should start with a verb, suggestive on the action,
and followed by one or two names, where the second name may start with a lower
case letter.  Some examples are :func:`.moveAtoms`, :func:`.wrapAtoms`,
:func:`.assignSecstr`, and :func:`.calcSubspaceOverlap`.

Abbreviations and/or truncated names should be used and obvious words
should be omitted to limit number of names to 20 characters.  For example,
:meth:`~.ANM.buildHessian` is preferred over :meth:`buildHessianMatrix`.
Another example is the change from using :meth:`getResidueNames` to
using :meth:`.AtomGroup.getResnames`.  In fact, this was part of a series of
major :ref:`changes` aimed at refining the library for interactive usage.

In addition, the following should be applied to enable grouping of methods and
functions based on their action and/or return value:

  * :meth:`buildSomething`: methods and functions that calculate a matrix
    should start with ``build``, e.g. :meth:`.GNM.buildKirchhoff` and
    :func:`.buildDistMatrix`
  * :meth:`calcSomething`: methods that calculate new data but does not
    necessarily return anything and especially those that take timely actions,
    should start with ``calc``, e.g. :meth:`.PCA.calcModes`
  * :meth:`getSomething`: methods, and sometimes functions, that return a copy
    of data should start with ``get``, such as :func:`.listReservedWords`
  * :meth:`setSomething`: methods, and sometimes functions, that alter internal
    data should start with ``set``


Variable Names
-------------------------------------------------------------------------------

Variable names in functions and methods should contain only lower case letters,
and may contain underscore characters to increase readability.
