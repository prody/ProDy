.. _testing:

.. currentmodule:: prody

Testing ProDy
=============

.. contents::
   :local:


Running Unittests
-----------------

The easiest way to run ProDy unit tests is using nose_. The following will
run all tests::

  $ nosetests prody

To skip tests that are slow, use the following::

  $ nosetests prody -a '!slow'

To run tests for a specific module do as follows::

  $ nosetests prody.tests.atomic prody.tests.sequence


Unittest Development
--------------------

Unit test development should follow these guidelines:

  #. For comparing Python numerical types and objects, e.g. int, list, tuple,
     use methods of :class:`unittest.TestCase`.

  #. For comparing Numpy arrays, use assertions available in
     :mod:`numpy.testing` module.

  #. All test files should be stored in :file:`tests` folder in the ProDy
     package directory, i.e. :file:`prody/tests/`

  #. All tests for functions and classes in a ProDy module should be in a
     single test file named after the module,
     e.g. :file:`test_atomic/test_select.py`.

  #. Data files for testing should be located in :file:`tests/test_datafiles`.

