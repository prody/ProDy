.. _unittest:

.. currentmodule:: prody

*******************************************************************************
Testing ProDy
*******************************************************************************

Running Unittest
===============================================================================

The easiest way to rung ProDy unit tests is using :ref:`prody-test` command::

  $ prody test
  $ prody test atomic.select -l full
  
See :ref:`prody-test` documentation for details.

Alternatively, you can run use :func:`prody.test` function in a Python 
session::

  prody.test()
  prody.test('atomic.select', label='full')

Unit Test Development
===============================================================================

Unit test development should follow these guidelines:

  #. For comparing Python numerical types and objects, e.g. int, list, tuple,
     use methods of :class:`unittest.Testcase`.

  #. For comparing Numpy arrays, use assertions available in 
    :mod:`numpy.testing` module.

  #. All tests for functions and classes in a ProDy module should be in a 
     single test file named after the module, e.g. :file:`test_proteins.py`.

  #. All test files should be stored in :file:`tests` folder in the ProDy 
     package directory, i.e. :file:`prody`

  #. If a test is parsing a file from :file:`tests/test_datafiles` folder, it 
     should be able to find those files when it is run from any folder.

