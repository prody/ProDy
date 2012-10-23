.. _unittest:

.. currentmodule:: prody

*******************************************************************************
Testing ProDy
*******************************************************************************

Running Unittests
===============================================================================

The easiest way to run ProDy unittests is using :ref:`prody-test` command::

  $ prody test
  $ prody test atomic.select -l full
  
See :ref:`prody-test` documentation for details.

Alternatively, you can use :func:`prody.test` function in a Python session::

  prody.test()
  prody.test('atomic.select', label='full')

Unittest Development
===============================================================================

Unit test development should follow these guidelines:

  #. For comparing Python numerical types and objects, e.g. int, list, tuple,
     use methods of :class:`unittest.TestCase`.

  #. For comparing Numpy arrays, use assertions available in 
     :mod:`numpy.testing` module.

  #. All tests for functions and classes in a ProDy module should be in a 
     single test file named after the module, 
     e.g. :file:`test_atomic/test_select.py`.

  #. All test files should be stored in :file:`tests` folder in the ProDy 
     package directory, i.e. :file:`prody/tests/`

  #. If a test is parsing a file from :file:`tests/test_datafiles` folder, it 
     should be able to find those files when it is run from any folder.

