.. _testing:

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


Pre-commit Testing 
===============================================================================

You can automatize testing of ProDy package using a git pre-commit hook.  
For example, the following script calls :file:`devel_test.sh` file that comes 
in the project directory ensures that the package imports and passes fast 
subset of tests:: 

  #!/bin/sh

  SRC="$(git diff --cached --name-only | grep -E "(\.py|\.c)$")"

  if [ "$SRC" ]
  then
      TEST="$(/home/abakan/Code/ProDy/devel_test.sh 3>&1 1>&2 2>&3)"
      echo "$TEST" >&2
      FAIL="$(echo $TEST | grep FAILED)"
      if [ "$FAIL" ]
      then
		  echo "ProDy unittests failed." >&2
		  exit 1
      fi
  fi


This script needs to be saved in :file:`.git/hooks/pre-commit` executable file.


Unittest Development
===============================================================================

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

