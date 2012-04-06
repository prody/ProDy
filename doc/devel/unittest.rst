.. _unittest:

.. currentmodule:: prody

*******************************************************************************
Unit Test Development
*******************************************************************************

ProDy releases have been tested using the code snippets in the documentation
of the package. Close to 1000 tests in the documentation files (:file:`.rst`)
and source code (:file:`.py`) are run all at once using Sphinx. Although
we have managed to make progress thus far with the advantage of having the 
test code as example for users, development of a comprehensive unit test suite
has become an urgent necessity. The long term aim is to develop tests
for every ProDy class and function. The ideal is to use the standard library 
module :mod:`unittest` for the development for now, but if limitations
arise other compatible options may be considered and incorporated. 

For now, unit test development should follow these guidelines:

  #. All tests for functions and classes in a ProDy module should be in a 
     single test file named after the module, e.g. :file:`test_proteins.py`.
     
  #. All test files should be stored in :file:`tests` folder in the ProDy 
     package directory, i.e. :file:`prody`
     
  #. If a test is parsing a file from :file:`tests/data` folder, it should
     be able to find those files when it is run from any folder.
