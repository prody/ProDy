.. _prody-test:

*******************************************************************************
prody test
*******************************************************************************

Usage
===============================================================================

Running :command:`prody test -h` displays::

  usage: prody test [-h] [--quiet] [--examples] [-l STR] [-v INT]
                    [mod [mod ...]]
  
  optional arguments:
    -h, --help            show this help message and exit
    --quiet               suppress info messages to stderr
    --examples            show usage examples and exit
  
  output options:
    -l STR, --label STR   label used when nose is available, fast (default) or
                          full
    -v INT, --verbose INT
                          verbosity level, default is 1
    mod                   ProDy module names

Examples
===============================================================================

Running :command:`prody test --examples` displays::

  ProDy unittests can be run as follows:
  
      $ prody test
  
  To run all tests, use -l/--label argument as follows:
  
      $ prody test -l full
  
  To increase verbosity, use -v/--verbose argument as follows:
  
      $ prody test -l fast -v 2
  
  To run tests for a specific module, pass module name:
  
      $ prody test proteins
      $ prody test atomic.select
