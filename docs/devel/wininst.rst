.. _wininst:

Making Windows Installers
=========================

`MinGW <http://www.mingw.org/>`_ can be used for compiling C modules when
making Windows installers.  Install MinGW and make :file:`distutils.cfg` file
in :file:`PythonXY\\Lib\\distutils` folder that contains::

  [build]
  compiler = mingw32
