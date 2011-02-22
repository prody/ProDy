.. _wininst:

.. currentmodule:: prody

*******************************************************************************
Creating Windows Installers
*******************************************************************************

When creating Windows installers http://www.mingw.org/ can be 
used for compiling C modules. In :file:`Lib\distutils` folder, make 
:file:`distutils.cfg` file that contains::

  [build]
  compiler = mingw32
