.. _wininst:

Making Windows Installers
=========================

`MinGW <http://www.mingw.org/>`_ (for 32-bit system) or `MinGW-w64 <http://mingw-w64.org/>`_ 
(for 64-bit system) can be used for compiling C modules when making Windows installers.  
Please follow the `instructions <https://wiki.python.org/moin/WindowsCompilers>`_ 
to install and configure them. 

`libpython <https://anaconda.org/anaconda/libpython>`_ is also required to be installed. 
If the compiler complains, such as, "'::hypot' has not been declared", please refer to this 
`link <https://stackoverflow.com/questions/10660524/error-building-boost-1-49-0-with-gcc-4-7-0>`_.
