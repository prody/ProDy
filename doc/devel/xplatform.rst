.. _xplatform:

.. currentmodule:: prody

*******************************************************************************
Cross-platform issues
*******************************************************************************

This section describes cross-platform issues that may emerge and provides 
possible solutions for them.

Numpy integer type
===============================================================================

Issues may arise when comparing Numpy integer types with Python :func:`int`.
Python :func:`int` equivalent Numpy integer type on Windows (Win7 64bit, 
Python 32bit) is :class:`~numpy.int32`, while on Linux (Ubuntu 64bit) it is 
:class:`~numpy.int64`.  For example, the statement 
``isinstance(np.array([1], np.int64), int)`` may return ``False`` resulting
in unexpected behavior in ProDy functions or methods.  If Numpy integer type 
needs to be specified, using :class:`~numpy.int` seems a safe option.  

Relative paths
===============================================================================

:func:`os.path.relpath` function raises exceptions when the working 
directory and the path of interest are on separate drives, e.g. trying 
to write a :file:`C:\\temp` while running tests on :file:`D:\\ProDy`.
Instead of this :func:`os.path.relpath`, ProDy function :func:`prody.relpath`
should be used to avoid problems. 
