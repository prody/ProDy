.. _getprody:

*******************************************************************************
Getting ProDy
*******************************************************************************

Required Software
===============================================================================

* `Python 2.6 or 2.7 <http://python.org/>`_ (for Windows, choose a 32bit 
  installer, which is supported by NumPy and ProDy packages)
* `NumPy 1.4+ <http://numpy.scipy.org/>`_

When compiling from source, on Linux for example, you will need a C compiler 
(e.g. :program:`gcc`) and Python developer libraries (i.e. :file:`python.h`).  
If you don't have Python developer libraries installed on your machine,
use your package manager to install :file:`python-dev` package. 

In addition, `Matplotlib <http://matplotlib.sourceforge.net/>`_ is required
for using plotting functions. ProDy and :ref:`commands` can be operated without
this package.   


Easy Install
===============================================================================

If you have `Easy Install <http://peak.telecommunity.com/DevCenter/EasyInstall>`_
installed, type the following::

  easy_install -U ProDy

If you don't have Easy Install, please download an installation file and 
follow the instructions.

Download Files
===============================================================================

Download a suitable ProDy installation file from http://python.org/pypi/ProDy 
or here:

.. csv-table::
   :file: _static/ProDy_downloads.csv
   :header-rows: 1
   :delim: ,
   :widths: 40, 20, 20, 10, 10

For details of ProDy releases see :ref:`changes` and :ref:`pypi-statistics`.

Windows installers are built using `MinGW <http://www.mingw.org/>`_.

Installation Instructions
===============================================================================

After installing the required packages, you will need to do one of the 
following:

Linux
-------------------------------------------------------------------------------

Download :file:`ProDy-0.{x}.{y}.tar.gz`. Extract tarball contents and run 
:file:`setup.py` as follows::

    $ tar -xzf ProDy-0.x.y.tar.gz
    $ cd ProDy-0.x.y
    $ python setup.py build
    $ python setup.py install

You may need root access for the last step.  

If you don't have root access, you can edit :envvar:`PYTHONPATH` system 
variable to specify the path to this package:
  
#. Move :file:`prody` directory from :file:`build/lib.linux-x86_{bb}-2.{z}` (or 
   similar) to a folder like :file:`/home/username/mypackages/`
#. Add a line to your :file:`.bashrc` (or similar) script as follows
   ``export PYTHONPATH=$PYTHONPATH:/home/username/mypackages/``

Mac OS
-------------------------------------------------------------------------------

For installing ProDy, please follow the Linux installation instructions.

Windows
-------------------------------------------------------------------------------

Download :file:`ProDy-0.{x}.{y}.win32-py2.{z}.exe` and run to install ProDy.

Plugins
===============================================================================

:ref:`nmwiz`, part of ProDy project, is a VMD plugin for visualization and 
animation of normal mode data. See :ref:`nmwiz` for installation instructions. 

Recommended Software
===============================================================================

* `Scipy <http://www.scipy.org/SciPy>`_, when installed, replaces
  linear algebra module of Numpy. The Scipy linear algebra module is more 
  flexible and can be faster depending on the situation.
* `MDAnalysis <http://code.google.com/p/mdanalysis/>`_ for analyzing molecular 
  dynamics trajectories.
* `IPython <http://ipython.scipy.org/>`_ for interactive ProDy sessions.
* `PyReadline <http://ipython.scipy.org/moin/PyReadline/Intro>`_ for 
  colorful interactive ProDy sessions on Windows.


Included in ProDy Package
===============================================================================
The following software is included in the ProDy installation packages:

* `Pyparsing 1.5.6 <http://pyparsing.wikispaces.com/>`_ 

  Pyparsing is used to define the atom selection grammar.

* `Biopython 1.56 <http://biopython.org/>`_ - KDTree package and pairwise2 
  module
   
  KDTree and pairwise2 components are used for distance based atom selections 
  and pairwise sequence alignment, respectively. 


Source Code
===============================================================================

The source code can be found at http://github.com/abakan/ProDy.

