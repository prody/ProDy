.. _getprody:

*******************************************************************************
Getting ProDy
*******************************************************************************

Required Software
===============================================================================

* `Python`_ 2.6 or 2.7 (for Windows, choose a 32bit installer that is supported
  by NumPy and ProDy packages)
* `NumPy`_ 1.4+

When compiling from source, on Linux for example, you will need a C compiler 
(e.g. :program:`gcc`) and Python developer libraries (i.e. :file:`python.h`).  
If you don't have Python developer libraries installed on your machine,
use your package manager to install :file:`python-dev` package. 

In addition, `matplotlib`_ is required for using plotting functions.  
ProDy and :ref:`commands` can be operated without this package.   

.. _Python: http://python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

Easy Install
===============================================================================

If you have `Easy Install`_ installed, type the following::

  easy_install -U ProDy

If you don't have Easy Install, please download an installation file and 
follow the instructions.

.. _Easy Install: "http://peak.telecommunity.com/DevCenter/EasyInstall

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

Windows installers are built using `MinGW`_.

.. _MinGW: http://www.mingw.org/

Installation Instructions
===============================================================================

After installing the required packages, you will need to do one of the 
following:

Linux
-------------------------------------------------------------------------------

Download :file:`ProDy-x.{y}.{z}.tar.gz`. Extract tarball contents and run 
:file:`setup.py` as follows::

    $ tar -xzf ProDy-x.y.z.tar.gz
    $ cd ProDy-x.y.z
    $ python setup.py build
    $ python setup.py install

You may need root access installation, i.e. ``sudo python setup.py install``.  

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

To be able use :ref:`commands`, append Python (e.g. :file:`C:\Python27`) and 
scripts (e.g. :file:`C:\Python27\Scripts`) folders to :envvar:`PATH` variable.

Plugins
===============================================================================

`NMWiz`_ is a `VMD`_ plugin for visualization and animation of normal mode data 
and comes with `VMD package`_.  See :ref:`nmwiz-tutorial` and :ref:`nmwiz` for
further details. 

.. _NMWiz: http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/
.. _VMD package: http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

Recommended Software
===============================================================================

* `Scipy`_, when installed, replaces linear algebra module of Numpy. 
  The Scipy linear algebra module is more flexible and can be faster.
* `IPython`_ is a must have for interactive ProDy sessions.
* `PyReadline`_ for colorful interactive ProDy sessions on Windows.
* `MDAnalysis`_ for analyzing molecular dynamics trajectories.

.. _Scipy: http://www.scipy.org/SciPy
.. _IPython: http://pypi.python.org/pypi/ipython
.. _PyReadline: http://pypi.python.org/pypi/pyreadline
.. _MDAnalysis: http://code.google.com/p/mdanalysis/

Included in ProDy Package
===============================================================================
The following software is included in the ProDy installation packages:

* `Pyparsing`_ is used to define the atom selection grammar.

* `Biopython`_ KDTree package and pairwise2 module are used for distance based
  atom selections and pairwise sequence alignment, respectively. 

* `argparse`_ is used to implement new :ref:`commands` and provided for 
  compatibility with Python 2.6.

.. _Pyparsing: http://pyparsing.wikispaces.com/
.. _Biopython: http://biopython.org/
.. _argparse: http://code.google.com/p/argparse/

Source Code
===============================================================================

The source code can be found at http://github.com/abakan/ProDy.

