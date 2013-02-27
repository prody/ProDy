.. _getprody:

*******************************************************************************
Getting ProDy
*******************************************************************************

Required Software
===============================================================================

* `Python`_ 2.7 or later
  
  *Windows*: Note that, NumPy and ProDy supports only **32-bit** Python 
  on Windows.
  
  *Python 2.6 and 3.1* can be used to install and run ProDy. However, note 
  that ProDy tests are performed using Python 2.7 and 3.2, so these are the 
  recommended versions. 
  
* `NumPy`_ 1.4+

When compiling from source, on Linux for example, you will need a C compiler 
(e.g. :program:`gcc`) and Python developer libraries (i.e. :file:`python.h`).  
If you don't have Python developer libraries installed on your machine,
use your package manager to install :file:`python-dev` package. 

In addition, `matplotlib`_ is required for using plotting functions.  
ProDy, :ref:`prody-apps`, and :ref:`evol-apps` can be operated without 
this package.   

.. _Python: http://www.python.org/download/
.. _NumPy: http://sourceforge.net/projects/numpy/files/NumPy/
.. _matplotlib: https://github.com/matplotlib/matplotlib/downloads


Easy Install
===============================================================================

If you have `Easy Install`_ installed, type the following::

  easy_install -U ProDy

If you don't have Easy Install, please download an installation file and 
follow the instructions.

.. _Easy Install: http://peak.telecommunity.com/DevCenter/EasyInstall


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

Remove all previously installed ProDy files.  You can find location of 
installation files as follows::

    $ python -c "import prody; print(prody.__path__)"

Download :file:`ProDy-x.{y}.{z}.tar.gz`.  Extract tarball contents and run 
:file:`setup.py` as follows::

    $ tar -xzf ProDy-x.y.z.tar.gz
    $ cd ProDy-x.y.z
    $ python setup.py build
    $ python setup.py install

You may need root access installation, i.e. ``sudo python setup.py install``.  

If you don't have root access, you can edit :envvar:`PYTHONPATH` environment 
variable to specify the path to this package:
  
#. Move :file:`prody` directory from :file:`build/lib.linux-x86_{bb}-2.{z}` (or 
   similar) to a folder like :file:`/home/username/mypackages/`
#. Add a line to your :file:`.bashrc` (or similar) script as follows
   ``export PYTHONPATH=$PYTHONPATH:/home/username/mypackages/``

You may also consult other alternate and custom installation schemes in
`Installing Python Modules <http://docs.python.org/install/index.html>`_.

Mac OS
-------------------------------------------------------------------------------

For installing ProDy, please follow the Linux installation instructions.

Windows
-------------------------------------------------------------------------------

Remove previously installed ProDy release from :program:`Uninstall a program` 
in :guilabel:`Control Panel`.
 
Download :file:`ProDy-0.{x}.{y}.win32-py2.{z}.exe` and run to install ProDy.

To be able use :ref:`prody-apps` and :ref:`evol-apps` in command prompt 
(:program:`cmd.exe`), append Python and scripts folders (e.g. 
:file:`C:\\Python27` and :file:`C:\\Python27\\Scripts`) to :envvar:`PATH` 
environment variable.

Testing
-------------------------------------------------------------------------------

You can test your ProDy installation using the following command::

    $ prody test

Note that :program:`prody` script/command must be reachable from your working
directory.  For more information on testing, see :ref:`testing`.


NMWiz
===============================================================================

:ref:`nmwiz` is a `VMD`_ plugin for comparative visual analysis of protein
dynamics modeled using theory or inferred from experimental structural
ensembles.  See :ref:`nmwiz` for available updates and installation 
instructions. 

.. _NMWiz: http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/
.. _VMD: http://www.ks.uiuc.edu/Research/vmd/


Recommended Software
===============================================================================

* `Scipy`_, when installed, replaces linear algebra module of Numpy. 
  Scipy linear algebra module is more flexible and can be faster.
* `IPython`_ is a must have for interactive ProDy sessions.
* `PyReadline`_ for colorful interactive ProDy sessions on Windows.
* `MDAnalysis`_ for analyzing molecular dynamics trajectories.

.. _Scipy: http://www.scipy.org/SciPy
.. _IPython: http://pypi.python.org/pypi/ipython
.. _PyReadline: http://pypi.python.org/pypi/pyreadline
.. _MDAnalysis: http://code.google.com/p/mdanalysis/


Included in ProDy Package
===============================================================================

Following software is included in the ProDy installation packages:

* `Pyparsing`_ is used to define the atom selection grammar.

* `Biopython`_ KDTree package and pairwise2 module are used for distance based
  atom selections and pairwise sequence alignment, respectively. 

* `argparse`_ is used to implement applications and provided for 
  compatibility with Python 2.6.

.. _Pyparsing: http://pyparsing.wikispaces.com/
.. _Biopython: http://biopython.org/
.. _argparse: http://code.google.com/p/argparse/


Source Code
===============================================================================

The source code can be found at https://bitbucket.org/abakan/prody.

