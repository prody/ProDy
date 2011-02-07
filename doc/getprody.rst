.. _getprody:

*******************************************************************************
Getting ProDy
*******************************************************************************

Easy Install
===============================================================================

If you have `Easy Install <http://peak.telecommunity.com/DevCenter/EasyInstall>`_
installed, following will work for you::

  easy_install -U ProDy

If you don't have, please download an installation file and follow the 
instructions.

Download Files
===============================================================================

Download a suitable ProDy installation file from here

.. csv-table::
   :file: _static/pypi_downloads.csv
   :header-rows: 1
   :delim: ;
   :widths: 40, 20, 20, 10, 10

or obtain them from 
`Python Package Index <http://pypi.python.org/pypi/ProDy#downloads>`_.

See :ref:`pypi-statistics`.

Required Software
===============================================================================

* `Python 2.6 or 2.7 <http://python.org/>`_ (for Windows, prefer 32 bit option 
  for compatibility)
* `Numpy 1.3+ <http://numpy.scipy.org/>`_
* `Biopython 1.54+ <http://biopython.org/>`_

Installation Instructions
===============================================================================

After installing required packages, you will need to do one of the following:

Linux
-------------------------------------------------------------------------------

Download :file:`ProDy-0.x.y.tar.gz`. Extract tarball contents and run 
:file:`setup.py` as follows::

    tar -xzf ProDy-0.x.y.tar.gz
    cd ProDy-0.x.y
    python setup.py install

You may need root access for the last step.
  
If you don't have root access, you can edit ``PYTHONPATH`` system variable to 
recognize the path to this package.
  
e.g. move :file:`pypdb` directory to a folder like :file:`/home/username/mypackages/`
then add a line to your .bashrc (or alike) script as follows
``export PYTHONPATH=$PYTHONPATH:/home/username/mypackages/``

Windows
-------------------------------------------------------------------------------

Download :file:`ProDy-0.x.y.Windows.exe` and run it to install ProDy.

Optional Software
===============================================================================

* `Matplotlib <http://matplotlib.sourceforge.net/>`_ is required for plotting 
  data.
* `Scipy 0.7+ <http://www.scipy.org/SciPy>`_, when installed, replaces
  linear algebra module of Numpy. Scipy linear algebra module is more flexible
  and can be faster depending on the situation.
* `NMWiz <http://code.google.com/p/nmwiz/>`_ for visualizing normal mode data 
  in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.
  
..
  * `Biopython 1.54+ <http://biopython.org/wiki/Main_Page>`_ required for 
    pairwise 
    sequence alignments and proximity based atom selections. Also, when 
    installed, Bio.KDTree is used in elastic network model calculations. It
    provides significant speed up when building Hessian (ANM) or Kirchoff (GNM) 
    matrices for large systems.

Recommended Software
===============================================================================


* `MDAnalysis <http://code.google.com/p/mdanalysis/>`_ for analyzing molecular 
  dynamics trajectories.
* `IPython <http://ipython.scipy.org/>`_ for interactive ProDy sessions.
* `PyReadline <http://ipython.scipy.org/moin/PyReadline/Intro>`_ for 
  colorful interactive ProDy sessions on Windows.

Source Code
===============================================================================

The source code can be found at http://github.com/abakan/ProDy.

.. Comes with ProDy
   ===============================================================================
   The following software is included in ProDy packages:
   * `Pyparsing <http://pyparsing.wikispaces.com/>`_ is used for the atom
     selection grammer.
   * `Biopython <http://biopython.org/>`_ KDTree, pairwise2, and Blast modules
     are included in ProDy packages.

