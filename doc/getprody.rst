.. _getprody:

*******************************************************************************
Getting ProDy
*******************************************************************************

ProDy is a Python package. Dependencies, which include a Python interpreter,
need to be installed before you can use it. In addition, installing the 
complementary/recommended software improves user's experience with ProDy.

Download
===============================================================================

ProDy is available on the `Python Package Index <http://pypi.python.org/pypi/ProDy>`_.

Dependencies
===============================================================================

**Required for all functionality:**

* `Python 2.6 <http://python.org/>`_ (for Windows, please choose 32 bit option)
* `Numpy 1.3+ <http://numpy.scipy.org/>`_

With NumPy, you can use the PDB parser and the powerful atom selector. For
dynamics analysis, you need to have the following as well.

**Required for certain functionality:**

* `Scipy 0.7+ <http://www.scipy.org/SciPy>`_ for NMA calculations.
* `Biopython 1.54+ <http://biopython.org/wiki/Main_Page>`_ for ENM calculations, 
  sequence alignment, and proximity based atom selections.
* `Matplotlib <http://matplotlib.sourceforge.net/>`_ for plotting normal mode 
  data.

Complementary Software
===============================================================================

* `NMWiz <http://code.google.com/p/nmwiz/>`_ for visualizing normal mode data 
  in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.

Recommended Software
===============================================================================

* `MDAnalysis <http://code.google.com/p/mdanalysis/>`_ for analyzing molecular 
  dynamics trajectories.
* `IPython <http://ipython.scipy.org/>`_ for interactive ProDy sessions.
* `PyReadline <http://ipython.scipy.org/moin/PyReadline/Intro>`_ for 
  colorful interactive ProDy sessions on Windows.

Installation
===============================================================================

After installing required packages, you will need to do one of the following:

**Easy Install**

If you have `Easy Install <http://peak.telecommunity.com/DevCenter/EasyInstall>`_
installed, following will work for you::

  easy_install -U ProDy

**Linux**

Download :file:`prody-0.x.y.tar.gz`. Extract tarball contents and run 
:file:`setup.py` as follows::

    tar -xzf pypdb-0.x.y.tar.gz
    cd pypdb-0.x.y
    python setup.py install

You may need root access for the last step.
  
If you don't have root access, you can edit ``PYTHONPATH`` system variable to 
recognize the path to this package.
  
e.g. move :file:`pypdb` directory to a folder like :file:`/home/username/mypackages/`
then add a line to your .bashrc (or alike) script as follows
``export PYTHONPATH=$PYTHONPATH:/home/username/mypackages/``

**Windows**

Download :file:`prody-0.x.y.windows.exe` and run it to install ProDy.

Source Code
===============================================================================

The source code can be found in a git repository, at 
http://github.com/abakan/ProDy/.
