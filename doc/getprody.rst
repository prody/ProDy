.. _getprody:

*******************************************************************************
Getting Prody
*******************************************************************************

ProDy is a Python package. Dependencies, which include a Python interpreter,
need to be installed before you can use it. In addition, installing the 
recommended software greatly improves user's experience with ProDy.

**Steps to get ProDy running**:

  1. Install dependencies.
  2. Install recommended software (optional)
  3. Install ProDy.

Dependencies
===============================================================================

**Required for all functionality:**

* `Python 2.6 <http://python.org/>`_ (On windows, please choose 32 bit option)
* `Numpy 1.3+ <http://numpy.scipy.org/>`_

**Required for certain functionality:**

* `Scipy 0.7+ <http://www.scipy.org/SciPy>`_ for NMA calculations.
* `Biopython 1.54+ <http://biopython.org/wiki/Main_Page>`_ for ENM calculations 
   and proximity based atom selections.
* `Matplotlib <http://matplotlib.sourceforge.net/>`_ for plotting normal mode 
   data.

Complementary Software
===============================================================================

* `NMWiz <http://code.google.com/p/nmwiz/>`_ for visualizing normal mode data 
  in `VMD <http://www.ks.uiuc.edu/Research/vmd/>`_.

Recommended Software
===============================================================================

* `IPython <http://ipython.scipy.org/>`_ for interactive ProDy sessions.
* `PyReadline <http://ipython.scipy.org/moin/PyReadline/Intro>`_ for 
   interactive ProDy sessions on Windows.

Related Software
===============================================================================

* `MDAnalysis <http://code.google.com/p/mdanalysis/>`_ for analyzing molecular dynamics trajectories.


Installation
===============================================================================

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
