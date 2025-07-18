.. _getprody:

Installation
============

Required Software
-----------------

* `Python`_ 3.10 or later. We recommend using `Anaconda`_, which provides the conda package and environment manager as well as many useful packages. 

.. _Anaconda: https://www.anaconda.com/products/individual

* `NumPy`_ 1.10 or later

* `SciPy` - we recommend that you use the latest version, but all versions should be supported.

.. _SciPy: https://sourceforge.net/projects/scipy/

* `Biopython` - we recommend that you use the latest version, but all versions should be supported.

.. _Biopython: http://biopython.org/wiki/Download/

When compiling from source, on Linux for example, you will need a C compiler
(e.g. :program:`gcc`) and Python developer libraries (i.e. :file:`python.h`).
If you don't have Python developer libraries installed on your machine,
use your package manager to install :file:`python-dev` package.

In addition, `matplotlib`_ is required for using plotting functions.
ProDy, :ref:`prody-apps`, and :ref:`evol-apps` can be operated without
this package.

Quick Install
-------------


We officially recommend installing through conda::

  conda install -c conda-forge prody


Installing From Source (not recommended)
----------------------------------------

**Linux**

Download :file:`ProDy-{x}.{y}.{z}.tar.gz`.  Extract tarball contents and run
:file:`setup.py` as follows::

    $ tar -xzf ProDy-x.y.z.tar.gz
    $ cd ProDy-x.y.z
    $ python setup.py build
    $ python setup.py install

If you need root access for installation, try ``sudo python setup.py install``.
If you don't have root access, please consult alternate and custom installation
schemes in `Installing Python Modules`_.

.. _Installing Python Modules: http://docs.python.org/install/index.html

**Mac OS**


For installing ProDy, please follow the Linux installation instructions.

Recommended Software
--------------------

* `IPython`_ is a must have for interactive ProDy sessions.
* `PyReadline`_ for colorful IPython sessions on Windows.
* `MDAnalysis`_ or `MDTraj`_ for reading molecular dynamics trajectories.


.. _PyReadline: http://ipython.org/pyreadline.html

Included in ProDy
-----------------

Following software is included in the ProDy installation packages:

* `pyparsing`_ is used to define the atom selection grammar.

* `Biopython`_ KDTree package and pairwise2 module are used for distance based
  atom selections and pairwise sequence alignment, respectively.

* `argparse`_ is used to implement applications and provided for
  compatibility with Python 2.6.

* `Scipy`_, when installed, replaces linear algebra module of Numpy.
  Scipy linear algebra module is more flexible and can be faster.

.. _argparse: http://code.google.com/p/argparse/


Source Code
-----------

Source code is available at https://github.com/prody/ProDy.