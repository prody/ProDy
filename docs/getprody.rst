.. _getprody:

Getting ProDy
=============

Required Software
-----------------

* `Python`_ 2.6, 2.7, 3.2 or later

  *Windows*: You need to use **32-bit** Python on Windows to be able to
  install NumPy and ProDy.

* `NumPy`_ 1.5+ (for Windows, select an installer built using a suitable
  version of NumPy)

When compiling from source, on Linux for example, you will need a C compiler
(e.g. :program:`gcc`) and Python developer libraries (i.e. :file:`python.h`).
If you don't have Python developer libraries installed on your machine,
use your package manager to install :file:`python-dev` package.

In addition, `matplotlib`_ is required for using plotting functions.
ProDy, :ref:`prody-apps`, and :ref:`evol-apps` can be operated without
this package.

Quick Install
-------------

If you have pip_ installed, type the following::

  pip install -U ProDy

If you don't have pip_, please download an installation file and
follow the instructions.


Download Files
--------------

Download a suitable ProDy installation file from http://python.org/pypi/ProDy
http://csb.pitt.edu/ProDy/#downloads. For details of ProDy releases see
:ref:`changes`.

Installation Instructions
-------------------------

After installing the required packages, you will need to do one of the
following:

**Linux**


Download :file:`ProDy-x.{y}.{z}.tar.gz`.  Extract tarball contents and run
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

**Windows**

Remove previously installed ProDy release from :program:`Uninstall a program`
in :guilabel:`Control Panel`.

Download :file:`ProDy-0.{x}.{y}.win32-py2.{z}.exe` and run to install ProDy.

To be able use :ref:`prody-apps` and :ref:`evol-apps` in command prompt
(:program:`cmd.exe`), append Python and scripts folders (e.g.
:file:`C:\\Python27` and :file:`C:\\Python27\\Scripts`) to :envvar:`PATH`
environment variable.

**Testing**

You can test your ProDy installation using the following command::

    $ prody test

Note that :program:`prody` script/command must be reachable from your working
directory.  For more information on testing, see :ref:`testing`.


NMWiz
-----

:ref:`nmwiz` is a `VMD`_ plugin for comparative visual analysis of protein
dynamics modeled using theory or inferred from experimental structural
ensembles.  See :ref:`nmwiz` for available updates and installation
instructions.


Recommended Software
--------------------

* `Scipy`_, when installed, replaces linear algebra module of Numpy.
  Scipy linear algebra module is more flexible and can be faster.
* `IPython`_ is a must have for interactive ProDy sessions.
* `PyReadline`_ for colorful IPython sessions on Windows.
* `MDAnalysis`_ for reading molecular dynamics trajectories.


.. _PyReadline: http://ipython.org/pyreadline.html

Included in ProDy
-----------------

Following software is included in the ProDy installation packages:

* `pyparsing`_ is used to define the atom selection grammar.

* `Biopython`_ KDTree package and pairwise2 module are used for distance based
  atom selections and pairwise sequence alignment, respectively.

* `argparse`_ is used to implement applications and provided for
  compatibility with Python 2.6.

.. _argparse: http://code.google.com/p/argparse/


Source Code
-----------

Source code can be found at https://github.com/abakan/ProDy.