.. _develop:

.. currentmodule:: prody

Contributing to ProDy
=====================

.. contents::
   :local:


Install Git and a GUI
---------------------

ProDy source code is managed using Git_ distributed revision controlling
system.  You need to install :program:`git`, and if you prefer a GUI for it,
on your computer to be able to contribute to development of ProDy.

On Debian/Ubuntu Linux, for example, you can run the following to install
:program:`git` and :program:`gitk`::

  $ sudo apt-get install git gitk

For other operating systems, you can obtain installation instructions and
files from Git_.

You will only need to use a few basic :program:`git` commands.  These commands
are provided below, but usually without an adequate description.  Please refer
to `Git book`_ and `Git docs`_ for usage details and examples.


Fork and Clone ProDy
--------------------

ProDy source code an issue tracker are hosted on Github_.  You need to create
an account on this service, if you do not have one already.


If you work on Mac OS or Windows, you may consider getting `GitHub Mac`_ or
`GitHub Windows`_ to help you manage a copy of the repository.

Once you have an account, you need to make a fork of ProDy, which is
creating a copy of the repository in your account.  You will see
a link for this on ProDy_ source code page.  You will have write access to
this fork and later will use it share your changes with others.

The next step is cloning the fork from your online account to your local
system.  If you are not using the GitHub software, you can do it as follows::

  $ git clone https://github.com/prody/ProDy.git
git

This will create :file:`ProDy` folder with a copy of the project files in it::

  $ cd ProDy
  $ ls
  bdist_wininst.bat  docs   INSTALL.rst  LICENSE.rst  Makefile
  MANIFEST.in        prody  README.rst   scripts      setup.py


Setup Working Environment
-------------------------

You can use ProDy directly from this clone by adding :file:`ProDy` folder
to your :envvar:`PYTHONPATH` environment variable, e.g.::

  export PYTHONPATH=$PYTHONPATH:$/home/USERNAME/path/to/ProDy

This will not be enough though, since you also need to compile C extensions.
You can run the following series of commands to build and copy C modules to where they need to be::

  $ cd ProDy
  $ python setup.py build_ext --inplace --force

or, on Linux you can::

  $ make build

You may also want to make sure that you can run :ref:`prody-apps` from anywhere
on your system.  One way to do this by adding :file:`ProDy/scripts` folder
to your :envvar:`PATH` environment variable, e.g.::

  export PATH=$PATH:$/home/USERNAME/path/to/ProDy/scripts


Modify, Test, and Commit
------------------------

When modifying ProDy files you may want to follow the :ref:`styleguide`.
Closely following the guidelines therein will allow for incorporation of your
changes to ProDy quickly.

If you changed :file:`.py` files, you should ensure to check the integrity
of the package.  To do this, you should at least run fast ProDy tests as
follows::

  $ cd ProDy
  $ nosetests

See :ref:`testing` for alternate and more comprehensive ways of testing.
ProDy unittest suit may not include a test for the function or the class
that you just changed, but running the tests will ensure that the ProDy
package can be imported and run without problems.

After ensuring that the package runs, you can commit your changes as follows::

  $ git commit modified_file_1.py modified_file_2.py

or::

  $ git commit -a

This command will open a text editor for you to describe the changes that
you just committed.


Push and Pull Request
---------------------

After you have committed your changes, you will need to push them to your
Bitbucket account::

  git push origin master

This step will ask for your account user name.  If you are going to push
to your GitHub/Bitbucket account frequently, you may add an SSH key for
automatic authentication.  To add an SSH key for your system, go to
:menuselection:`Edit Your Profile --> SSH keys` page on GitHub or
:menuselection:`Manage Account --> SSH keys` page on Bitbucket.

After pushing your changes, you will need to make a pull request from your
to notify ProDy developers of the changes you made and facilitate their
incorporation to ProDy.


Update Local Copy
-----------------

You can also keep an up-to-date copy of ProDy by pulling changes from the
master ProDy_ repository on a regular basis.  You need add to the master
repository as a remote to your local copy.  You can do this running the
following command from the ProDy project folder::

  $ cd prody
  $ git remote add prodymaster git@github.com:abakan/ProDy.git

or::

  $ cd prody
  $ git remote add prodymaster git@bitbucket.org:abakan/prody.git


You may use any name other than `prodymaster`, but `origin`, which points to
the ProDy fork in your account.

After setting up this remote, calling :program:`git pull` command will
fetch latest changes from ProDy_ master repository and merge them to your
local copy::

  $ git pull prodymaster master

Note that when there are changes in C modules, you need to run the following
commands again to update the binary module files::

  $ python setup.py build_ext --inplace --force


.. _Git: http://git-scm.com/downloads
.. _Git book: http://git-scm.com/book
.. _Git docs: http://git-scm.com/docs
.. _GitHub Mac: http://mac.github.com
.. _GitHub Windows: http://windows.github.com
