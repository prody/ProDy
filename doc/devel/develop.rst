.. _develop:

.. currentmodule:: prody

*******************************************************************************
Contributing to ProDy
*******************************************************************************

Install Git and a GUI
===============================================================================

ProDy source code is managed using Git_ distributed revision controlling 
system.  You need to install :program:`git`, and if you prefer a GUI for it,
on your computer to be able to contribute to development of ProDy. 

On Debian/Ubuntu Linux, for example, you can run the following to install 
:program:`git` and :program:`gitk`::

  $ sudo apt-get install git gitk
 
For other operating systems, you can obtain installation instructions and 
files from Git_.

Also, if you are using Mac OS, you may consider installing SourceTree_, which
is a GUI client that integrates well with Bitbucket_ that hosts ProDy_ source 
code.

You will only need to use a few basic :program:`git` commands.  These commands
are provided below, but usually without an adequate description.  You can refer
to `Git book`_ and `Git docs`_ for detailed explanations of and examples on 
these commands.


Fork and Clone ProDy
===============================================================================

After installing :program:`git` and maybe a GUI, you need a Bitbucket_ 
account.  Note that by using an educational email address you can keep
unlimited private repositories on Bitbucket_ for free.
     
Once you have an account, you need to make a fork of ProDy, which is 
creating an exact copy of the repository in your account.  You will see
a link for this on ProDy_ source code page.  You will have write access to
this fork and later will use it share your changes with others.

The next step is cloning the fork from your Bitbucket account to your local 
system.  You can do it as follows::
   
  $ git clone https://UNAME@bitbucket.org/UNAME/prody.git

This will create :file:`prody` folder with a copy of the project files in it::

  $ cd prody
  $ ls
  CHANGES.rst  lib          logo.svg     README.rst  ubiquitin.png
  doc          LICENSE.rst  MANIFEST.in  scripts
  INSTALL.rst  logo.png     plugins      setup.py


Setup Working Environment
===============================================================================

You can use ProDy directly from this clone by adding :file:`prody/lib` folder
to your :envvar:`PYTHONPATH` environment variable, e.g.::

  export PYTHONPATH=$PYTHONPATH:$/home/abakan/Code/prody/lib

This will not be enough though, since you also need compiled C modules to be 
in the :file:`lib` folder.  You can run the following series of commands to 
build and copy C modules to where they need to be::

  $ # go to ProDy folder
  $ cd prody
  $ python setup.py build
  $ python setup.py copy

You may also want to make sure that you can run :ref:`prody-apps` from anywhere
on your system.  One way to do this by adding :file:`prody/scripts` folder
to your :envvar:`PATH` environment variable, e.g.::

  export PATH=$PATH:$/home/abakan/Code/prody/scripts


Modify, Test, and Commit
===============================================================================

When modifying ProDy files you may want to follow the :ref:`styleguide`. 
Closely following the guidelines therein will allow for incorporation of your
changes to ProDy quickly.

If you changed :file:`.py` files, you should ensure to check the integrity
of the package.  To do this, you should at least run fast ProDy tests as 
follows::

  $ prody test

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
===============================================================================

After you have committed your changes, you will need to push them to your 
Bitbucket account::
  
  git push origin master

This step will ask for your account user name.  If you are going to push
to your Bitbucket account frequently, you may add an SSH key for automatic
authentication.  To add an SSH key for your system, go to 
:menuselection:`Manage Account --> SSH keys` page on Bitbucket_.

After pushing your changes, you will need to make a pull request from your 
to notify ProDy developers of the changes you made and facilitate their
incorporation to ProDy.


Update Local Copy
===============================================================================

You can also keep an up-to-date copy of ProDy by pulling changes from the
master ProDy_ repository on a regular basis.  You need add to the master 
repository as a remote to your local copy.  You can do this running the 
following command from the ProDy project folder::

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

  $ python setup.py build
  $ python setup.py copy


.. _Git: http://git-scm.com/downloads
.. _Git book: http://git-scm.com/book
.. _Git docs: http://git-scm.com/docs
.. _SourceTree: http://www.sourcetreeapp.com
.. _Bitbucket: https://bitbucket.org/
.. _ProDy: https://bitbucket.org/abakan/prody
