.. _build-website:

.. currentmodule:: prody

Building the Website
=================

.. contents::
   :local:


This is a short guide for building the ProDy website.

Environment Setup
--------------

First log in to your ProDy webserver and then run the following::

  $ conda deactivate

This will then bring you into an environment with sphinx-build 1.3.5 and python 2.7, 
which is necessary for building the website.

Next change directory to the website root dir::

  $ cd /var/www/html/prody

In this directory, you will find a number of directories, one of which will be 
attached to a symbolic link to ProDy-website. That directory will contain the 
current website and should not be changed!

Instead navigate to one of the others and build the website in there, such as 
ProDy-website-workdir. You can then copy files back over afterwards. 

It's recommended to have the symbolic link called test_prody pointing to 
your build directory instead and then you can monitor changes by going to 
http://yourdomainname/test_prody/_build/html/ in your web browser.


Updating from GitHub
----------------------

Changes to the website code are made using reStructuredText, which 
is stored in plain-text files with :file:`.rst` extension,
and converted to HTML and PDF pages using `Sphinx`_.

.. _reStructuredText: http://docutils.sourceforge.net/rst.html
.. _reStructuredText Primer: http://sphinx-doc.org/rest.html
.. _Sphinx: http://sphinx-doc.org/

More information about the style etc. can be found in the instructions 
on making tutorials. 

Any changes should be made on your local computer and added to the 
ProDy-Website GitHub repository via pull requests. These can then be 
pulled onto the ProDy webserver in your working directory::

  $ make pull

You may also need to install the ProDy in that directory again to 
make it get used during the building of the website. This can be 
done as follows::

  $ cd ProDy
  $ pip install -U . --user
  $ cd ..

DruGUI and its tutorial are handled by their own GitHub repo. The 
following command can be used to update them::

  $ make drugui

Publishing Changes
-------------------

For generally making :file:`.rst` files convert to HTML format, use the following
command::

  $ make html

You will find HTML files in :file:`_build/html` folder.

It may also help to run the following first::

  $ make clean

For tutorials, once a tutorial is complete and looks good in HTML (no code execution
problems), the following commands can be used to generate a PDF file and
tutorial file achieves::

  $ make pdf
  $ make files

ProDy online documentation will contain these files as well as tutorial pages
in HTML format.

You can then copy files over to the main ProDy-Website directory to have them 
incorporated into the main prody website.
