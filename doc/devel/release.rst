.. _release:

.. currentmodule:: prody

*******************************************************************************
How to Make a New Release
*******************************************************************************

#. Make sure ProDy imports and passes unit tests::
     
     $ python
     >>> import prody
     >>> prody.test()

#. Clean the contents of the :file:`Documentation` folder::

     $ cd doc
     $ make clean
  
#. Make sure ProDy passes doctests::

     $ make doctest
  

#. Update the version numbers in: 
    
   * :file:`doc/conf.py`.
   * :file:`prody/__init__.py` 


#. Update the latest release date in:
   
   * :file:`doc/changes.rst`.

#. Make sure the following files are  file is up-to-date.
    
   * :file:`README.txt`
   * :file:`MANIFEST.in`
   * :file:`setup.py`

#. Generate the source distributions::

     $ cd ..
     $ python setup.py sdist --formats=gztar,zip
  
#. Generate Windows installers using Python 2.6 and 2.7
   (see :ref:`wininst`)::

     $ C:\python26\python setup.py bdist_wininst
     $ C:\python27\python setup.py bdist_wininst

#. Test the installers.

#. Register new release to PyPI::

     $ python setup.py register

#. Upload the new release files to the 
   `PyPI <http://pypi.python.org/pypi/ProDy/>`_.

#. Update the release download file names::
   
     $ cd doc
     $ ./pypi.py

#. Update the links to NMWiz installers in:

   * :file:`doc/plugins/getnmwiz.rst`.

#. Copy the documentation files packed with source release::

     $ make copytxt

#. Generate HTML and PDF documentation::

     $ make html
     $ make latexpdf
  

#. Commit the final changes::

     $ cd ..
     $ git commit -a
  
#. Tag the repository with the current version number::

     $ git tag vX.Y
  
#. Push the changes with tags to `github <https://github.com/abakan/prody>`_::

     $ git push --tags origin master

