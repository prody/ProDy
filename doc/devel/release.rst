.. _release:

.. currentmodule:: prody

*******************************************************************************
How to Make a Release
*******************************************************************************

#. Make sure ProDy imports and passes unit tests (Python 2.7)::
     
     $ python
     >>> import prody
     >>> prody.test()

#. Clean the contents of the :file:`Documentation` folder::

     $ cd doc
     $ make clean
  
#. Make sure ProDy passes doctests::

     $ make doctest
  

#. Update the version numbers in: 
    
   * :file:`doc/conf.py`
   * :file:`prody/__init__.py` 


#. Update the latest release date in:
   
   * :file:`doc/changes.rst`.

#. Make sure the following files are  file is up-to-date.
    
   * :file:`README.txt`
   * :file:`MANIFEST.in`
   * :file:`setup.py`

#. Copy the documentation files packed with source release::

     $ make copytxt

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

#. Update new release file links::
   
     $ cd doc
     $ make stats

#. Generate HTML and PDF documentation::

     $ make html
     $ make latexpdf

#. Commit the final changes::

     $ cd ..
     $ git commit -a
  
#. Tag the repository with the current version number::

     $ git tag vX.Y
  
#. Push the changes with tags to `https://bitbucket.org/abakan/prody>`_::

     $ git push --tags origin master

