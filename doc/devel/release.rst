.. _release:

.. currentmodule:: prody

*******************************************************************************
How to Make a Release
*******************************************************************************

#. Make sure ProDy imports and passes all unit tests::
     
     $ prody test -l full
     
   This will run tests using the default Python installation.  The tests 
   should be executed for both Python 2.7 and Python 3.2, and using nose
   via :program:`nosetests` command may help::
   
     $ cd prody/tests
     $ nosetests
     $ nosetests3
     

   See :ref:`testing` for alternative testing methods. 


#. Make sure ProDy passes doctests::

     $ cd doc
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
  
  
#. Prepare and test Windows installers (see :ref:`wininst`)::

     $ C:\Python27\python setup.py bdist_wininst
     $ C:\Python32\python setup.py bdist_wininst


#. Register new release to PyPI::

     $ python setup.py register


#. Upload the new release files to the 
   `PyPI <http://pypi.python.org/pypi/ProDy/>`_.


#. Update new release file links::
   
     $ cd doc
     $ make stats


#. Generate HTML and PDF documentation::

     $ make clean
     $ make html
     $ make latexpdf


#. Commit the final changes::

     $ cd ..
     $ git commit -a
  
  
#. Tag the repository with the current version number and push the changes 
   with tags to `Bitbucket <https://bitbucket.org/abakan/prody>`_::

     $ git tag vX.Y
  
