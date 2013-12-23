.. _release:

.. currentmodule:: prody

How to Make a Release
=====================

#. Make sure ProDy imports and passes all unit tests both Python 2 and
   Python 3, and using nose :program:`nosetests` command::

     $ cd ProDy
     $ nosetests
     $ nosetests3


   See :ref:`testing` for more on testing.


#. Update the version number in:

   * :file:`prody/__init__.py`


#. Update the latest release date in:

   * :file:`docs/release/vX.Y_series.rst`.


#. Make sure the following files are file is up-to-date.

   * :file:`README.txt`
   * :file:`MANIFEST.in`
   * :file:`setup.py`


   If there is a new folder with source code, that is a folder in
   :file:`prody`, it should be listed in :file:`MANIFEST.in`.


#. Generate the source distributions::

     $ cd ..
     $ python setup.py sdist --formats=gztar,zip


#. Prepare and test Windows installers (see :ref:`wininst`)::

     $ C:\Python26\python setup.py bdist_wininst
     $ C:\Python27\python setup.py bdist_wininst
     $ C:\Python32\python setup.py bdist_wininst
     $ C:\Python33\python setup.py bdist_wininst


   Alternatively, use :program:`bdist_wininst.bat` to run these commands.
   When there is a newer Python major release, it should be added to this
   list.

#. Register new release to PyPI::

     $ python setup.py register


#. Upload the new release files to the PyPI_.


#. Commit final changes, if there are any::

     $ cd ..
     $ git commit -a


#. Tag the repository with the current version number::

     $ git tag vX.Y


#. Rebase ``devel`` branch to ``master``::

     $ git checkout master
     $ git rebase devel


#. Push the changes with the new tag::

     $ git checkout master
     $ git push --tags
     $ git checkout devel
     $ git push --tags

#. Finally, update the documentation on ProDy_ website.  See :ref:`document`.