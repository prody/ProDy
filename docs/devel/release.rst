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

   Also, commend ``+ '-dev'`` out, so that documentation will build
   for a stable release.


#. Update the most recent changes and the latest release date in:

   * :file:`docs/release/vX.Y_series.rst`.

   If there is a new incremental release, start a new file.


#. Make sure the following files are up-to-date.

   * :file:`README.txt`
   * :file:`MANIFEST.in`
   * :file:`setup.py`


   If there is a new file format, that is a new extensions not captured in
   :file:`MANIFEST.in`, it should be included.

   If there is a new C extension, it should be listed in :file:`setup.py`.

   After checking these files, commit change and push them to GitHub_.


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
   list. Don't forget to pull most recent changes to your Windows machine.

   A good practice is installing ProDy using all newly created installers
   and checking that it works. ProDy script can be used to check that, e.g.::


      C:\> C:\Python33\Scripts\prody.bat anm 1ubi


   If this command runs for all supported Python versions, release is good
   to go.

#. Register new release to PyPI::

     $ python setup.py register

   This will offer a number of options. ProDy on PyPI is owned by user
   ``prody.devel``.


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

#. Update the documentation on ProDy_ website.  See :ref:`document`.

#. Now that you made a release, you can go back to development.
   You may start with appending ``'-dev'`` to ``__release__`` in
   :file:`prody/__init__.py`.