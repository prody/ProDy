# ProDy: A Python Package for Protein Dynamics Analysis
#
# Copyright (C) 2010-2012 Ahmet Bakan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines functions for handling files and paths."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
import sys
from os import sep as pathsep
from glob import glob as pyglob
import pickle as pypickle
import zipfile
import platform
import os.path
from os.path import isfile, isdir, join, split, splitext
from os.path import getsize, isabs, exists, abspath
from shutil import copy

PLATFORM = platform.system()
USERHOME = os.getenv('USERPROFILE') or os.getenv('HOME')

__all__ = ['gunzip', 'backupFile', 'openFile',
           'openDB', 'openSQLite', 'openURL', 'copyFile',
           'isExecutable', 'isReadable', 'isWritable',
           'makePath', 'relpath', 'sympath', 'which',
           'pickle', 'unpickle', 'glob', 'addext',
           'PLATFORM', 'USERHOME']

major, minor = sys.version_info[:2]
if major > 2 and minor < 3:
    import gzip
    from gzip import GzipFile
    import io

    class TextIOWrapper(io.TextIOWrapper):

        def _getlines(self):

            try:
                lines = self._lines
            except AttributeError:
                self._lines = None

            if self._lines is None:
                self._lines = self.read().split('\n')
            return self._lines

        def readline(self, *args):

            lines = self._getlines()
            if lines:
                return lines.pop(0)
            else:
                return ''

        def readlines(self, size=None):

            lines = self._getlines()
            if size is None:
                self._lines = []
                return lines
            else:
                self._lines = lines[size:]
                return lines[:size]


    def gzip_open(filename, mode="rb", compresslevel=9,
             encoding=None, errors=None, newline=None):
        """Open a gzip-compressed file in binary or text mode.

        The filename argument can be an actual filename (a str or bytes object), or
        an existing file object to read from or write to.

        The mode argument can be "r", "rb", "w", "wb", "a" or "ab" for binary mode,
        or "rt", "wt" or "at" for text mode. The default mode is "rb", and the
        default compresslevel is 9.

        For binary mode, this function is equivalent to the GzipFile constructor:
        GzipFile(filename, mode, compresslevel). In this case, the encoding, errors
        and newline arguments must not be provided.

        For text mode, a GzipFile object is created, and wrapped in an
        io.TextIOWrapper instance with the specified encoding, error handling
        behavior, and line ending(s).

        """
        if "t" in mode:
            if "b" in mode:
                raise ValueError("Invalid mode: %r" % (mode,))
        else:
            if encoding is not None:
                raise ValueError("Argument 'encoding' not supported in binary mode")
            if errors is not None:
                raise ValueError("Argument 'errors' not supported in binary mode")
            if newline is not None:
                raise ValueError("Argument 'newline' not supported in binary mode")

        gz_mode = mode.replace("t", "")
        if isinstance(filename, (str, bytes)):
            binary_file = GzipFile(filename, gz_mode, compresslevel)
        elif hasattr(filename, "read") or hasattr(filename, "write"):
            binary_file = GzipFile(None, gz_mode, compresslevel, filename)
        else:
            raise TypeError("filename must be a str or bytes object, or a file")

        if "t" in mode:
            return TextIOWrapper(binary_file, encoding, errors, newline)
        else:
            return binary_file
else:
    import gzip
    def gzip_open(filename, *args, **kwargs):
        if args and "t" in args[0]:
            args = (args[0].replace("t", ""), ) + args[1:]
        if isinstance(filename, str):
            return gzip.open(filename, *args, **kwargs)
        else:
            return gzip.GzipFile(filename, *args, **kwargs)

if (major, minor) >= (3, 2):
    from gzip import compress as gzip_compress
    from gzip import decompress as gzip_decompress

OPEN = {
    '.gz': gzip_open,
    '.zip': zipfile.ZipFile,
}


def backupFile(filename, backup=None, backup_ext='.BAK', **kwargs):
    """Rename *filename* with *backup_ext* appended to its name for backup
    purposes, if *backup* is **True** or if automatic backups is turned on
    using :func:`.confProDy`.  Default extension :file:`.BAK` is used when
    one is not set using :func:`.confProDy`.  If *filename* does not exist,
    no action will be taken and *filename* will be returned.  If file is
    successfully renamed, new filename will be returned."""

    try:
        exists = isfile(filename)
    except Exception as err:
        raise TypeError('filename must be a string ({0})'.format(str(err)))

    from prody import SETTINGS
    if exists and (backup or SETTINGS.get('backup', False)):
        if backup_ext == '.BAK':
            backup_ext = SETTINGS.get('backup_ext', '.BAK')
        bak = filename + backup_ext
        if isfile(bak):
            try:
                os.remove(bak)
            except Exception as err:
                pass
        try:
            os.rename(filename, bak)
        except Exception as err:
            pass
        return bak
    else:
        return filename


def openFile(filename, *args, **kwargs):
    """Open *filename* for reading, writing, or appending.  First argument in
    *args* is treated as the mode.  Opening :file:`.gz` and :file:`.zip` files
    for reading and writing is handled automatically.

    :arg backup: backup existing file using :func:`.backupFile` when opening
        in append or write modes, default is obtained from package settings
    :type backup: bool

    :arg backup_ext: extension for backup file, default is :file:`.BAK`
    :type backup_ext: str"""

    from prody import SETTINGS
    try:
        exists = isfile(filename)
    except Exception as err:
        raise TypeError('filename must be a string ({0})'.format(str(err)))

    folder = kwargs.pop('folder', None)
    if folder:
        filename = join(folder, filename)

    backup = kwargs.pop('backup', None)
    if backup is not None and backup and args and args[0][0] in ('a', 'w'):
        backupFile(filename, backup=backup,
                   backup_ext=kwargs.pop('backup_ext', None))

    ext = splitext(filename)[1]
    return OPEN.get(ext.lower(), open)(filename, *args, **kwargs)


def gunzip(filename, outname=None):
    """Return output name that contains decompressed contents of *filename*.
    When no *outname* is given, *filename* is used as the output name as it
    is or after :file:`.gz` extension is removed.  *filename* may also be a
    string buffer, in which case decompressed string buffer or *outname* that
    contains buffer will be returned."""

    if len(filename) < 1000:
        try:
            afile = isfile(filename)
        except TypeError:
            afile = False
    else:
        afile = False

    if afile:
        if outname is None:
            if filename.endswith('.gz'):
                outname = filename[:-3]
            elif filename.endswith('.tgz'):
                outname = filename[:-4] + '.tar'
            elif filename.endswith('.gzip'):
                outname = filename[:-5]
            else:
                outname = filename

        inp = gzip_open(filename, 'rb')
        data = inp.read()
        inp.close()
        out = open(outname, 'wb')
        out.write(data)
        out.close()
        return outname
    else:
        result = None
        try:
            from StringIO import StringIO
        except ImportError:
            from io import BytesIO
            buff = gzip_open(BytesIO(filename))
            if outname is None:
                try:
                    result = buff.read()
                except IOError:
                    pass
            else:
                buff = buff.read()
                if isinstance(buff, bytes):
                    with open(outname, 'wb') as out:
                        out.write(buff)
                else:
                    with open(outname, 'wb') as out:
                        out.write(buff)
                return outname
        else:
            from StringIO import StringIO
            buff = gzip.GzipFile(fileobj=StringIO(filename))
            try:
                result = buff.read()
            except IOError:
                pass

        if result is not None:
            if outname is None:
                return result
            else:
                with open(outname, 'w') as out:
                    out.write(result)
                return outname
        raise ValueError('filename is not a valid path or a compressed'
                         ' string buffer')


def isExecutable(path):
    """Return true if *path* is an executable."""

    return (isinstance(path, str) and exists(path) and
        os.access(path, os.X_OK))


def isReadable(path):
    """Return true if *path* is readable by the user."""

    return (isinstance(path, str) and exists(path) and
        os.access(path, os.R_OK))


def isWritable(path):
    """Return true if *path* is writable by the user."""

    return (isinstance(path, str) and exists(path) and
        os.access(path, os.W_OK))


def relpath(path):
    """Return *path* on Windows, and relative path elsewhere."""

    if PLATFORM == 'Windows':
        return path
    else:
        return os.path.relpath(path)


def sympath(path, beg=2, end=1, ellipsis='...'):
    """Return a symbolic path for a long *path*, by replacing folder names
    in the middle with *ellipsis*.  *beg* and *end* specified how many folder
    (or file) names to include from the beginning and end of the path."""

    abs_items = abspath(path).split(pathsep)
    rel_items = relpath(path).split(pathsep)
    if len(abs_items) <= len(rel_items):
        items = abs_items
    else:
        items = rel_items
    if len(items) <= beg + end:
        return pathsep.join(items)
    else:
        return pathsep.join(items[:beg+1] + [ellipsis] + items[-end:])


def makePath(path):
    """Make all directories that does not exist in a given *path*."""

    if not isdir(path):
        dirs = path.split(pathsep)
        for i, dirname in enumerate(dirs):
            if not dirname:
                continue
            dirname = pathsep.join(dirs[:i+1])
            try:
                if not isdir(dirname):
                    os.mkdir(dirname)
            except OSError:
                raise OSError('{0} could not be created, please '
                            'specify another path'.format(path))
    return path


def which(program):
    """This function is based on the example in:
    http://stackoverflow.com/questions/377017/"""

    fpath, fname = os.path.split(program)
    if fpath and isExecutable(program):
        return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = os.path.join(path, program)
            if isExecutable(path):
                return path
    return None


def pickle(obj, filename, protocol=2, **kwargs):
    """Pickle *obj* using :func:`pickle.dump` in *filename*.  *protocol* is set
    to 2 for compatibility between Python 2 and 3."""

    out = openFile(filename, 'wb', **kwargs)
    pypickle.dump(obj, out, protocol=2)
    out.close()
    return filename


def unpickle(filename, **kwargs):
    """Unpickle object in *filename* using :func:`pickle.load`."""

    inf = openFile(filename, 'rb', **kwargs)
    obj = pypickle.load(inf)
    inf.close()
    return obj


def openDB(filename, *args):
    """Open a database with given *filename*."""

    try:
        import anydbm as dbm
    except ImportError:
        import dbm
    return anydbm.open(filename, *args)


def openSQLite(filename, *args):
    """Return a connection to SQLite database *filename*.  If ``'n'`` argument
    is passed, remove any existing databases with the same name and return
    connection to a new empty database."""

    if 'n' in args and isfile(filename):
        os.remove(filename)
    import sqlite3
    return sqlite3.connect(filename)


def openURL(url, timeout=5, **kwargs):
    """Open *url* for reading. Raise an :exc:`IOError` if *url* cannot be
    reached.  Small *timeout* values are suitable if *url* is an ip address.
    *kwargs* will be used to make :class:`urllib.request.Request` instance
    for opening the *url*."""

    try:
        from urllib2 import urlopen, URLError, Request
    except ImportError:
        from urllib.request import urlopen, Request
        from urllib.error import URLError

    if kwargs:
        request = Request(url, **kwargs)
    else:
        request = str(url)

    try:
        return urlopen(request, timeout=int(timeout))
    except URLError:
        raise IOError('{0} could not be opened for reading, invalid URL or '
                      'no internet connection'.format(repr(request)))


def glob(*pathnames):
    """Return concatenation of ordered lists of paths matching patterns in
    *pathnames*."""

    paths = []
    for pathname in pathnames:
        matches = pyglob(pathname)
        matches.sort()
        paths.extend(matches)
    return paths


def copyFile(src, dst):
    """Return *dst*, a copy of *src*."""

    copy(src, dst)
    return dst


def addext(filename, extension):
    """Return *filename*, with *extension* if it does not have one."""

    return filename + ('' if splitext(filename)[1] else extension)
