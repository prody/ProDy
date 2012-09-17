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
import gzip
from glob import glob as pyglob
import pickle as pypickle
import zipfile
import platform
import os.path
from os.path import isfile, isdir, join, split, splitext
from os.path import getsize, isabs, exists

import prody as pkg

PLATFORM = platform.system()
USERHOME = os.getenv('USERPROFILE') or os.getenv('HOME')

__all__ = ['gunzip', 'openFile', 'openDB', 
           'isExecutable', 'isReadable', 'isWritable', 
           'makePath', 'relpath', 'which', 
           'pickle', 'unpickle', 'glob',
           'PLATFORM', 'USERHOME']
        

OPEN = {
    '.gz': gzip.open,
    '.GZ': gzip.open,
    '.zip': zipfile.ZipFile,
    '.ZIP': zipfile.ZipFile,
}

def openFile(filename, *args, **kwargs):
    """Open *filename* for reading, writing, or appending.  First argument in
    *args* is treated as the mode.
    
    :arg backup: backup existing file when opening for appending or writing,
        default is obtained from package settings
    :type backup: bool
    :arg backup_ext: extension for backup file, default is :file:`.BAK`
    :type backup_ext: str"""

    if not isinstance(filename, str):
        raise TypeError('filename must be a string')
    folder = kwargs.pop('folder', None)
    if folder:
        filename = join(folder, filename)
    ext = splitext(filename)[1]
    backup = kwargs.pop('backup', pkg.SETTINGS.get('backup', False))
    backup_ext = kwargs.pop('backup_ext', 
                            pkg.SETTINGS.get('backup_ext', '.BAK'))
    if args and args[0][0] in ('a', 'w'):
        if isfile(filename) and backup:
            bak = filename + backup_ext
            if isfile(bak):
                os.remove(bak)
            os.rename(filename, bak)
    return OPEN.get(ext, open)(filename, *args, **kwargs)
    
    
def gunzip(filename, outname=None):
    """Decompress *filename* and save as *outname*.  When ``outname=None``, 
    *filename* is used as the output name.  Returns output filename upon 
    successful completion."""

    if not isinstance(filename, str):
        raise TypeError('filename must be a string')
    if not isfile(filename):
        raise ValueError('{0:s} does not exist'.format(filename))
    if outname is None: 
        outname = filename
    inp = gzip.open(filename, 'rb')
    data = inp.read()
    inp.close()
    out = open(outname, 'w')
    out.write(data)
    out.close()
    return outname


def isExecutable(path):
    """Return true if *path* is an executable."""
    
    return isinstance(path, str) and exists(path) and \
        os.access(path, os.X_OK)


def isReadable(path):
    """Return true if *path* is readable by the user."""
    
    return isinstance(path, str) and exists(path) and \
        os.access(path, os.R_OK)


def isWritable(path):
    """Return true if *path* is writable by the user."""
    
    return isinstance(path, str) and exists(path) and \
        os.access(path, os.W_OK)


def relpath(path):
    """Return *path* on Windows, and relative path elsewhere."""
    
    if PLATFORM == 'Windows':
        return path
    else:
        return os.path.relpath(path)


def makePath(path):
    """Make all directories that does not exist in a given path."""
    
    if isabs(path):
        path = relpath(path)
    if not isdir(path):
        dirs = path.split(os.sep)
        for i in range(len(dirs)):
            dirname = os.sep.join(dirs[:i+1])
            try:
                if not isdir(dirname): 
                    os.mkdir(dirname)
            except OSError:
                raise OSError('{0:s} could not be created, please '
                            'specify another path'.format(path))
                return os.getcwd()
    return join(os.getcwd(), path)


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


def pickle(obj, filename, **kwargs):
    """Pickle *obj* using :func:`pickle.dump` in *filename*."""
    
    out = openFile(filename, 'wb', **kwargs)
    pypickle.dump(obj, out)
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

    if filename.endswith('sql'):    
        if 'n' in args and isfile(filename):
            os.remove(filename)
        import sqlite3
        return sqlite3.connect(filename)
    else:
        import anydbm
        return anydbm.open(filename, *args)


def glob(*pathnames):
    """Return concatenation of ordered lists of paths matching patterns in 
    *pathnames*."""
    
    paths = []
    for pathname in pathnames:
        matches = pyglob(pathname)
        matches.sort()
        paths.extend(matches)
    return paths
        
        
