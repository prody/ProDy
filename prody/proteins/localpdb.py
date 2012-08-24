# -*- coding: utf-8 -*-
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

"""This module defines functions for handling local PDB folders."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from glob import glob
from os.path import abspath, isdir, join

from prody import LOGGER, SETTINGS

__all__ = ['getPDBLocalFolder', 'getPDBMirrorPath', 
           'setPDBLocalFolder', 'setPDBMirrorPath',
           'iterPDBFilenames']
           

def getPDBLocalFolder():
    """Return the path to a local PDB folder and folder structure specifier. 
    If a local folder is not set, **None** will be returned."""

    folder = SETTINGS.get('pdb_local_folder')
    if folder is not None:
        if isinstance(folder, str) and isdir(folder):
            return folder, SETTINGS.get('pdb_local_divided', True)
        else:
            LOGGER.warning('PDB local folder {0:s} is not a accessible.'
                           .format(repr(folder)))


def setPDBLocalFolder(folder, divided=False):
    """Set a local PDB folder.  Setting a local PDB folder will make 
    :func:`fetchPDB` function to seek that folder for presence of requested
    PDB files.  Also, files downloaded from `wwPDB <http://www.wwpdb.org/>`_ 
    FTP servers will be saved in this folder.  This may help users to store 
    PDB files in a single place and have access to them in different working 
    directories.
    
    If *divided* is **True**, the divided folder structure of wwPDB servers 
    will be assumed when reading from and writing to the local folder.  For 
    example, a structure with identifier **1XYZ** will be present as 
    :file:`pdblocalfolder/yz/pdb1xyz.pdb.gz`. 
    
    If *divided* is **False**, a plain folder structure will be expected and 
    adopted when saving files.  For example, the same structure will be 
    present as :file:`pdblocalfolder/1xyz.pdb.gz`.
    
    Finally, in either case, lower case letters will be used and compressed
    files will be stored."""
    
    if not isinstance(folder, str):
        raise TypeError('folder must be a string')
    assert isinstance(divided, bool), 'divided must be a boolean'
    if isdir(folder):
        folder = abspath(folder)
        LOGGER.info('Local PDB folder is set: {0:s}'.format(repr(folder)))
        if divided:
            LOGGER.info('When using local PDB folder, wwPDB divided '
                        'folder structure will be assumed.')
        else:
            LOGGER.info('When using local PDB folder, a plain folder structure '
                        'will be assumed.')
        SETTINGS['pdb_local_folder'] = folder
        SETTINGS['pdb_local_divided'] = divided
        SETTINGS.save()
    else:
        raise IOError('No such directory: {0:s}'.format(repr(folder)))


def getPDBMirrorPath():
    """Return the path to a local PDB mirror, or **None** if a mirror path is 
    not set."""

    path = SETTINGS.get('pdb_mirror_path')
    if isinstance(path, str):
        if isdir(path):
            return path
        else:
            LOGGER.warning('PDB mirror path {0:s} is not a accessible.'
                           .format(repr(path)))


def setPDBMirrorPath(path):
    """Set the path to a local PDB mirror."""
    
    if not isinstance(path, str):
        raise TypeError('path must be a string')
    if isdir(path):
        path = abspath(path)
        LOGGER.info('Local PDB mirror path is set: {0:s}'.format(repr(path)))
        SETTINGS['pdb_mirror_path'] = path
        SETTINGS.save()
    else:
        raise IOError('No such directory: {0:s}'.format(repr(path)))


def iterPDBFilenames(path=None, sort=False):
    """Yield PDB filenames in local PDB mirror (see :func:`.getPDBMirrorPath`)
    or in *path* specified by the user."""

    if path is None:
        path = getPDBMirrorPath()
        if path is None:
            raise ValueError('path must be specified or PDB mirror path '
                               'must be set')
        path = join(path, 'data/structures/divided/pdb/', '*/*.ent.gz')
        pdbs = glob(path)
    else:
        pdbs = []
        for ext in ['.pdb', 'pdb.gz', '.ent', '.ent.gz']:
            pdbs.extend(glob(join(path, '*' + ext)))
    if sort:
        pdbs.sort()
    for fn in pdbs:
        yield fn
