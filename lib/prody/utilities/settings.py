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

"""This module defines class for handling and storing package settings."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
import cPickle

from pathtools import pickle, unpickle, isWritable, USERHOME
import prody as pkg

__all__ = ['PackageSettings', 'getPackagePath', 'setPackagePath',]


class PackageSettings(object):
    
    """A class for managing package settings.  Settings are saved in user's 
    home director.  When settings are changed by the users, the changes are 
    automatically saved.  Settings are stored in a :class:`dict` instance. 
    The dictionary is pickled in user's home directory for permanent storage.
    """
    
    def __init__(self, pkg, rcfile=None, logger=None):
        """*rcfile* is the filename for pickled settings dictionary, and by 
        default is set to :file:`.pkgrc`."""
        
        self._package = pkg
        if rcfile is None:
            self._rcfile = os.path.join(USERHOME, '.' + pkg + 'rc')
        else:
            self._rcfile = rcfile

        self._logger = logger
        self._settings = {}
        
    def __getitem__(self, key):
        
        return self._settings[key]
        
    def __setitem__(self, key, value):
        
        self._settings[key] = value
        
    def get(self, key, default=None):
        """Return value corresponding to specified *key*, or *default* if *key*
        is not found."""
        
        return self._settings.get(key, default)
        
    def pop(self, key, default=None):
        """Remove specified *key* and return corresponding value.  If *key* is 
        not found, *default* is returned."""
        
        return self._settings.pop(key, default)
    
    def update(self, *args, **kwargs):
        """Update settings dictionary. """
        
        for arg in args:
            self._settings.update(arg)
        if kwargs:
            self._settings.update(kwargs)
        
    def load(self):
        """Load settings by unpickling the settings dictionary."""

        if os.path.isfile(self._rcfile):
            try:
                settings = unpickle(self._rcfile)
            except Exception as err:
                if self._logger:
                    self._logger.warning("{0:s} configuration file '{1:s}' "
                                "could not be loaded ({2:s})."
                                .format(self._package, self._rcfile, err))
            else:                    
                if isinstance(settings, dict):
                    self._settings.update(settings)

    def save(self):
        """Save settings by pickling the settings dictionary."""
        
        if isWritable(USERHOME):
            try:
                pickle(self._settings, self._rcfile, backup=False)
            except Exception as err:
                if self._logger:
                    self._logger.warning("{0:s} cannot write configuration "
                                "file '{1:s}' ({2:s})."
                                .format(self._package, self._rcfile, err))
        elif self._logger:
            self._logger.warning("{0:s} cannot write configuration file to "
                                "'{1:s}', user does not have write access."
                                .format(self._package, USERHOME))


def setPackagePath(path):
    """Set package path."""
    
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except Exception as err:
            pkg.LOGGER.warning('Failed to make folder "{0:s}": {1:s}'
                           .format(path, err.strerror))
            return False
    pkg.SETTINGS['package_path'] = path
    pkg.SETTINGS.save()
    return path    


def getPackagePath():
    """Return package path."""
    
    path = pkg.SETTINGS.get('package_path', None)
    
    update = False
    if path is None:
        pkg.LOGGER.warning('{0:s} package path is not yet set by the user.'
                       .format(__package__))
        update = True
    elif not os.path.isdir(path):
        pkg.LOGGER.warning("{0:s} package path '{1:s}' does not exist."
                       .format(__package__, path))
        update = True
    elif not os.access(path, os.W_OK):
        pkg.LOGGER.warning("User does not have write access to {0:s} package "
                           "path '{1:s}'.".format(__package__, path))
        update = True
    if update:
        default = os.path.join(USERHOME, '.' + __package__)
        path = raw_input('Please specify a folder for storing {0:s} data '
                         '(press enter for "{1:s}"):'
                         .format(__package__, default)) or default
        while not setPackagePath(path):
            path = raw_input('Please specify a valid folder name with write ' 
                             'access:')
    return path
