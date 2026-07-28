# Druggability: Python Package and VMD GUI for Druggability Index Analysis
# Copyright (C) 2010  Ahmet Bakan <ahb12@pitt.edu>
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""This module defines a base class for other classes. Important features of 
this base class are its logging and self pickling functionalities.

Classes:

* :class:`ABase`
    
Functions:

* :func:`get_logger`
    
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'
__version__ = '0.5.2'

import pickle
import gzip
import logging
import logging.handlers
import os
import os.path

LOGGING_LEVELS = {'debug': logging.DEBUG,
                  'info': logging.INFO,
                  'warning': logging.WARNING,
                  'error': logging.ERROR,
                  'critical': logging.CRITICAL}
LOGGING_LEVELS.setdefault(logging.INFO)

SIGNATURE = '@>'

__all__ = ['ABase']

def _set_workdir(workdir):
    """Set a working directory, by creating if it doesn't exist."""
    if os.path.isabs(workdir):
        workdir = os.path.relpath(workdir)
    if not os.path.isdir(workdir):
        dirs = workdir.split(os.sep)
        for i in range(len(dirs)):
            dirname = os.sep.join(dirs[:i+1])
            try:
                if not os.path.isdir(dirname): 
                    os.mkdir(dirname)
            except OSError:
                return os.getcwd()
    return os.path.join(os.getcwd(), workdir)

def get_logger(name, **kwargs):
    """Return a logger.
    
    :arg name: name of the logger instance
    :type name: str
    
    :keyword verbose: loglevel for console verbosity
    :type verbose: str, default is "info" 
    
    :keyword writelog: control logging in a file
    :type writelog: bool, default is True
    
    :keyword workdir: location of logfile
    :type workdir: str, default is "."
    
    :keyword loglevel: loglevel for logfile verbosity
    :type loglevel: str, default is "debug"
    
    :keyword logfilemode: mode in which logfile will be opened 
    :type logfilemode: str, default is "w"
    
    :keyword backupcount: number of old *name.log* files to save
    :type backupcount: int, default is 3

    ======== ==================================================================
    Loglevel Description
    ======== ==================================================================
    debug    Eveything will be printed on the colsole or written into logfile.
    info     Only brief information will be printed or written.
    warning  Only critical information will be printed or written.
    error    This loglevel is equivalent to *warning* in package.
    critical This loglevel is equivalent to *warning* in package.
    ======== ==================================================================

        
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    
    if logger.handlers: 
        return logger

    
    console = logging.StreamHandler()
    console.setLevel(LOGGING_LEVELS[kwargs.get('verbose', 'info')])
    console.setFormatter(logging.Formatter(SIGNATURE + ' %(message)s'))
    logger.addHandler(console)
    
    if not ('writelog' in kwargs and not kwargs['writelog']):
        logfilename = os.path.join(kwargs.get('workdir', '.'), name+'.log')
        rollover = False 
        # if filemode='a' is provided, rollover is not performed
        if os.path.isfile(logfilename) and kwargs.get('filemode', None) != 'a':
            rollover = True
        logfile = logging.handlers.RotatingFileHandler(logfilename, 
                    mode=kwargs.get('filemode', 'a'), maxBytes=0,
                    backupCount=kwargs.get('backupcount', 3))
        logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
        logfile.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(logfile)
        if rollover:
            logger.debug('Saving this log file and starting a new one.')
            logfile.doRollover()
    return logger

class ABase(object):
    
    """A base class that provides logging and self pickling functionality.


    .. attribute:: name

       Name of the instance, which is frequently used as a prefix to output 
       names.

    .. attribute:: workdir

       Working directory for the instance, into which outputs are written.
    
    .. attribute:: logger

       A Python logger for the instance. By default a log file is started for 
       the instance with the name :file:`workdir/name.log`.
    
    """
    
    
    def __init__(self, name, **kwargs):
        """Instantiate class using an instance name.

        :arg name: name of the class instance
        :type name: str
        
        :keyword workdir: location of all outputs and logfile
        :type workdir: str, default is "."
        
        :keyword logger: logger instance to log actions and method calls
        :type logger: logging.Logger, default None
        
        Unless an existing logger is passed as a keyword argument, a new logger
        is started for the object. All keyword arguments are passed to the 
        :meth:`get_logger` method.
        
        """
        self.name = name
        
        if 'wordkdir' in kwargs:
            workdir = kwargs['workdir']
        else: 
            workdir = '.'
        
        kwargs['workdir'] = _set_workdir(workdir)
                        
        if 'logger' in kwargs:
            self.logger = kwargs['logger']
        else:
            self.logger = get_logger(name, **kwargs)
        self.logger.info('{0:s} is initialized.'.format(self.name))
        self.set_workdir(workdir)
        
    def set_workdir(self, workdir):
        """Change working directory.
        
        If *workdir* does not exist, it is created.
        
        :arg workdir: new working directory
        :type workdir: str
        
        """
        self.workdir = _set_workdir(workdir)
        self.logger.info('{0:s} working directory is set to "{1:s}".'
                         .format(self.name, os.path.relpath(workdir)))
        
    def set_logger(self, **kwargs):
        """Setup a logger.
        
        This method can be used to reset current logger or to restart logger
        after the object is unpickled.
        
        All keyword arguments are passed to the :func:`get_logger` method.
        
        """
        if not kwargs.has_key('filemode'):
            kwargs['filemode'] = 'a'
        self.logger = get_logger(self.name, workdir=self.workdir, **kwargs)
        
    def pickle(self, filename=None, compress=True):
        """cPickle the object into a file with .dso(.gz) extension.

        Handler objects of logger attribute prevents pickling of this
        object. Hence, this method is defined. It temporarily removes 
        handlers from the logger, pickles the object, and restores the 
        handlers.
        
        To restore the object, use pickler method in Functions module. 
                
        :arg filename: name of the file to dump object (without an extension)
        :type filename: str or None, default is :attr:`ABase.name`
        
        :arg compress: gzip output
        :type compress: bool, default is True
            
        """
        if filename is None: 
            filename = os.path.join(self.workdir, self.name)
        filename += '.dso'
        if compress: 
            filename += '.gz'
            out = gzip.open(filename, 'w')
        else: 
            out = open(filename, 'w')
            
        # spare logger
        logger = self.logger
        self.logger = None
            
        pickle.dump(self, out)
        out.close()
        
        # restore logger and kdtree
        self.logger = logger
        
        self.logger.info('{0:s} is cPickled into file {1:s}.'
                         .format(self.name, os.path.relpath(filename)))
