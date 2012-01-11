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
"""This module contains tools for handling files, logging, type checking."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
import sys
import gzip
import time
import cPickle
import logging
import os.path
import zipfile
import datetime
import platform
import logging.handlers

today = datetime.date.today
now = datetime.datetime.now
getsize = os.path.getsize

import numpy as np

pkg = __import__(__package__)

__all__ = ['PackageLogger', 'PackageSettings',
        'getPackagePath', 'setPackagePath',
        'checkCoordsArray', 
        'gunzip', 'openFile', 'openDB',
        'isExecutable', 'isReadable', 'isWritable', 
        'makePath', 'relpath', 'which', 
        'pickle', 'unpickle',
        'rangeString',
        'today', 'now', 'getsize',
        'PLATFORM', 'USERHOME',
        'alnum']

USERHOME = os.getenv('USERPROFILE') or os.getenv('HOME')

PLATFORM = platform.system()
LOGGING_LEVELS = {'debug': logging.DEBUG,
                'info': logging.INFO,
                'warning': logging.WARNING,
                'error': logging.ERROR,
                'critical': logging.CRITICAL,
                'none': logging.CRITICAL}
for key, value in LOGGING_LEVELS.items():
    LOGGING_LEVELS[value] = key
LOGGING_LEVELS.setdefault(logging.INFO)

class PackageLogger(object):
    
    def __init__(self, name, **kwargs):
        """Start logger for the package. Returns a logger instance."""
    
        self._level = logging.DEBUG
        logger = logging.getLogger(name)
        logger.setLevel(self._level)

        prefix = kwargs.get('prefix', '@> ')
        assert isinstance(prefix, str), 'prefix must be as string'
        self._prefix = prefix
        
        for handler in logger.handlers: 
            handler.close()
        logger.handlers = []
        
        console = logging.StreamHandler()
        console.setLevel(LOGGING_LEVELS[kwargs.get('console', 'debug')])
        console.setFormatter(logging.Formatter(self._prefix + '%(message)s'))
        logger.addHandler(console)
        self._logger = logger
        
        self._info = kwargs.get('info', '')
        self._warning = kwargs.get('warning', 'WARNING ')
        self._error = kwargs.get('error', 'ERROR ')
        
        self._n = None
        self._last = None
        self._start = None
        self._barlen = None
        self._prev = None
        self._line = None

    def getVerbosity(self):
        """Return verbosity *level* of the logger."""
        
        return LOGGING_LEVELS.get(self._logger.handlers[0].level)
    
    def setVerbosity(self, level):
        """Change verbosity *level* of the logger for the current session.  
        Default verbosity level **debug**.  Log messages are written to 
        ``sys.stderr``.  This function accepts one of the following as 
        *level* argument:
        
        ========  ===========================================
        Level     Description
        ========  ===========================================
        debug     Everything will be printed on the console.
        info      Only brief information will be printed.
        warning   Only warning information will be printed.
        none      ProDy will not log any messages.
        ========  ==========================================="""
        
        lvl = LOGGING_LEVELS.get(level, None)
        if lvl is None: 
            self.warning('{0:s} is not a valid log level.'.format(level))
        else:
            self._logger.handlers[0].level = lvl
            self._level = lvl 
            
    def getPrefix(self):
        """Return string prefixed to console messages."""
        
        return self._prefix

    def info(self, msg):
        """Log *msg* with severity 'INFO'."""

        self.clear()
        self._logger.info(msg)

    def critical(self, msg):
        """Log *msg* with severity 'CRITICAL'."""
        
        self.clear()
        self._logger.critical(msg)

    def debug(self, msg):
        """Log *msg* with severity 'DEBUG'."""

        self.clear()
        self._logger.debug(msg)
        
    def warning(self, msg):
        """Log *msg* with severity 'WARNING'."""
        
        self.clear()
        self._logger.warning(self._warning + msg)

    def error(self, msg):
        """Log *msg* with severity 'ERROR'."""
        
        self.clear()
        self._logger.error(self._error + msg)
    
    def addHandler(self, hdlr):
        """Add the specified handler to this logger."""
        
        self._logger.addHandler(hdlr)
        
    def getHandlers(self):
        """Return handlers."""
        
        return self._logger.handlers

    def delHandler(self, index):
        """Remove handler at given *index* from the logger instance."""
        
        self._logger.handlers.pop(index)

    def progress(self, msg, steps, **kwargs):
        """Instantiate with message number of steps."""
        
        assert isinstance(steps, int) and steps > 0, \
            'steps must be a positive integer'
        self._steps = steps
        self._last = 0
        self._start = time.time()
        self._prev = (0, 0)
        self._msg = msg
        self._line = ''
    
    def update(self, step):
        """Update progress status to current line in the console."""
        
        assert isinstance(step, int), 'step must be a positive integer'
        n = self._steps
        i = step
        if self._level < logging.WARNING and n > 0 and i <= n and \
            i > self._last:
            self._last = i
            percent = 100 * i / n
            if percent > 3:
                seconds = int(np.ceil((time.time()-self._start) * (n-i)/i))
                prev = (percent, seconds)
            else:
                prev = (percent, 0)
            if self._prev == prev:
                return
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
            if percent > 3:
                line = self._prefix + self._msg + \
                    ' [%3d%%] %ds' % (percent, seconds)
            else:
                line = self._prefix + self._msg + ' [%3d%%]' % percent
            sys.stderr.write(line)
            sys.stderr.flush()
            self._prev = prev
            self._line = line
    
    def clear(self):
        """Clear sys.stderr."""
        
        if self._line and self._level < logging.WARNING:
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
            self._line = ''

    def write(self, line):
        """Write *line* to sys.stderr."""
        
        self._line = str(line)
        if self._level < logging.WARNING:
            sys.stderr.write(self._line)
            sys.stderr.flush()
            
    def sleep(self, seconds, msg=''):
        """Sleep for seconds while updating screen message every second. 
        Message will start with ``"Waiting for XXs"``"""
        
        assert isinstance(seconds, int), 'seconds must be an integer'
        assert isinstance(msg, str), 'msg must be a string'
        for second in range(seconds, 0, -1):
            self.write('Waiting for {0:d}s{1:s}'.format(second, msg))
            time.sleep(1)
            self.clear()
            
    def startLogfile(self, filename, **kwargs):
        """Start a file to save logs.
        
        :keyword filename: name of the logfile
        :keyword mode: mode in which logfile will be opened, default is "w" 
        :keyword backupcount: number of existing *filename.log* files to 
            backup, default is 1
        """

        assert isinstance(filename, str), 'filename must be a string'
        logfilename = filename
        if not logfilename.endswith('.log'):
            logfilename += '.log'
        rollover = False 
        # if filemode='a' is provided, rollover is not performed
        if os.path.isfile(logfilename) and kwargs.get('filemode', None) != 'a':
            rollover = True
        logfile = logging.handlers.RotatingFileHandler(logfilename, 
                    mode=kwargs.get('mode', 'a'), maxBytes=0,
                    backupCount=kwargs.get('backupcount', 1))
        logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
        logfile.setFormatter(logging.Formatter('%(message)s'))
        self._logger.addHandler(logfile)
        if rollover:
            logfile.doRollover()
        self.info("Logging into '{0:s}'.".format(logfilename))

    def closeLogfile(self, filename):
        """Close log file *filename*."""
        
        filename = str(filename)
        if not filename.endswith('.log'):
            filename += '.log'
        for index, handler in enumerate(self.getHandlers()):
            if isinstance(handler, logging.handlers.RotatingFileHandler):
                if handler.stream.name in (filename,os.path.abspath(filename)):
                    self.info("Closing logfile '{0:s}'".format(filename))
                    handler.close()
                    self.delHandler(index)
                    return
        self.warning("Logfile '{0:s}' was not found.".format(filename))

    def timeit(self):
        """Start timing a process.  Use :meth:`timing` to report time."""
        
        self._start = time.time()
        
    def timing(self, msg=None):
        """If *msg* is none, return time passes since timing started. If 
        a message is given, e.g. ``"Completed in %.2fs."``, report the time 
        it took to complete the process at *debug* log level."""
        
        if msg is None:
            return time.time() - self._start
        else:
            self.debug(msg % (time.time() - self._start))
        

class PackageSettings(object):
    
    """A class for managing package settings.  Settings are saved in user's 
    home director.  When settings are changed by the users, the changes are 
    automatically saved.  Settings are stored in a :class:`dict` instance. 
    The dictionary is pickled in user's home directory for permanent storage.
    """
    
    def __init__(self, pkg=__package__, rcfile=None, logger=None):
        """*rcfile* is the filename for pickled settings dictionary, and by 
        default is set to :file:`.pkgrc`."""
        
        self._package = pkg
        if rcfile is None:
            self._rcfile = os.path.join(USERHOME, '.' + pkg + 'rc')
        else:
            self._rcfile = rcfile
        if isinstance(logger, PackageLogger):
            self._logger = logger
        else:
            self._logger = None
        
        self._settings = {}
        
    def __getitem__(self, key):
        
        return self._settings[key]
        
    def __setitem__(self, key, value):
        
        self._settings[key] = value
        
    def get(self, key, default=None):
        
        return self._settings.get(key, default)
        
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

def checkCoordsArray(array, arg='array', cset=False, n_atoms=None, 
                    reshape=None, dtype=(float,)):
    """Return array if checks pass, otherwise raise an exception."""

    assert isinstance(arg, str), 'arg must be a string'
    assert isinstance(cset, bool), 'cset must be a boolean'
    assert n_atoms is None or isinstance(n_atoms, int) and n_atoms >= 0, \
        'n_atoms must be a positive integer'
    assert reshape is None or isinstance(reshape, bool), \
        'reshape must be a boolean'
    if not isinstance(dtype, tuple):
        dtype = (dtype, )

    if not isinstance(array, np.ndarray):
        raise TypeError(arg + ' must be a Numpy array')
    elif cset and array.ndim not in (2,3): 
        raise ValueError(arg + '.ndim must be  2 or 3'.format(ndim))
    elif not cset and array.ndim != 2:
        raise ValueError(arg + '.ndim must be 2')
    elif array.shape[-1] != 3:
        raise ValueError(arg + '.shape[-1] of 3, i.e. ([n_csets,]n_atoms,3)')
    if n_atoms is not None and n_atoms != 0 and array.shape[-2] != n_atoms:
        raise ValueError(arg + ' size do not match number of atoms')
    if array.dtype not in dtype:
        try:
            array = array.astype(dtype[0])
        except ValueError:
            raise ValueError(arg + '.astype(' + str(dtype[0]) + ') fails, '
                            'float type could not be assigned')
    if cset and reshape and array.ndim == 2:
        array = array.reshape([1, array.shape[0], 3])
    return array

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
        filename = os.path.join(folder, filename)
    ext = os.path.splitext(filename)[1]
    backup = kwargs.pop('backup', pkg.SETTINGS.get('backup', False))
    backup_ext = kwargs.pop('backup_ext', 
                            pkg.SETTINGS.get('backup_ext', '.BAK'))
    if args and args[0][0] in ('a', 'w'):
        if os.path.isfile(filename) and backup:
                os.rename(filename, filename + backup_ext)
    return OPEN.get(ext, open)(filename, *args, **kwargs)
    
    
def gunzip(filename, outname=None):
    """Decompress *filename* and save as *outname*.  When ``outname=None``, 
    *filename* is used as the output name.  Returns output filename upon 
    successful completion."""

    if not isinstance(filename, str):
        raise TypeError('filename must be a string')
    if not os.path.isfile(filename):
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
    
    return isinstance(path, str) and os.path.exists(path) and \
        os.access(path, os.X_OK)

def isReadable(path):
    """Return true if *path* is readable by the user."""
    
    return isinstance(path, str) and os.path.exists(path) and \
        os.access(path, os.R_OK)


def isWritable(path):
    """Return true if *path* is writable by the user."""
    
    return isinstance(path, str) and os.path.exists(path) and \
        os.access(path, os.W_OK)


def relpath(path):
    """Return *path* on Windows, and relative path elsewhere."""
    
    if PLATFORM == 'Windows':
        return path
    else:
        return os.path.relpath(path)

def makePath(path):
    """Make all directories that does not exist in a given path."""
    
    if os.path.isabs(path):
        path = relpath(path)
    if not os.path.isdir(path):
        dirs = path.split(os.sep)
        for i in range(len(dirs)):
            dirname = os.sep.join(dirs[:i+1])
            try:
                if not os.path.isdir(dirname): 
                    os.mkdir(dirname)
            except OSError:
                raise OSError('{0:s} could not be created, please '
                            'specify another path'.format(path))
                return os.getcwd()
    return os.path.join(os.getcwd(), path)


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
    """Pickle *obj* using :mod:`cPickle` and dump in *filename*."""
    
    out = openFile(filename, 'wb')
    cPickle.dump(obj, out)
    out.close()
    return filename

def unpickle(filename, **kwargs):
    """Unpickle object in *filename* using :mod:`cPickle`."""
    
    inf = openFile(filename, 'rb')
    obj = cPickle.load(inf)
    inf.close()
    return obj

def rangeString(lint, sep=' ', rng=' to '):
    """Return a structured string for a given list of integers.
    
    :arg lint: integer list or array
    :arg sep: range or number separator         
    :arg rng: inclusive range symbol

    E.g. for ``sep=' '`` and ``rng=' to '``: 
        ``[1, 2, 3, 4, 10, 15, 16, 17]`` -> ``"1 to 4 10 15 to 17"``
    for ``sep=','`` and ``rng='-'``:
        ``[1, 2, 3, 4, 10, 15, 16, 17]`` -> ``"1-4,10,15-17"``
    """
    lint = np.unique(lint)
    strint = ''
    i = -1
    for j in lint:
        if j < 0:
            continue
        if i < 0:
            i = j
        diff = j - i
        if diff == 0:
            strint += str(j)
        elif diff > 1: 
            strint += rng + str(i) + sep + str(j)
            k = j 
        i = j
    if diff == 1: 
        strint += rng + str(i)
    elif diff > 1 and k != j: 
        strint += rng + str(i) + sep + str(j)
    return strint

def openDB(filename, *args):
    import anydbm
    return anydbm.open(filename, *args)

def alnum(string, alt='_'):
    """Replace non alpha numeric characters with *alt*."""
    
    result = ''
    for char in string:
        if char.isalnum():
            result += char
        else:
            result += alt
    return result
