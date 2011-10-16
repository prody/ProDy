# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'
__version__ = '0.8.2'

import logging
import logging.handlers
import os
import os.path
import cPickle
import sys
import time
import math
import platform

PLATFORM = platform.system()
_PY3K = sys.version_info[0] > 2

USERHOME = os.getenv('USERPROFILE') or os.getenv('HOME')
PACKAGE_PATH = os.path.join(USERHOME, '.' + __package__)
PACKAGE_CONF =  os.path.join(USERHOME, '.' + __package__ + 'rc')
if not os.path.isfile(PACKAGE_CONF) and os.path.isfile(PACKAGE_CONF[:-2]):
    os.rename(PACKAGE_CONF[:-2], PACKAGE_CONF)

def isExecutable(path):
    return os.path.exists(path) and os.access(path, os.X_OK)

def relpath(path):
    """Return *path* on Windows, and relative path elsewhere."""
    
    if PLATFORM == 'Windows':
        return path
    else:
        return os.path.relpath(path)

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

def importLA():
    try:
        import scipy.linalg as linalg
        dynamics.scipyla = True
    except ImportError:
        dynamics.scipyla = False
        try:
            import numpy.linalg as linalg
        except:
            raise ImportError('scipy.linalg or numpy.linalg is required for '
                              'NMA calculations and aligning structures.')
    dynamics.linalg = linalg
    measure.linalg = linalg

def importScipySparse():
    try:
        from scipy import sparse
    except ImportError:    
        raise ImportError('scipy.sparse is required, but could not be '
                          'imported.')
    dynamics.scipy_sparse = sparse
    
def importScipySparseLA():
    try:
        from scipy.sparse import linalg
    except ImportError:    
        raise ImportError('scipy.sparse is required, but could not be '
                          'imported.')
    dynamics.scipy_sparse_la = linalg

def importPyPlot():
    try:
        import matplotlib.pyplot as pyplot
        from mpl_toolkits.mplot3d import Axes3D
        dynamics.plt = pyplot
        ensemble.plt = pyplot
        dynamics.Axes3D = Axes3D 
    except ImportError:
        dynamics.plt = False
        ensemble.plt = False
        LOGGER.warning('Matplotlib is required for plotting.')

def importBioKDTree():
    try:
        from prody.KDTree import KDTree
    except ImportError:
        try:
            from Bio.KDTree import KDTree
        except ImportError:
            dynamics.KDTree = False
            select.KDTree = False
            LOGGER.warning('KDTree module could not be imported. '
                           'Reinstalling ProDy or installing BioPython '
                           'can resolve the problem.')
        else:
            dynamics.KDTree = KDTree
            select.KDTree = KDTree
    else:
        dynamics.KDTree = KDTree
        select.KDTree = KDTree

def importBioPairwise2():
    try:
        import pairwise2
    except ImportError:
        try:
            from Bio import pairwise2
        except ImportError:
            compare.pairwise2 = False
            LOGGER.warning('pairwise2 module could not be imported. '
                           'Reinstalling ProDy or installing BioPython '
                           'can resolve the problem.')
        else:
            compare.pairwise2 = pairwise2
    else:
        compare.pairwise2 = pairwise2

_ProDySettings = None 

def _loadProDySettings():
    settings = None
    if os.path.isfile(PACKAGE_CONF):
        inf = open(PACKAGE_CONF)
        settings = cPickle.load(inf)
        inf.close()
        if not isinstance(settings, dict):
            settings = None
    if settings is None:
        settings = {
        'loglevel': 'debug',
        }
    global _ProDySettings
    _ProDySettings = settings

def _saveProDySettings():
    out = open(PACKAGE_CONF, 'wb')
    cPickle.dump(_ProDySettings, out)
    out.close()
    
_loadProDySettings()

def setPackagePath(path):
    if not os.path.isdir(path):
        try:
            os.mkdir(path)
        except Exception as err:
            LOGGER.warning('Failed to make folder "{0:s}": {1:s}'
                           .format(path, err.strerror))
            return False
    _ProDySettings['package_path'] = path
    _saveProDySettings()
    return path    

def getPackagePath():
    
    path = _ProDySettings.get('package_path', None)
    
    update = False
    if path is None:
        LOGGER.warning('{0:s} package path is not yet set by the user.'
                       .format(__package__))
        update = True
    elif not os.path.isdir(path):
        LOGGER.warning('{0:s} package path "{1:s}" does not exist.'
                       .format(__package__, path))
        update = True
    elif not os.access(path, os.W_OK):
        LOGGER.warning('User does not have write access to {0:s} package path '
                       ' "{1:s}".'
                       .format(__package__, path))
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
    

#def ProDyPrintSettings():
#    print _ProDySettings

#ProDySettings = property(ProDySettings) 

LOGGING_LEVELS = {'debug': logging.DEBUG,
                  'info': logging.INFO,
                  'warning': logging.WARNING,
                  'error': logging.ERROR,
                  'critical': logging.CRITICAL,
                  'none': logging.CRITICAL}
LOGGING_LEVELS.setdefault(logging.INFO)
PackageSignature = '@>'

class PackageLogger(object):
    
    def __init__(self, name='.'+__package__, **kwargs):
        """Start logger for the package. Returns a logger instance."""
    
        self._level = logging.DEBUG
        logger = logging.getLogger(name)
        logger.setLevel(self._level)
        
        for handler in logger.handlers: 
            handler.close()
        logger.handlers = []
        
        console = logging.StreamHandler()
        console.setLevel(LOGGING_LEVELS[kwargs.get('console', 'debug')])
        console.setFormatter(logging.Formatter(PackageSignature + 
                                               ' %(message)s'))
        logger.addHandler(console)
        self._logger = logger

    def getVerbosityLevel(self):
        """Return verbosity level of the logger."""
        
        return self._logger.handlers[0].level
    
    def setVerbosityLevel(self, level):
        """Set verbosity level of the logger."""
        
        lvl = LOGGING_LEVELS.get(level, None)
        if lvl is None: 
            self.warning('{0:s} is not a valid log level.'.format(level))
        else:
            self._logger.handlers[0].level = lvl
            self._level = lvl 

    def info(self, msg):
        """Log *msg* with severity 'INFO'."""
        
        self._logger.info(msg)

    def critical(self, msg):
        """Log *msg* with severity 'CRITICAL'."""
        
        self._logger.critical(msg)

    def debug(self, msg):
        """Log *msg* with severity 'DEBUG'."""
        
        self._logger.debug(msg)
        
    def warning(self, msg):
        """Log *msg* with severity 'WARNING'."""
        
        self._logger.warning(msg)

    def error(self, msg):
        """Log *msg* with severity 'ERROR'."""
        
        self._logger.error(msg)
    
    def addHandler(self, hdlr):
        """Add the specified handler to this logger."""
        
        self._logger.addHandler(hdlr)
        
    def getHandlers(self):
        """Return handlers."""
        
        return self._logger.handlers

    def delHandler(self, index):
        """Remove handler at given *index* from the logger instance."""
        
        self._logger.handlers.pop(index)

    def progress(self, n, **kwargs):
        """Instantiate with number of steps."""
        
        self._n = n
        self._level = getVerbosityLevel()
        self._start = time.time()
        self._prefix = str(kwargs.get('prefix', ''))
        self._barlen = int(kwargs.get('barlen', 30))
        self._prev = (0, 0)
        self._line = ''
    
    def report(self, i):
        """Print status to current line in the console."""
        
        if self._level < logging.WARNING:   
            n = self._n
            percent = 100 * i / n
            if percent > 3:
                seconds = int(math.ceil((time.time()-self._start) * (n-i)/i))
                prev = (percent, seconds)
            else:
                prev = (percent, 0)
            if self._prev == prev:
                return
            prefix = '\r' + PackageSignature + self._prefix
            bar = ''
            barlen = self._barlen
            if barlen > 10:
                bar = int(round(percent*barlen/100.))
                bar = '='*bar + '>' + ' '*(barlen-bar)
                bar = ' [' + bar  + '] '
            
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
            if percent > 3:
                line = (prefix + ' %2d%%' + bar + '%ds')%(percent, seconds)
            else:
                line = (prefix + ' %2d%%' + bar) % percent
            sys.stderr.write(line)
            sys.stderr.flush()
            self._prev = prev
            self._line = line
    
    def clear(self):
        """Clear sys.stderr."""
        
        if self._level < logging.WARNING:
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')

    def write(self, line):
        """Write a line to sys.stderr."""
        
        self._line = str(line)
        if self._level < logging.WARNING:
            sys.stderr.write(self._line)
            sys.stderr.flush()

LOGGER = PackageLogger()

def plog(*text):
    """Log *text* using ProDy logger at log level info.
    
    .. versionadded:: 0.6.2
    
    .. versionchanged:: 0.7
       Multiple arguments are accepted. They will converted to string and
       joined using a white space as delimiter. 
    
    """
    
    LOGGER.info(' '.join([str(s) for s in text]))

def startLogfile(filename, **kwargs):
    """Start a file to save ProDy logs.
    
    :keyword filename: name of the logfile
    
    :keyword mode: mode in which logfile will be opened, default is "w" 
    
    :keyword backupcount: number of old *filename.log* files to save, default is 1

    """
    
    """
    :keyword loglevel: loglevel for logfile verbosity
    :type loglevel: str, default is *debug*
    """
    
    logfilename = filename
    if not logfilename.endswith('.log'):
        logfilename += '.log'
    if isinstance(logfilename, str):
        rollover = False 
        # if filemode='a' is provided, rollover is not performed
        if os.path.isfile(logfilename) and kwargs.get('filemode', None) != 'a':
            rollover = True
        logfile = logging.handlers.RotatingFileHandler(logfilename, 
                    mode=kwargs.get('mode', 'a'), maxBytes=0,
                    backupCount=kwargs.get('backupcount', 1))
        logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
        logfile.setFormatter(logging.Formatter('%(message)s'))
        LOGGER.addHandler(logfile)
        if rollover:
            LOGGER.info('Saving existing logfile "{0:s}" and starting a new one.'
                         .format(filename))
            logfile.doRollover()
        else:
            LOGGER.info('Logfile "{0:s}" has been started.'.format(filename))


def closeLogfile(filename):
    """Close logfile with *filename*."""
    filename = str(filename)
    if not filename.endswith('.log'):
        filename += '.log'
    for index, handler in enumerate(LOGGER.getHandlers()):
        if isinstance(handler, logging.handlers.RotatingFileHandler):
            if handler.stream.name in (filename, os.path.abspath(filename)):
                LOGGER.info('Closing logfile {0:s}'.format(filename))
                handler.close()
                LOGGER.delHandler(index)
                return
    LOGGER.warning('Logfile "{0:s}" was not found.'.format(filename))


def changeVerbosity(level):
    """Set ProDy console verbosity *level*.
    
    By default, console verbosity *level* is debug. This function accepts
    one of the following:
    
    ======== ===========
    Level    Description
    ======== ===========
    debug    Eveything will be printed on the colsole or written into logfile.
    info     Only brief information will be printed or written.
    warning  Only warning information will be printed or written.
    none     ProDy will not log any messages.
    ======== ===========
    
    >>> from prody import *
    >>> changeVerbosity('none')
    >>> plog('test')
    >>> changeVerbosity('debug')
    >>> plog('test')
    
    """

    LOGGER.setVerbosityLevel(level)

def getVerbosityLevel():
    """Return ProDy console verbosity level."""
    
    return LOGGER.getVerbosityLevel()

def checkUpdates():
    """Check latest ProDy release and compare with the installed one.
    
    .. versionadded:: 0.5.3
    
    User is informed whether or not latest version is installed.
    
    """
    
    import xmlrpclib
    pypi = xmlrpclib.Server('http://pypi.python.org/pypi')
    releases = pypi.package_releases('ProDy')
    if releases[0] == __version__:
        LOGGER.info('You are using the latest ProDy release ({0:s}).'
                         .format(__version__))
    else:
        LOGGER.info('ProDy {0:s} has been released. '
                    'You are using release {1:s}.'.format(releases[0],
                                                          __version__))

try:
    import numpy as np
except ImportError, err:
    raise ImportError('numpy not found, it is a required package')

class ProDyException(Exception):
    pass

__all__ = ['startLogfile', 'closeLogfile', 'changeVerbosity',
           'checkUpdates', 'plog']

def test(verbosity=2, descriptions=True, stream=sys.stderr):
    
    if sys.version_info[:2] > (2,6):    
        try:
            import prody.tests
        except ImportError:
            LOGGER.warning('Could not import test modules, you might be '
                           'using a Windows machine.')
        else:
            tests.test(verbosity, descriptions, stream)
    else:
        LOGGER.warning('Unit tests are compatible with Python 2.7 and later.')
        

from . import atomic 
from atomic import *
__all__.extend(atomic.__all__)
__all__.append('atomic')

from . import select
from select import *
__all__.extend(select.__all__)
__all__.append('select')
ProDyAtomSelect = Select()

import proteins
from proteins import *  
__all__.extend(proteins.__all__)

from . import measure
from measure import *
__all__.extend(measure.__all__)

from . import compare
from compare import *
__all__.extend(compare.__all__)

from . import dynamics
from .dynamics import *
__all__.extend(dynamics.__all__)

from . import ensemble
from .ensemble import *
__all__.extend(ensemble.__all__)

#from . import volume
#from .volume import *
#__all__.extend(volume.__all__)


import prody

__all__.append('prody')
