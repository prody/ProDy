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
__version__ = '0.8.3'

release = tuple([int(x) for x in __version__.split('.')])

import logging
import logging.handlers
import os
import os.path
import cPickle
import sys
import time
import math
import platform

import warnings

DEPRECATION_WARNINGS = False
if release >= (0,9):
    warnings.filterwarnings('default', category=DeprecationWarning)

def deprecate(dep, alt, rel=(0,9)):
    """Issue a deprecation warning for *dep* and recommend using *alt*, if 
    *rel* is greater or equal to the current release."""

    if DEPRECATION_WARNINGS or release >= rel:
        # Assume that the deprecated method or function will be removed
        # when the next major release is made
        rel = float('.'.join((str(x) for x in rel[:2]))) + 0.1
        warnings.warn('`{0:s}` is deprecated and will be removed in v{1:.1f}, '
                      'use `{2:s}`.'.format(dep, rel, alt), 
                      DeprecationWarning, stacklevel=2)

def turnonDepracationWarnings():
    global DEPRECATION_WARNINGS
    DEPRECATION_WARNINGS = True
    warnings.filterwarnings('always', category=DeprecationWarning)
    
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

class PackageLogger(object):
    
    def __init__(self, name='.'+__package__, **kwargs):
        """Start logger for the package. Returns a logger instance."""
    
        self._level = logging.DEBUG
        logger = logging.getLogger(name)
        logger.setLevel(self._level)

        prefix = kwargs.get('prefix', '@>')
        assert isinstance(prefix, str), 'prefix must be as string'
        self._prefix = prefix
        
        for handler in logger.handlers: 
            handler.close()
        logger.handlers = []
        
        console = logging.StreamHandler()
        console.setLevel(LOGGING_LEVELS[kwargs.get('console', 'debug')])
        console.setFormatter(logging.Formatter(self._prefix + ' %(message)s'))
        logger.addHandler(console)
        self._logger = logger
        
        self._n = None
        self._last = None
        self._start = None
        self._barlen = None
        self._prev = None
        self._line = None

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
            
    def getPrefix(self):
        """Return string prefixed to console messages."""
        
        return self._prefix
        

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
        
        assert isinstance(n, int) and n > 0, 'n must be a positive integer'
        self._n = n
        self._last = 0
        self._start = time.time()
        self._barlen = int(kwargs.get('barlen', 30))
        self._prev = (0, 0)
        self._line = ''
    
    def report(self, i):
        """Print status to current line in the console."""
        
        n = self._n
        if self._level < logging.WARNING and n > 0 and i <= n and \
            i > self._last:
            self._last = i
            head = '>'
            if i == n:
                head = ''
            percent = 100 * i / n
            if percent > 3:
                seconds = int(math.ceil((time.time()-self._start) * (n-i)/i))
                prev = (percent, seconds)
            else:
                prev = (percent, 0)
            if self._prev == prev:
                return
            prefix = '\r' + self._prefix
            bar = ''
            barlen = self._barlen
            if barlen > 10:
                bar = int(round(percent*barlen/100.))
                bar = '='*bar + head + ' '*(barlen-bar)
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
            
    def sleep(self, seconds, msg=''):
        """Sleep for seconds while updating screen message every second. 
        Message will start with ``"Waiting for XXs"``"""
        
        assert isinstance(seconds, int), 'seconds must be an integer'
        assert isinstance(msg, str), 'msg must be a string'
        for second in range(seconds, 0, -1):
            self.write('Waiting for {0:d}s{1:s}'.format(second, msg))
            time.sleep(1)
            self.clear()
            
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
    
    :keyword backupcount: number of old *filename.log* files to save, 
        default is 1

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
            LOGGER.info("Saving existing logfile '{0:s}' and starting a new "
                        "one.".format(filename))
            logfile.doRollover()
        else:
            LOGGER.info("Logfile '{0:s}' has been started.".format(filename))


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
    debug    Everything will be printed on the console or written into logfile.
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

def checkCoordsArray(array, arg='array', cset=False, n_atoms=None, 
                     reshape=None):
    """Return array if checks pass, otherwise raise an exception."""

    assert isinstance(arg, str), 'arg must be a string'
    assert isinstance(cset, bool), 'cset must be a boolean'
    assert n_atoms is None or isinstance(n_atoms, int) and n_atoms >= 0, \
        'n_atoms must be a positive integer'
    assert reshape is None or isinstance(reshape, bool), \
        'reshape must be a boolean'

    if not isinstance(array, np.ndarray):
        raise TypeError(arg + ' must be a Numpy array')
    elif cset and array.ndim not in (2,3): 
        raise ValueError(arg + '.ndim must be  2 or 3'.format(ndim))
    elif not cset and array.ndim != 2:
        raise ValueError(arg + '.ndim must be 2')
    elif array.shape[-1] != 3:
        raise ValueError(arg + '.shape[-1] of 3, i.e. ([n_csets,]n_atoms,3)')
    if n_atoms is not None and n_atoms != 0 and array.shape[-2] != n_atoms:
        raise ValueError(arg + ' size do nut match number of atoms')
    if array.dtype != float:
        try:
            array = array.astype(float)
        except ValueError:
            raise ValueError(arg + '.astype(float) fails, float type could '
                             'not be assigned')
    if cset and reshape and array.ndim == 2:
        array = array.reshape([1, array.shape[0], 3])
    return array

__all__ = ['startLogfile', 'closeLogfile', 'changeVerbosity',
           'checkUpdates', 'plog']

def test(**kwargs):
    """Run ProDy tests. See :mod:`prody.testing` documentation for more 
    details, i.e. ``import prody.tests; help(prody.tests)``"""
    
    try:
        import prody.tests
    except ImportError:
        LOGGER.warning('Could not import test modules.')
    else:
        tests.test(**kwargs)

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
__all__.append('proteins')

from . import measure
from measure import *
__all__.extend(measure.__all__)
__all__.append('measure')

from . import compare
from compare import *
__all__.extend(compare.__all__)
__all__.append('compare')

from . import dynamics
from .dynamics import *
__all__.extend(dynamics.__all__)
__all__.append('dynamics')

from . import ensemble
from .ensemble import *
__all__.extend(ensemble.__all__)
__all__.append('ensemble')

#from . import volume
#from .volume import *
#__all__.extend(volume.__all__)

import prody
__all__.append('prody')

