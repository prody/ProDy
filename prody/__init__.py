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
__version__ = '0.7.2'

import logging
import logging.handlers
import os
import os.path
import cPickle
import sys
import time
import math

_PY3K = sys.version_info[0] > 2

def isExecutable(path):
    return os.path.exists(path) and os.access(path, os.X_OK)

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
        ProDyLogger.warning('Matplotlib is required for plotting.')

def importBioKDTree():
    try:
        from prody.KDTree import KDTree
    except ImportError:
        try:
            from Bio.KDTree import KDTree
        except ImportError:
            dynamics.KDTree = False
            select.KDTree = False
            ProDyLogger.warning('KDTree module could not be imported. '
                                'Reinstalling ProDy or installing BioPython '
                                'can resolve the problem.')
        else:
            dynamics.KDTree = KDTree
            select.KDTree = KDTree
    else:
        dynamics.KDTree = KDTree
        select.KDTree = KDTree

def importBioBlast():
    import bioblast
    proteins.BioBlast = bioblast

def importBioPairwise2():
    try:
        import pairwise2
    except ImportError:
        try:
            from Bio import pairwise2
        except ImportError:
            compare.pairwise2 = False
            ProDyLogger.warning('pairwise2 module could not be imported. '
                                'Reinstalling ProDy or installing BioPython '
                                'can resolve the problem.')
        else:
            compare.pairwise2 = pairwise2
    else:
        compare.pairwise2 = pairwise2

_ProDySettings = None 

def _loadProDySettings():
    fn = os.path.join(os.path.expanduser('~'), '.prody')
    settings = None
    if os.path.isfile(fn):
        inf = open(fn)
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
    fn = os.path.join(os.path.expanduser('~'), '.prody')
    out = open(fn, 'w')
    cPickle.dump(_ProDySettings, out)
    out.close()
    
_loadProDySettings()

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
ProDySignature = '@>'

def _startLogger(**kwargs):
    """Set a global logger for ProDy."""
    
    name = '.prody'
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    
    for handler in logger.handlers: 
        handler.close()
    logger.handlers = []
    
    console = logging.StreamHandler()
    console.setLevel(LOGGING_LEVELS[kwargs.get('console', 'debug')])
    console.setFormatter(logging.Formatter(ProDySignature + ' %(message)s'))
    logger.addHandler(console)
    return logger

ProDyLogger = _startLogger()

def plog(*text):
    """Log *text* using ProDy logger at log level info.
    
    .. versionadded:: 0.6.2
    
    .. versionchanged:: 0.7
       Multiple arguments are accepted. They will converted to string and
       joined using a white space as delimiter. 
    
    """
    
    ProDyLogger.info(' '.join([str(s) for s in text]))

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
    logger = ProDyLogger
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
        logger.addHandler(logfile)
        if rollover:
            logger.info('Saving existing logfile "{0:s}" and starting a new one.'
                         .format(filename))
            logfile.doRollover()
        else:
            logger.info('Logfile "{0:s}" has been started.'.format(filename))


def closeLogfile(filename):
    """Close logfile with *filename*."""
    filename = str(filename)
    if not filename.endswith('.log'):
        filename += '.log'
    logger = ProDyLogger
    for index, handler in enumerate(logger.handlers):
        if isinstance(handler, logging.handlers.RotatingFileHandler):
            if handler.stream.name in (filename, os.path.abspath(filename)):
                logger.info('Closing logfile {0:s}'.format(filename))
                handler.close()
                logger.handlers.pop(index)
                return
    logger.warning('Logfile "{0:s}" was not found.'.format(filename))

class ProDyProgress(object):
    
    """A class to print progress to sys.stderr."""
    
    def __init__(self, n, **kwargs):
        """Instantiate with number of steps."""
        
        self._n = n
        self._level = getVerbosityLevel()
        self._start = time.time()
        self._prefix = str(kwargs.get('prefix', ''))
        self._barlen = int(kwargs.get('barlen', 30))
        self._prev = (0, 0)
    
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
            prefix = ProDySignature + self._prefix
            bar = ''
            barlen = self._barlen
            if barlen > 10:
                barlen = int(round(0.2*percent))
                bar = ' [' + '=' * (barlen-1) + '>' + ' ' * (20-barlen) + '] '
            if percent > 3:
                sys.stderr.write(('\r' + prefix + ' %2d%%' + bar + '%ds    ') % 
                 (percent, seconds))
            else:
                sys.stderr.write(('\r' + prefix + ' %2d%%' + bar) % percent) 
            sys.stderr.flush()
            self._prev = prev
    
    def clean(self):
        """Clean sys.stderr."""
        
        if self._level < logging.WARNING:
            sys.stderr.write('\r' + ' ' * (self._barlen + 20) + '\r')
    

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
    lvl = LOGGING_LEVELS.get(level, None)
    if lvl is None: 
        ProDyLogger.warning('{0:s} is not a valid log level.'.format(level))
    else:
        ProDyLogger.handlers[0].level = lvl

def getVerbosityLevel():
    """Return ProDy console verbosity level."""
    
    return ProDyLogger.handlers[0].level

def checkUpdates():
    """Check latest ProDy release and compare with the installed one.
    
    .. versionadded:: 0.5.3
    
    User is informed whether or not latest version is installed.
    
    """
    
    import xmlrpclib
    pypi = xmlrpclib.Server('http://pypi.python.org/pypi')
    releases = pypi.package_releases('ProDy')
    if releases[0] == __version__:
        ProDyLogger.info('You are using the latest ProDy release ({0:s}).'
                         .format(__version__))
    else:
        ProDyLogger.info('ProDy {0:s} has been released. '
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

from . import volume
from .volume import *
__all__.extend(volume.__all__)


import prody

__all__.append('prody')
