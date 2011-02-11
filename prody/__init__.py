# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010  Ahmet Bakan
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
__copyright__ = 'Copyright (C) 2010 Ahmet Bakan'
__version__ = '0.5.3'

import logging
import logging.handlers
import os
import os.path


    

re = None
def importRE():
    import re as RE
    global re
    re = RE

def importLA():
    try:
        import scipy.linalg as linalg
        dynamics.scipyla = True
    except ImportError:
        dynamics.scipyla = False
        try:
            import numpy.linalg as linalg
            ProDyLogger.warning('scipy.linalg is not found, numpy.linalg is used instead.')
        except:
            raise ImportError('scipy.linalg or numpy.linalg is required for NMA calculations and superpositioning.')
    dynamics.linalg = linalg
    measure.linalg = linalg
    
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
        from Bio.KDTree import KDTree
    except ImportError:
        dynamics.KDTree = False
        select.KDTree = False
        ProDyLogger.warning('Bio.KDTree is not found.')
    else:
        dynamics.KDTree = KDTree
        select.KDTree = KDTree

NCBIWWW = None
NCBIXML = None
def importBioBlast():
    try:
        from Bio.Blast import NCBIWWW as ncbiwww 
        from Bio.Blast import NCBIXML as ncbixml 
    except ImportError as err:
        raise ImportError('BioPython is required for Blast searching PDB.org.\n' + str(err))
    global NCBIWWW, NCBIXML
    NCBIWWW = ncbiwww
    NCBIXML = ncbixml

def importBioPairwise2():
    try:
        from Bio import pairwise2
        compare.pairwise2 = pairwise2
    except ImportError:
        compare.pairwise2 = False
        ProDyLogger.warning('BioPython is required for pairwise alignments.')

_ProDySettings = {
    'loglevel': 'debug',
    'local_pdb_folder': None}


def ProDyPrintSettings():
    print _ProDySettings

#ProDySettings = property(ProDySettings) 

LOGGING_LEVELS = {'debug': logging.DEBUG,
                  'info': logging.INFO,
                  'warning': logging.WARNING,
                  'error': logging.ERROR,
                  'critical': logging.CRITICAL,
                  'none': logging.CRITICAL}
LOGGING_LEVELS.setdefault(logging.INFO)
ProDySignature = '@>'

def _ProDyStartLogger(**kwargs):
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

ProDyLogger = _ProDyStartLogger()

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
            logger.info('Backing up existing logfile "{0:s}" and starting a new one.'
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
    >>> changeVerbosity('debug')
    
    """
    lvl = LOGGING_LEVELS.get(level, None)
    if lvl is None: 
        ProDyLogger.warning('{0:s} is not a valid log level.'.format(level))
    else:
        ProDyLogger.handlers[0].level = lvl 

def checkUpdates():
    """Check latest ProDy release and compare with the installed one.
    
    .. versionadded:: 0.6
    
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
           'checkUpdates']

from . import atomic 
from atomic import *
__all__ += atomic.__all__

from . import select
from select import *
__all__ += select.__all__
ProDyAtomSelect = Select()

import proteins
from proteins import *  
__all__ += proteins.__all__

from . import measure
from measure import *
__all__ += measure.__all__

from . import compare
from compare import *
__all__ += compare.__all__

from . import dynamics
from .dynamics import *
__all__ += dynamics.__all__

from . import ensemble
from .ensemble import *
__all__ += ensemble.__all__

try:
    from . import future
except ImportError:
    pass
else:
    from .future import *
    __all__ += future.__all__

import prody

__all__.append('prody')


def ProDySetLocalPDBFolder(path):
    pass
