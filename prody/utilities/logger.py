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

"""This module defines class that can be used a package wide logger."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import sys
import math
import time
import os.path
import logging
import datetime
import logging.handlers

__all__ = ['PackageLogger']

LOGGING_LEVELS = {'debug': logging.DEBUG,
                'info': logging.INFO,
                'warning': logging.WARNING,
                'error': logging.ERROR,
                'critical': logging.CRITICAL,
                'none': logging.CRITICAL}
for key, value in LOGGING_LEVELS.items():
    LOGGING_LEVELS[value] = key
LOGGING_LEVELS.setdefault(logging.INFO)

now = datetime.datetime.now

class PackageLogger(object):
    
    """A class for package wide logging functionality."""
    
    def __init__(self, name, **kwargs):
        """Start logger for the package. Returns a logger instance.
        
        :arg prefix: prefix to console log messages, default is ``'@> '``
        :arg console: log level for console (``sys.stderr``) messages,
            default is ``'debug'``
        :arg info: prefix to log messages at *info* level
        :arg warning: prefix to log messages at *warning* level, default is 
            ``'WARNING '``
        :arg error: prefix to log messages at *error* level, default is 
            ``'ERROR '``
        """
    
        self._level = logging.DEBUG
        self._logger = logger = logging.getLogger(name)
        logger.setLevel(self._level)

        
        for handler in logger.handlers: 
            handler.close()
        logger.handlers = []
        
        console = logging.StreamHandler()
        console.setLevel(LOGGING_LEVELS[kwargs.get('console', 'debug')])
        logger.addHandler(console)
        self.prefix = kwargs.get('prefix', '@> ')
        
        self._info = kwargs.get('info', '')
        self._warning = kwargs.get('warning', 'WARNING ')
        self._error = kwargs.get('error', 'ERROR ')
        
        self._n = None
        self._last = None
        self._start = None
        self._barlen = None
        self._prev = None
        self._line = None
        self._timer = None

    def getVerbosity(self):
        """Deprecated for removal in v1.3, get :attr:`verbosity` directly."""

        from prody import deprecate
        deprecate('getVerbosity', 'verbosity')
        
        return LOGGING_LEVELS.get(self._logger.handlers[0].level)
    
    def setVerbosity(self, level):
        """Deprecated for removal in v1.3, set :attr:`verbosity` directly."""

        from prody import deprecate
        deprecate('setVerbosity', 'verbosity')
        
        lvl = LOGGING_LEVELS.get(str(level).lower(), None)
        if lvl is None: 
            self.warning('{0:s} is not a valid log level.'.format(level))
        else:
            self._logger.handlers[0].level = lvl
            self._level = lvl 
           
    verbosity = property(getVerbosity, setVerbosity, doc=
        """Verbosity *level* of the logger, default level is **debug**.  Log 
        messages are written to ``sys.stderr``.  Following logging levers are
        recognized:
        
        ========  =============================================
        Level     Description
        ========  =============================================
        debug     Everything will be printed to the sys.stderr.
        info      Only brief information will be printed.
        warning   Only warning messages will be printed.
        none      Nothing will be printed.
        ========  =============================================""")
            
    def getPrefix(self):
        """Deprecated for removal in v1.3, get :attr:`prefix` directly."""

        from prody import deprecate
        deprecate('getPrefix', 'prefix')
        return self._prefix

    def _getprefix(self):
        
        return self._prefix

    def _setprefix(self, prefix):
        
        self._prefix = str(prefix)
        prefix += '%(message)s'
        self._logger.handlers[0].setFormatter(logging.Formatter(prefix))

    prefix = property(_getprefix, _setprefix, doc='String prepended to console'
                      ' log messages.')

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
    
    warn = warning

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
                seconds = int(math.ceil((time.time()-self._start) * (n-i)/i))
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
        """Clear current line in ``sys.stderr``."""
        
        if self._line and self._level < logging.WARNING:
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
            self._line = ''

    def write(self, line):
        """Write *line* into ``sys.stderr``."""
        
        self._line = str(line)
        if self._level < logging.WARNING:
            sys.stderr.write(self._line)
            sys.stderr.flush()
            
    def sleep(self, seconds, msg=''):
        """Sleep for seconds while updating screen message every second. 
        Message will start with ``'Waiting for Xs '``"""
        
        msg = str(msg)
        seconds = int(seconds)
        for second in range(seconds, 0, -1):
            self.write('Waiting for {0:d}s {1:s}'.format(second, msg))
            time.sleep(1)
            self.clear()
            
    def startLogfile(self, filename, **kwargs):
        """Deprecated, use :meth:`start` instead."""
        
        from prody import deprecate
        deprecate('startLogfile', 'start')
        self.start(filename, **kwargs)
        
    def start(self, filename, **kwargs):
        """Start a logfile.  If *filename* does not have an extension. 
        :file:`.log` will be appended to it. 
        
        :arg filename: name of the logfile
        :arg mode: mode in which logfile will be opened, default is "w" 
        :arg backupcount: number of existing *filename.log* files to 
            backup, default is 1"""

        filename = str(filename)
        if os.path.splitext(filename)[1] == '':
            filename += '.log'
        rollover = False 
        # if filemode='a' is provided, rollover is not performed
        if os.path.isfile(filename) and kwargs.get('filemode', None) != 'a':
            rollover = True
        logfile = logging.handlers.RotatingFileHandler(filename, 
                    mode=kwargs.get('mode', 'a'), maxBytes=0,
                    backupCount=kwargs.get('backupcount', 1))
        logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
        logfile.setFormatter(logging.Formatter('%(message)s'))
        self.info("Logging into file: {0:s}".format(filename))
        self._logger.addHandler(logfile)
        if rollover:
            logfile.doRollover()
        self.info("Logging started at {0:s}".format(str(now())))

    def closeLogfile(self, filename):
        """Deprecated, use :meth:`close` instead."""
        
        from prody import deprecate
        deprecate('closeLogfile', 'close')
        self.start(filename)
        
    def close(self, filename):
        """Close logfile *filename*."""
        
        filename = str(filename)
        if os.path.splitext(filename)[1] == '':
            filename += '.log'
        for index, handler in enumerate(self.getHandlers()):
            if isinstance(handler, logging.handlers.RotatingFileHandler):
                if handler.stream.name in (filename,os.path.abspath(filename)):
                    self.info("Logging stopped at {0:s}".format(str(now())))
                    handler.close()
                    self.delHandler(index)
                    self.info("Closing logfile: {0:s}".format(filename))
                    return
        self.warning("Logfile '{0:s}' was not found.".format(filename))

    def timeit(self):
        """Start timing a process.  Use :meth:`timing` to report time."""
        
        self._timer = time.time()
        
    def timing(self, msg='Completed in %.2fs.'):
        """Write *msg* with timing information at *debug* logging level."""
        
        self.debug(msg % (time.time() - max(self._timer, self._start)))
