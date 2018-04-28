"""This module defines class that can be used a package wide logger."""

import sys
import math
import time
import os.path
import logging
import datetime
import logging.handlers
import numbers

__all__ = ['PackageLogger', 'LOGGING_LEVELS']

LOGGING_PROGRESS = logging.INFO + 5

LOGGING_LEVELS = {'debug': logging.DEBUG,
                'info': logging.INFO,
                'progress': LOGGING_PROGRESS,
                'warning': logging.WARNING,
                'error': logging.ERROR,
                'critical': logging.CRITICAL,
                'none': logging.CRITICAL}
LOGGING_INVERSE = {}
for key, value in LOGGING_LEVELS.items(): # PY3K: OK
    LOGGING_INVERSE[value] = key

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
        self._barlen = None
        self._prev = None
        self._line = None
        self._times = {}

        self._n_progress = 0

    # ====================
    # Attributes
    # ====================

    def _getverbosity(self):

        return LOGGING_INVERSE.get(self._logger.handlers[0].level)

    def _setverbosity(self, level):
        lvl = LOGGING_LEVELS.get(str(level).lower(), None)
        if lvl is None:
            self.warn('{0} is not a valid log level.'.format(level))
        else:
            self._logger.handlers[0].level = lvl
            self._level = lvl

    verbosity = property(_getverbosity, _setverbosity, doc=
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

    def _getprefix(self):

        return self._prefix

    def _setprefix(self, prefix):

        self._prefix = str(prefix)
        prefix += '%(message)s'
        self._logger.handlers[0].setFormatter(logging.Formatter(prefix))

    prefix = property(_getprefix, _setprefix, doc='String prepended to console'
                      ' log messages.')

    # ====================
    # Logging methods
    # ====================

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
        """Log *msg* with severity 'ERROR' and terminate with status 2."""

        self.clear()
        self._logger.error(self._error + msg)
        self.exit(2)

    def write(self, line):
        """Write *line* into ``sys.stderr``."""

        self._line = str(line)
        if self._level < logging.WARNING:
            sys.stderr.write(self._line)
            sys.stderr.flush()

    def clear(self):
        """Clear current line in ``sys.stderr``."""

        if self._level != LOGGING_PROGRESS: 
            if self._line and self._level < logging.WARNING:
                sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
                self._line = ''

    def exit(self, status=0):
        """Exit the interpreter."""

        sys.exit(status)

    # ====================
    # Handlers & logfiles
    # ====================

    def addHandler(self, hdlr):
        """Add the specified handler to this logger."""

        self._logger.addHandler(hdlr)

    def getHandlers(self):
        """Returns handlers."""

        return self._logger.handlers

    def delHandler(self, index):
        """Remove handler at given *index* from the logger instance."""

        self._logger.handlers.pop(index)

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
        # if mode='a' is provided, rollover is not performed
        if os.path.isfile(filename) and kwargs.get('mode', None) != 'a':
            rollover = True
        logfile = logging.handlers.RotatingFileHandler(filename,
                    mode=kwargs.get('mode', 'a'), maxBytes=0,
                    backupCount=kwargs.get('backupcount', 1))
        logfile.setLevel(LOGGING_LEVELS[kwargs.get('loglevel', 'debug')])
        logfile.setFormatter(logging.Formatter('%(message)s'))
        self.info("Logging into file: {0}".format(filename))
        self._logger.addHandler(logfile)
        if rollover:
            logfile.doRollover()
        self.info("Logging started at {0}".format(str(now())))

    def close(self, filename):
        """Close logfile *filename*."""

        filename = str(filename)
        if os.path.splitext(filename)[1] == '':
            filename += '.log'
        for index, handler in enumerate(self.getHandlers()):
            if isinstance(handler, logging.handlers.RotatingFileHandler):
                if handler.stream.name in (filename,os.path.abspath(filename)):
                    self.info("Logging stopped at {0}".format(str(now())))
                    handler.close()
                    self.delHandler(index)
                    self.info("Closing logfile: {0}".format(filename))
                    return
        self.warning("Logfile '{0}' was not found.".format(filename))

    # ====================
    # Progress and timing
    # ====================

    def progress(self, msg, steps, label=None, **kwargs):
        """Instantiate a labeled process with message and number of steps."""

        assert isinstance(steps, numbers.Integral) and steps > 0, \
            'steps must be a positive integer'
        self._steps = steps
        self._last = 0
        self._times[label] = time.time()
        self._prev = (-1, 0)
        self._msg = msg
        self._line = ''

        if not hasattr(self, '_verb'):
            self._verb = self._getverbosity()
            self._setverbosity('progress')
        self._n_progress += 1

    def update(self, step, msg=None, label=None):
        """Update progress status to current line in the console."""

        assert isinstance(step, numbers.Integral), 'step must be a positive integer'
        if msg is not None:
            self._msg = msg
        n = self._steps
        i = step
        if self._level < logging.WARNING and n > 0 and i <= n and \
            i > self._last:
            start = self._times[label]
            self._last = i
            percent = 100 * i / n
            #if percent > 3:
            seconds = int(math.ceil((time.time()-start) * (n-i)/i))
            prev = (percent, seconds)
            #else:
                #prev = (percent, 0)
            #if self._prev == prev:
            #    return
            sys.stderr.write('\r' + ' ' * (len(self._line)) + '\r')
            #if percent > 3:
            line = self._prefix + self._msg + ' [%3d%%] %ds' % (percent, seconds)
            #else:
            #    line = self._prefix + self._msg + ' [%3d%%]' % percent
            sys.stderr.write(line)
            sys.stderr.flush()
            self._prev = prev
            self._line = line

    def finish(self):
        self._n_progress -= 1
        if self._n_progress < 0:
            self._n_progress = 0
        if self._n_progress == 0:
            if hasattr(self, '_verb'):
                self._setverbosity(self._verb)
                del self._verb
                self.clear()

    def sleep(self, seconds, msg=''):
        """Sleep for seconds while updating screen message every second.
        Message will start with ``'Waiting for Xs '`` followed by *msg*."""

        msg = str(msg)
        for second in range(int(seconds), 0, -1):
            self.write('Waiting for {0}s {1}'.format(second, msg))
            time.sleep(1)
            self.clear()

    def timeit(self, label=None):
        """Start timing a process.  Use :meth:`timing` and :meth:`report` to
        learn and report timing, respectively."""

        self._times[label] = time.time()

    def timing(self, label=None):
        """Returns timing for a labeled or default (**None**) process."""

        return time.time() - self._times.get(label, 0)

    def report(self, msg='Completed in %.2fs.', label=None):
        """Write *msg* with timing information for a labeled or default process
        at *debug* logging level."""

        self.debug(msg % (time.time() - self._times[label]))
