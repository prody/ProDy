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

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'
__version__ = '1.0.4'

release = [int(x) for x in __version__.split('.')]

import os
import sys
import os.path
import platform

if not (2,6) <= sys.version_info[:2] <= (2,7):
    raise Exception("prody is compatible with Python 2.6 and 2.7, you are "
                    "using Python " + platform.python_version())

try:
    import numpy as np
except ImportError:
    raise ImportError('numpy not found, it is a required package')

from utilities import *

import warnings

DEPRECATION_WARNINGS = False

CONFIGURATION = {
    'backup': False,
    'backup_ext': '.BAK',
    'ligand_xml_save': False,
    'typo_warnings': True,
    'check_updates': 0,
    'auto_secondary': False,
}

from datetime import date
today = date.today()

if today.year == 2012:
    warnings.filterwarnings('default', category=DeprecationWarning)

def deprecate(dep, alt, ver='1.2', sl=3):
    """Issue a deprecation warning for *dep* and recommend using *alt*."""

    warnings.warn('`{0:s}` is deprecated for removal in v{1:s}, use `{2:s}`.'
                  .format(dep, ver, alt), DeprecationWarning, stacklevel=sl)

def turnonDepracationWarnings(action='always'):
    """Turn on deprecation warnings for the current session.  By default
     (``action='always'``), deprecation warnings will be printed every time
     a function is called to help identification of multiple occurrence 
     of deprecated function and method names.  When ``action='default'``
     is passed, warning will be issued at the first call of a function.
     The latter behavior will automatically kick in when v0.9 is released.
     Until v0.9 is released, restarting the session will turn of warnings.
     This function must be called as ``prody.turnonDepracationWarnings``."""
    
    global DEPRECATION_WARNINGS
    DEPRECATION_WARNINGS = True
    warnings.filterwarnings(action, category=DeprecationWarning)

    
_PY3K = sys.version_info[0] > 2

PACKAGECONF =  os.path.join(USERHOME, '.' + __package__ + 'rc')
if not os.path.isfile(PACKAGECONF) and os.path.isfile(PACKAGECONF[:-2]):
    os.rename(PACKAGECONF[:-2], PACKAGECONF)

LOGGER = PackageLogger('.prody')

SETTINGS = PackageSettings('prody', logger=LOGGER) 
SETTINGS.load()


_max_key_len = max([len(_key) for _key in CONFIGURATION.iterkeys()])
_max_def_len = 68 - _max_key_len
_table_lines = '=' * _max_key_len +  '  ' + '=' * _max_def_len

docstring = """

    {0:s}
    {1:s}  {2:s}         
    {0:s}""".format(_table_lines, 'Option'.ljust(_max_key_len), 
                    'Default setting (type)')

_ = {}

from textwrap import wrap
_keys = CONFIGURATION.keys()
_keys.sort()
for _key in _keys:
    _value = CONFIGURATION[_key]
    if isinstance(_value, str):
        _val = '"' + _value + '"'
    else: 
        _val = str(_value)
    _lines = wrap(_val + ' (' + type(_value).__name__  + ')', _max_def_len)
    docstring += """
    {0:s}  {1:s}""".format(_key.ljust(_max_key_len), _lines[0])
    for _line in _lines[1:]:
        docstring += """
        {0:s}  {1:52s}""".format(''.ljust(_max_key_len), _line)
    if SETTINGS.get(_key) is None:
        _[_key] = _value
docstring += """
    {0:s}""".format(_table_lines)

if _:
    SETTINGS.update(_)

def confProDy(*args, **kwargs):
    """Configure ProDy."""

    settings = None
    if args:
        if len(args) == 1:
            settings = SETTINGS.get(args[0])
        else:
            settings = [SETTINGS.get(option) for option in args]

    for option, value in kwargs.iteritems():    
        if isinstance(option, str):
            if option in CONFIGURATION:
                type_ = type(CONFIGURATION[option])
                if type(value) == type_:
                    SETTINGS[option] = value
                    SETTINGS.save()
                    LOGGER.info('ProDy is configured: {0:s}={1:s}'
                                .format(option, repr(value)))
                else:
                    raise TypeError('value must be a ' + type_.__name__)
            else:
                raise KeyError("'{0:s}' is not a valid option".format(option))
        else:
            raise TypeError('option must be a string')
    return settings

confProDy.__doc__ += docstring

confProDy.__doc__ += """

    Usage example::
            
      confProDy('backup')
      confProDy('backup', 'backup_ext')
      confProDy(backup=True, backup_ext='.bak')
      confProDy(backup_ext='.BAK')"""

def plog(*text):
    """Log *text* using ProDy logger at log level info.  Multiple arguments 
    are accepted.  Each argument will be converted to string and joined using 
    a white space as delimiter."""
    
    LOGGER.info(' '.join([str(s) for s in text]))

def startLogfile(filename, **kwargs):
    
    LOGGER.startLogfile(filename, **kwargs)
    
startLogfile.__doc__ = LOGGER.startLogfile.__doc__
    
def closeLogfile(filename):
    """Close logfile with *filename*."""

    LOGGER.closeLogfile(filename)

def setVerbosity(level):
    """    

    >>> from prody import *
    >>> setVerbosity('none')
    >>> plog('test')
    >>> setVerbosity('debug')
    >>> plog('test')"""

    LOGGER.setVerbosity(level)

setVerbosity.__doc__ = LOGGER.setVerbosity.__doc__.replace('\n    ', '\n') + \
                       setVerbosity.__doc__ 

def getVerbosity():
    """Return ProDy console verbosity level."""
    
    return LOGGER.getVerbosity()

getVerbosity.__doc__ = LOGGER.getVerbosity.__doc__

def checkUpdates():
    """Check PyPI to see if there is a newer ProDy version available.  Setting
    ProDy configuration parameter *check_updates* to a positive integer will 
    make ProDy automatically check updates, e.g.::
       
      confProDy(check_updates=7) # check at most once a week
      confProDy(check_updates=0) # do not auto check updates
      confProDy(check_updates=-1) # check at the start of every session"""
    
    import xmlrpclib
    pypi = xmlrpclib.Server('http://pypi.python.org/pypi')
    releases = pypi.package_releases('ProDy')
    if releases[0] == __version__:
        LOGGER.info('You are using the latest ProDy release (v{0:s}).'
                    .format(__version__))
    else:
        LOGGER.info('ProDy v{0:s} is available, you are using {1:s}.'
                    .format(releases[0], __version__))
    if SETTINGS['check_updates']:
        import time
        SETTINGS['last_check'] = time.time()
        SETTINGS.save()

if SETTINGS['check_updates']: 
    
    if SETTINGS.get('last_check') is None:
        SETTINGS['last_check'] = 0
    import time
    if ((time.time() - SETTINGS.get('last_check')) / 3600 / 24 > 
        SETTINGS['check_updates']):
        LOGGER.info('Checking PyPI for ProDy updates:')
        checkUpdates()

def test(**kwargs):
    """Run ProDy tests, ``prody.test()``. See :mod:`prody.tests` 
    documentation for more details."""
    
    if sys.version_info[:2] > (2,6):
        try:
            import prody.tests
        except ImportError:
            LOGGER.warning('Could not import ProDy unit tests, '
                           'please check your installation.')
        else:
            prody.tests.test(**kwargs)
    else:
        LOGGER.warning('ProDy tests are available for Python 2.7') 

__all__ = ['checkUpdates', 'confProDy', 'getVerbosity', 'setVerbosity',
           'startLogfile', 'closeLogfile', 
           'plog']

import kdtree
from kdtree import *
__all__.extend(kdtree.__all__)
__all__.append('kdtree')

import atomic 
from atomic import *
__all__.extend(atomic.__all__)
__all__.append('atomic')

from atomic import SELECT

import proteins
from proteins import *  
__all__.extend(proteins.__all__)
__all__.append('proteins')

import measure
from measure import *
__all__.extend(measure.__all__)
__all__.append('measure')

import dynamics
from dynamics import *
__all__.extend(dynamics.__all__)
__all__.append('dynamics')

import ensemble
from ensemble import *
__all__.extend(ensemble.__all__)
__all__.append('ensemble')

import trajectory
from trajectory import *
__all__.extend(trajectory.__all__)
__all__.append('trajectory')

import prody
__all__.append('prody')
