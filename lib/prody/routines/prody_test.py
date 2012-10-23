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

"""Run ProDy unit tests from command line."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from .actions import *

__all__ = ['prody_test']

def prody_test(*mods, **kwargs):
    """Run ProDy unittests."""

    import prody
    prody.test(*mods,
               label=kwargs.get('label', 'fast'), 
               verbose=kwargs.get('verbose', 1))


def addCommand(commands):

    subparser = commands.add_parser('test', 
    help='run ProDy unit tests')

    subparser.add_argument('--quiet', help="suppress info messages to stderr",
        action=Quiet, nargs=0)
        
    subparser.add_argument('--examples', action=UsageExample, nargs=0,
        help='show usage examples and exit')

    subparser.set_defaults(usage_example=
    """ProDy unittests can be run as follows:
    
    $ prody test
      
To run all tests, use -l/--label argument as follows:
        
    $ prody test -l full 

To increase verbosity, use -v/--verbose argument as follows:
        
    $ prody test -l fast -v 2
    
To run tests for a specific module, pass module name:
    
    $ prody test proteins
    $ prody test atomic.select""", test_examples=[])


    group = subparser.add_argument_group('output options')
    
    group.add_argument('-l', '--label', dest='label', metavar='STR', type=str,
        help='label used when nose is available, fast (default) or full',
        default='fast')

    group.add_argument('-v', '--verbose', dest='verbose', metavar='INT', 
        type=int, help='verbosity level, default is 1',
        default=1)
        

    group.add_argument('mod', nargs='*',
        help='ProDy module names')
        
    subparser.set_defaults(func=lambda ns: prody_test(*ns.mod, **ns.__dict__))
    subparser.set_defaults(subparser=subparser)
