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

"""This module define tools for application development."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from copy import copy
import argparse
import textwrap

__all__ = ['Quiet', 'UsageExample', 'ProDyCitation', 'ProDyVersion', 
           'DevelApp']

class Quiet(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import prody
        prody.LOGGER.verbosity = 'warning'

class UsageExample(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        tw = textwrap.TextWrapper()
        for line in namespace.usage_example.splitlines():
            if line.lstrip().startswith('$'):
                print(line)
            else:
                print('\n'.join(tw.wrap(line)))
        parser.exit()

class ProDyCitation(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print("Bakan A, Meireles LM, Bahar I "
              "ProDy: Protein Dynamics Inferred from Theory and Experiments "
              "Bioinformatics 2011 27(11):1575-1577.")
        parser.exit()

class ProDyVersion(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import prody
        print("ProDy " + prody.__version__)
        parser.exit()

class DevelApp(object):
    
    def __init__(self, name, help):
        
        self._name = name
        self._help = help
        self._groups = []
        self._group_desc = {}
        self._group_args = {'positional': [], 'ungrouped': [], }
        self._args = {}
        self._example = None
        self._egtests = None
        self._function = None
    
    def _getKwargs(self, arg):
        """Return keyword arguments."""
        
        kwargs = copy(self._args[arg])
        default = kwargs.get('default')
        if default is not None and 'action' not in kwargs:
            kwargs['help'] += ' (default: {0:s})'.format(repr(default))
        return kwargs
    
    def addGroup(self, name, description):

        if name in self._group_args:
            raise ValueError('group {0:s} is already defined'.format(name))
                 
        self._groups.append(name)
        self._group_desc[name] = description
        self._group_args[name] = []

    def addArgument(self, *args, **kwargs):
        
        if args in self._args:
            raise ValueError('argument {0:s} is already defined'
                             .format(str(args)))
        if len(args) == 1 and args[0][0] != '-':
            group = 'positional'
        else:
            group = kwargs.pop('group', 'ungrouped')
        self._group_args[group].append(args)
        self._args[args] = kwargs 

    def setExample(self, example, tests=None):
        """Set usage *example* string and list of examples for *tests*."""
        
        self._example = example
        self._egtests = tests or []

    def addApplication(self, parser):
        
        if self._function is None:
            raise ValueError('an application function is not set')
        
        sub = parser.add_parser(self._name, help=self._help)
        
        sub.add_argument('--quiet', help="suppress info messages to stderr",
            action=Quiet, nargs=0)

        if self._example:
            sub.add_argument('--examples', action=UsageExample, nargs=0,
                help='show usage examples and exit')
            sub.set_defaults(usage_example=self._example, 
                             test_examples=self._egtests)
        
        for arg in self._group_args['positional']:
            sub.add_argument(*arg, **self._args[arg])
        for arg in self._group_args['ungrouped']:
            sub.add_argument(*arg, **self._getKwargs(arg))
        for name in self._groups:
            group = sub.add_argument_group(self._group_desc[name])
            for arg in self._group_args[name]:
                group.add_argument(*arg, **self._getKwargs(arg))
        
        positional = [(arg, self._args[arg].get('nargs', None)) 
                      for arg in self._group_args['positional']]
        
        def callback(ns, func=self._function, pos=positional):
            args = [ns.__dict__.pop(arg) for arg in pos]
            func(*args, **ns.__dict__)
            
        sub.set_defaults(func=callback)
        sub.set_defaults(subparser=sub)
    
    def addDocstring(self, app):
        """Add documentation string to *app* function."""
    
        pass
    
    def setFunction(self, function):
        """Set *function*."""
        
        self._function = function
