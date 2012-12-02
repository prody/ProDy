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
try:
    import argparse
except ImportError:
    from . import argparse
import textwrap

from prody.utilities import wrapText

__all__ = ['Quiet', 'UsageExample', 'ProDyCitation', 'ProDyVersion', 
           'DevelApp']
           
FIGARGS = {
    ('-F', '--figure-format'): { 
        'dest': 'figformat', 
        'type': str, 
        'default': 'pdf', 
        'metavar': 'STR',
        'help': 'figure file format', 
        'choices': set('eps pdf png ps raw rgba svg svgz'.split())},
    ('-D', '--dpi'): { 
        'dest': 'figdpi', 
        'type': int, 
        'default': 300, 
        'metavar': 'INT', 
        'help': 'figure resolution (dpi)'},
    ('-W', '--width'): { 
        'dest': 'figwidth', 
        'type': float, 
        'default': 8, 
        'metavar': 'FLOAT', 
        'help': 'figure width (inch)'},        
    ('-H', '--height'): {
        'dest': 'figheight', 
        'type': float, 
        'default': 6, 
        'metavar': 'FLOAT', 
        'help': 'figure height (inch)'},
}

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
        self._figures = []
        self._figargs = []
    
    def _getKwargs(self, arg):
        """Return keyword arguments."""
        
        kwargs = copy(self._args[arg])
        default = kwargs.get('default')
        help = '' 
        if 'choices' in kwargs: 
            help = ', one of ' + ', '.join([str(ch) 
                                       for ch in kwargs['choices']])


        if default is not None and 'action' not in kwargs:
            help += ' (default: %(default)s)'.format()
        kwargs['help'] += help
        return kwargs

    def _docArg(self, doc, arg):
        """Add documentation string for *arg* to *doc*."""
        
        kwargs = self._args[arg]
        dest = kwargs.get('dest', arg[0])
        desc = ':arg {0}: {1}'.format(dest, kwargs['help'])
        choices = kwargs.get('choices')
        if choices:
            desc += ', one of ' + ', '.join(['``' + repr(ch) + '``' 
                                             for ch in choices])
        default = kwargs.get('default')
        if default:
            desc += ', default is ``{0}``'.format(repr(default))
        doc.extend(wrapText(desc, join=False, subsequent_indent='    '))
        try:
            type = kwargs['type']
        except KeyError: 
            try:
                action = kwargs['action']
            except KeyError:
                type = None
            else:
                if action.startswith('store') and action.endswith('e'):
                    type = bool
                else:
                    type = None
        if type is not None:
            doc.append(':type {0}: {1}'.format(dest, type.__name__))
        doc.append('')
    
    def addGroup(self, name, description, mutexc=False, required=False):
        """If *mutexc* is **True** add a mutually exclusive group. """

        if name in self._group_args:
            raise ValueError('group {0} is already defined'.format(name))
                 
        self._groups.append((name, mutexc, required))
        self._group_desc[name] = description
        self._group_args[name] = []

    def addArgument(self, *args, **kwargs):
        
        if args in self._args:
            raise ValueError('argument {0} is already defined'
                             .format(str(args)))
        if len(args) == 1 and args[0][0] != '-':
            group = 'positional'
        else:
            group = kwargs.pop('group', 'ungrouped')
            #if 'choices' in kwargs and 'default' not in kwargs:
            #    raise ValueError('argument has multiple choices, '
            #                     'but no default value')
        self._group_args[group].append(args)
        self._args[args] = kwargs 

    def addFigure(self, *args, **kwargs):
        """Add figure argument.  Figure options will be displayed at the end.
        When multiple figures are added, ``-A, --all-figures`` option will
        be added automatically."""
        
        if args in self._args:
            raise ValueError('argument {0} is already defined'
                             .format(str(args)))
        if args[0][0] != '-':
            raise ValueError('figure argument cannot be a positional argument')
        
        if not self._figures:
            self._args.update(FIGARGS)

        self._figures.append(args)
        self._args[args] = kwargs

    def addFigarg(self, *args, **kwargs):
        
        if args in self._args:
            raise ValueError('argument {0} is already defined'
                             .format(str(args)))
        if args[0][0] != '-':
            raise ValueError('figure argument cannot be a positional argument')

        self._figargs.append(args)
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
        #sub.add_argument('--debug', help="print error traceback to screen",
        #    action='store_true', dest='debug')

        if self._example:
            sub.add_argument('--examples', action=UsageExample, nargs=0,
                help='show usage examples and exit')
            sub.set_defaults(usage_example=self._example, 
                             test_examples=self._egtests)
        
        for arg in self._group_args['positional']:
            sub.add_argument(*arg, **self._args[arg])
        for arg in self._group_args['ungrouped']:
            sub.add_argument(*arg, **self._getKwargs(arg))
        for name, mutexc, required in self._groups:
            if mutexc:
                group = sub.add_argument_group(self._group_desc[name])
                group = group.add_mutually_exclusive_group(required=required)
            else:
                group = sub.add_argument_group(self._group_desc[name])
            for arg in self._group_args[name]:
                group.add_argument(*arg, **self._getKwargs(arg))
        
        if self._figures:
            if len(self._figures) > 1:
                group = sub.add_argument_group('figures')
                group.add_argument('-A', '--all-figures', dest='figall', 
                    action='store_true', default=False)
            else:
                group = sub.add_argument_group('figure options')

            for arg in self._figures:
                group.add_argument(*arg, **self._getKwargs(arg))

            if len(self._figures) > 1:
                group = sub.add_argument_group('figure options')


            for arg in self._figargs:            
                group.add_argument(*arg, **self._getKwargs(arg))
                
            for arg in FIGARGS.keys():            
                group.add_argument(*arg, **self._getKwargs(arg))
        
        
        positional = [(arg[0], self._args[arg].get('nargs', None)) 
                      for arg in self._group_args['positional']]
        
        def callback(ns, func=self._function, pos=positional):
            args = []
            for arg, nargs in pos:
                if nargs is None:
                    args.append(ns.__dict__.pop(arg))
                else:
                    args.extend(ns.__dict__.pop(arg))
            try:
                func(*args, **ns.__dict__)
            except Exception as err:
                import sys, traceback, tempfile
                traceback.print_exc(file=sys.stderr)
                sys.stderr.write('An exception occurred when executing '
                 'your command.  If this is not a user error, please '
                 'report it to ProDy developers at: '
                 'https://bitbucket.org/abakan/prody/issues\n')
                sys.stderr.write('Error: ' + str(err) + '\n')
            
        sub.set_defaults(func=callback)
        sub.set_defaults(subparser=sub)
    
    def addDocstring(self, function, help=True):
        """Add documentation string to *function* function.  If *help* is 
        **True** help text will also be added."""
    
        doc = []
        if help:
            help = self._help
            for arg in self._group_args['positional']:
                arg = arg[0]
                help = help.replace(arg, '*' + arg + '*')
            doc.append(help[0].upper() + help[1:] + '.')
            doc.append('')
    
        for arg in self._group_args['positional']:
            self._docArg(doc, arg)
        for arg in self._group_args['ungrouped']:
            self._docArg(doc, arg)
        for name, _, _ in self._groups:
            doc.append('*' + self._group_desc[name].title() + '*')
            doc.append('')
            for arg in self._group_args[name]:
                self._docArg(doc, arg)

        function.__doc__ = '\n'.join(doc)
    
    def setFunction(self, function, docstring=True):
        """Set *function*, and add *docstring* based on argument settings.
        Function should be set after all arguments are added."""
        
        self._function = function
        if docstring:
            self.addDocstring(function)
