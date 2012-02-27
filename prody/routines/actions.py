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

"""This module defines argument parser actions."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import argparse
import textwrap

__all__ = ['Quiet', 'UsageExample', 'ProDyCitation', 'ProDyVersion']

class Quiet(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        import prody
        prody.setVerbosity('warning')

class UsageExample(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        tw = textwrap.TextWrapper()
        for line in namespace.usage_example.splitlines():
            print("\n".join(tw.wrap(line)))
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
