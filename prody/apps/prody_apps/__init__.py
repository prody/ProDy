"""This module defines structure and dynamics analysis applications."""

import importlib
import sys

try:
    import argparse
except ImportError:
    from .. import argparse

from prody.apps.apptools import *
from prody.utilities.misctools import impLoadModule

if sys.version_info[0] == 2:
    import imp
    path_prody = imp.find_module('prody')[1]
else:
    path_prody = importlib.util.find_spec("prody").submodule_search_locations[0]

try:
    import imp
    path_apps = imp.find_module('apps', [path_prody])[1]
    path_apps = imp.find_module('prody_apps', [path_apps])[1]
except ModuleNotFoundError:
    path_apps = importlib.util.find_spec("prody.apps").submodule_search_locations[0]
    path_apps += '/prody_apps/'

PRODY_APPS = ['anm', 'gnm', 'pca', 'eda', 'align', 'blast', 'biomol',
                  'catdcd', 'contacts', 'fetch', 'select', 'energy', 
                  'clustenm']

__all__ = ['prody_main']

prody_parser = argparse.ArgumentParser(
    description="ProDy: A Python Package for Protein Dynamics Analysis",
    epilog="See 'prody <command> -h' for more information on a specific "
           "command."
    )

prody_parser.add_argument('-c', '--cite',
    help="print citation info and exit",
    action=ProDyCitation, nargs=0)

prody_parser.add_argument('-v', '--version',
    help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

prody_commands = prody_parser.add_subparsers(
    title='subcommands')


for cmd in PRODY_APPS:
    cmd = 'prody_' + cmd
    mod = impLoadModule('prody.apps.prody_apps.', cmd, path_apps)
    mod.addCommand(prody_commands)

def prody_main():

    if len(sys.argv) == 1:
        prody_parser.print_help()
    else:
        namespace = prody_parser.parse_args()
        namespace.func(namespace)


if __name__ == '__main__':
    prody_main()
