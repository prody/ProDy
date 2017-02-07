"""This module defines structure and dynamics analysis applications."""

import imp
import sys

try:
    import argparse
except ImportError:
    from .. import argparse

from ..apptools import *

path_prody = imp.find_module('prody')[1]
path_apps = imp.find_module('apps', [path_prody])[1]
path_apps = imp.find_module('prody_apps', [path_apps])[1]

PRODY_APPS = ['anm', 'gnm', 'pca', 'eda', 'align', 'blast', 'biomol',
                  'catdcd', 'contacts', 'fetch', 'select', 'saxs',]

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
    mod = imp.load_module('prody.apps.prody_apps.' + cmd,
                          *imp.find_module(cmd, [path_apps]))
    mod.addCommand(prody_commands)


def prody_main():

    if len(sys.argv) == 1:
        prody_parser.print_help()
    else:
        namespace = prody_parser.parse_args()
        namespace.func(namespace)


if __name__ == '__main__':
    prody_main()
