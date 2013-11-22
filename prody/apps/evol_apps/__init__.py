"""This module defines some sequence evolution applications."""

import imp
import sys

try:
    import argparse
except ImportError:
    from .. import argparse

from ..apptools import *

path_prody = imp.find_module('prody')[1]
path_apps = imp.find_module('apps', [path_prody])[1]
path_apps = imp.find_module('evol_apps', [path_apps])[1]

EVOL_APPS = ['search', 'fetch', 'filter', 'refine', 'merge', 'occupancy',
             'conserv', 'coevol', 'rankorder']

__all__ = ['evol_main']

evol_parser = argparse.ArgumentParser(
    description="Evol: Sequence Evolution and Dynamics Analysis",
    epilog="See 'evol <command> -h' for more information on a specific "
           "command."
    )

evol_parser.add_argument('-c', '--cite',
    help="print citation info and exit",
    action=ProDyCitation, nargs=0)

evol_parser.add_argument('-v', '--version',
    help="print ProDy version and exit",
    action=ProDyVersion, nargs=0)

evol_parser.add_argument('-e', '--examples', action=UsageExample, nargs=0,
    help='show usage examples and exit')


evol_parser.set_defaults(usage_example=
"""Sequence coevolution analysis involves several steps that including
retrieving data and refining it for calculations.  These steps are
illustrated below for RnaseA protein family.

Search Pfam database:

  $  evol search 2w5i

Download Pfam MSA file:

  $  evol fetch RnaseA

Refine MSA file:

  $ evol refine RnaseA_full.slx -l RNAS1_BOVIN --seqid 0.98 --rowocc 0.8

Checking occupancy:

  $ evol occupancy RnaseA_full.slx -l RNAS1_BOVIN -o col -S

Conservation analysis:

  $ evol conserv RnaseA_full_refined.slx

Coevolution analysis:

  $ evol coevol RnaseA_full_refined.slx -S -c apc

Rank order analysis:


  $ evol rankorder RnaseA_full_refined_mutinfo_corr_apc.txt -p \
2w5i_1-121.pdb --seq-sep 3 """, test_examples=[(0,1,2)])

evol_commands = evol_parser.add_subparsers(
    title='subcommands')


for cmd in EVOL_APPS:
    cmd = 'evol_' + cmd
    mod = imp.load_module('prody.apps.evol_apps.' + cmd,
                          *imp.find_module(cmd, [path_apps]))
    mod.APP.addApplication(evol_commands)


def evol_main():

    if len(sys.argv) == 1:
        evol_parser.print_help()
    else:
        namespace = evol_parser.parse_args()
        namespace.func(namespace)


if __name__ == '__main__':
    evol_main()
