"""Perform EDA calculations and output the results in plain text, NMD, and
graphical formats.

Download example :download:`MDM2 trajectory files </doctest/mdm2.tar.gz>`."""

import os.path

from ..apptools import *
from .nmaoptions import *

def addCommand(commands):

    subparser = commands.add_parser('eda',
        parents=[commands.choices.get('pca')],
        help='perform essential dynamics analysis calculations',
        add_help=False)

