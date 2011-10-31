#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

from os import system
from glob import glob
from os import path

from prody import routines

for cmd in routines.PRODY_COMMANDS: 
    system('prody ' + cmd + ' -h > ' + path.join('scripts', 'prody_' + cmd + '.txt'))
    system('prody ' + cmd + ' --examples > ' + path.join('scripts', 'prody_' + cmd + '_eg.txt'))
