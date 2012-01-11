#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

from os import system
from glob import glob
from os import path

from prody import routines

system('prody > prody.txt')
for cmd in routines.commands.choices.keys(): 
    system('prody ' + cmd + ' -h > ' + 'prody_' + cmd + '.txt')
    system('prody ' + cmd + ' --examples > ' + 'prody_' + cmd + '_eg.txt')
