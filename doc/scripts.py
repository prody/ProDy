#!/usr/bin/python
__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010  Ahmet Bakan'

from os import system
from glob import glob
from os import path

for script in glob('../scripts/*py'): 
    name = path.splitext(path.split(script)[1])[0]
    system(script + ' -h > ' + path.join('scripts', name + '.txt'))
    system(script + ' --examples > ' + path.join('scripts', name + '_eg.txt'))
