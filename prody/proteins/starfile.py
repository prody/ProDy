# -*- coding: utf-8 -*-
"""
This module defines functions for parsing `STAR files`_.

.. _STAR files: https://www2.mrc-lmb.cam.ac.uk/relion/index.php/Conventions_%26_File_formats#The_STAR_format
"""


from collections import defaultdict
import os.path

import numpy as np

from prody.utilities import openFile
from prody import LOGGER, SETTINGS

__all__ = ['parseSTAR','writeSTAR']

def parseSTAR(filename):
    """Returns a dictionary containing data
    parsed from a Relion STAR file.

    :arg filename: a filename
        The .star extension can be omitted.
    """

    if not os.path.isfile(filename) and not os.path.isfile(filename + '.star'):
        raise IOError('There is no file with that name.')

    starfile = open(filename, 'r')
    lines = starfile.readlines()
    starfile.close()

    finalDictionary = {}
    currentLoop = -1
    fieldCounter = 0
    dataItemsCounter = 0

    for line in lines:
        if line.startswith('data_'):
            currentDataBlock = line[5:].strip()
            finalDictionary[currentDataBlock] = {}
            currentLoop = -1

        elif line.startswith('loop_'):
            currentLoop += 1
            finalDictionary[currentDataBlock][currentLoop] = {}
            finalDictionary[currentDataBlock][currentLoop]['fields'] = {}
            finalDictionary[currentDataBlock][currentLoop]['data'] = {}
            fieldCounter = 0

        elif line.startswith('_'):
            currentField = line.strip() 
            finalDictionary[currentDataBlock][currentLoop]['fields'][fieldCounter] = currentField
            fieldCounter += 0
            dataItemsCounter = 0

        elif line.strip() == '':
            pass

        elif len(line.split()) == fieldCounter:
            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter] = {}
            fieldCounter = 0
            for fieldEntry in line.strip().split():
                currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][fieldCounter]
                finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                fieldCounter += 1
            dataItemsCounter += 1

        else:
            raise TypeError('This file does not conform to the STAR file format.')

    return finalDictionary

def writeSTAR(filename, starDict):
    """Writes a STAR file from a dictionary containing data
    such as that parsed from a Relion STAR file.

    :arg filename: a filename
        The .star extension can be omitted.

    :arg dictionary: a dictionary in STAR format
        This should have nested entries starting with data blocks then loops/tables then
        field names and finally data.
    """

    star = open(filename,'w')

    for dataBlockKey in starDict:
        star.write('\ndata_' + dataBlockKey + '\n')
        for loopNumber in starDict[dataBlockKey]:
            star.write('\nloop_\n')
            for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                star.write('_' + starDict[dataBlockKey][loopNumber]['fields'][fieldNumber] + '\n')
            for dataItemNumber in starDict[dataBlockKey][loopNumber]['data']:
                for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                    currentField = starDict[dataBlockKey][loopNumber]['fields'][fieldNumber]
                    star.write(starDict[dataBlockKey][loopNumber]['data'][dataItemNumber][currentField] + ' ')
                star.write('\n')

    star.close()
    return

