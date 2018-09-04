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

from .emdfile import parseEMD

__all__ = ['parseSTAR', 'writeSTAR', 'parseImagesFromSTAR']


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
    lineNumber = 0
    for line in lines:
        if line.startswith('data_'):
            currentDataBlock = line[5:].strip()
            finalDictionary[currentDataBlock] = {}
            currentLoop = -1
            inLoop = False
            startingBlock = True
            fieldCounter = 0

        elif line.startswith('loop_'):
            currentLoop += 1
            inLoop = True
            finalDictionary[currentDataBlock][currentLoop] = {}
            finalDictionary[currentDataBlock][currentLoop]['fields'] = {}
            finalDictionary[currentDataBlock][currentLoop]['data'] = {}
            fieldCounter = 0

        elif line.startswith('_'):
            currentField = line.strip().split()[0]

            if inLoop:
                finalDictionary[currentDataBlock][currentLoop]['fields'][fieldCounter + 1] = currentField
                dataItemsCounter = 0
            else:
                if startingBlock:
                    finalDictionary[currentDataBlock]['fields'] = {}
                    finalDictionary[currentDataBlock]['data'] = {}
                    startingBlock = False
                    dataItemsCounter = 0
                finalDictionary[currentDataBlock]['fields'][fieldCounter + 1] = currentField
                finalDictionary[currentDataBlock]['data'][dataItemsCounter] = {}
                finalDictionary[currentDataBlock]['data'][dataItemsCounter][currentField] = line.strip().split()[1]
                dataItemsCounter += 1

            fieldCounter += 1

        elif line.strip() == '':
            pass

        elif len(line.split()) == fieldCounter:
            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter] = {}
            fieldCounter = 0
            for fieldEntry in line.strip().split():
                currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][fieldCounter + 1]
                finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                fieldCounter += 1
            dataItemsCounter += 1

        elif line.startswith('#'):
            pass

        else:
            raise TypeError('This file does not conform to the STAR file format.'
                            'There is a problem with line {0}:\n {1}'.format(lineNumber, line))

        lineNumber += 1

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

    star = open(filename, 'w')

    for dataBlockKey in starDict:
        star.write('\ndata_' + dataBlockKey + '\n')
        for loopNumber in starDict[dataBlockKey]:
            star.write('\nloop_\n')
            for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                star.write('_' + starDict[dataBlockKey]
                           [loopNumber]['fields'][fieldNumber] + '\n')
            for dataItemNumber in starDict[dataBlockKey][loopNumber]['data']:
                for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                    currentField = starDict[dataBlockKey][loopNumber]['fields'][fieldNumber]
                    star.write(starDict[dataBlockKey][loopNumber]
                               ['data'][dataItemNumber][currentField] + ' ')
                star.write('\n')

    star.close()
    return


def parseImagesFromSTAR(particlesSTAR, indices, **kwargs):
    from skimage.transform import rotate

    saveImageArrays = kwargs.get('saveImages', False)
    saveDirectory = kwargs.get('saveDirectory', None)
    rotateImages = kwargs.get('rotateImages', True)

    image_stacks = {}
    images = []

    if not isinstance(particlesSTAR, dict):
        try:
            particlesSTAR = parseSTAR(particlesSTAR)
        except:
            raise ValueError('particlesSTAR should be a dictionary parsed from a STAR file ' 
                             'or a filename corresponding to one')

    if indices is None:
        indices = list(particlesSTAR.keys())

    for i in indices:
        particle = particlesSTAR[i]
        image_index = int(particle['_rlnImageName'].split('@')[0])-1
        filename = particle['_rlnImageName'].split('@')[1]
        if not filename in list(image_stacks.keys()):
            image_stacks[filename] = parseEMD(filename).density

        image = image_stacks[filename][image_index]

        if saveImageArrays:
            if saveDirectory is not None:
                np.save('{0}/{1}'.format(saveDirectory, i), image)
            else:
                np.save('{1}'.format(i), image)

        if rotateImages:
            images.append(rotate(image, float(particle['_rlnAnglePsi']),
                                 center=(180-float(particle['_rlnOriginX']),
                                         180-float(particle['_rlnOriginY']))))
        else:
            images.append(image)

    return images
