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


class StarDict:
    def __init__(self, parsingDict, title='unnamed'):
        self._title = title
        self._dict = parsingDict
        self.dataBlocks = [StarDataBlock(self._dict, key)
                           for key in list(self._dict.keys())]
        self.numDataBlocks = len(self.dataBlocks)

    def __getitem__(self, key):
        try:
            return np.array(self.dataBlocks)[key]
        except:
            try:
                key = np.where(np.array(list(self._dict.keys())) == key)[0][0]
                return self.dataBlocks[key]
            except:
                raise ValueError('The key for getting items should '
                                 'be the name or number of a data block')

    def getTitle(self):
        return self._title

    def setTitle(self, value):
        self._title = value

    def getDict(self):
        return self._dict

    def __repr__(self):
        if self.numDataBlocks == 1:
            return '<StarDict: {0} (1 data block)>'.format(self._title)
        return '<StarDict: {0} ({1} data blocks)>'.format(self._title, self.numDataBlocks)

    def __iter__(self):
        """Yield StarDataBlock instances."""
        for key in list(self._dict.keys()):
            yield StarDataBlock(self, key)


class StarDataBlock:
    def __init__(self, starDict, key):
        self._title = key
        self._dict = starDict[key]
        self.loops = [StarLoop(self, index) for index in list(self._dict.keys())]
        self.numLoops = len(self.loops)

    def getLoop(self, index):
        try:
            return self.loops[index]
        except:
            raise ValueError('There is no loop with that index')

    def getTitle(self):
        return self._title

    def setTitle(self, title):
        self._title = title

    def __getitem__(self, key):
        try:
            return np.array(self.loops)[key]
        except:
            try:
                key = np.where(np.array(list(self._dict.keys())) == key)[0][0]
                return self.loops[key]
            except:
                raise ValueError(
                    'The key for getting items should be the name or number of a loop')

    def __repr__(self):
        if self.numLoops == 1:
            return '<StarDataBlock: {0} ({1} loop)>'.format(self._title, self.numLoops)
        return '<StarDataBlock: {0} ({1} loops)>'.format(self._title, self.numLoops)

    def __iter__(self):
        """Yield StarLoop instances."""
        for key in list(self._dict.keys()):
            yield StarLoop(self, key)


class StarLoop:
    def __init__(self, dataBlock, key):
        self._dict = dataBlock._dict[key]
        self.fields = list(self._dict['fields'].values())
        self.data = list(self._dict['data'].values())
        self.numFields = len(self.fields)
        self.numRows = len(self.data)
        self._title = dataBlock._title + ' loop ' + str(key)

    def getData(self, key):
        if key in self.fields:
            return [row[key] for row in self.data]
        else:
            raise ValueError('That field is not present in this loop')

    def getTitle(self):
        return self._title

    def setTitle(self, title):
        self._title = title

    def getDict(self):
        return self._dict

    def __getitem__(self, key):
        try:
            return np.array(self.data)[key]
        except:
            try:
                key = np.where(np.array(list(self._dict.keys())) == key)[0][0]
                return self.data[key]
            except:
                try:
                    return self.getData(key)
                except:
                    raise ValueError('The key for getting items should be fields, data, '
                                     'or a field name or number corresponding to a '
                                     'row or column of data')

    def __repr__(self):
        if self.numFields == 1 and self.numRows != 1:
            return '<StarLoop: {0} (1 column and {2} rows)>'.format(self._title, self.numRows)
        elif self.numFields != 1 and self.numRows == 1:
            return '<StarLoop: {0} ({1} columns and 1 row)>'.format(self._title, self.numFields)
        elif self.numFields == 1 and self.numRows == 1:
            return '<StarLoop: {0} (1 column and 1 row)>'.format(self._title)
        else:
            return '<StarLoop: {0} ({1} columns and {2} rows)>'.format(self._title, self.numFields, self.numRows)


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

    parsingDict = parseSTARStream(lines)

    return StarDict(parsingDict, filename)


def parseSTARStream(stream):
    finalDictionary = {}
    currentLoop = -1
    fieldCounter = 0
    dataItemsCounter = 0
    lineNumber = 0
    for line in stream:
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

        elif line.startswith('_') or line.startswith(' _'):
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
    '''
    Parses particle images using data from a STAR file containing information about them.


    arg particlesSTAR: a dictionary containing STAR file data about particles or
        a filename for a STAR file from which such data can be parsed.
    type particlesSTAR: str, StarDict, StarDataBlock, StarLoop

    arg indices: row indices of images to be parsed from data loop tables
        These can be integers or list-like entries with top level keys as well as the integers.
    type indices: list

    arg saveImageArrays: whether to save the numpy array for each image to file
        default is False
    type saveImageArrays: bool

    arg saveDirectory: directory where numpy image arrays are saved
        default is None, which means save to the current working directory
    type saveDirectory: str, None

    arg rotateImages: whether to apply in plane translations and rotations using 
        provided psi and origin data, default is True
    type rotateImages: bool 
    '''
    from skimage.transform import rotate

    saveImageArrays = kwargs.get('saveImageArrays', False)
    saveDirectory = kwargs.get('saveDirectory', None)
    rotateImages = kwargs.get('rotateImages', True)

    if not isinstance(particlesSTAR, (StarDict, StarDataBlock, StarLoop)):
        try:
            particlesSTAR = parseSTAR(particlesSTAR)
        except:
            raise ValueError('particlesSTAR should be a dictionary parsed from a STAR file '
                             'or a filename corresponding to one')

    loops = []
    if indices is None:
        if isinstance(particlesSTAR, StarDict):
            maxLoops = 0
            for dataBlock in particlesSTAR:
                if dataBlock.numLoops > maxLoops:
                    maxLoops = dataBlock.numLoops

            totalRows = 0
            for dataBlock in particlesSTAR:
                for loop in dataBlock:
                    totalRows += loop.numRows

            indices = np.array((particlesSTAR.numDataBlocks,maxLoops,))
            for i, entry in enumerate(reversed(indices)):          
                for j, loop in enumerate(particlesSTAR[entry[0]][:]):
                    foundImageField = False
                    if ('_image' in loop.fields) or ('_rlnImageName' in loop.fields):
                        indices[particlesSTAR.numDataBlocks-1-i].extend([j, list(loop.getDict()['data'].keys())])
                        foundImageField = True
                        loops.append(loop)
                if not foundImageField:
                    indices.pop(particlesSTAR.numDataBlocks-1-i)

        elif isinstance(particlesSTAR, StarDataBlock):
            indices = [[loop.getTitle().split(' ')[-1]] for loop in particlesSTAR.loops]
            for i, entry in enumerate(reversed(indices)):
                loop = particlesSTAR[entry[0]]
                foundImageField = False
                if ('_image' in loop.fields) or ('_rlnImageName' in loop.fields):
                    indices[particlesSTAR.numLoops-1-i].extend([list(loop.getDict()['data'].keys())])
                    foundImageField = True
                    loops.append(loop)
            if not foundImageField:
                indices.pop(particlesSTAR.numLoops-1-i)

        elif isinstance(particlesSTAR, StarLoop):
            indices = list(loop.getDict()['data'].keys())
            loops.append(particlesSTAR)

    elif len(indices[0]) > 1:
        if isinstance(particlesSTAR, StarDict):
            kw_indices = indices
            indices = [[dataBlock.getTitle()] for dataBlock in particlesSTAR.dataBlocks]
            for i, entry in enumerate(reversed(indices)):          
                for j, loop in enumerate(particlesSTAR[entry[0]][:]):
                    foundImageField = False
                    if ('_image' in loop.fields) or ('_rlnImageName' in loop.fields):
                        indices[particlesSTAR.numDataBlocks-1-i].extend([j, kw_indices])
                        foundImageField = True
                        loops.append(loop)
                if not foundImageField:
                    indices.pop(particlesSTAR.numDataBlocks-1-i)

        

    if indices == []:
        raise ValueError('particlesSTAR does not contain any data loop tables with image fields')

    image_stacks = {}
    images = []
    for i in indices:
        if isinstance(particlesSTAR, StarDict):
            particle = particlesSTAR[i[0]][i[1]][i[2]]
        elif isinstance(particlesSTAR, StarDataBlock):
            particle = particlesSTAR[i[0]][i[1]]
        elif isinstance(particlesSTAR, StarLoop):
            particle = particlesSTAR[i]

        try:
            image_index = int(particle['_rlnImageName'].split('@')[0])-1
            filename = particle['_rlnImageName'].split('@')[1]
        except:
            try:
                image_index = int(particle['_image'].split('@')[0])-1
                filename = particle['_image'].split('@')[1]
            except:
                raise ValueError('particlesSTAR does not contain data about particle image '
                                 'location in either RELION or XMIPP format')

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
