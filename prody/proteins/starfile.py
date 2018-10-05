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

__all__ = ['parseSTAR', 'writeSTAR', 'parseImagesFromSTAR', 
           'StarDict', 'StarDataBlock', 'StarLoop',]


class StarDict:
    def __init__(self, parsingDict, prog, title='unnamed'):
        self._title = title
        self._dict = parsingDict
        self._prog = prog
        self.dataBlocks = [StarDataBlock(self, key)
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
                                 'be the names or numbers of data blocks')

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

    def pop(self, index):
        """Pop dataBlock with the given index from the list of dataBlocks in StarDict"""
        self.dataBlocks.pop(index)
        self.numDataBlocks -= 1

class StarDataBlock:
    def __init__(self, starDict, key):
        self._title = key
        self._dict = starDict._dict[key]
        self._prog = starDict._prog

        if list(self._dict.keys()) == ['fields','data']:
            self.loops = []
            self.numLoops = 0
            self.data = list(self._dict['data'].values())
            self.fields = list(self._dict['fields'].values())
            self.numEntries = len(self.data)
            self.numFields = len(self.fields)
        else:
            self.loops = [StarLoop(self, index)
                          for index in list(self._dict.keys())]
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

    def getDict(self):
        return self._dict

    def __getitem__(self, key):
        if self.loops == []:
            try:
                return np.array(self._dict['data'][key])

            except:
                try:
                    return np.array(self.data)[key]
                except:
                    raise ValueError('The key for getting items should be the data entry number')

        else:
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
        if self.numLoops == 0:
            if self.numEntries == 1:
                return '<StarDataBlock: {0} ({1} entry)>'.format(self._title, self.numEntries)
            else:
                return '<StarDataBlock: {0} ({1} entries)>'.format(self._title, self.numEntries)
        elif self.numLoops == 1:
            return '<StarDataBlock: {0} ({1} loop containing ' \
                    '{2} columns and {3} rows)>'.format(self._title, self.numLoops, \
                                                        self.loops[0].numFields, self.loops[0].numRows)
        return '<StarDataBlock: {0} ({1} loops)>'.format(self._title, self.numLoops)

    def __iter__(self):
        """Yield StarLoop instances."""
        for key in list(self._dict.keys()):
            yield StarLoop(self, key)

    def pop(self, index):
        """Pop loop with the given index from the list of loops in dataBlock"""
        self.loops.pop(index)
        self.numLoops -= 1


class StarLoop:
    def __init__(self, dataBlock, key):
        self._dict = dataBlock._dict[key]
        self._prog = dataBlock._prog
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

    parsingDict, prog = parseSTARStream(lines)

    return StarDict(parsingDict, prog, filename)


def parseSTARStream(stream):
    prog = 'RELION'
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
            inLoop = False

        elif len(line.split()) == fieldCounter:
            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter] = {}
            fieldCounter = 0
            for fieldEntry in line.strip().split():
                currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][fieldCounter + 1]
                finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                fieldCounter += 1
            dataItemsCounter += 1

        elif line.startswith('#'):
            if line.startswith('# XMIPP'):
                prog = 'XMIPP'

        else:
            raise TypeError('This file does not conform to the STAR file format.'
                            'There is a problem with line {0}:\n {1}'.format(lineNumber, line))

        lineNumber += 1

    return finalDictionary, prog


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


def parseImagesFromSTAR(particlesSTAR, **kwargs):
    '''
    Parses particle images using data from a STAR file 
    containing information about them.

    arg particlesSTAR: a filename for a STAR file.
    type particlesSTAR: str

    arg block_indices: indices for data blocks containing rows 
        corresponding to images of interest
        The indexing scheme is similar to that for numpy arrays.
        Default behavior is use all data blocks about images
    type block_indices: list, :class:`~numpy.ndarray`

    arg row_indices: indices for rows corresponding to images of interest
        The indexing scheme is similar to that for numpy arrays. 
        row_indices should be a 1D or 2D array-like.
        2D row_indices should contain an entry for each relevant loop. 
        If a 1D array-like is given the same row indices 
        will be applied to all loops.
        Default behavior is to use all rows about images
    type row_indices: list, :class:`~numpy.ndarray`

    arg particle_indices: indices for particles regardless of STAR structure
        default is take all particles
        Please note: this acts after block_indices and row_indices
    type particle_indices: list, :class"`~numpy.ndarray`

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

    block_indices = kwargs.get('block_indices', None)
    # No loop_indices because data blocks about particle images contain 1 loop 
    row_indices = kwargs.get('row_indices', None)
    particle_indices = kwargs.get('particle_indices',None)

    saveImageArrays = kwargs.get('saveImageArrays', False)
    saveDirectory = kwargs.get('saveDirectory', None)
    rotateImages = kwargs.get('rotateImages', True)

    try:
        particlesSTAR = parseSTAR(particlesSTAR)
    except:
        raise ValueError('particlesSTAR should be a filename for a STAR file')

    # Check dimensions/contents of particlesSTAR and generate full indices
    dataBlocks = []
    loops = []
    maxLoops = 0
    maxRows = 0
    dataBlock_goodness = []
    for dataBlock in particlesSTAR:

        foundImageField = False
        for loop in dataBlock:
            if ('_image' in loop.fields) or ('_rlnImageName' in loop.fields):
                foundImageField = True
                loops.append(loop)
                if loop.numRows > maxRows:
                    maxRows = loop.numRows
            else:
                dataBlock.pop(int(loop.getTitle().split(' ')[-1]))

        if dataBlock.numLoops > maxLoops:
            maxLoops = dataBlock.numLoops

        if foundImageField:
            dataBlocks.append(dataBlock)
            dataBlock_goodness.append(True)
        else:
            dataBlock_goodness.append(False)

    indices = np.zeros((len(dataBlocks),maxLoops,maxRows,3),dtype=int)
    i = -1
    for n, dataBlock in enumerate(particlesSTAR):
        if dataBlock_goodness[n]:
            i += 1
            for j, loop in enumerate(dataBlock):
                for k in range(maxRows):
                    if k < loop.numRows:
                        indices[i,j,k] = np.array([n,j,k])
                    else:
                        indices[i,j,k] = np.array([0,0,0])

    dataBlocks = np.array(dataBlocks)
    loops = np.array(loops)

    # Convert keyword indices to valid indices if possible
    if block_indices is not None:
        if np.array_equal(dataBlocks, np.array([])):
            raise TypeError('particlesSTAR must have data blocks to use block_indices')

        try:
            block_indices = np.array(block_indices)
        except:
            raise TypeError('block_indices should be array-like')

        if block_indices.ndim != 1:
            raise ValueError('block_indices should be a 1-dimensional array-like')

        for i, index in enumerate(list(reversed(block_indices))):
            try:
                block = particlesSTAR[index]
                if not isinstance(block, StarDataBlock):
                    LOGGER.warn('There is no block corresponding to block_index {0}. '
                                'This index has been removed.'.format(block_indices.shape[0]-i-1))
                    block_indices = np.delete(block_indices, i, 0)
            except:
                LOGGER.warn('There is no block corresponding to block_index {0}. '
                            'This index has been removed.'.format(block_indices.shape[0]-i-1))
                block_indices = np.delete(block_indices, i, 0)

        if not np.array_equal(block_indices, np.array([])):
            indices = np.concatenate(([indices[np.where(indices[:,0,0,0] == item)] 
                                       for item in block_indices]),axis=0)
        else:
            LOGGER.warn('None of the block_indices corresponded to dataBlocks. '
                        'Default block indices corresponding to all dataBlocks '
                        'will be used instead.')

    dataBlocks = particlesSTAR[block_indices]

    if row_indices is not None:
        try:
            row_indices = np.array(row_indices)
        except:
            raise TypeError('row_indices should be array-like')

        if row_indices.ndim == 1:
            if isinstance(row_indices[0], int):
                # row_indices provided was truly 1D so 
                # we will use same row indices for all data blocks 
                # and warn the user we are doing so
                if len(dataBlocks) != 1:
                    LOGGER.warn('row_indices is 1D but there are multiple data blocks '
                                'so the same row indices will be used for each')
                
                row_indices = np.array([row_indices for i in range(len(dataBlocks))])
                # This also works if len(dataBlocks) == 1

            elif isinstance(row_indices[0], (list, tuple)):
                # A list-like of list-likes of different sizes was provided
                # We turn it into a proper 2D array by filling the short 
                # list likes with zeros

                if len(row_indices) != len(dataBlocks):
                    raise ValueError('There should be an entry in row indices for '
                                     'each data block')

                max_len = 0
                for entry in row_indices:
                    if not np.isscalar(entry):
                        if len(entry) > max_len:
                            max_len = len(entry)
                            
                row_indices_list_entries = []
                for entry in row_indices:
                    if isinstance(entry, int):
                        list_entry = [entry]
                    else:
                        list_entry = list(entry)

                    while len(list_entry) < max_len:
                        list_entry.append(0)

                    row_indices_list_entries.append(list_entry)

                row_indices = np.array(row_indices_list_entries)

        elif row_indices.ndim == 2:
            # A list-like of list-likes of the same size was provided
            if row_indices.shape[0] != len(dataBlocks):
                if len(row_indices) == 1:
                    # we will use same row indices for all data blocks 
                    # and warn the user we are doing so
                    if len(dataBlocks) != 1:
                        LOGGER.warn('row_indices has one entry but there are multiple data blocks '
                                    'so the same row indices will be used for each')
                    
                    row_indices = np.array([row_indices[0] for i in range(len(dataBlocks))])
                    # This also works if len(dataBlocks) == 1
                else:
                    raise ValueError('There should be an entry in row indices for '
                                    'each data block')
        
        else:
            raise ValueError('row_indices should be 1D or 2D array-like objects')

        # indices need updating
        good_indices_list = []
        for i, index_i in enumerate(indices):
            good_indices_list.append([])
            for j, index_j in enumerate(index_i):
                good_indices_list[i].append([])
                for r, index_r in enumerate(row_indices[i]):
                    for k, index_k in enumerate(index_j):
                        if k == index_r:
                            if not (r != 0 and index_r == 0):
                                good_indices_list[i][j].append(index_k)
                            else:
                                good_indices_list[i][j].append(np.array([0,0,0]))


        indices = np.array(good_indices_list)

    if indices is np.array([]):
        raise ValueError('selection does not contain any rows with image fields')

    # Use indices to collect particle data dictionaries
    particles = []

    for i, index_i in enumerate(indices):
        for j, index_j in enumerate(index_i):
            for k, index_k in enumerate(index_j):
                if not (np.array_equal(index_k, np.array([0,0,0])) 
                and not (i == 0 and j == 0 and k == 0)):
                    particles.append(particlesSTAR[index_k[0]][index_k[1]][index_k[2]])

    if particle_indices is None:
        particle_indices = list(range(len(particles)))

    # Parse images using particle dictionaries
    image_stacks = {}
    images = []
    parsed_images_data = []
    stk_images = []
    if particlesSTAR._prog == 'RELION':
        imageFieldKey = '_rlnImageName'
    else:
        imageFieldKey = '_image'
        
    for i in particle_indices:
        particle = particles[i]

        try:
            image_field = particle[imageFieldKey]
            image_index = int(image_field.split('@')[0])-1
            filename = image_field.split('@')[1]
        except:
            raise ValueError('particlesSTAR does not contain data about particle image '
                             '{0} location in either RELION or XMIPP format'.format(i))

        if filename.endswith('.stk'):
            stk_images.append(str(i))
            continue

        if not filename in list(image_stacks.keys()):
            image_stacks[filename] = parseEMD(filename).density

        image = image_stacks[filename][image_index]
        parsed_images_data.append(image_field)

        if saveImageArrays:
            if saveDirectory is not None:
                np.save('{0}/{1}'.format(saveDirectory, i), image)
            else:
                np.save('{1}'.format(i), image)

        if rotateImages:
            if particlesSTAR._prog == 'RELION':
                anglePsi = float(particle['_rlnAnglePsi'])
                originX = float(particle['_rlnOriginX'])
                originY = float(particle['_rlnOriginY'])
            elif particlesSTAR._prog == 'XMIPP':
                anglePsi = float(particle['_anglePsi'])
                originX = float(particle['_shiftX'])
                originY = float(particle['_shiftY'])
            images.append(rotate(image, anglePsi,
                                 center=(float(image.shape[0])-originX,
                                         float(image.shape[1])-originY)))
        else:
            images.append(image)

    if len(stk_images) > 0:
        LOGGER.warn('ProDy currently cannot parse images from XMIPP .stk files. '
                    'Please be aware that images {0} and {1} will be missing '
                    'from the final array.'.format(', '.join(stk_images[:-1]),stk_images[-1]))

    return np.array(images), parsed_images_data
