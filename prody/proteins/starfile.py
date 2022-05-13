# -*- coding: utf-8 -*-
"""
This module defines functions for parsing `STAR`_ (Self-defining Text Archiving and Retrieval) files.
This includes metadata files from cryo-EM image analysis programs including RELION and XMIPP as well as
the Crystallographic Information File (CIF) format used for much of the PDB.

.. _STAR: https://www2.mrc-lmb.cam.ac.uk/relion/index.php/Conventions_%26_File_formats#The_STAR_format
"""

from collections import OrderedDict
import os.path
from numbers import Integral
import numpy as np
import sys

from prody.utilities import openFile, split, pystr, isListLike
from prody import LOGGER

from .emdfile import parseEMD

__all__ = ['parseSTAR', 'writeSTAR', 'parseImagesFromSTAR',
           'StarDict', 'StarDataBlock', 'StarLoop', 
           'parseSTARSection']


class StarDict:
    def __init__(self, parsingDict, prog, title='unnamed', indices=None):
        self._title = title
        self._dict = parsingDict
        self._prog = prog
        self._indices = indices

        if indices is None:
            self.dataBlocks = [StarDataBlock(self, key)
                               for key in self._dict.keys()]
        else:
            self.dataBlocks = []
            for idx in indices:
                if isListLike(idx):
                    self.dataBlocks.append(StarDataBlock(self, idx[0], idx[1]))
                else:
                    self.dataBlocks.append(StarDataBlock(self, idx))

            self._dict = OrderedDict()
            for i, idx in enumerate(indices):
                self._dict[idx[0]] = self.dataBlocks[i]._dict

        self._n_data_blocks = len(self.dataBlocks)

    def numDataBlocks(self):
        return self._n_data_blocks
    
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
        if self._n_data_blocks == 1:
            return '<StarDict: {0} (1 data block)>'.format(self._title)
        return '<StarDict: {0} ({1} data blocks)>'.format(self._title, self.numDataBlocks())

    def __iter__(self):
        """Yield StarDataBlock instances."""
        for key in list(self._dict.keys()):
            yield StarDataBlock(self, key)

    def pop(self, index):
        """Pop dataBlock with the given index from the list of dataBlocks in StarDict"""
        self.dataBlocks.pop(index)
        self._n_data_blocks -= 1

    def search(self, substr, return_indices=False):
        indices = []
        for data_block in self:
            if data_block.search(substr).numEntries() != 0 or data_block.search(substr).numLoops() != 0:
                indices.append((data_block._title,
                                data_block.search(substr, return_indices=True)[0]))

        if return_indices:
            return indices, StarDict(self._dict, self._prog, self._title, indices=indices)

        return StarDict(self._dict, self._prog, self._title, indices=indices)

    def printData(self):
        for data_block in self:
            sys.stdout.write('data_' + data_block._title + '\n')
            data_block.printData()
            sys.stdout.write('\n')


class StarDataBlock:
    def __init__(self, starDict, key, indices=None):
        self._title = key
        self._prog = starDict._prog
        self._starDict = starDict

        if indices is None:
            try:
                self._dict = starDict._dict[key]
            except:
                self._dict = list(starDict._dict)[key]

            keys = list(self._dict.keys())
        else:
            keys = [idx[0] for idx in indices]
            self._dict = OrderedDict()
            self._dict['data'] = OrderedDict()
            self._dict['fields'] = OrderedDict()
            for idx in indices:
                if idx[0] == 'data':
                    self._dict[idx[0]][idx[1]] = starDict._dict[self._title][idx[0]][idx[1]]
                    if not 'fields' in keys:
                        for k, v in self._starDict._dict[self._title]['fields'].items():
                            if v == idx[1]:
                                self._dict['fields'][k] = v
                else:
                    self._dict[idx[0]] = OrderedDict()
                    self._dict[idx[0]]['fields'] = starDict._dict[self._title][idx[0]]['fields']
                    self._dict[idx[0]]['data'] = OrderedDict()
                    for id1 in idx[1]:
                        self._dict[idx[0]]['data'][id1] = starDict._dict[self._title][idx[0]]['data'][id1]

        if set(keys) == set(['data', 'fields']):
            self.loops = []
            self._n_loops = 0

            self.data = np.array(list(self._dict['data'].values()))
            self.fields = np.array(list(self._dict['fields'].values()))

            if not isListLike(self.data):
                self.data = [self.data]

            if not isListLike(self.fields):
                self.fields = [self.fields]

            self._n_entries = len(self.data)
            self._n_fields = len(self.fields)

        elif 'data' in keys and 'fields' in keys:
            if indices is not None:
                self.loops = [StarLoop(self, key, idx)
                              for (key, idx) in indices
                              if key not in ['data', 'fields']]
            else:
                self.loops = [StarLoop(self, key) for key in keys
                              if key not in ['data', 'fields']]

            self.data = np.array(list(self._dict['data'].values()))
            self.fields = np.array(list(self._dict['fields'].values()))

            if not isListLike(self.data):
                self.data = [self.data]

            if not isListLike(self.fields):
                self.fields = [self.fields]

            self._n_entries = len(self.data)
            self._n_fields = len(self.fields)
            self._n_loops = len(self.loops)

        elif 'data' in keys:
            if indices is not None:
                self.loops = [StarLoop(self, key, idx)
                              for (key, idx) in indices
                              if key != 'data']
            else:
                self.loops = [StarLoop(self, key)
                              for key in keys
                              if key != 'data']

            self.data = np.array(list(self._dict['data'].values()))
            self.fields = np.array(list(self._dict['fields'].values()))

            if not isListLike(self.data):
                self.data = [self.data]

            if not isListLike(self.fields):
                self.fields = [self.fields]

            self._n_loops = len(self.loops)
            self._n_entries = len(self.data)
            self._n_fields = 0

        elif 'fields' in keys:
            if indices is not None:
                self.loops = [StarLoop(self, key, idx)
                              for (key, idx) in indices
                              if key != 'fields']
            else:
                self.loops = [StarLoop(self, key)
                              for key in keys
                              if key != 'fields']

            self.data = np.array(list(self._dict['data'].values()))
            self.fields = np.array(list(self._dict['fields'].values()))

            if not isListLike(self.data):
                self.data = [self.data]

            if not isListLike(self.fields):
                self.fields = [self.fields]

            self._n_loops = len(self.loops)
            self._n_entries = len(self.data)
            self._n_fields = 0

        else:
            if indices is not None:
                self.loops = [StarLoop(self, key, idx)
                              for (key, idx) in indices]
            else:
                self.loops = [StarLoop(self, key) for key in keys]

            self._n_loops = len(self.loops)
            self._n_entries = 0
            self._n_fields = 0

    def numLoops(self):
        return self._n_loops
    
    def numEntries(self):
        return self._n_entries
    
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
        if key == 'data':
            try:
                return self._dict[key]
            except:
                raise ValueError('This StarDataBlock has no non-loop data')
        if key == 'fields':
            try:
                return self._dict[key]
            except:
                raise ValueError(
                    'This StarDataBlock has no non-loop data fields')
        else:
            try:
                return StarLoop(self, key)
            except:
                try:
                    return self.loops[key]
                except:
                    raise ValueError('The key for getting items should be data, fields, '
                                     'or the name or number of a loop')

    def __repr__(self):
        if self.numLoops() == 0:
            if self.numEntries() == 0:
                return '<StarDataBlock: {0} (no entries)>'.format(self._title)
            if self.numEntries() == 1:
                return '<StarDataBlock: {0} (1 entry)>'.format(self._title)
            else:
                return '<StarDataBlock: {0} ({1} entries)>'.format(self._title, self.numEntries())
        elif self.numEntries() == 0:
            if self.numLoops() == 1:
                return '<StarDataBlock: {0} (1 loop containing ' \
                    '{1} columns and {2} rows)>'.format(self._title,
                                                        self.loops[0].numFields(), self.loops[0].numRows())
            return '<StarDataBlock: {0} ({1} loops)>'.format(self._title, self.numLoops())
        else:
            if self.numLoops() == 1:
                if self.numEntries() == 1:
                    return '<StarDataBlock: {0} (1 entry and 1 loop)>'.format(self._title)
                else:
                    return '<StarDataBlock: {0} ({1} entries and 1 loop)>'.format(self._title, self.numEntries())
            else:
                if self.numEntries() == 1:
                    return '<StarDataBlock: {0} (1 entry and {1} loops)>'.format(self._title, self.numLoops())
                else:
                    return '<StarDataBlock: {0} ({1} entries and {2} loops)>'.format(self._title, self.numEntries(),
                                                                                     self.numLoops())

    def __iter__(self):
        """Yield StarLoop instances."""
        for key in list(self._dict.keys()):
            yield StarLoop(self, key)

    def pop(self, index):
        """Pop loop with the given index from the list of loops in dataBlock"""
        self.loops.pop(index)
        self._n_loops -= 1

    def search(self, substr, return_indices=False):
        indices = []
        for key, value in self._dict.items():
            if key == 'fields':
                pass
            elif key == 'data':
                for field, val in value.items():
                    if field.find(substr) != -1 or val.find(substr) != -1:
                        indices.append(('data', field))
            else:
                idx, loop = self[key].search(substr, return_indices=True)
                if loop.numRows() != 0:
                    indices.append((key, idx))

        if return_indices:
            return indices, StarDataBlock(self._starDict, self._title, indices)

        return StarDataBlock(self._starDict, self._title, indices)

    def printData(self):
        for key, value in self._dict.items():
            if key == 'fields':
                longest_len = 0
                for field in self[key].values():
                    if len(field) > longest_len:
                        longest_len = len(field)
            elif key == 'data':
                longest_len = 0
                for field, val in value.items():
                    sys.stdout.write(field + ' '*(longest_len - len(field)))
                    sys.stdout.write('\t' + val + '\n')
                sys.stdout.write('\n')
            else:
                sys.stdout.write('_loop\n')
                self[key].printData()
                sys.stdout.write('\n')


class StarLoop:
    def __init__(self, dataBlock, key, indices=None):
        self._key = key
        self._dataBlock = dataBlock

        if indices is None:
            self._dict = dataBlock._dict[self._key]
        else:
            self._dict = OrderedDict()
            self._dict['fields'] = dataBlock._dict[self._key]['fields']
            self._dict['data'] = OrderedDict()

            indices_are_keys = True
            for index in indices:
                try:
                    self._dict['data'][index] = dataBlock._dict[self._key]['data'][index]
                except:
                    indices_are_keys = False
                    break

            if not indices_are_keys:
                self._dict['data'] = OrderedDict()
                for index in indices:
                    self._dict['data'][index] = list(self._dataBlock._dict[self._key]['data'].values())[index]

        self._prog = dataBlock._prog
        self.fields = list(self._dict['fields'].values())
        self.data = list(self._dict['data'].values())
        self._n_fields = len(self.fields)
        self._n_rows = len(self.data)
        self._title = str(dataBlock._title) + ' loop ' + str(key)

    def numRows(self):
        return self._n_rows
    
    def numFields(self):
        return self._n_fields
    
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
                try:
                    return self._dict[key]
                except KeyError:
                    key = list(self._dict.keys()).index(key)
                    return self.data[key]
            except:
                try:
                    return self.getData(key)
                except:
                    raise ValueError('The key for getting items should be fields, data, '
                                     'or a field name or number corresponding to a '
                                     'row or column of data')

    def __repr__(self):
        if self.numFields() == 1 and self.numRows() != 1:
            return '<StarLoop: {0} (1 column and {1} rows)>'.format(self._title, self.numRows())
        elif self.numFields() != 1 and self.numRows() == 1:
            return '<StarLoop: {0} ({1} columns and 1 row)>'.format(self._title, self.numFields())
        elif self.numFields() == 1 and self.numRows() == 1:
            return '<StarLoop: {0} (1 column and 1 row)>'.format(self._title)
        else:
            return '<StarLoop: {0} ({1} columns and {2} rows)>'.format(self._title, self.numFields(), self.numRows())

    def search(self, substr, return_indices=False):
        indices = []
        for j, row in self._dict['data'].items():
            found_it = False
            for entry in row.items():
                field, value = entry
                if field.find(substr) != -1 or value.find(substr) != -1:
                    found_it = True
                    break
            if found_it:
                indices.append(j)

        if return_indices:
            return indices, StarLoop(self._dataBlock, self._key, indices)

        return StarLoop(self._dataBlock, self._key, indices)

    def printData(self):
        for field in self.fields:
            sys.stdout.write(field + '\n')
        for row in self.data:
            for entry in row.values():
                sys.stdout.write(entry + '\t')
            sys.stdout.write('\n')


def parseSTAR(filename, **kwargs):
    """Returns a dictionary containing data parsed from a STAR file.

    :arg filename: a filename
        The .star extension can be omitted.
    :type filename: str

    :arg start: line number for starting
        Default is **None**, meaning start at the beginning
    :type start: int, None

    :arg stop: line number for stopping
        Default is **None**, meaning don't stop.
    :type stop: int, None

    :arg shlex: whether to use shlex for splitting lines so as to preserve quoted substrings
        Default is **False**
    :type shlex: bool
    """
    if not os.path.isfile(filename) and not os.path.isfile(filename + '.star'):
        raise IOError('There is no file called {0}.'.format(filename))

    start = kwargs.get('start', None)
    if start is not None and not isinstance(start, Integral):
        raise TypeError('start should be an integer or None')

    stop = kwargs.get('stop', None)
    if stop is not None and not isinstance(stop, Integral):
        raise TypeError('stop should be an integer or None')

    shlex = kwargs.get('shlex', False)
    if not isinstance(shlex, bool):
        raise TypeError('shlex should be a boolean')

    starfile = openFile(filename, 'r')
    lines = [pystr(line) for line in starfile.readlines()]
    starfile.close()

    parsingDict, prog = parseSTARLines(lines, **kwargs)

    return StarDict(parsingDict, prog, filename)


def parseSTARLines(lines, **kwargs):
    start = kwargs.get('start', None)
    if start is None:
        start = 0

    stop = kwargs.get('stop', None)
    if stop is None:
        stop = len(lines)

    prog = kwargs.get('prog', None)
    shlex = kwargs.get('shlex', False)

    finalDictionary = OrderedDict()
    currentLoop = -1
    block_fieldCounter = 0
    loop_fieldCounter = 0
    active_fieldCounter = 0
    dataItemsCounter = 0
    lineNumber = 0
    inLoop = False
    inShortBlock = False
    for line in lines[start:stop]:
        if line.startswith('data_'):
            currentDataBlock = line[5:].strip()
            finalDictionary[currentDataBlock] = OrderedDict()
            currentLoop = -1
            inLoop = False
            inShortBlock = False
            startingBlock = True
            block_fieldCounter = 0

        elif line.startswith('loop_'):
            currentLoop += 1
            inLoop = True
            inShortBlock = False
            finalDictionary[currentDataBlock][currentLoop] = OrderedDict()
            finalDictionary[currentDataBlock][currentLoop]['fields'] = OrderedDict()
            finalDictionary[currentDataBlock][currentLoop]['data'] = OrderedDict()
            loop_fieldCounter = 0

        elif line.startswith('_') or line.startswith(' _'):
            # This marks a field identifier
            currentField = split(line.strip(), shlex=shlex)[0]

            if inLoop:
                # We expect to only have the field identifier and no data until after

                if len(split(line.strip(), shlex=shlex)) == 1:
                    # This is what we expect for a data loop
                    finalDictionary[currentDataBlock][currentLoop]['fields'][loop_fieldCounter] = currentField
                    dataItemsCounter = 0
                    loop_fieldCounter += 1

                else:
                    # This is contrary to that so we leave the loop
                    inLoop = False

                    # We populate fields and data together, continuing the regular data block
                    finalDictionary[currentDataBlock]['fields'][block_fieldCounter] = currentField
                    finalDictionary[currentDataBlock]['data'][currentField] = split(line.strip(),
                                                                                    shlex=shlex)[1]
                    block_fieldCounter += 1

            else:
                # Outside a loop, populate fields and data together in a regular data block
                if startingBlock:
                    # Initialise the data block first
                    finalDictionary[currentDataBlock]['fields'] = OrderedDict()
                    finalDictionary[currentDataBlock]['data'] = OrderedDict()
                    startingBlock = False
                    block_fieldCounter = 0

                finalDictionary[currentDataBlock]['fields'][block_fieldCounter] = currentField

                if len(split(line.strip(), shlex=shlex)) > 1:
                    # This is the usual behaviour so we can fill in the data from the rest of the line
                    finalDictionary[currentDataBlock]['data'][currentField] = split(line.strip(),
                                                                                    shlex=shlex)[1]
                else:
                    # In this case, we will look for the data in a short block over the following line(s).
                    # If a single field takes multiple lines, these lines start and end with a semi-colon.
                    # We'll handle that in the data section.
                    finalDictionary[currentDataBlock]['data'][currentField] = ''
                    inShortBlock = True
                    startingShortBlock = True

                block_fieldCounter += 1

        elif line.strip() == '#':
            inLoop = False
            inShortBlock = False

        elif line.strip() == '':
            pass

        elif inLoop:
            # Here we handle the data part of the loop.
            # Data outside a loop is handled in line with the fields above or in shortDataBlocks below.
            if not inShortBlock and len(split(line, shlex=shlex)) == loop_fieldCounter:
                # This is the usual case where each entry in the line corresponds to a field
                finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter] = OrderedDict()
                active_fieldCounter = 0
                for fieldEntry in split(line.strip(), shlex=shlex):
                    currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][active_fieldCounter]
                    finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                    active_fieldCounter += 1
                dataItemsCounter += 1
            else:
                # The data is now being broken across lines.
                if not inShortBlock:
                    inShortBlock = True
                    finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter] = OrderedDict()
                    active_fieldCounter = 0
                    if not line.startswith(';'):
                        # Then we haven't got a split field and can treat fields as normal
                        inSplitField = False
                        for fieldEntry in split(line.strip(), shlex=shlex):
                            currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][active_fieldCounter]
                            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                            active_fieldCounter += 1
                    else:
                        # We have a single field split over many lines
                        inSplitField = True
                        currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][active_fieldCounter]
                        finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = line.strip() + ' '
                else:
                    if not inSplitField:
                        # check if we are entering one
                        if line.startswith(';'):
                            inSplitField = True
                            currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][active_fieldCounter]
                            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = line.strip() + ' '
                        else:
                            # continue as normal
                            for fieldEntry in split(line.strip(), shlex=shlex):
                                currentField = finalDictionary[currentDataBlock][currentLoop]['fields'][active_fieldCounter]
                                finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] = fieldEntry
                                active_fieldCounter += 1
                    else:
                        finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] += line.strip()
                        if line.strip() == ';':
                            # This marks the end of the split field
                            inSplitField = False
                            active_fieldCounter += 1
                        else:
                            # Prepare for the next line
                            finalDictionary[currentDataBlock][currentLoop]['data'][dataItemsCounter][currentField] += ' '

                    if active_fieldCounter == loop_fieldCounter:
                        inShortBlock = False
                        dataItemsCounter += 1

        elif inShortBlock:
            # We can now append the data in the lines here.
            finalDictionary[currentDataBlock]['data'][currentField] += line.strip()
            if startingShortBlock:
                startingShortBlock = False
                if not line.startswith(';'):
                    # We only expect one line if there's no semi-colon
                    inShortBlock = False
                else:
                    # Prepare for the next line
                    finalDictionary[currentDataBlock]['data'][currentField] += ' '
            else:
                if line.strip() == ';':
                    # This marks the end of the field so we've filled it
                    inShortBlock = False
                else:
                    # Prepare for the next line
                    finalDictionary[currentDataBlock]['data'][currentField] += ' '

        elif line.startswith('#'):
            if line.startswith('# XMIPP'):
                prog = 'XMIPP'

        else:
            raise TypeError('This file does not conform to the STAR file format. '
                            'There is a problem with line {0}:\n {1}'.format(lineNumber, line))

        lineNumber += 1

    return finalDictionary, prog


def writeSTAR(filename, starDict, **kwargs):
    """Writes a STAR file from a dictionary containing data
    such as that parsed from a Relion STAR file.

    :arg filename: a filename
        The .star extension can be omitted.
    :type filename: str

    :arg starDict: a dictionary in STAR format
        This should have nested entries starting with data blocks then loops/tables then
        field names and finally data.
    :type starDict: dict

    kwargs can be given including the program style to follow (*prog*)
    """
    prog=kwargs.get('prog', 'XMIPP')

    star = open(filename, 'w')

    for dataBlockKey in starDict:
        star.write('\ndata_' + dataBlockKey + '\n')
        for loopNumber in starDict[dataBlockKey]:
            star.write('\nloop_\n')
            for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                if prog == 'XMIPP':
                    star.write(' ')
                star.write(starDict[dataBlockKey][loopNumber]['fields'][fieldNumber] + '\n')
            for dataItemNumber in starDict[dataBlockKey][loopNumber]['data']:
                if prog == 'XMIPP':
                    star.write('\t')
                for fieldNumber in starDict[dataBlockKey][loopNumber]['fields']:
                    currentField = starDict[dataBlockKey][loopNumber]['fields'][fieldNumber]
                    star.write(starDict[dataBlockKey][loopNumber]['data'][dataItemNumber][currentField] + '\t')
                star.write('\n')

    star.close()
    return


def parseImagesFromSTAR(particlesSTAR, **kwargs):
    """
    Parses particle images using data from a STAR file 
    containing information about them.

    :arg particlesSTAR: a filename for a STAR file.
    :type particlesSTAR: str

    :arg block_indices: indices for data blocks containing rows 
        corresponding to images of interest
        The indexing scheme is similar to that for numpy arrays.
        Default behavior is use all data blocks about images
    :type block_indices: list, :class:`~numpy.ndarray`

    :arg row_indices: indices for rows corresponding to images of interest
        The indexing scheme is similar to that for numpy arrays. 
        row_indices should be a 1D or 2D array-like.
        2D row_indices should contain an entry for each relevant loop. 
        If a 1D array-like is given the same row indices 
        will be applied to all loops.
        Default behavior is to use all rows about images
    :type row_indices: list, :class:`~numpy.ndarray`

    :arg particle_indices: indices for particles regardless of STAR structure
        default is take all particles
        Please note: this acts after block_indices and row_indices
    :type particle_indices: list, :class:`~numpy.ndarray`

    :arg saveImageArrays: whether to save the numpy array for each image to file
        default is False
    :type saveImageArrays: bool

    :arg saveDirectory: directory where numpy image arrays are saved
        default is **None**, which means save to the current working directory
    :type saveDirectory: str

    :arg rotateImages: whether to apply in plane translations and rotations using 
        provided psi and origin data, default is True
    :type rotateImages: bool 
    
    """

    try:
        from skimage.transform import rotate
    except ImportError:
        raise ImportError('This function requires scikit-image.')

    block_indices = kwargs.get('block_indices', None)
    # No loop_indices because data blocks about particle images contain 1 loop
    row_indices = kwargs.get('row_indices', None)
    particle_indices = kwargs.get('particle_indices', None)

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
                if loop.numRows() > maxRows:
                    maxRows = loop.numRows()
            else:
                dataBlock.pop(int(loop.getTitle().split(' ')[-1]))

        if dataBlock.numLoops() > maxLoops:
            maxLoops = dataBlock.numLoops()

        if foundImageField:
            dataBlocks.append(dataBlock)
            dataBlock_goodness.append(True)
        else:
            dataBlock_goodness.append(False)

    indices = np.zeros((len(dataBlocks), maxLoops, maxRows, 3), dtype=int)
    i = -1
    for n, dataBlock in enumerate(particlesSTAR):
        if dataBlock_goodness[n]:
            i += 1
            for j, loop in enumerate(dataBlock):
                for k in range(maxRows):
                    if k < loop.n_rows:
                        indices[i, j, k] = np.array([n, j, k])
                    else:
                        indices[i, j, k] = np.array([0, 0, 0])

    dataBlocks = np.array(dataBlocks)
    loops = np.array(loops)

    # Convert keyword indices to valid indices if possible
    if block_indices is not None:
        if np.array_equal(dataBlocks, np.array([])):
            raise TypeError(
                'particlesSTAR must have data blocks to use block_indices')

        try:
            block_indices = np.array(block_indices)
        except:
            raise TypeError('block_indices should be array-like')

        if block_indices.ndim != 1:
            raise ValueError(
                'block_indices should be a 1-dimensional array-like')

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
            indices = np.concatenate(([indices[np.where(indices[:, 0, 0, 0] == item)]
                                       for item in block_indices]), axis=0)
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

                row_indices = np.array(
                    [row_indices for i in range(len(dataBlocks))])
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

                    row_indices = np.array([row_indices[0]
                                            for i in range(len(dataBlocks))])
                    # This also works if len(dataBlocks) == 1
                else:
                    raise ValueError('There should be an entry in row indices for '
                                     'each data block')

        else:
            raise ValueError(
                'row_indices should be 1D or 2D array-like objects')

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
                                good_indices_list[i][j].append(
                                    np.array([0, 0, 0]))

        indices = np.array(good_indices_list)

    if indices is np.array([]):
        raise ValueError(
            'selection does not contain any rows with image fields')

    # Use indices to collect particle data dictionaries
    particles = []

    for i, index_i in enumerate(indices):
        for j, index_j in enumerate(index_i):
            for k, index_k in enumerate(index_j):
                if not (np.array_equal(index_k, np.array([0, 0, 0]))
                        and not (i == 0 and j == 0 and k == 0)):
                    particles.append(
                        particlesSTAR[index_k[0]][index_k[1]][index_k[2]])

    if particle_indices is None:
        particle_indices = list(range(len(particles)))

    # Parse images using particle dictionaries
    image_stacks = OrderedDict()
    images = []
    parsed_images_data = []
    stk_images = []
    if particlesSTAR._prog == 'XMIPP':
        imageFieldKey = '_image'
    else:
        imageFieldKey = '_rlnImageName'

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
                    'from the final array.'.format(', '.join(stk_images[:-1]), stk_images[-1]))

    return np.array(images), parsed_images_data


def parseSTARSection(lines, key):
    """Parse a section of data from *lines* from a STAR file 
    corresponding to a *key* (part before the dot). 
    This can be a loop or data block.
    
    Returns data encapulated in a list and the associated fields."""

    if not isinstance(key, str):
        raise TypeError("key should be a string")

    if not key.startswith("_"):
        key = "_" + key

    i = 0
    fields = OrderedDict()
    fieldCounter = -1
    foundBlock = False
    foundBlockData = False
    doneBlock = False
    start = 0
    stop = 0

    while not doneBlock and i < len(lines):
        line = lines[i]
        if line.split(".")[0] == key:
            fieldCounter += 1
            fields[line.split(".")[1].strip()] = fieldCounter
            if not foundBlock:
                foundBlock = True

        if foundBlock:
            if not line.startswith("#"):
                if not foundBlockData:
                    start = i
                    foundBlockData = True
            else:
                if foundBlockData:
                    doneBlock = True
                    stop = i

        i += 1

    if i < len(lines):
        star_dict, _ = parseSTARLines(lines[:2] + lines[start-1: stop], shlex=True)
        loop_dict = list(star_dict.values())[0]

        if lines[start - 1].strip() == "loop_":
            data = list(loop_dict[0]["data"].values())
        else:
            data = [loop_dict["data"]]
    else:
        LOGGER.warn("Could not find {0} in lines.".format(key))
        return []

    return data