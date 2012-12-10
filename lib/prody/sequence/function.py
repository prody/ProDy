# -*- coding: utf-8 -*-
# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2012 Ahmet Bakan
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

"""This module defines MSA analysis functions."""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import arange

from .analysis import *
from prody import LOGGER
from prody.utilities import openFile
import numpy as np

__all__ = ['parseHeatmap', 'writeHeatmap', 'showHeatmap']

HM_HEADERLENGTH = 8
hmheaders = {}
for title, val in [
    ('-min', None),
    ('-xlabel', ''),
    ('-xorigin', '1'),   
    ('-xstep', '1'),
    ('-ylabel', ''),
    ('-max', None),
    ('-title', 'Heatmap'),
    ('-numbering', '')]:
    hmheaders[title] = val

def showHeatmap(heatmap, clim=None, *args, **kwargs):
    """This shows a heatmap plot. Input can be an 2D matrix
    or an outfile file from VMD's "Heatmapper" plugin. 
    
    HeatMap is plotted using :func:`~matplotlib.pyplot.imshow`
    function."""

    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except AttributeError:
        headers, heatmap = parseHeatmap(heatmap)
        ndim, shape = heatmap.ndim, heatmap.shape
    
    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')
    
    if not kwargs.has_key('interpolation'):
        kwargs['interpolation'] = 'nearest'
    if not kwargs.has_key('origin'):
        kwargs['origin'] = 'lower'
    
    y, x = shape
    indices = kwargs.pop('indices', None)
    if indices is not None:
        start = indices[0] - 0.5
        end = start + y
        extent = [-0.5, x - 0.5, start, end]
    else:
        extent = [-0.5, x - 0.5, -0.5, y - 0.5]
    
    xlabel = kwargs.pop('xlabel', None)
    ylabel = kwargs.pop('ylabel', None)
    title = kwargs.pop('title', 'HeatMap')
    format = kwargs.pop('format', True)
    
    import matplotlib.pyplot as plt
    show = plt.imshow(heatmap, extent=extent, *args, **kwargs), plt.colorbar()
    
    if clim is not None:
        try:
            cmin, cmax = clim
        except:
            try:
                cmin, cmax = clim[0], clim[1]
            except:
                LOGGER.warn('clim should be a tuple or list containing min and '
                        'max values. Could not set limits')
        else:
            if cmin < cmax:
                show[0].set_clim(cmin, cmax)
            else:
                LOGGER.warn('first element of clim should be smaller than the'
                            ' second. Could not set limits')
            
    if format:
        if xlabel is not None:
            plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)
        plt.title(title)
    return show


def parseHeatmap(hmfile, **kwargs):
    """This parses a HeatMap outfile file from VMD's "Heatmapper" plugin.
    Input is ``".hm"`` file. Returns a tuple containing a dictionary of
    headers and a numpy array containing the vales plotted in Heatmapper."""    
    
    from os.path import isfile

    if isfile(str(hmfile)):
        with openFile(hmfile) as file:
            head = [file.next() for x in xrange(HM_HEADERLENGTH)]
    else:
        raise ValueError('Input file not a valid file')
    
    headers = {}
    for i in xrange(8):
        temp = head[i].split()
        if len(temp) < 2:
            raise TypeError('Cannot parse header of hm, check file format')
        headers[temp[0]] = ' '.join(temp[1:])
               
    if set.difference(set(headers.keys()), set(hmheaders.keys())):
        raise TypeError('Cannot parse header of hm, check file format')
    
    try:
        hmArray = np.loadtxt(hmfile, delimiter=';',
                             skiprows=HM_HEADERLENGTH, dtype='|S')
    except:
        raise TypeError('Cannot parse heatmap file {0:s}, check file format'
                        .format(hmfile))
    if np.all(hmArray[:,-1] == ' '):
        hmArray = np.delete(hmArray,hmArray.shape[1]-1, -1)
    y, x = hmArray.shape
    
    for i in xrange(y):
        temp1, temp2, firstelem = hmArray[i,0].split(':')
        hmArray[i,0] = firstelem
    hmArray = hmArray.astype(np.float)
        
    return headers, hmArray
            
        
def writeHeatmap(heatmap, prefix=None, **kwargs):
    """This writes the input 2D array into a ``"".hm""`` file in the same
    format as VMD's Heapmapper plugin. Input arguments may contain the files`
    header information like ``"-max"``, ``"-min""`` , ``"xorigin"``,
    ``""xstep""``, ``""xlabel""``,  ``""ylabel""``, ``"title"`` or
    ``"numbering"``. y axis can also be labeled based on input indices."""
    
    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except:
        raise TypeError('heatmap should be an array object')
    if ndim!= 2:
        raise TypeError('heatmap should be a 2D array')
    
    if not prefix:
        prefix = 'HeatMapFile'
    outname = prefix + '.hm'
    file = openFile(outname, 'wb')
    
    for key in hmheaders:
        item = kwargs.get(key[1:], hmheaders[key])
        if key == '-min':
            if item is None:
                item = '{0:.2f}'.format(heatmap.min())
        if key == '-max':
            if item is None:
                item = '{0:.2f}'.format(heatmap.max())
        line = key + ' "' + str(item) + '"\n'
        file.write(line)
    
    col1 = kwargs.get('indices', None)
    if col1 == None or len(col1) != shape[0]:
        col1 = np.arange(1, shape[0]+1)
    col2 = np.arange(0, shape[0])
    
    firstelem = map(lambda x,y,z: x+':'+y+':'+z, map(str, col1),
                    map(str, col2), map(str, heatmap[:,0]))
    heatmap = np.array(map(str, heatmap.reshape(heatmap.size))).reshape(
        (shape[0], shape[1]))
    heatmap[:,0] = firstelem
    for i in xrange(shape[0]):
        file.write((';'.join(heatmap[i,:])+'\n'))
    file.close()
    LOGGER.info('Heatmap file is written as {0:s}'.format(outname))
        
        
    