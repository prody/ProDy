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

"""This module defines functions for supporting VMD Heatmapper_ plugin format. 

. _Heatmapper: http://www.ks.uiuc.edu/Research/vmd/plugins/heatmapper/"""

__author__ = 'Ahmet Bakan, Anindita Dutta'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from numpy import arange, fromstring, array

from .analysis import *

from prody.utilities import openFile, intorfloat

from prody import LOGGER


__all__ = ['parseHeatmap', 'writeHeatmap', 'showHeatmap']

HMTYPES = {
    'max': float,
    'min': float,
    'numbering': lambda line: line.split(':'),
    'title': str,
    'xlabel': str,
    'ylabel': str,
    'xorigin': intorfloat,
    'xstep': intorfloat,
}



def showHeatmap(heatmap, clim=None, *args, **kwargs):
    """Show *heatmap*, which can be an 2D matrix or an :file:`.hm` file in 
    VMD Heatmapper plugin  format. 
    
    Heatmap is plotted using :func:`~matplotlib.pyplot.imshow`
    function."""

    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except AttributeError:
        headers, heatmap = parseHeatmap(heatmap)
        ndim, shape = heatmap.ndim, heatmap.shape
    
    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')
    
    if not 'interpolation' in kwargs:
        kwargs['interpolation'] = 'nearest'
    
    if not 'origin' in kwargs:
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


def parseHeatmap(heatmap, **kwargs):
    """Return a two dimensional array and a dictionary with information parsed
    from *heatmap*, which may be an input stream or an :file:`.hm` file in VMD
    plugin Heatmapper format."""    
    
    try:
        readline, close = heatmap.readline, lambda: None
    except AttributeError:
        heatmap = openFile(heatmap)
        readline, close = heatmap.readline, heatmap.close 
    
    meta = {}
    arrs = []
    
    line = readline()
    while line:
        if line.startswith('-'):
            label, data = line[1:].split(None, 1)
            data = data.strip()
            if data[0] == data[-1] == '"':
                data = data[1:-1]
            label = label.strip()
            try:
                meta[label] = HMTYPES[label](data)
            except KeyError: 
                LOGGER.warn('Unrecognized label encountered: {0}'
                            .format(repr(label)))
                meta[label] = HMTYPES[label](data)
            except TypeError:
                LOGGER.warn('Could not parse data with label {0}.'
                            .format(repr(label)))
        else:
            arrs.append(line.rstrip())
        line = readline()
    close()
    nnums = len(meta.get('numbering', '')) 
    heatmap = []
    numbers = []

    for arr in arrs:
        if nnums:
            items = arr.split(':', nnums + 1)
            numbers.append(items[:nnums])
        else:
            items = [arr]
        heatmap.append(fromstring(items[-1], float, sep=';'))
        
    heatmap = array(heatmap)
    if nnums:
        numbering = meta['numbering']
        try:
            numbers = array(numbers, int)
        except ValueError:
            try:
                numbers = array(numbers, float)
            except ValueError:            
                LOGGER.warn('Numbering for y-axis could not be parsed.')
                numbering = []
        for i, label in enumerate(numbering):
            meta[label] = numbers[:, i].copy()
        
    return heatmap, meta
            
        
def writeHeatmap(filename, heatmap, **kwargs):
    """Return *filename* that contains *heatmap* in :file:`.hm` format of VMD 
    plugin Heatmapper.  *filename* may also be an output stream.
    
    :arg title: title of the heatmap
    :type title: str
    
    :arg xlabel: x-axis lab, default is ``'unknown'``
    :type xlabel: str

    :arg ylabel: y-axis lab, default is ``'unknown'``
    :type ylabel: str

    :arg xorigin: x-axis origin, default is 0
    :type xorigin: float

    :arg xstep: x-axis step, default is 1
    :type xstep: float

    :arg min: minimum value, default is minimum in *heatmap*
    :type min: float

    :arg max: maximum value, default is maximum in *heatmap*
    :type max: float

    Other keyword arguments that are arrays with length equal to the y-axis
    (second dimension of heatmap) will be considered as *numbering*."""
    
    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except:
        raise TypeError('heatmap must be an array object')
    if ndim!= 2:
        raise TypeError('heatmap must be a 2D array')
    
    try:
        write, close = filename.write, lambda: None
    except AttributeError: 
        out = openFile(filename, 'wb')
        write, close = out.write, out.close
    
    write('-min "{0}"\n'.format(kwargs.pop('min', heatmap.min())))
    write('-max "{0}"\n'.format(kwargs.pop('max', heatmap.max())))
    for label, default in [
        ('title', 'unknown'),
        ('xlabel', 'unknown'),
        ('xorigin', 0),
        ('xstep', 1),
        ('ylabel', 'unknown'),
    ]:
        write('-{0} "{1}"\n'.format(label, kwargs.pop(label, default)))
    
    numbering = []
    numlabels = []
    for key, val in kwargs.items():
        try:
            length = len(val)
        except TypeError:
            LOGGER.warn('Keyword argument {0} is not used.'.format(key))
        if length == shape[0]:
            numlabels.append(key)
            numbering.append(val)
    if not numbering:
        numlabels.append('unknown')
        numbering.append(arange(1, shape[0] + 1))
    
    write('-numbering "{0}"\n'.format(':'.join(numlabels)))
    
    for i, row in enumerate(heatmap):
        write(':'.join(str(nums[i])for nums in numbering) + ':')
        write(';'.join(row.astype('|S10')))
        write(';\n')
    
    close()
    return filename
