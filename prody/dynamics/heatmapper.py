# -*- coding: utf-8 -*-
"""This module defines functions for supporting VMD plugin `Heatmapper`_
format files.

.. _Heatmapper: http://www.ks.uiuc.edu/Research/vmd/plugins/heatmapper/"""

__author__ = 'Ahmet Bakan, Anindita Dutta'

from numpy import arange, fromstring, frombuffer, array

from prody.utilities import openFile, intorfloat, startswith, addext

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


def showHeatmap(heatmap, *args, **kwargs):
    """Show *heatmap*, which can be an two dimensional array or a Heat Mapper
    :file:`.hm` file.

    Heatmap is plotted using :func:`~matplotlib.pyplot.imshow` function.
    Default values passed to this function are ``interpolation='nearest'``,
    ``aspect='auto'``, and ``origin='lower'``."""

    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except AttributeError:
        heatmap, headers = parseHeatmap(heatmap)
        ndim, shape = heatmap.ndim, heatmap.shape

        xorigin = headers.pop('xorigin', 0)
        xextent = headers.pop('xstep', 1) * shape[0]

        ylabel = kwargs.get('ylabel', '').lower()
        indices = None
        if ylabel:
            for key in headers.get('numbering', []):
                if startswith(ylabel, key.lower()):
                    indices = headers.get(key)
        if indices is not None:
            extent = [indices[0] - .5, indices[0] + len(indices) - .5,
                      xorigin - .5, xextent - .5]
        else:
            extent = [-.5, shape[1] * 2 - .5, xorigin - .5, xextent - .5]
        kwargs.setdefault('extent', extent)

        for key in ['numbering', 'min', 'max'] + headers.get('numbering', []):
            headers.pop(key, None)
        headers.update(kwargs)
        kwargs = headers


    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('origin', 'lower')
    kwargs.setdefault('aspect', 'auto')

    xlabel = kwargs.pop('xlabel', None)
    ylabel = kwargs.pop('ylabel', None)
    title = kwargs.pop('title', None)

    import matplotlib.pyplot as plt
    show = plt.imshow(heatmap, *args, **kwargs), plt.colorbar()

    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    return show


def parseHeatmap(heatmap, **kwargs):
    """Returns a two dimensional array and a dictionary with information parsed
    from *heatmap*, which may be an input stream or an :file:`.hm` file in VMD
    plugin Heat Mapper format."""

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
            
        try:
            heatmap.append(fromstring(items[-1], float, sep=';'))
        except:
            heatmap.append(frombuffer(items[-1], float, sep=';'))

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
    """Returns *filename* that contains *heatmap* in Heat Mapper :file:`.hm`
    file (extension is automatically added when not found).  *filename* may
    also be an output stream.

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

    :arg format: number format, default is ``'%f'``
    :type format: str

    Other keyword arguments that are arrays with length equal to the y-axis
    (second dimension of heatmap) will be considered as *numbering*."""

    try:
        ndim, shape = heatmap.ndim, heatmap.shape
    except:
        raise TypeError('heatmap must be an array object')
    if ndim!= 2:
        raise TypeError('heatmap must be a 2D array')

    try:
        write, close, stream = filename.write, lambda: None, filename
    except AttributeError:
        out = openFile(addext(filename, '.hm'), 'w')
        write, close, stream = out.write, out.close, out

    format = kwargs.pop('format', '%f')
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
        else:
            if length == shape[0]:
                numlabels.append(key)
                numbering.append(val)
    if not numbering:
        numlabels.append('unknown')
        numbering.append(arange(1, shape[0] + 1))

    write('-numbering "{0}"\n'.format(':'.join(numlabels)))

    for i, row in enumerate(heatmap):
        write(':'.join(str(nums[i]) for nums in numbering) + ':')
        row.tofile(stream, sep=';', format=format)
        write(';\n')

    close()
    return filename
