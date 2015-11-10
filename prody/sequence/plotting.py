# -*- coding: utf-8 -*-
"""This module defines MSA analysis functions."""

__author__ = 'Ahmet Bakan, Chakra Chennubhotla'

from numpy import arange

from .analysis import *
from prody import LOGGER

__all__ = ['showMSAOccupancy', 'showShannonEntropy', 'showMutinfoMatrix', 
           'showDirectInfoMatrix', 'showSCAMatrix']


def pickSequence(msa):
    """Pick a sequence without gaps and deletions and return its residue
    numbers and labels to be used as indices and X-axis label, or a pair
    of **None** at failure."""

    try:
        counts = calcMSAOccupancy(msa, 'row', count=True)
    except TypeError:
        return None, None
    else:
        length = msa.numResidues()
        split, msa.split = msa.split, True
        rows = (counts == length).nonzero()[0]
        for row in rows:
            try:
                label, indices = msa[row].getLabel(), msa[row].getResnums()
            except:
                break
            else:
                return (indices, 'Residue number ({0})'.format(label))
        return None, None


def showMSAOccupancy(msa, occ='res', indices=None, count=False, **kwargs):
    """Show a bar plot of occupancy calculated for :class:`.MSA` instance *msa*
    using :func:`.calcMSAOccupancy`.  *occ* may be ``'res'`` or ``'col'``, or a
    a pre-calculated occupancy array.  If x-axis *indices* are not specified,
    they will be inferred from *msa* or given *label* that may correspond to
    a sequence in the msa.

    Occupancy is plotted using :func:`~matplotlib.pyplot.bar` function."""

    try:
        numseq, lenseq = msa.numSequences(), msa.numResidues()
    except AttributeError:
        raise TypeError('msa must be an MSA instance')

    try:
        length = len(occ)
    except TypeError:
        raise TypeError("occ must be 'res', 'col', or an occupancy array")

    try:
        sw = occ.startswith
    except (TypeError, AttributeError):
        try:
            ndim = occ.ndim
        except AttributeError:
            raise TypeError("occ must be 'res', 'col', or an occupancy array")
        else:
            if ndim != 1:
                raise ValueError('occ must be a 1-dimensional array')
    else:
        occ = calcMSAOccupancy(msa, occ, count)
        length = len(occ)

    if length == numseq:
        xlabel = kwargs.pop('xlabel', None) or 'MSA sequence index'

    if length == lenseq:
        label = kwargs.pop('label', None)
        if label is not None:
            try:
                indices = msa[label].getResnums()
            except:
                LOGGER.info('Specified label not in msa.')
        xlabel = kwargs.pop('xlabel', None) or label or 'MSA column index'
        if xlabel is None and indices is None:
            indices, xlabel = pickSequence(msa)
    if indices is None:
        indices = arange(1, length + 1)

    ylabel = kwargs.pop('ylabel', 'Count' if count else 'Occupancy')
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)
    import matplotlib.pyplot as plt
    show = plt.bar(indices, occ, **kwargs)
    if format:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title is None:
            title = 'Occupancy: ' + msa.getTitle()
        plt.title(title)
    return show


def showShannonEntropy(entropy, indices=None, **kwargs):
    """Show a bar plot of Shannon *entropy* array.  :class:`.MSA` instances
    or Numpy character arrays storing sequence alignments are also accepted
    as *entropy* argument, in which case :func:`.calcShannonEntropy` will
    be used for calculations.  *indices* may be residue numbers, when not
    specified they will be inferred from *msa* or indices starting from 1
    will be used.

    Entropy is plotted using :func:`~matplotlib.pyplot.bar` function."""

    msa = None
    try:
        ndim = entropy.ndim
    except AttributeError:
        msa = entropy
        entropy = calcShannonEntropy(msa)
        ndim = entropy.ndim

    if ndim != 1:
        raise ValueError('entropy must be a 1D array')

    msa = kwargs.pop('msa', msa)
    xlabel = kwargs.pop('xlabel', None)
    if indices is None:
        length = len(entropy)
        if msa is not None:
            indices, xlabel = pickSequence(msa)
        if indices is None:
            indices = arange(1, length + 1)
        xlabel = xlabel or 'MSA column index'

    ylabel = kwargs.pop('ylabel', 'Shannon entropy')
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)
    import matplotlib.pyplot as plt
    show = plt.bar(indices, entropy, **kwargs)
    if format:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if title is None:
            if msa is None:
                title = 'Entropy'
            else:
                title = 'Entropy: ' + str(msa)
        plt.title(title)
    return show


def showMutinfoMatrix(mutinfo, *args, **kwargs):
    """Show a heatmap of mutual information array.  :class:`.MSA` instances
    or Numpy character arrays storing sequence alignment are also accepted
    as *mutinfo* argument, in which case :func:`.buildMutinfoMatrix` will
    be used for calculations.  Note that x, y axes contain indices of the
    matrix starting from 1.

    Mutual information is plotted using :func:`~matplotlib.pyplot.imshow`
    function. vmin and vmax values can be set by user to achieve better
    signals using this function."""

    msa = None
    try:
        ndim, shape = mutinfo.ndim, mutinfo.shape
    except AttributeError:
        msa = mutinfo
        mutinfo = buildMutinfoMatrix(mutinfo)
        ndim, shape = mutinfo.ndim, mutinfo.shape

    msa = kwargs.pop('msa', msa)
    if ndim != 2:
        raise ValueError('mutinfo must be a 2D matrix')
    y, x = shape
    if x != y:
        raise ValueError('mutinfo matrix must be a square matrix')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('origin', 'lower')

    if msa is not None:
        indices, msalabel = pickSequence(msa)
        if indices is not None:
            start = indices[0] + 0.5
            end = start + x
            extent = [start, end, start, end]
        else:
            extent = [0.5, x + 0.5, 0.5, y + 0.5]
    else:
        msalabel = None
        extent = [0.5, x + 0.5, 0.5, y + 0.5]

    xlabel = kwargs.pop('xlabel', None)
    if xlabel is None:
        xlabel = msalabel or 'MSA column index'
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)

    import matplotlib.pyplot as plt
    show = plt.imshow(mutinfo, extent=extent, *args, **kwargs), plt.colorbar()

    if format:
        plt.xlabel(xlabel)
        plt.ylabel(xlabel)
        if title is None:
            if msa is None:
                title = 'Mutual Information'
            else:
                title = 'Mutual Information: ' + str(msa)
        plt.title(title)
    return show

def showDirectInfoMatrix(dirinfo, *args, **kwargs):
    """Show a heatmap of direct information array.  :class:`.MSA` instances
    or Numpy character arrays storing sequence alignment are also accepted
    as *dirinfo* argument, in which case :func:`.buildDirectInfoMatrix` will
    be used for calculations.  Note that x, y axes contain indices of the
    matrix starting from 1.

    Direct information is plotted using :func:`~matplotlib.pyplot.imshow`
    function. vmin and vmax values can be set by user to achieve better
    signals using this function."""

    msa = None
    try:
        ndim, shape = dirinfo.ndim, dirinfo.shape
    except AttributeError:
        msa = dirinfo
        dirinfo = buildDirectInfoMatrix(dirinfo)
        ndim, shape = dirinfo.ndim, dirinfo.shape

    msa = kwargs.pop('msa', msa)
    if ndim != 2:
        raise ValueError('dirinfo must be a 2D matrix')
    y, x = shape
    if x != y:
        raise ValueError('dirinfo matrix must be a square matrix')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('origin', 'lower')

    if msa is not None:
        indices, msalabel = pickSequence(msa)
        if indices is not None:
            start = indices[0] + 0.5
            end = start + x
            extent = [start, end, start, end]
        else:
            extent = [0.5, x + 0.5, 0.5, y + 0.5]
    else:
        msalabel = None
        extent = [0.5, x + 0.5, 0.5, y + 0.5]

    xlabel = kwargs.pop('xlabel', None)
    if xlabel is None:
        xlabel = msalabel or 'MSA column index'
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)

    import matplotlib.pyplot as plt
    show = plt.imshow(dirinfo, extent=extent, *args, **kwargs), plt.colorbar()

    if format:
        plt.xlabel(xlabel)
        plt.ylabel(xlabel)
        if title is None:
            if msa is None:
                title = 'Direct Information'
            else:
                title = 'Direct Information: ' + str(msa)
        plt.title(title)
    return show


def showSCAMatrix(scainfo, *args, **kwargs):
    """Show a heatmap of SCA (statistical coupling analysis) array.  
    :class:`.MSA` instances. blah

    or Numpy character arrays storing sequence alignment are also accepted
    as *scainfo* argument, in which case :func:`.buildSCAMatrix` will
    be used for calculations.  Note that x, y axes contain indices of the
    matrix starting from 1.

    SCA information is plotted using :func:`~matplotlib.pyplot.imshow`
    function. vmin and vmax values can be set by user to achieve better
    signals using this function."""

    msa = None
    try:
        ndim, shape = scainfo.ndim, scainfo.shape
    except AttributeError:
        msa = scainfo
        scainfo = buildSCAMatrix(scainfo)
        ndim, shape = scainfo.ndim, scainfo.shape
    

    msa = kwargs.pop('msa', msa)
    if ndim != 2:
        raise ValueError('scainfo must be a 2D matrix')
    y, x = shape
    if x != y:
        raise ValueError('scainfo matrix must be a square matrix')

    kwargs.setdefault('interpolation', 'nearest')
    kwargs.setdefault('origin', 'lower')

    if msa is not None:
        indices, msalabel = pickSequence(msa)
        if indices is not None:
            start = indices[0] + 0.5
            end = start + x
            extent = [start, end, start, end]
        else:
            extent = [0.5, x + 0.5, 0.5, y + 0.5]
    else:
        msalabel = None
        extent = [0.5, x + 0.5, 0.5, y + 0.5]

    xlabel = kwargs.pop('xlabel', None)
    if xlabel is None:
        xlabel = msalabel or 'MSA column index'
    title = kwargs.pop('title', None)
    format = kwargs.pop('format', True)

    import matplotlib.pyplot as plt
    show = plt.imshow(scainfo, extent=extent, *args, **kwargs), plt.colorbar()

    if format:
        plt.xlabel(xlabel)
        plt.ylabel(xlabel)
        if title is None:
            if msa is None:
                title = 'SCA Information'
            else:
                title = 'SCA Information: ' + str(msa)
        plt.title(title)
    return show



if __name__ == '__main__':
    from prody import *
    msa = parseMSA('piwi_seed.sth')
    print(repr(msa))
    msa = refineMSA(msa, label=msa[0][0])
    print(repr(msa))
    print(calcMSAOccupancy(msa, 'row', count=True))
    print(pickSequence(msa))
