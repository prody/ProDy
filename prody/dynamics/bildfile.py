"""This module defines input and output functions for BILD format for ChimeraXs.

.. _bild-format: https://rbvi.ucsf.edu/chimerax/docs/user/formats/bild.html

"""


__all__ = ['writeBILD']

import os
from os.path import join, split, splitext

import numpy as np

from prody import LOGGER
from prody.utilities import openFile, getCoords

from .nma import NMA
from .anm import ANM
from .gnm import GNM, ZERO
from .pca import PCA
from .mode import Vector, Mode
from .modeset import ModeSet

def writeBILD(filename, modes, atoms, scale=25.):
    """Returns *filename* that contains *modes* and *atoms* data in BILD format
    described at :ref:`bild-format`.  :file:`.bild` extension is appended to
    filename, if it does not have an extension. 
    
    Arrows will be scaled by the square root of the variance times *scale*.

    .. note::
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`.Vector` instance is given, it will be normalized
          before it is written. It's length before normalization will be
          written as the scaling factor of the vector.
       #. It is recommended to only use one mode to avoid having lots of arrows.
    """

    if not filename.lower().endswith(".bild"):
        filename += '.bild'

    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0}'.format(type(modes)))
    if modes.numAtoms() != len(atoms):
        raise Exception('number of atoms do not match')
    out = openFile(filename, 'w')

    name = modes.getTitle()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = splitext(split(filename)[1])[0]
    out.write('.comment name {0}\n'.format(name))
    try:
        coords = getCoords(atoms)
    except:
        raise ValueError('coordinates could not be retrieved from atoms')
    if coords is None:
        raise ValueError('atom coordinates are not set')

    try:
        data = atoms.getNames()
        if data is not None:
            out.write('.comment atomnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnames()
        if data is not None:
            out.write('.comment resnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnums()
        if data is not None:
            out.write('.comment resids ')
            data.tofile(out, ' ')
            out.write('\n')
    except:
        pass
    try:
        data = atoms.getChids()
        if data is not None:
            out.write('.comment chainids {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getSegnames()
        if data is not None:
            out.write('.comment segnames {0}\n'.format(' '.join(data)))
    except:
        pass

    try:
        data = atoms.getBetas()
        if data is not None:
            out.write('.comment bfactors ')
            data.tofile(out, ' ', '%.2f')
            out.write('\n')
    except:
        pass

    format = '{0:.3f}'.format
    out.write('.comment coordinates ')
    coords.tofile(out, ' ', '%.3f')
    out.write('\n')
    count = 0
    if isinstance(modes, Vector):
        out.write('.comment mode 1 {0:.2f} \n'.format(abs(modes)))
        normed = modes.getNormed()._getArray() * scale
        for i in range(modes.numAtoms()):
            out.write('.arrow ')
            coords[3*i:3*i+1].tofile(out, ' ', '%.3f')
            normed[3*i:3*i+1].tofile(out, ' ', '%.3f')
            out.write('\n')
        out.write('\n')
        count += 1
    else:
        if isinstance(modes, Mode):
            modes = [modes]
        else:
            LOGGER.warning('It is recommended to only provide one mode '
                'to writeBILD as otherwise there will be a mess of arrows.')
        for mode in modes:
            if mode.getEigval() < ZERO:
                continue
            out.write('.comment mode {0} {1:.2f} \n'.format(
                       mode.getIndex()+1, mode.getVariance()**0.5))
            
            arr = mode._getArray() * mode.getVariance()**0.5 * scale
            for i in range(mode.numAtoms()):
                out.write('.arrow ')
                coords[i].tofile(out, ' ', '%.3f')
                out.write(' ')
                
                arr_out = arr[3*i:3*i+3] + coords[i]
                arr_out.tofile(out, ' ', '%.3f')
                out.write('\n')
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. '
                       'Given modes might have 0 eigenvalues.')
    out.close()
    return filename
