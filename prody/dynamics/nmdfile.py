"""This module defines input and output functions for NMD format.

.. _nmd-format:

NMD Format
-------------------------------------------------------------------------------

Description
^^^^^^^^^^^

NMD files (extension :file:`.nmd`) are plain text files that contain at
least normal mode and system coordinate data.

NMD files can be visualized using :ref:`nmwiz`.  ProDy functions
:func:`.writeNMD` and :func:`.parseNMD` can be used to read and write NMD
files.

Data fields
^^^^^^^^^^^

Data fields in bold face are required. All data arrays and lists must be in a
single line and items must be separated by one or more space characters.

**coordinates**: system coordinates as a list of decimal numbers
  Coordinate array is the most important line in an NMD file. All mode array
  lengths must match the length of the coordinate array. Also, number of atoms
  in the system is deduced from the length of the coordinate array.

::

  coordinates 27.552 4.354 23.629 24.179 4.807 21.907 ...

**mode**: normal mode array as a list of decimal numbers
  Optionally, mode index and a scaling factor may be provided
  in the same line as a mode array. Both of these must precede the mode array.
  Providing a scaling factor enables relative scaling of the mode arrows and
  the amplitude of the fluctuations in animations. For NMA, scaling factors
  may be chosen to be the square-root of the inverse-eigenvalue associated
  with the mode. Analogously, for PCA data, scaling factor would be the
  square-root of the eigenvalue.

  If a mode line contains numbers preceding the mode array, they are evaluated
  based on their type. If an integer is encountered, it is considered the mode
  index. If a decimal number is encountered, it is considered the scaling
  factor. Scaling factor may be the square-root of the inverse eigenvalue
  if data is from an elastic network model, or the square-root of the
  eigenvalue if data is from an essential dynamics (or principal component)
  analysis.

  For example, all of the following lines are valid. The first line contains
  mode index and scaling factor. Second and third lines contain mode index or
  scaling factor. Last line contains only the mode array.

::

  mode 1 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 1    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 2.37    0.039 0.009 0.058 0.038 -0.011 0.052  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 0.043  ...

*name*: name of the model

The length of all following data fields must be equal to the number of atoms in
the system. NMWiz uses such data when writing a temporary PDB files for
loading coordinate data into VMD.

*atomnames*: list of atom names
  If not provided, all atom names are set to "CA".

*resnames*: list of residue names
  If not provided, all residue names are set to "GLY".

*chainids*: list of chain identifiers
  If not provided, all chain identifiers are set to "A".

*resids*: list of residue numbers
  If not provided, residue numbers are started from 1 and incremented by one
  for each atom.

*bfactors*: list of experimental beta-factors
  If not provided, all beta-factors are set to zero.
  Beta-factors can be used to color the protein representation.

NMD files may contain additional lines. Only lines that start with one of the
above field names are evaluated by NMWiz.


Autoload Trick
^^^^^^^^^^^^^^

By adding a special line in an NMD file, file content can be automatically
loaded into VMD at startup. The first line calls a NMWiz function to load the
file itself (:file:`xyzeros.nmd`).

::

  nmwiz_load xyzeros.nmd
  coordinates 0 0 0 0 0 0  ...
  mode 0.039 0.009 0.058 0.038 -0.011 0.052 ...
  mode -0.045 -0.096 -0.009 -0.040 -0.076 -0.010 ...
  mode 0.007 -0.044 0.080 0.015 -0.037 0.062 ...


In this case, VMD must be started from the command line by typing
:program:`vmd -e xyzeros.nmd`."""


__all__ = ['parseNMD', 'writeNMD', 'pathVMD', 'getVMDpath', 'setVMDpath',
           'viewNMDinVMD']

import os
from os.path import abspath, join, split, splitext

import numpy as np

from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup
from prody.utilities import openFile, isExecutable, which, PLATFORM, addext

from .nma import NMA
from .anm import ANM
from .gnm import GNM, ZERO
from .pca import PCA
from .mode import Vector, Mode
from .modeset import ModeSet


def pathVMD(*path):
    """Return VMD path, or set it to be a user specified *path*."""

    if not path:
        path = SETTINGS.get('vmd', None)
        if isExecutable(path):
            return path
        else:
            LOGGER.warning('VMD path is not set by user, looking for it.')

            vmdbin = None
            vmddir = None
            if PLATFORM == 'Windows':
                if PY3K:
                    import winreg
                else:
                    import _winreg as winreg  # PY3K: OK
                for vmdversion in ('1.8.7', '1.9', '1.9.1'):
                    try:
                        key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                                'Software\\University of Illinois\\VMD\\' +
                                vmdversion)
                        vmddir = winreg.QueryValueEx(key, 'VMDDIR')[0]
                        vmdbin = join(vmddir, 'vmd.exe')
                    except:
                        pass
                    try:
                        key = winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                    'Software\\WOW6432node\\University of Illinois\\VMD\\' +
                    vmdversion)
                        vmddir = winreg.QueryValueEx(key, 'VMDDIR')[0]
                        vmdbin = join(vmddir, 'vmd.exe')
                    except:
                        pass
            else:
                vmdbin = which('vmd')
                if False:
                    pipe = os.popen('which vmd')
                    vmdbin = pipe.next().strip()
                    vmdfile = open(vmdbin)
                    for line in vmdfile:
                        if line.startswith('defaultvmddir='):
                            vmddir = line.split('=')[1].replace('"', '')
                            break
                    vmdfile.close()
            if isExecutable(vmdbin):
                setVMDpath(vmdbin)
                return vmdbin
    elif len(path) == 1:
        path = path[0]
        if isExecutable(path):
            SETTINGS['vmd'] = path
            SETTINGS.save()
            LOGGER.info("VMD path is set to '{0}'.".format(path))
        else:
            raise OSError('{0} is not executable.'.format(str(path)))
    else:
        raise ValueError('specify a single path string')


def getVMDpath():
    """Deprecated for removal in v1.5, use :func:`pathVMD` instead."""

    return pathVMD()


def setVMDpath(path):
    """Deprecated for removal in v1.5, use :func:`pathVMD` instead."""

    return pathVMD(path)


NMD_LABEL_MAP = {
    'atomnames': 'name',
    'resnames': 'resname',
    'resnums': 'resnum',
    'resids': 'resnum',
    'chainids': 'chain',
    'bfactors': 'beta',
    'segnames': 'segment',
    'segments': 'segment',
}


def parseNMD(filename, type=None):
    """Return :class:`.NMA` and :class:`.AtomGroup` instances storing data
    parsed from *filename* in :file:`.nmd` format.  Type of :class:`.NMA`
    instance, e.g. :class:`.PCA`, :class:`.ANM`, or :class:`.GNM` will
    be determined based on mode data."""

    assert not isinstance(type, NMA), 'type must be NMA, ANM, GNM, or PCA'

    atomic = {}
    atomic.update([(label, None) for label in NMD_LABEL_MAP])
    atomic['coordinates'] = None
    atomic['name'] = None
    modes = []

    with open(filename) as nmd:
        for i, line in enumerate(nmd):
            try:
                label, data = line.split(None, 1)
            except ValueError:
                pass

            if label == 'mode':
                modes.append((i + 1, data))
            elif label in atomic:
                if atomic[label] is None:
                    atomic[label] = (i + 1, data)
                else:
                    LOGGER.warn('Data label {0} is found more than once in '
                                '{1}.'.format(repr(label), repr(filename)))

    name = atomic.pop('name', '')[1].strip() or splitext(split(filename)[1])[0]
    ag = AtomGroup(name)
    dof = None
    n_atoms = None

    line, coords = atomic.pop('coordinates', None)
    if coords is not None:
        coords = np.fromstring(coords, dtype=float, sep=' ')
        dof = coords.shape[0]
        if dof % 3 != 0:
            LOGGER.warn('Coordinate data in {0} at line {1} is corrupt '
                        'and will be omitted.'.format(repr(filename), line))
        else:
            n_atoms = dof / 3
            coords = coords.reshape((n_atoms, 3))
            ag.setCoords(coords)

    from prody.atomic import ATOMIC_FIELDS

    for label, data in atomic.items():  # PY3K: OK
        if data is None:
            continue
        line, data = data
        data = data.split()
        if n_atoms is None:
            n_atoms = len(data)
            dof = n_atoms * 3
        elif len(data) != n_atoms:
            LOGGER.warn('Data with label {0} in {1} at line {2} is '
                        'corrupt, expected {2} values, parsed {3}.'.format(
                        repr(label), repr(filename), line, n_atoms, len(data)))
            continue
        label = NMD_LABEL_MAP[label]
        data = np.array(data, dtype=ATOMIC_FIELDS[label].dtype)
        ag.setData(label, data)

    if not modes:
        return None, ag

    length = len(modes[0][1].split())
    is3d = length > n_atoms + 2
    if dof is None:
        dof = length - (length % 3)
    elif not is3d:  # GNM
        dof = n_atoms

    array = np.zeros((dof, len(modes)))
    less = 0
    eigvals = []
    count = 0
    for i, (line, mode) in enumerate(modes):
        mode = np.fromstring(mode, dtype=float, sep=' ')
        diff = len(mode) - dof
        if diff < 0 or diff > 2:
            LOGGER.warn('Mode data in {0} at line {1} is corrupt.'
                        .format(repr(filename), line))
            continue
        array[:, i - less] = mode[diff:]
        count += 1
        eigvals.append(mode[:diff])

    if count == 0:
        return None, ag

    try:
        eigvals = np.array(eigvals, dtype=float)
    except TypeError:
        LOGGER.warn('Failed to parse eigenvalues from {0}.'
                    .format(repr(filename)))

    if eigvals.shape[1] > 2:
        LOGGER.warn('Failed to parse eigenvalues from {0}.'
                    .format(repr(filename)))
        eigvals = None
    elif eigvals.shape[1] == 1:
        if np.all(eigvals % 1 == 0):
            LOGGER.warn('Failed to parse eigenvalues from {0}.'
                        .format(repr(filename)))
            eigvals = None
        else:
            eigvals = eigvals.flatten() ** 2
    else:
        eigvals = eigvals[:, 1] ** 2

    if is3d:
        if eigvals is not None and np.all(eigvals[:-1] >= eigvals[1:]):
            nma = PCA(name)
        else:
            nma = ANM(name)
    else:
        nma = GNM(name)
    if count != array.shape[1]:
        array = array[:, :count].copy()

    nma.setEigens(array, eigvals)
    return nma, ag


def writeNMD(filename, modes, atoms, zeros=False):
    """Return *filename* that contains *modes* and *atoms* data in NMD format
    described in :ref:`nmd-format`.  :file:`.nmd` extension is appended to
    filename, if it does not have an extension.

    .. note::
       #. If zeros is **False** (by default), this function skips modes 
          with zero eigenvalues. If zeros is **True**, modes with zero 
          eigenvalues are written out, their scaling factor being the 
          square-root of the inverse of the mode number times 0.0001.
          This provides descending factors consistent with the NMA modes.
       #. If a :class:`.Vector` instance is given, it will be normalized
          before it is written. It's length before normalization will be
          written as the scaling factor of the vector."""

    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0}'.format(type(modes)))
    if modes.numAtoms() != atoms.numAtoms():
        raise Exception('number of atoms do not match')
    out = openFile(addext(filename, '.nmd'), 'w')

    #out.write('#!{0} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0}\n'.format(abspath(filename)))
    name = modes.getTitle()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = splitext(split(filename)[1])[0]
    out.write('name {0}\n'.format(name))
    try:
        coords = atoms.getCoords()
    except:
        raise ValueError('coordinates could not be retrieved from atoms')
    if coords is None:
        raise ValueError('atom coordinates are not set')

    try:
        data = atoms.getNames()
        if data is not None:
            out.write('atomnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnames()
        if data is not None:
            out.write('resnames {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnums()
        if data is not None:
            out.write('resids ')
            data.tofile(out, ' ')
            out.write('\n')
    except:
        pass
    try:
        data = atoms.getChids()
        if data is not None:
            out.write('chainids {0}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getSegnames()
        if data is not None:
            out.write('segnames {0}\n'.format(' '.join(data)))
    except:
        pass

    try:
        data = atoms.getBetas()
        if data is not None:
            out.write('bfactors ')
            data.tofile(out, ' ', '%.2f')
            out.write('\n')
    except:
        pass

    format = '{0:.3f}'.format
    out.write('coordinates ')
    coords.tofile(out, ' ', '%.3f')
    out.write('\n')
    count = 0
    if isinstance(modes, Vector):
        out.write('mode 1 {0:.2f} '.format(abs(modes)))
        modes.getNormed()._getArray().tofile(out, ' ', '%.3f')
        out.write('\n')
        count += 1
    else:
        if isinstance(modes, Mode):
            modes = [modes]
        for mode in modes:
            if (mode.getEigval() < ZERO) and not zeros:
                continue
            elif (mode.getEigval() < ZERO) and zeros:
                out.write('mode {0} {1:.2f} '.format(
                mode.getIndex()+1, np.sqrt(1/(0.0001*(mode.getIndex()+1)))))
            else:
                out.write('mode {0} {1:.2f} '.format(
                mode.getIndex()+1, mode.getVariance()**0.5))
            arr = mode._getArray().tofile(out, ' ', '%.3f')
            out.write('\n')
            count += 1
    if count == 0:
        LOGGER.warning('No normal mode data was written. '
                       'Given modes might have 0 eigenvalues.')
    out.close()
    return filename


def viewNMDinVMD(filename):
    """Start VMD in the current Python session and load NMD data."""

    vmd = pathVMD()
    if vmd:
        os.system('{0} -e {1}'.format(vmd, abspath(filename)))
