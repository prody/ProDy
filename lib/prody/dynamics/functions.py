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

"""This module defines input and output functions.""" 

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import os
from os.path import abspath, join, isfile, isdir, split, splitext

import numpy as np

from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import AtomGroup
from prody.utilities import openFile, isExecutable, which, PLATFORM, addext

from .nma import NMA
from .anm import ANM
from .gnm import GNM, GNMBase, ZERO
from .pca import PCA, EDA
from .mode import Vector, Mode
from .modeset import ModeSet

__all__ = ['parseArray', 'parseModes', 'parseNMD',
           'writeArray', 'writeModes', 'parseSparseMatrix',
           'writeNMD', 
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD',]

           
def saveModel(nma, filename=None, matrices=False, **kwargs):
    """Save *nma* model data as :file:`filename.nma.npz`.  By default, 
    eigenvalues, eigenvectors, variances, trace of covariance matrix, 
    and name of the model will be saved.  If *matrices* is ``True``, 
    covariance, Hessian or Kirchhoff matrices are saved too, whichever 
    are available.  If *filename* is ``None``, name of the NMA instance 
    will be used as the filename, after ``" "`` (white spaces) in the name 
    are replaced with ``"_"`` (underscores).  Extension may differ based 
    on the type of the NMA model.  For ANM models, it is :file:`.anm.npz`.
    Upon successful completion of saving, filename is returned. This 
    function makes use of :func:`numpy.savez` function."""
    
    if not isinstance(nma, NMA):
        raise TypeError('invalid type for nma, {0}'.format(type(nma)))
    if len(nma) == 0:
        raise ValueError('nma instance does not contain data')
    
    dict_ = nma.__dict__
    attr_list = ['_title', '_trace', '_array', '_eigvals', '_vars', '_n_atoms',
                 '_dof', '_n_modes']
    if filename is None:
        filename = nma.getTitle().replace(' ', '_')
    if isinstance(nma, GNMBase):
        attr_list.append('_cutoff')
        attr_list.append('_gamma')
        if matrices:
            attr_list.append('_kirchhoff')
            if isinstance(nma, ANM):
                attr_list.append('_hessian')
        if isinstance(nma, ANM):
            type_ = 'ANM'
        else:
            type_ = 'GNM'
    elif isinstance(nma, EDA):
        type_ = 'EDA'
    elif isinstance(nma, PCA):
        type_ = 'PCA'
    else:
        type_ = 'NMA'  
    
    if matrices:
        attr_list.append('_cov')
    attr_dict = {'type': type_}
    for attr in attr_list:
        value = dict_[attr]
        if value is not None:
            attr_dict[attr] = value
    filename += '.' + type_.lower() + '.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadModel(filename):
    """Return NMA instance after loading it from file (*filename*).  
    This function makes use of :func:`numpy.load` function.  See 
    also :func:`saveModel`."""
    
    attr_dict = np.load(filename)
    try:
        type_ = attr_dict['type']
    except KeyError:
        raise IOError('{0} is not a valid NMA model file'.format(filename))
    try:
        title = str(attr_dict['_title'])
    except KeyError: 
        title = str(attr_dict['_name'])
    if type_ == 'ANM':
        nma = ANM(title)
    elif type_ == 'PCA':
        nma = PCA(title)
    elif type_ == 'EDA':
        nma = EDA(title)
    elif type_ == 'GNM':
        nma = GNM(title)
    elif type_ == 'NMA':
        nma = NMA(title)
    else:
        raise IOError('NMA model type is not recognized'.format(type_))
    dict_ = nma.__dict__ 
    for attr in attr_dict.files:
        if attr in ('type', '_name', '_title'): 
            continue
        elif attr in ('_trace', '_cutoff', '_gamma'):
            dict_[attr] = float(attr_dict[attr])
        elif attr in ('_dof', '_n_atoms', '_n_modes'):
            dict_[attr] = int(attr_dict[attr])
        else:
            dict_[attr] = attr_dict[attr]
    return nma


def saveVector(vector, filename, **kwargs):
    """Save *vector* data as :file:`filename.vec.npz`.  Upon successful 
    completion of saving, filename is returned.  This function makes use 
    of :func:`numpy.savez` function."""
    
    if not isinstance(vector, Vector):
        raise TypeError('invalid type for vector, {0}'.format(type(vector)))
    attr_dict = {}
    attr_dict['title'] = vector.getTitle()
    attr_dict['array'] = vector._getArray()
    attr_dict['is3d'] = vector.is3d()
    filename += '.vec.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadVector(filename):
    """Return :class:`.Vector` instance after loading it from *filename* using
    :func:`numpy.load`.  See also :func:`saveVector`."""
    
    attr_dict = np.load(filename)
    try:
        title = str(attr_dict['title'])
    except KeyError:
        title = str(attr_dict['name'])
    return Vector(attr_dict['array'], title, bool(attr_dict['is3d']))


def getVMDpath():
    """Return VMD path set by user or one identified automatically."""
    
    path = SETTINGS.get('vmd', None)
    if isExecutable(path):
        return path   
    else:
        LOGGER.warning('VMD path is not set by user, looking for it.')    

        from types import StringType, UnicodeType
        vmdbin = None
        vmddir = None
        if PLATFORM == 'Windows': 
            if PY3K:
                import winreg
            else:
                import _winreg as winreg # PY3K: OK
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
        

def setVMDpath(path):
    """Set path to a VMD executable."""
    
    if isExecutable(path):
        SETTINGS['vmd'] = path
        SETTINGS.save()
        LOGGER.info("VMD path is set to '{0}'.".format(path))
    else:
        raise OSError('{0} is not executable.'.format(str(path)))


NMD_LABEL_MAP = {
    'atomnames': 'name', 
    'resnames': 'resname',
    'resnums': 'resnum', 
    'resids': 'resnum',
    'chainids': 'chain', 
    'bfactors': 'beta'
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

    for label, data in atomic.items(): # PY3K: OK
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
    is3d =  length > n_atoms + 2
    if dof is None: 
        dof = length - (length % 3)
    elif not is3d: # GNM
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
        array = array[:,:count].copy()
        
    nma.setEigens(array, eigvals)
    return nma, ag
    

def writeNMD(filename, modes, atoms):
    """Return *filename* that contains *modes* and *atoms* data in NMD format
    described in :ref:`nmd-format`.  :file:`.nmd` extension is appended to
    filename, if it does not have an extension.
    
    .. note:: 
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`~.Vector` instance is given, it will be normalized 
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
            if mode.getEigval() < ZERO:
                continue
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
    
    vmd = SETTINGS.get('vmd')
    if vmd:
        os.system('{0} -e {1}'.format(vmd, abspath(filename)))
        
def writeModes(filename, modes, format='%.18e', delimiter=' '):
    """Write *modes* (eigenvectors) into a plain text file with name 
    *filename*. See also :func:`writeArray`."""
    
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0}'
                        .format(type(modes)))
    return writeArray(filename, modes._getArray(), format=format, 
                      delimiter=delimiter)

def parseModes(normalmodes, eigenvalues=None, nm_delimiter=None, 
               nm_skiprows=0, nm_usecols=None, ev_delimiter=None, 
               ev_skiprows=0, ev_usecols=None, ev_usevalues=None):
    """Return :class:`~.NMA` instance with normal modes parsed from 
    *normalmodes*.
    
    In normal mode file *normalmodes*, columns must correspond to modes 
    (eigenvectors).  Optionally, *eigenvalues* can be parsed from a separate 
    file. If eigenvalues are not provided, they will all be set to 1.
    
    :arg normalmodes: File or filename that contains normal modes. 
        If the filename extension is :file:`.gz` or :file:`.bz2`, the file is 
        first decompressed.
    :type normalmodes: str or file
    
    :arg eigenvalues: Optional, file or filename that contains eigenvalues. 
        If the filename extension is :file:`.gz` or :file:`.bz2`, 
        the file is first decompressed.
    :type eigenvalues: str or file

    :arg nm_delimiter: The string used to separate values in *normalmodes*. 
        By default, this is any whitespace.
    :type nm_delimiter: str

    :arg nm_skiprows: Skip the first *skiprows* lines in *normalmodes*. 
        Default is ``0``.
    :type nm_skiprows: 0

    :arg nm_usecols: Which columns to read from *normalmodes*, with 0 being the 
        first. For example, ``usecols = (1,4,5)`` will extract the 2nd, 5th and 
        6th columns. The default, ``None``, results in all columns being read.
    :type nm_usecols: list

    :arg ev_delimiter: The string used to separate values in *eigenvalues*. 
        By default, this is any whitespace.
    :type ev_delimiter: str

    :arg ev_skiprows: Skip the first *skiprows* lines in *eigenvalues*. 
        Default is ``0``.
    :type ev_skiprows: 0

    :arg ev_usecols: Which columns to read from *eigenvalues*, with 0 being the 
        first. For example, ``usecols = (1,4,5)`` will extract the 2nd, 5th and 
        6th columns. The default, ``None``, results in all columns being read.
    :type ev_usecols: list

    :arg ev_usevalues: Which columns to use after the eigenvalue column is
        parsed from *eigenvalues*, with 0 being the first. 
        This can be used if *eigenvalues* contains more values than the
        number of modes in *normalmodes*.
    :type ev_usevalues: list
    
    See :func:`parseArray` for details of parsing arrays from files."""
    
    modes = parseArray(normalmodes, delimiter=nm_delimiter, 
                       skiprows=nm_skiprows, usecols=nm_usecols)
    if eigenvalues is not None:
        values = parseArray(eigenvalues, delimiter=ev_delimiter, 
                            skiprows=ev_skiprows, usecols=ev_usecols)
        values = values.flatten()
        if ev_usevalues is not None:
            values = values[ev_usevalues]
    nma = NMA(splitext(split(normalmodes)[1])[0])
    nma.setEigens(modes, values)
    return nma
    
    
def writeArray(filename, array, format='%d', delimiter=' '):
    """Write 1-d or 2-d array data into a delimited text file.
    
    This function is using :func:`numpy.savetxt` to write the file, after 
    making some type and value checks.  Default *format* argument is ``"%d"``.
    Default *delimiter* argument is white space, ``" "``.
    
    *filename* will be returned upon successful writing."""
    
    if not isinstance(array, np.ndarray):
        raise TypeError('array must be a Numpy ndarray, not {0}'
                        .format(type(array)))
    elif not array.ndim in (1, 2):
        raise ValueError('array must be a 1 or 2-dimensional Numpy ndarray, '
                         'not {0}-d'.format(type(array.ndim)))
    np.savetxt(filename, array, format, delimiter)
    return filename

def parseArray(filename, delimiter=None, skiprows=0, usecols=None, 
               dtype=float):
    """Parse array data from a file.
    
    This function is using :func:`numpy.loadtxt` to parse the file.  Each row 
    in the text file must have the same number of values.
    
    :arg filename: File or filename to read. If the filename extension is 
        :file:`.gz` or :file:`.bz2`, the file is first decompressed.
    :type filename: str or file
    
    :arg delimiter: The string used to separate values. By default, 
        this is any whitespace.
    :type delimiter: str
    
    :arg skiprows: Skip the first *skiprows* lines, default is ``0``.
    :type skiprows: int
     
    :arg usecols: Which columns to read, with 0 being the first. For example, 
        ``usecols = (1,4,5)`` will extract the 2nd, 5th and 6th columns. 
        The default, ``None``, results in all columns being read.
    :type usecols: list
    
    :arg dtype: Data-type of the resulting array, default is :func:`float`. 
    :type dtype: :class:`numpy.dtype`."""

    array = np.loadtxt(filename, dtype=dtype, delimiter=delimiter, 
                       skiprows=skiprows, usecols=usecols)
    return array
        
def parseSparseMatrix(filename, symmetric=False, delimiter=None, skiprows=0,
                      irow=0, icol=1, first=1):
    """Parse sparse matrix data from a file.
    
    This function is using :func:`parseArray` to parse the file.
    Input must have the following format::
        
       1       1    9.958948135375977e+00
       1       2   -3.788214445114136e+00
       1       3    6.236155629158020e-01
       1       4   -7.820609807968140e-01
    
    Each row in the text file must have the same number of values.
    
    :arg filename: File or filename to read. If the filename extension is 
        :file:`.gz` or :file:`.bz2`, the file is first decompressed.
    :type filename: str or file
    
    :arg symmetric: Set ``True`` if the file contains triangular part of a 
        symmetric matrix, default is ``False``.
    :type symmetric: bool
    
    :arg delimiter: The string used to separate values. By default, 
        this is any whitespace.
    :type delimiter: str
    
    :arg skiprows: Skip the first *skiprows* lines, default is ``0``.
    :type skiprows: int
    
    :arg irow: Index of the column in data file corresponding to row indices,
        default is ``0``. 
    :type irow: int 
        
    :arg icol: Index of the column in data file corresponding to row indices,
        default is ``0``. 
    :type icol: int
    
    :arg first: First index in the data file (0 or 1), default is ``1``. 
    :type first: int

    Data-type of the resulting array, default is :func:`float`."""
    
    irow = int(irow)
    icol = int(icol)
    first = int(first)
    assert 0 <= irow <= 2 and 0 <= icol <= 2, 'irow/icol may be 0, 1, or 2'
    assert icol != irow, 'irow and icol must not be equal' 
    idata = [0, 1, 2]
    idata.pop(idata.index(irow))
    idata.pop(idata.index(icol))
    idata = idata[0]
    sparse = parseArray(filename, delimiter, skiprows)
    dof = sparse[:,[irow, icol]].max() 
    matrix = np.zeros((dof,dof))
    irow = (sparse[:,irow] - first).astype(int)
    icol = (sparse[:,icol] - first).astype(int)
    matrix[irow, icol] = sparse[:,idata]
    if symmetric:
        matrix[icol, irow] = sparse[:,idata]
    return matrix
