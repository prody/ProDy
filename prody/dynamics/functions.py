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

import os.path

import numpy as np

from prody.atomic import AtomGroup
from prody.tools import openFile, isExecutable, which, PLATFORM

from nma import NMA
from anm import ANM
from gnm import GNM, GNMBase, ZERO
from pca import PCA, EDA
from mode import Vector, Mode
from modeset import ModeSet

__all__ = ['parseArray', 'parseModes', 'parseNMD',
           'writeArray', 'writeModes', 'parseSparseMatrix',
           'writeNMD', 
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           'getVMDpath', 'setVMDpath', 'viewNMDinVMD',]

pkg = __import__(__package__)
LOGGER = pkg.LOGGER
SETTINGS = pkg.SETTINGS
           
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
        raise TypeError('invalid type for nma, {0:s}'.format(type(nma)))
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
        raise IOError('{0:s} is not a valid NMA model file'.format(filename))
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
        raise TypeError('invalid type for vector, {0:s}'.format(type(vector)))
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
    """Return :class:`~.Vector` instance after loading it from file 
    (*filename*). This function makes use of :func:`numpy.load` function.  
    See also :func:`saveVector`."""
    
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
            import _winreg
            for vmdversion in ('1.8.7', '1.9', '1.9.1'): 
                try:
                    key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                            'Software\\University of Illinois\\VMD\\' + 
                            vmdversion)
                    vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                    vmdbin = os.path.join(vmddir, 'vmd.exe') 
                except:    
                    pass
                try:
                    key = _winreg.OpenKey(_winreg.HKEY_LOCAL_MACHINE, 
                'Software\\WOW6432node\\University of Illinois\\VMD\\' + 
                vmdversion)
                    vmddir = _winreg.QueryValueEx(key, 'VMDDIR')[0]
                    vmdbin = os.path.join(vmddir, 'vmd.exe') 
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
        if False and \
           isinstance(vmdbin, (StringType, UnicodeType)) and \
           isinstance(vmddir, (StringType, UnicodeType)) and \
           os.path.isfile(vmdbin) and os.path.isdir(vmddir): 
            pass#return vmdbin, vmddir
        if isExecutable(vmdbin):
            setVMDpath(vmdbin)
            return vmdbin
        

def setVMDpath(path):
    """Set path to a VMD executable."""
    
    if isExecutable(path):
        SETTINGS['vmd'] = path
        SETTINGS.save()
        LOGGER.info("VMD path is set to '{0:s}'.".format(path))
    else:
        raise OSError('{0:s} is not executable.'.format(str(path)))

def parseNMD(filename, type=NMA):
    """Returns normal mode and atomic data parsed from an NMD file.
    Normal mode data is returned in an :class:`~.NMA` instance. Atomic
    data is returned in an :class:`~.AtomGroup` instance."""

    assert not isinstance(type, NMA), 'type must be NMA, ANM, GNM, or PCA'
    atomic = dict()
    modes = []
    nmd = open(filename)
    for line in nmd:
        split = line.find(' ')
        if line[:split] == 'mode':
            modes.append(line[split:].strip())
        elif line[:split] in ('coordinates', 'atomnames', 'resnames', 
                              'resnums', 'resids', 'chainids', 'bfactors',
                              'name'):
            atomic[line[:split]] = line[split:].strip()
    nmd.close()
    
    name = atomic.pop('name', os.path.splitext(os.path.split(filename)[1])[0])
    coords = atomic.pop('coordinates', None)
    dof = None
    if coords is not None:
        coords = np.fromstring( coords, dtype=float, sep=' ')
        dof = coords.shape[0]
        ag = None
        n_atoms = dof / 3
        coords = coords.reshape((n_atoms, 3))
        ag = AtomGroup(name)
        ag.setCoords(coords)
        data = atomic.pop('atomnames', None)
        if data is not None:
            ag.setNames(data.split())
        data = atomic.pop('resnames', None)
        if data is not None:
            ag.setResnames(data.split())
        data = atomic.pop('chainids', None)
        if data is not None:
            ag.setChids(data.split())
        data = atomic.pop('resnums', None)
        if data is not None:
            ag.setResnums(np.fromstring(data, int, sep=' '))
        data = atomic.pop('resids', None)
        if data is not None:
            ag.setResnums(np.fromstring(data, int, sep=' '))
        data = atomic.pop('bfactors', None)
        if data is not None:
            ag.setBetas(np.fromstring(data, float, sep=' '))
    nma = type(name)
    for mode in modes:
        
        items = mode.split()
        diff = len(items) - dof
        mode = np.array(items[diff:]).astype(float)
        if len(mode) != dof:
            pass
        if diff == 1 and not items[0].isdigit():
            value = float(items[0])
        else:
            if not items[0].isdigit():
                value = float(items[0])
            elif not items[1].isdigit():
                value = float(items[1])
            else:
                value = 1.0
        nma.addEigenpair(mode, value)
    return nma, ag
    

def writeNMD(filename, modes, atoms):
    """Writes an NMD file for given *modes* and includes applicable data from 
    *atoms*.  Returns *filename*, if file is successfully written.  NMD file 
    format is described at :ref:`nmd-format`.
    
    .. note:: 
       #. This function skips modes with zero eigenvalues.
       #. If a :class:`~.Vector` instance is given, it will be normalized 
          before it is written. It's length before normalization will be 
          written as the scaling factor of the vector."""
    
    if not isinstance(modes, (NMA, ModeSet, Mode, Vector)):
        raise TypeError('modes must be NMA, ModeSet, Mode, or Vector, '
                        'not {0:s}'.format(type(modes)))
    if modes.numAtoms() != atoms.numAtoms():
        raise Exception('number of atoms do not match')
    out = openFile(filename, 'w')
    
    #out.write('#!{0:s} -e\n'.format(VMDPATH))
    out.write('nmwiz_load {0:s}\n'.format(os.path.abspath(filename)))
    name = modes.getTitle()
    name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = str(atoms)
        name = name.replace(' ', '_').replace('.', '_')
    if not name.replace('_', '').isalnum() or len(name) > 30:
        name = os.path.splitext(os.path.split(filename)[1])[0]
    out.write('name {0:s}\n'.format(name))
    try:
        coords = atoms.getCoords()
    except:
        raise ValueError('coordinates could not be retrieved from atoms')
    if coords is None:
        raise ValueError('atom coordinates are not set')
    
    try:
        data = atoms.getNames()
        if data is not None:
            out.write('atomnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnames()
        if data is not None:
            out.write('resnames {0:s}\n'.format(' '.join(data)))
    except:
        pass
    try:
        data = atoms.getResnums()
        if data is not None:
            out.write('resids {0:s}\n'.format(' '.join(data.astype('|S5'))))
    except:
        pass
    try:
        data = atoms.getChids()
        if data is not None:
            out.write('chainids {0:s}\n'.format(' '.join(data)))
    except:
        pass
    
    try:
        data = atoms.getBetas()
        if data is not None:
            out.write('bfactors {0:s}\n'.format(' '.join(
                            ['{0:.3f}'.format(x) for x in data.flatten()])))
    except:
        pass
    
    out.write('coordinates {0:s}\n'.format(
                    ' '.join(['{0:.3f}'.format(x) for x in coords.flatten()])))
    
    count = 0
    if isinstance(modes, Vector):
        out.write('mode 1 {0:.2f} {1:s}\n'.format(abs(modes), ' '.join(
                ['{0:.3f}'.format(x) for x in modes.getNormed()._getArray()])))
        count += 1
    else:
        if isinstance(modes, Mode):
            modes = [modes]
        for mode in modes:
            if mode.getEigenvalue() < ZERO:
                continue
            out.write('mode {0:d} {1:.2f} {2:s}\n'.format(
                       mode.getIndex()+1, mode.getVariance()**0.5, 
                       ' '.join(
                            ['{0:.3f}'.format(x) for x in mode._getArray()])))
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
        os.system('{0:s} -e {1:s}'.format(vmd, os.path.abspath(filename)))
        
def writeModes(filename, modes, format='%.18e', delimiter=' '):
    """Write *modes* (eigenvectors) into a plain text file with name 
    *filename*. See also :func:`writeArray`."""
    
    if not isinstance(modes, (NMA, ModeSet, Mode)):
        raise TypeError('modes must be NMA, ModeSet, or Mode, not {0:s}'
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
    nma = NMA(os.path.splitext(os.path.split(normalmodes)[1])[0])
    nma.setEigens(modes, values)
    return nma
    
    
def writeArray(filename, array, format='%d', delimiter=' '):
    """Write 1-d or 2-d array data into a delimited text file.
    
    This function is using :func:`numpy.savetxt` to write the file, after 
    making some type and value checks.  Default *format* argument is ``"%d"``.
    Default *delimiter* argument is white space, ``" "``.
    
    *filename* will be returned upon successful writing."""
    
    if not isinstance(array, np.ndarray):
        raise TypeError('array must be a Numpy ndarray, not {0:s}'
                        .format(type(array)))
    elif not array.ndim in (1, 2):
        raise ValueError('array must be a 1 or 2-dimensional Numpy ndarray, '
                         'not {0:d}-d'.format(type(array.ndim)))
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
