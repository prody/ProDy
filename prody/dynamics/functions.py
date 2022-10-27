# -*- coding: utf-8 -*-
"""This module defines input and output functions."""

from collections import OrderedDict
import datetime

import os
from os.path import abspath, join, isfile, isdir, split, splitext

import numpy as np

from prody import LOGGER, SETTINGS, PY3K
from prody.atomic import Atomic, AtomSubset
from prody.utilities import openFile, openSQLite, isExecutable, which, PLATFORM, addext, wrapModes
from prody.proteins import parseSTAR, writeSTAR, StarDict, alignChains, parsePDB
from prody.ensemble import PDBEnsemble

from .nma import NMA, MaskedNMA
from .anm import ANM, ANMBase, MaskedANM
from .analysis import calcCollectivity, calcScipionScore
from .analysis import calcProjection
from .analysis import calcCollectivity
from .gnm import GNM, GNMBase, ZERO, MaskedGNM
from .exanm import exANM, MaskedExANM
from .rtb import RTB
from .pca import PCA, EDA
from .imanm import imANM
from .exanm import exANM
from .mode import Vector, Mode, VectorBase
from .modeset import ModeSet
from .editing import sliceModel, reduceModel, trimModel
from .editing import sliceModelByMask, reduceModelByMask, trimModelByMask

__all__ = ['parseArray', 'parseModes', 'parseSparseMatrix',
           'parseGromacsModes', 'parseScipionModes',
           'writeArray', 'writeModes', 'writeScipionModes',
           'saveModel', 'loadModel', 'saveVector', 'loadVector',
           'calcENM', 'realignModes']


def saveModel(nma, filename=None, matrices=False, **kwargs):
    """Save *nma* model data as :file:`filename.nma.npz`.  By default,
    eigenvalues, eigenvectors, variances, trace of covariance matrix,
    and name of the model will be saved.  If *matrices* is **True**,
    covariance, Hessian or Kirchhoff matrices are saved too, whichever
    are available.  If *filename* is **None**, name of the NMA instance
    will be used as the filename, after ``" "`` (white spaces) in the name
    are replaced with ``"_"`` (underscores).  Extension may differ based
    on the type of the NMA model.  For ANM models, it is :file:`.anm.npz`.
    Upon successful completion of saving, filename is returned. This
    function makes use of :func:`~numpy.savez` function."""

    if not isinstance(nma, NMA):
        raise TypeError('invalid type for nma, {0}'.format(type(nma)))
    #if len(nma) == 0:
    #    raise ValueError('nma instance does not contain data')

    add_attr = kwargs.pop('attr', [])

    dict_ = nma.__dict__
    attr_list = ['_title', '_trace', '_array', '_eigvals', '_vars', '_n_atoms',
                 '_dof', '_n_modes']

    if add_attr:
        for attr in add_attr:
            if attr not in attr_list:
                attr_list.append(attr)
    if filename is None:
        filename = nma.getTitle().replace(' ', '_')
    if isinstance(nma, GNMBase):
        attr_list.append('_cutoff')
        attr_list.append('_gamma')
        if matrices:
            attr_list.append('_kirchhoff')
            if isinstance(nma, ANMBase):
                attr_list.append('_hessian')
        if isinstance(nma, ANMBase):
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

    if isinstance(nma, MaskedNMA):
        if isinstance(nma, MaskedGNM):
            attr_dict['type'] = 'mGNM'
        elif isinstance(nma, MaskedANM):
            attr_dict['type'] = 'mANM'
        else:
            raise TypeError('invalid MaskedNMA type: %s'%(str(type(nma))))

        attr_dict['mask'] = nma.mask
        attr_dict['masked'] = nma.masked
    
    if isinstance(nma, RTB):
        attr_dict['type'] = 'RTB'
        if matrices:
            attr_dict['_project'] = nma._project

    if isinstance(nma, imANM):
        attr_dict['type'] = 'imANM'

    if isinstance(nma, exANM):
        attr_dict['type'] = 'exANM'

    suffix = '.' + attr_dict['type'].lower()
    if not filename.lower().endswith('.npz'):
        if not filename.lower().endswith(suffix):
            filename += suffix + '.npz'
        else:
            filename += '.npz'
    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadModel(filename, **kwargs):
    """Returns NMA instance after loading it from file (*filename*).
    This function makes use of :func:`~numpy.load` function.  See
    also :func:`saveModel`."""

    if not 'encoding' in kwargs:
        kwargs['encoding'] = 'latin1'

    if not 'allow_pickle' in kwargs:
        kwargs['allow_pickle'] = True 

    with np.load(filename, **kwargs) as attr_dict:
        try:
            type_ = attr_dict['type']
        except KeyError:
            raise IOError('{0} is not a valid NMA model file'.format(filename))

        if isinstance(type_, np.ndarray):
            type_ = np.asarray(type_, dtype=str)

        type_ = str(type_)

        try:
            title = attr_dict['_title']
        except KeyError:
            title = attr_dict['_name']

        if isinstance(title, np.ndarray):
            title = np.asarray(title, dtype=str)
        title = str(title)
        if type_ == 'ANM':
            nma = ANM(title)
        elif type_ == 'PCA':
            nma = PCA(title)
        elif type_ == 'EDA':
            nma = EDA(title)
        elif type_ == 'GNM':
            nma = GNM(title)
        elif type_ == 'mGNM':
            nma = MaskedGNM(title)
        elif type_ == 'mANM':
            nma = MaskedANM(title)
        elif type_ == 'exANM':
            nma = exANM(title)
        elif type_ == 'imANM':
            nma = imANM(title)
        elif type_ == 'NMA':
            nma = NMA(title)
        elif type_ == 'RTB':
            nma = RTB(title)
        else:
            raise IOError('NMA model type is not recognized: {0}'.format(type_))

        dict_ = nma.__dict__
        for attr in attr_dict.files:
            if attr in ('type', '_name', '_title'):
                continue
            elif attr in ('_trace', '_cutoff', '_gamma'):
                dict_[attr] = attr_dict[attr][()]
            elif attr in ('_dof', '_n_atoms', '_n_modes'):
                dict_[attr] = int(attr_dict[attr])
            elif attr in ('masked', ):
                dict_[attr] = bool(attr_dict[attr])
            elif attr in ('mask', ):
                if not attr_dict[attr].shape:
                    dict_[attr] = bool(attr_dict[attr])
                else:
                    dict_[attr] = attr_dict[attr]
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

    if not filename.lower().endswith('.npz'):
        if not filename.lower().endswith('.vec'):
            filename += '.vec.npz'
        else:
            filename += '.npz'

    ostream = openFile(filename, 'wb', **kwargs)
    np.savez(ostream, **attr_dict)
    ostream.close()
    return filename


def loadVector(filename):
    """Returns :class:`.Vector` instance after loading it from *filename* using
    :func:`numpy.load`.  See also :func:`saveVector`."""

    attr_dict = np.load(filename)
    try:
        title = str(attr_dict['title'])
    except KeyError:
        title = str(attr_dict['name'])
    return Vector(attr_dict['array'], title, bool(attr_dict['is3d']))


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
    """Returns :class:`.NMA` instance with normal modes parsed from
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
        6th columns. The default, **None**, results in all columns being read.
    :type nm_usecols: list

    :arg ev_delimiter: The string used to separate values in *eigenvalues*.
        By default, this is any whitespace.
    :type ev_delimiter: str

    :arg ev_skiprows: Skip the first *skiprows* lines in *eigenvalues*.
        Default is ``0``.
    :type ev_skiprows: 0

    :arg ev_usecols: Which columns to read from *eigenvalues*, with 0 being the
        first. For example, ``usecols = (1,4,5)`` will extract the 2nd, 5th and
        6th columns. The default, **None**, results in all columns being read.
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


def parseScipionModes(metadata_file, title=None, pdb=None):
    """Returns :class:`.NMA` containing eigenvectors and eigenvalues 
    parsed from a ContinuousFlex FlexProtNMA Run directory.

    :arg run_path: path to the Run directory
    :type run_path: str
    
    :arg title: title for :class:`.NMA` object
    :type title: str
    """
    run_path = os.path.split(metadata_file)[0]
    top_dirs = os.path.split(run_path)[0][:-4]
    run_name = os.path.split(run_path)[-1]

    if metadata_file.endswith('.xmd'):
        star_data = parseSTAR(metadata_file)
        
    elif metadata_file.endswith('.sqlite'):
        # reconstruct star data from sqlite

        sql_con = openSQLite(metadata_file)
        cursor = sql_con.cursor()

        star_dict = OrderedDict()
        star_block_dict = OrderedDict()
        
        star_loop_dict = OrderedDict()
        star_loop_dict["fields"] = OrderedDict([(0, '_enabled'),
                                            (1, '_nmaCollectivity'),
                                            (2, '_nmaModefile'),
                                            (3, '_nmaScore'),
                                            (4, '_nmaEigenval'),
                                            (5, '_order_')])
        star_loop_dict["data"] = OrderedDict()

        for row in cursor.execute("SELECT * FROM Objects;"):
        
            id_ = row[0]
            key = id_ - 1 # sqlite ids start from 1 not 0

            star_loop_dict["data"][key] = OrderedDict()

            star_loop_dict["data"][key]['_order_'] = id_

            star_loop_dict["data"][key]['_enabled'] = row[1]

            star_loop_dict["data"][key]['_nmaCollectivity'] = row[6]

            star_loop_dict["data"][key]['_nmaModefile'] = row[5]

            star_loop_dict["data"][key]['_nmaScore'] = row[7]

            if len(row) > 8:
                star_loop_dict["data"][key]['_nmaEigenval'] = row[8]

        star_block_dict[0] = star_loop_dict
        star_dict[0] = star_block_dict

        star_data = StarDict(star_dict, prog='XMIPP')

    else:
        raise ValueError("Metadata file should be an xmd or sqlite file")

    star_loop = star_data[0][0]
    
    n_modes = star_loop.numRows()
    
    row1 = star_loop[0]
    mode1 = parseArray(top_dirs + row1['_nmaModefile']).reshape(-1)
    dof = mode1.shape[0]

    if pdb is not None:
        atoms = parsePDB(pdb)
        n_atoms = atoms.numAtoms()
    else:
        # assume standard NMA
        n_atoms = dof//3

    vectors = np.zeros((dof, n_modes))
    vectors[:, 0] = mode1

    eigvals = np.zeros(n_modes)

    try:
        eigvals[0] = float(row1['_nmaEigenval'])
        found_eigvals = True
    except:
        found_eigvals = False

    for i, row in enumerate(star_loop[1:]):
        vectors[:, i+1] = parseArray(top_dirs + row['_nmaModefile']).reshape(-1)
        if found_eigvals:
            eigvals[i+1] = float(row['_nmaEigenval'])
    
    if not found_eigvals:
        log_fname = run_path + '/logs/run.stdout'
        fi = open(log_fname, 'r')
        lines = fi.readlines()
        fi.close()

        for line in lines:
            if line.find('Eigenvector number') != -1:
                j = int(line.strip().split()[-1]) - 1
            if line.find('Corresponding eigenvalue') != -1:
                eigvals[j] = float(line.strip().split()[-1])
                if not found_eigvals:
                    found_eigvals = True
        
    if title is None:
        title = run_name

    if not found_eigvals:
        LOGGER.warn('No eigenvalues found')
        eigvals=None

    if dof == n_atoms * 3:
        nma = NMA(title)
    else:
        nma = GNM(title)

    nma.setEigens(vectors, eigvals)
    return nma


def writeScipionModes(output_path, modes, write_star=False, scores=None,
                      only_sqlite=False, collectivityThreshold=0.):
    """Writes *modes* to a set of files that can be recognised by Scipion.
    A directory called **"modes"** will be created if it doesn't already exist. 
    Filenames inside will start with **"vec"** and have the mode number as the extension.
    
    :arg output_path: path to the directory where the modes directory will be
    :type output_path: str

    :arg modes: modes to be written to files
    :type modes: :class:`.Mode`, :class:`.ModeSet`, :class:`.NMA`

    :arg write_star: whether to write modes.xmd STAR file.
        Default is **False** as qualifyModesStep writes it with scores.
    :type write_star: bool

    :arg scores: scores from qualifyModesStep for re-writing sqlite
        Default is **None** and then it uses :func:`.calcScipionScore`
    :type scores: list

    :arg only_sqlite: whether to write only the sqlite file instead of everything.
        Default is **False** but it can be useful to set it to **True** for updating the sqlite file.
    :type only_sqlite: bool    

    :arg collectivityThreshold: collectivity threshold below which modes are not enabled
        Default is 0.
    :type collectivityThreshold: float
    """
    if not isinstance(output_path, str):
        raise TypeError('output_path should be a string, not {0}'
                        .format(type(output_path)))

    if not isdir(output_path):
        raise ValueError('output_path should be a working path')

    if not isinstance(modes, (NMA, ModeSet, VectorBase)):
        raise TypeError('rows must be NMA, ModeSet, or Mode, not {0}'
                        .format(type(modes)))

    if not isinstance(write_star, bool):
        raise TypeError('write_star should be boolean, not {0}'
                        .format(type(write_star)))

    if scores is not None:
        if not isinstance(scores, list):
            raise TypeError('scores should be a list or None, not {0}'
                            .format(type(scores)))

    if not isinstance(only_sqlite, bool):
        raise TypeError('only_sqlite should be boolean, not {0}'
                        .format(type(only_sqlite)))

    if not isinstance(collectivityThreshold, float):
        raise TypeError('collectivityThreshold should be float, not {0}'
                        .format(type(collectivityThreshold)))

    if modes.numModes() == 1 and not isinstance(modes, NMA):
        old_modes = modes
        modes = NMA(old_modes)
        modes.setEigens(old_modes.getArray().reshape(-1, 1))

    modes_dir = output_path + '/modes/'
    if not isdir(modes_dir):
        os.mkdir(modes_dir)

    modefiles = []
    for mode in modes:
        mode_num = mode.getIndex() + 1
        if mode.is3d():
            modefiles.append(writeArray(modes_dir + 'vec.{0}'.format(mode_num),
                                        mode.getArrayNx3(), '%12.4e', ''))
        else:
            modefiles.append(writeArray(modes_dir + 'vec.{0}'.format(mode_num),
                                        mode.getArray(), '%12.4e', ''))

    if modes.numModes() > 1:
        order = modes.getIndices()
        collectivities = list(calcCollectivity(modes))
        eigvals = modes.getEigvals()
        enabled = [1 if eigval > ZERO and collectivities[i] > collectivityThreshold else -1
                   for i, eigval in enumerate(eigvals)]
        if scores is None:
            scores = list(calcScipionScore(modes))
    else:
        mode = modes[0]
        eigvals = np.array([mode.getEigval()])
        collectivities = [calcCollectivity(mode)]
        order = [mode.getIndex()]
        enabled = [1 if mode.getEigval() > ZERO and collectivities[0] > collectivityThreshold else -1]
        if scores is None:
            scores = [calcScipionScore(mode)[0]]

    modes_sqlite_fn = output_path + '/modes.sqlite'
    sql_con = openSQLite(modes_sqlite_fn, 'n')
    cursor = sql_con.cursor()
    
    cursor.execute('''CREATE TABLE Properties(key,value)''')
    properties = [('self', 'SetOfNormalModes'),
                  ('_size', str(modes.numModes())),
                  ('_streamState', '2'),
                  ('_mapperPath', '{0}, '.format(modes_sqlite_fn))]
    cursor.executemany('''INSERT INTO Properties VALUES(?,?);''', properties);
    
    cursor.execute('''CREATE TABLE Classes(id primary key, label_property, column_name, class_name)''')
    classes = [(1, 'self', 'c00', 'NormalMode'),
               (2, '_modeFile', 'c01', 'String'),
               (3, '_collectivity', 'c02', 'Float'),
               (4, '_score', 'c03', 'Float'),
               (5, '_eigenval', 'c04', 'Float')]
    cursor.executemany('''INSERT INTO Classes VALUES(?,?,?,?);''', classes);
    
    cursor.execute('''CREATE TABLE Objects(id primary key, enabled, label, comment, creation, c01, c02, c03, c04)''')
    
    star_dict = OrderedDict()

    star_dict['noname'] = OrderedDict() # Data Block with title noname
    loop_dict = star_dict['noname'][0] = OrderedDict() # Loop 0

    loop_dict['fields'] = OrderedDict()
    fields = ['_enabled', '_nmaCollectivity', '_nmaModefile', '_nmaScore',
              '_nmaEigenval', '_order_']
    for j, field in enumerate(fields):
        loop_dict['fields'][j] = field

    loop_dict['data'] = OrderedDict()

    now = datetime.datetime.now()    
    creation = now.strftime("%Y-%m-%d %H:%M:%S")
    for i, mode in enumerate(modes):
        loop_dict['data'][i] = OrderedDict()
        
        id = mode.getIndex() + 1
        loop_dict['data'][i]['_order_'] = str(id)
        
        enab = enabled[i]
        loop_dict['data'][i]['_enabled'] = '%2i' % enab
        if enab != 1:
            enab = 0
            
        label = ''
        comment = ''
        
        c01 = loop_dict['data'][i]['_nmaModefile'] = modefiles[i]
        
        collec = collectivities[i]
        loop_dict['data'][i]['_nmaCollectivity'] = '%8.6f' % collec
        c02 = float('%6.4f' % collec)
        
        c03 = scores[i]
        loop_dict['data'][i]['_nmaScore'] = '%8.6f' % c03

        c04 = eigvals[i]
        if float('%9.6f' % c04) > 0:
            loop_dict['data'][i]['_nmaEigenval'] = '%9.6f' % c04
        else:
            loop_dict['data'][i]['_nmaEigenval'] = '%9.6e' % c04
        
        cursor.execute('''INSERT INTO Objects VALUES(?,?,?,?,?,?,?,?,?)''',
                       (id, enab, label, comment, creation, c01, c02, c03, c04))

    if write_star:
        writeSTAR(output_path + '/modes.xmd', star_dict)

    sql_con.commit()
    sql_con.close()

    return modes_dir


def writeArray(filename, array, format='%3.2f', delimiter=' '):
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
        The default, **None**, results in all columns being read.
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

    :arg symmetric: Set **True** if the file contains triangular part of a
        symmetric matrix, default is **True**.
    :type symmetric: bool

    :arg delimiter: The string used to separate values. By default,
        this is any whitespace.
    :type delimiter: str

    :arg skiprows: Skip the first *skiprows* lines, default is ``0``.
    :type skiprows: int

    :arg irow: Index of the column in data file corresponding to row indices,
        default is ``0``.
    :type irow: int

    :arg icol: Index of the column in data file corresponding to column indices,
        default is ``1``.
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
    if symmetric:
        dim1 = dim2 = int(sparse[:, [irow, icol]].max())
    else:
        dim1, dim2 = sparse[:, [irow, icol]].max(0).astype(int)
    matrix = np.zeros((dim1, dim2))
    irow = (sparse[:, irow] - first).astype(int)
    icol = (sparse[:, icol] - first).astype(int)
    matrix[irow, icol] = sparse[:, idata]
    if symmetric:
        matrix[icol, irow] = sparse[:, idata]
    return matrix

def calcENM(atoms, select=None, model='anm', trim='trim', gamma=1.0, 
            title=None, n_modes=None, **kwargs):
    """Returns an :class:`.ANM` or :class:`.GNM` instance and *atoms* used for the 
    calculations. The model can be trimmed, sliced, or reduced based on 
    the selection.

    :arg atoms: atoms on which the ENM is performed. It can be any :class:`Atomic` 
        class that supports selection or a :class:`~numpy.ndarray`.
    :type atoms: :class:`.Atomic`, :class:`.AtomGroup`, :class:`.Selection`, :class:`~numpy.ndarray`

    :arg select: part of the atoms that is considered as the system. 
        If set to **None**, then all atoms will be considered as the system
    :type select: str, :class:`.Selection`, :class:`~numpy.ndarray`

    :arg model: type of ENM that will be performed. It can be either ``"anm"`` 
        or ``"gnm"`` or ``"exanm"``
    :type model: str

    :arg trim: type of method that will be used to trim the model. It can 
        be either ``"trim"`` , ``"slice"``, or ``"reduce"``. If set to ``"trim"``, 
        the parts that is not in the selection will simply be removed
    :type trim: str
    """
    
    if isinstance(select, (str, AtomSubset)):
        if not isinstance(atoms, Atomic):
            raise TypeError('atoms must be a Atomic instance in order to be selected')
    try:
        if title is None:
            title = atoms.getTitle()
    except AttributeError:
        title = 'Unknown'

    mask = kwargs.pop('mask', None)
    zeros = kwargs.pop('zeros', False)
    turbo = kwargs.pop('turbo', True)

    if model is GNM:
        model = 'gnm'
    elif model is ANM:
        model = 'anm'
    else:
        model = str(model).lower().strip() 

    if trim is reduceModel:
        trim = 'reduce'
    elif trim is sliceModel:
        trim = 'slice'
    elif trim is None:
        trim = 'trim'
    else:
        trim = str(trim).lower().strip()
    
    enm = None
    MaskedModel = None
    if model == 'anm':
        anm = ANM(title)
        anm.buildHessian(atoms, gamma=gamma, **kwargs)
        enm = anm
        MaskedModel = MaskedANM
    elif model == 'gnm':
        gnm = GNM(title)
        gnm.buildKirchhoff(atoms, gamma=gamma, **kwargs)
        enm = gnm
        MaskedModel = MaskedGNM
    elif model == 'exanm':
        exanm = exANM(title)
        exanm.buildHessian(atoms, gamma=gamma, **kwargs)
        enm = exanm
        MaskedModel = MaskedExANM
    else:
        raise TypeError('model should be either ANM or GNM instead of {0}'.format(model))
    
    if select is None:
        enm.calcModes(n_modes=n_modes, zeros=zeros, turbo=turbo)
    else:
        if trim == 'slice':
            enm.calcModes(n_modes=n_modes, zeros=zeros, turbo=turbo)
            if isinstance(select, np.ndarray):
                enm = sliceModelByMask(enm, select)
                atoms = select
            else:
                enm, atoms = sliceModel(enm, atoms, select)
        elif trim == 'reduce':
            if isinstance(select, np.ndarray):
                enm = reduceModelByMask(enm, select)
                atoms = select
            else:
                enm, atoms = reduceModel(enm, atoms, select)
            enm.calcModes(n_modes=n_modes, zeros=zeros, turbo=turbo)
        elif trim == 'trim':
            if isinstance(select, np.ndarray):
                enm = trimModelByMask(enm, select)
                atoms = select
            else:
                enm, atoms = trimModel(enm, atoms, select)
            enm.calcModes(n_modes=n_modes, zeros=zeros, turbo=turbo)
        else:
            raise ValueError('trim can only be "trim", "reduce", or "slice"')
    
    if mask is not None:
        enm = MaskedModel(enm, mask)
    return enm, atoms


def parseGromacsModes(run_path, title="", model='nma', **kwargs):
    """Returns :class:`.NMA` containing eigenvectors and eigenvalues parsed from a run directory 
    containing results from gmx covar or gmx nmeig followed by gmx anaeig 
    including eigenvalues in an xvg file and eigenvectors in pdb files
    (see http://www.strodel.info/index_files/lecture/html/analysis-9.html).

    :arg run_path: path to the run directory
    :type run_path: str
    
    :arg title: title for resulting object
        Default is ``""``
    :type title: str

    :arg model: type of calculated that was performed. It can be either ``"nma"`` 
        or ``"pca"``. If it is not changed to ``"pca"`` then ``"nma"`` will be assumed.
    :type model: str

    :arg eigval_fname: filename or path for xvg file containing eigenvalues
        Default is ``"eigenval.xvg"`` as this is the default from Gromacs
    :type eigval_fname: str

    :arg eigvec_fname: filename or path for trr file containing eigenvectors
        Default is ``"eigenvec.trr"`` as this is the default from Gromacs
    :type eigvec_fname: str

    :arg pdb_fname: filename or path for pdb file containing the reference structure
        Default is ``"average.pdb"`` although this is probably suboptimal
    :type pdb_fname: str
    """ 
    try:
        from mdtraj import load_trr
    except ImportError:
        raise ImportError('Please install mdtraj in order to use parseGromacsModes.')

    if not isinstance(run_path, str):
        raise TypeError('run_path should be a string')

    if not run_path.endswith('/'):
        run_path += '/'

    if not isinstance(title, str):
        raise TypeError('title should be a string')

    if model == 'pca':
        result = PCA(title)
    else:
        if model != 'nma':
            LOGGER.warn('model not recognised so using NMA')
        result = NMA(title)


    eigval_fname = kwargs.get('eigval_fname', 'eigenval.xvg')
    if not isinstance(eigval_fname, str):
        raise TypeError('eigval_fname should be a string')

    if isfile(eigval_fname):
        vals_fname = eigval_fname
    elif isfile(run_path + eigval_fname):
        vals_fname = run_path + eigval_fname
    else:
        raise ValueError('eigval_fname should point be a path to a file '
                         'either relative to run_path or an absolute one')


    eigvec_fname = kwargs.get('eigvec_fname', 'eigenvec.trr')
    if not isinstance(eigvec_fname, str):
        raise TypeError('eigvec_fname should be a string')

    if isfile(eigvec_fname):
        vecs_fname = eigval_fname
    elif isfile(run_path + eigvec_fname):
        vecs_fname = run_path + eigvec_fname
    else:
        raise ValueError('eigvec_fname should point be a path to a file '
                         'either relative to run_path or an absolute one')


    pdb_fname = kwargs.get('pdb_fname', 'average.pdb')
    if not isinstance(pdb_fname, str):
        raise TypeError('pdb_fname should be a string')

    if isfile(pdb_fname):
        pdb = eigval_fname
    elif isfile(run_path + pdb_fname):
        pdb = run_path + pdb_fname
    else:
        raise ValueError('pdb_fname should point be a path to a file '
                         'either relative to run_path or an absolute one')
    
    
    fi = open(vals_fname, 'r')
    lines = fi.readlines()
    fi.close()
    
    eigvals = []
    for line in lines:
        if not (line.startswith('@') or line.startswith('#')):
            eigvals.append(float(line.strip().split()[-1])*100) # convert to A**2 from nm**2

    eigvals = np.array(eigvals)

    # Parse eigenvectors trr with mdtraj, which uses nm so doesn't rescale
    vecs_traj = load_trr(vecs_fname, top=pdb)

    # format vectors appropriately, skipping initial and average structures
    vectors = np.array([frame.xyz.flatten() for frame in vecs_traj[2:]]).T

    result.setEigens(vectors, eigvals)
    return result

def realignModes(modes, atoms, ref):
    """Align *modes* in the original frame based on *atoms*
    onto another frame based on *ref* using the transformation 
    from alignment of *atoms* to *ref*
    
    :arg modes: multiple 3D modes
    :type modes: :class:`.ModeSet`, :class:`.ANM`, :class:`.PCA`

    :arg atoms: central structure related to *modes* to map onto *ref*
        Inserting *atoms* into an ensemble and projecting onto *modes*
        should give all zeros
    :type atoms: :class:`.Atomic`
    
    :arg ref: reference structure for mapping
    :type ref: :class:`.Atomic`
    """
    if not isinstance(modes, (ModeSet, NMA)):
        raise TypeError('modes should be a ModeSet of NMA instance')

    if not modes.is3d():
        raise ValueError('modes should be 3D for this function to work')

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms should be an Atomic instance')

    if not isinstance(ref, Atomic):
        raise TypeError('ref should be an Atomic instance')

    n_atoms = modes.numAtoms()

    if atoms.numAtoms() != n_atoms:
        raise ValueError('atoms and modes should have the same number of atoms')

    def_coords = np.array([atoms.getCoords() + mode.getArrayNx3()
                           for mode in modes])

    def_ens = PDBEnsemble('applied eigvecs')
    def_ens.setCoords(atoms)
    def_ens.setAtoms(atoms)
    def_ens.addCoordset(atoms)
    def_ens.addCoordset(def_coords)

    if not np.allclose(calcProjection(def_ens[0], modes),
                       np.zeros(modes.numModes())):
        raise ValueError('projection of atoms onto modes (via an ensemble) '
                         'is not all zeros so atoms is not appropriate')

    if ref.numAtoms() != n_atoms:
        ref = alignChains(ref, atoms)[0]
    
    def_ens.setCoords(ref)
    def_ens.superpose()

    new_vectors = np.array([np.array(coords - def_ens.getCoordsets()[0]).flatten()
                            for coords in def_ens.getCoordsets()[1:]]).T

    # initialise a new modes object with the same type
    result = type(modes)()

    result.setEigens(new_vectors, modes.getEigvals())
    return result
