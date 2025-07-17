# -*- coding: utf-8 -*-
"""This module defines functions for interfacing BioExcel-CV19 simulation database."""

__author__ = 'James Krieger'

from os.path import join, isfile, basename, splitext
from numbers import Number
import numpy as np

from prody import LOGGER, PY3K
from prody.utilities import makePath

from prody.atomic.atomgroup import AtomGroup
from prody.proteins.pdbfile import parsePDB
from prody.trajectory.psffile import parsePSF, writePSF
from prody.trajectory.dcdfile import parseDCD

__all__ = ['fetchBioexcelPDB', 'parseBioexcelPDB', 'convertXtcToDcd',
           'fetchBioexcelTrajectory', 'parseBioexcelTrajectory',
           'fetchBioexcelTopology', 'parseBioexcelTopology']

cv19_prefix = 'https://bioexcel-cv19.bsc.es/api/rest/v1/projects/'
mddb_prefix = 'https://irb.mddbr.eu/api/rest/v1/projects/'
mddb_dev_prefix = 'https://irb-dev.mddbr.eu/api/rest/v1/projects/'
dot_json_str = '.json'

def fetchBioexcelPDB(acc, **kwargs):
    """Returns a path to the downloaded PDB file corresponding
    to the starting structure from a BioExcel-CV19 trajectory.
    
    
    :arg acc: BioExcel-CV19 project accession or ID
    :type acc: str

    :arg timeout: timeout for blocking connection attempt in seconds, 
        default is 60
    :type timeout: int

    :arg folder: path to folder for storing the results
        default is the current working directory
    :type folder: str

    :arg outname: out filename, default is input ``'acc.pdb'``
    :type outname: str

    :arg selection: atom selection for download
        options are ``'_C'`` for all carbon atoms, 
        ``'backbone'`` for backbone atoms,
        or ``'backbone and _C'`` for both
    type selection: str

    :arg db: database to use, options are ``'mddb'``, ``'mddb-dev'`` and ``'cv19'``
        default is ``'mddb-dev'``
    :type db: str

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info
    """
    acc, _, selection, filepath, timeout, _ = checkInputs(acc, **kwargs)
    if not filepath.endswith('.pdb'):
        filepath += '.pdb'

    db = kwargs.get('db', 'mddb-dev')
    if db == 'cv19':
        prefix = cv19_prefix
    elif db == 'mddb-dev':
        prefix = mddb_dev_prefix
    else:
        prefix = mddb_prefix

    url = prefix + acc + "/structure"
    if selection is not None:
        url += '?selection=' + selection.replace(" ","%20")
    
    filepath = requestFromUrl(url, timeout, filepath, source='pdb', **kwargs)

    return filepath


def fetchBioexcelTrajectory(acc, **kwargs):
    """Returns a path to the downloaded BioExcel-CV19 trajectory
    in xtc format or in dcd format if mdtraj is installed and 
    *convert* is **True**.

    :arg acc: BioExcel-CV19 project accession or ID
    :type acc: str

    :arg timeout: timeout for blocking connection attempt in seconds, 
        default is 60
    :type timeout: int

    :arg folder: path to folder for storing the results
        default is the current working directory
    :type folder: str

    :arg outname: out filename, default is input ``'acc.pdb'``
    :type outname: str

    :arg frames: which frames to select in 
        e.g. ``'1-5,11-15'`` or ``'10:20:2'``
        default is to not specify and get all
    :type frames: str

    :arg selection: atom selection for download
        options are ``'_C'`` for all carbon atoms, 
        ``'backbone'`` for backbone atoms,
        or ``'backbone and _C'`` for both
    type selection: str

    :arg db: database to use, options are ``'mddb'``, ``'mddb-dev'`` and ``'cv19'``
        default is ``'mddb-dev'``
    :type db: str

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info

    :arg convert: convert to dcd if mdtraj is installed
        default is True
    type convert: bool
    """
    acc, convert, selection, filepath, timeout, frames = checkInputs(acc, **kwargs)
    if not filepath.endswith('.xtc'):
        filepath += '.xtc'

    db = kwargs.get('db', 'mddb-dev')
    if db == 'cv19':
        prefix = cv19_prefix
    elif db == 'mddb-dev':
        prefix = mddb_dev_prefix
    else:
        prefix = mddb_prefix

    url = prefix + acc + "/trajectory?format=xtc"

    if frames is not None:
        url += '&frames=' + frames

    if selection is not None:
        url += '&selection=' + selection.replace(" ","%20")

    filepath = requestFromUrl(url, timeout, filepath, source='xtc', **kwargs)

    if convert:
        filepath = convertXtcToDcd(filepath, **kwargs)

    return filepath


def fetchBioexcelTopology(acc, **kwargs):
    """Returns a path to the downloaded BioExcel-CV19 topology
    in json format or in psf format if *convert* is **True**.
    
    :arg acc: BioExcel-CV19 project accession or ID
    :type acc: str

    :arg timeout: timeout for blocking connection attempt in seconds, 
        default is 60
    :type timeout: int

    :arg folder: path to folder for storing the results
        default is the current working directory
    :type folder: str

    :arg outname: out filename, default is input ``'acc.pdb'``
    :type outname: str

    :arg db: database to use, options are ``'mddb'``, ``'mddb-dev'`` and ``'cv19'``
        default is ``'mddb-dev'``
    :type db: str

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info
    """
    if isfile(acc):
        filepath = acc
    else:
        acc, convert, _, filepath, timeout, _ = checkInputs(acc, **kwargs)
        if not filepath.endswith(dot_json_str) and not filepath.endswith('.psf'):
            filepath += dot_json_str

    if filepath.endswith('.psf'):
        convert = False

    db = kwargs.get('db', 'mddb-dev')
    if db == 'cv19':
        prefix = cv19_prefix
    elif db == 'mddb-dev':
        prefix = mddb_dev_prefix
    else:
        prefix = mddb_prefix

    if not isfile(filepath):
        url = prefix + acc + "/topology"
        filepath = requestFromUrl(url, timeout, filepath, source='json', **kwargs)

    if convert:
        ag = parseBioexcelTopology(filepath, **kwargs)
        filepath = filepath.replace(dot_json_str, '.psf')
        writePSF(filepath, ag)

    return filepath


def parseBioexcelTopology(query, **kwargs):
    """Parse a BioExcel-CV19 topology json into an :class:`.AtomGroup`,
    fetching it if needed using **kwargs
    """
    query = checkQuery(query)

    kwargs['convert'] = False
    if not isfile(query):
        filename = fetchBioexcelTopology(query, **kwargs)
    else:
        filename = query

    if filename.endswith(dot_json_str):
        import json

        fp = open(filename, 'r')
        data = json.load(fp)
        fp.close()

        title = basename(splitext(filename)[0])
        ag = AtomGroup(title)

        ag.setNames(data['atom_names'])
        ag.setElements(data['atom_elements'])
        ag.setCharges(data['atom_charges'])
        ag.setResnums(data['atom_residue_indices'])

        # set false n_csets and acsi to allow nodes in ag
        ag._n_csets = 1
        ag._acsi = 0

        indices = np.ix_(*[np.array(data['atom_residue_indices'])])

        chids = np.array([data['chain_names'][chain_index]
                          for chain_index in data['residue_chain_indices']])
        ag.setChids(chids[indices])

        ag.setResnames(np.array(data['residue_names'])[indices])
        ag.setResnums(np.array(data['residue_numbers'])[indices])

        if data['residue_icodes'] is not None:
            ag.setIcodes(np.array(data['residue_icodes'])[indices])

        # restore acsi and n_csets to defaults
        ag._acsi = None
        ag._n_csets = 0
    else:
        ag = parsePSF(filename)

    selection = checkSelection(**kwargs)

    if selection == '_C':
        ag = ag.select('element C').copy()
    elif selection == 'backbone':
        ag = ag.select('backbone').copy()
    elif selection == 'backbone and _C':
        ag = ag.select('backbone and element C')

    return ag

def parseBioexcelTrajectory(query, **kwargs):
    """Parse a BioExcel-CV19 topology json into an :class:`.Ensemble`,
    fetching it if needed using **kwargs

    :arg top: topology filename
    :type top: str
    """
    kwargs['convert'] = True
    if isfile(query) and query.endswith('.dcd'):
        filename = query
    elif isfile(query + '.dcd'):
        filename = query + '.dcd'
    elif isfile(query) and query.endswith('.xtc'):
        filename = convertXtcToDcd(query, **kwargs)
    elif isfile(query + '.xtc'):
        filename = convertXtcToDcd(query + '.xtc', **kwargs)
    else:
        filename = fetchBioexcelTrajectory(query, **kwargs)

    return parseDCD(filename)

def parseBioexcelPDB(query, **kwargs):
    """Parse a BioExcel-CV19 topology json into an :class:`.Ensemble`,
    fetching it if needed using **kwargs
    """
    kwargs['convert'] = True
    if isfile(query):
        filename = query
    elif isfile(query + '.pdb'):
        filename = query + '.pdb'
    else:
        filename = fetchBioexcelPDB(query, **kwargs)

    ag = parsePDB(filename)
    if ag is None:
        filename = fetchBioexcelPDB(query, **kwargs)
        ag = parsePDB(filename)

    acc = basename(splitext(filename)[0])
    ag2 = parseBioexcelTopology(acc, **kwargs)

    ag.setElements(ag2.getElements())
    return ag

def convertXtcToDcd(filepath, **kwargs):
    """Convert xtc trajectories to dcd files using mdtraj.
    Returns path to output dcd file.

    :arg top: topology filename
    :type top: str    
    """
    topFile = kwargs.get('top', None)
    if topFile is not None:
        acc = topFile
    else:
        acc = basename(splitext(filepath)[0])

    if not isfile(acc):
        acc = fetchBioexcelTopology(acc, **kwargs)

    try:
        import mdtraj
    except ImportError:
        raise ImportError('Please install mdtraj to convert to dcd.')
    else:
        top = mdtraj.load_topology(acc)
        traj = mdtraj.load_xtc(filepath, top=top)
        filepath = filepath.replace('xtc', 'dcd')
        traj.save_dcd(filepath)

    return filepath

def requestFromUrl(url, timeout, filepath, source=None, **kwargs):
    """Helper function to make a request from a url and return the response"""
    import requests
    import json
    import mdtraj
    import tempfile

    db = kwargs.get('db', 'mddb-dev')
    if db == 'cv19':
        prefix = cv19_prefix
    elif db == 'mddb-dev':
        prefix = mddb_dev_prefix
    else:
        prefix = mddb_prefix

    acc = url.split(prefix)[1].split('/')[0]

    LOGGER.timeit('_bioexcel')
    response = None
    sleep = 2
    while LOGGER.timing('_bioexcel') < timeout:
        try:
            response = requests.get(url).content

            if source == 'json':
                json.loads(response)

                if PY3K:
                    response = response.decode()

                fo = open(filepath, 'w')
                fo.write(response)
                fo.close()

            elif source == 'xtc':
                fo = open(filepath, 'wb')
                fo.write(response)
                fo.close()
                
                top = mdtraj.load_psf(fetchBioexcelTopology(acc, timeout=timeout))
                mdtraj.load_xtc(filepath, top=top)

            elif source == 'pdb':
                if PY3K:
                    response = response.decode()

                fo = open(filepath, 'w')
                fo.write(response)
                fo.close()

                ag = parsePDB(filepath)
                numAtoms = ag.numAtoms()

        except Exception:
            pass
        else:
            break
        
        sleep = 100 if int(sleep * 1.5) >= 100 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    return filepath

def checkSelection(**kwargs):
    """Helper function to check selection"""
    selection = kwargs.get('selection', None)
    if selection is not None:
        if not isinstance(selection, str):
            raise TypeError('selection should be a string')

        if selection not in ['_C', 'backbone', 'backbone and _C']:
            raise ValueError("selection should be '_C', 'backbone' or 'backbone and _C'")

    return selection

def checkQuery(input):
    """Check query or acc argument, which should be a string as follows:
    
    :arg acc: BioExcel-CV19 project accession or ID
    :type acc: str
    """
    if not isinstance(input, str):
        raise TypeError('query should be string')
    return input

def checkConvert(**kwargs):
    convert = kwargs.get('convert', True)
    if not isinstance(convert, (bool, type(None))):
        raise TypeError('convert should be bool')
    return convert

def checkTimeout(**kwargs):
    timeout = kwargs.get('timeout', 200)
    if not isinstance(timeout, (Number, type(None))):
        raise TypeError('timeout should be number')
    return timeout


def checkFilePath(query, **kwargs):
    """Check folder and outname and path, making it if needed"""

    folder = kwargs.get('folder', '.')
    if not isinstance(folder, str):
        try:
            folder = str(folder)
        except Exception:
            raise TypeError('folder should be a string')
    
    outname = kwargs.get('outname', None)
    if not outname:
        outname = query
    elif not isinstance(outname, str):
        raise TypeError('outname should be a string')

    return join(makePath(folder), outname)


def checkFrames(**kwargs):
    """Check frames kwarg matches the following:

    :arg frames: which frames to select in 
        e.g. ``'1-5,11-15'`` or ``'10:20:2'``
        default is to not specify and get all
    :type frames: str
    """
    frames = kwargs.get('frames', None)
    if frames is None:
        return frames
    
    if not isinstance(frames, str):
        raise TypeError('frames should be a string')
 
    for framesRange in frames.split(','):
        if framesRange.find('-') != -1 and framesRange.find(':') != -1:
            raise ValueError('Frames should have a set of comma-separated ranges containing either a hyphen or colons, not both')
        
        if framesRange.count('-') > 0:
            if not np.all([item.isnumeric() for item in framesRange.split('-')]):
                raise ValueError('Each frames range should only have numbers and hyphens (or colons), not spaces or anything else')
            if framesRange.count('-') > 1:
                raise ValueError('Each frames range can only have one hyphen')
        
        if framesRange.count(':') not in [0,1,2]:
            raise ValueError('Each frames range can have no more than 2 colons')

        if framesRange.count(':') > 0:
            if not np.all([item.isnumeric() for item in framesRange.split(':')]):
                raise ValueError('Each frames range should only have numbers and colons (or hyphens), not spaces or anything else')

        if (framesRange.find('-') == -1 and framesRange.find(':') == -1 and 
            not framesRange.isnumeric()):
            raise ValueError('Each frames range should only have numbers (and colons or hyphens), not spaces or anything else')

    return frames

def checkInputs(query, **kwargs):
    query = checkQuery(query)
    timeout = checkTimeout(**kwargs)
    convert = checkConvert(**kwargs)
    selection = checkSelection(**kwargs)
    filepath = checkFilePath(query, **kwargs)
    frames = checkFrames(**kwargs)
    return query, convert, selection, filepath, timeout, frames
