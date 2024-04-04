# -*- coding: utf-8 -*-
"""This module defines functions for interfacing BioExcel-CV19 simulation database."""

__author__ = 'James Krieger'

from os.path import join, isfile, basename, splitext

from prody import LOGGER, PY3K
from prody.utilities import makePath

from prody.atomic.atomgroup import AtomGroup
from prody.atomic.functions import extendAtomicData
from prody.proteins.pdbfile import parsePDB
from prody.trajectory.psffile import parsePSF, writePSF
from prody.trajectory.dcdfile import parseDCD

__all__ = ['fetchBioexcelPDB', 'parseBioexcelPDB', 'convertXtcToDcd',
           'fetchBioexcelTrajectory', 'parseBioexcelTrajectory',
           'fetchBioexcelTopology', 'parseBioexcelTopology']

prefix = 'https://bioexcel-cv19.bsc.es/api/rest/v1/projects/'
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

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info
    """

    url = prefix + acc + "/structure"

    selection = checkSelection(**kwargs)
    if selection is not None:
        url += '?selection=' + selection.replace(" ","%20")

    LOGGER.timeit('_bioexcel')
    timeout = kwargs.get('timeout', 60)
    response = requestFromUrl(url, timeout)

    if PY3K:
        response = response.decode()

    folder = str(kwargs.get('folder', '.'))
    outname = kwargs.get('outname', None)
    if not outname:
        outname = acc
    if not outname.endswith('.pdb'):
        outname += '.pdb'
    filepath = join(makePath(folder), outname)
    fo = open(filepath, 'w')
    fo.write(response)
    fo.close()

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

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info

    :arg convert: convert to dcd if mdtraj is installed
        default is True
    type convert: bool
    """

    url = prefix + acc + "/trajectory?format=xtc"

    convert = kwargs.get('convert', True)
    if not isinstance(convert, bool):
        raise TypeError('convert should be a bool')

    frames = kwargs.get('frames', None)
    if frames is not None:
        if not isinstance(frames, str):
            raise TypeError('frames should be a string')
        
        url += '&frames=' + frames

    selection = checkSelection(**kwargs)
    if selection is not None:
        url += '&selection=' + selection.replace(" ","%20")

    LOGGER.timeit('_bioexcel')
    timeout = kwargs.get('timeout', 60)
    response = requestFromUrl(url, timeout)

    folder = str(kwargs.get('folder', '.'))
    outname = kwargs.get('outname', None)
    if not outname:
        outname = acc
    if not outname.endswith('.xtc'):
        outname += '.xtc'
    filepath = join(makePath(folder), outname)
    fo = open(filepath, 'wb')
    fo.write(response)
    fo.close()

    if convert:
        filepath = convertXtcToDcd(filepath)

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

    See https://bioexcel-cv19.bsc.es/api/rest/docs for more info
    """

    url = prefix + acc + "/topology"

    convert = kwargs.get('convert', True)
    if not isinstance(convert, bool):
        raise TypeError('convert should be a bool')

    LOGGER.timeit('_bioexcel')
    timeout = kwargs.get('timeout', 60)
    response = requestFromUrl(url, timeout)

    if PY3K:
        response = response.decode()

    folder = str(kwargs.get('folder', '.'))
    outname = kwargs.get('outname', None)
    if not outname:
        outname = acc
    if not outname.endswith(dot_json_str):
        outname += dot_json_str
    filepath = join(makePath(folder), outname)
    fo = open(filepath, 'w')
    fo.write(response)
    fo.close()

    if convert:
        ag = parseBioexcelTopology(filepath, **kwargs)
        filepath = filepath.replace(dot_json_str, '.psf')
        writePSF(filepath, ag)

    return filepath


def parseBioexcelTopology(query, **kwargs):
    """Parse a BioExcel-CV19 topology json into an :class:`.AtomGroup`,
    fetching it if needed using **kwargs
    """
    kwargs.pop('convert', True)
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

        nodes = ag.select('name N')

        residue_chids = [data['chain_names'][chain_index] for chain_index in data['residue_chain_indices']]
        chids, _ = extendAtomicData(residue_chids, nodes, ag)
        ag.setChids(chids)

        resnames, _ = extendAtomicData(data['residue_names'], nodes, ag)
        ag.setResnames(resnames)

        resnums, _ = extendAtomicData(data['residue_numbers'], nodes, ag)
        ag.setResnums(resnums)

        if data['residue_icodes'] is not None:
            icodes, _ = extendAtomicData(data['residue_icodes'], nodes, ag)
            ag.setIcodes(icodes)

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
    """
    kwargs.pop('convert', True)
    kwargs['convert'] = True
    if isfile(query):
        filename = query
    elif isfile(query + '.xtc'):
        filename = convertXtcToDcd(query + '.xtc')
    else:
        filename = fetchBioexcelTrajectory(query, **kwargs)

    return parseDCD(filename)

def parseBioexcelPDB(query, **kwargs):
    """Parse a BioExcel-CV19 topology json into an :class:`.Ensemble`,
    fetching it if needed using **kwargs
    """
    kwargs.pop('convert', True)
    kwargs['convert'] = True
    if not isfile(query):
        filename = fetchBioexcelPDB(query, **kwargs)
    else:
        filename = query

    return parsePDB(filename)

def convertXtcToDcd(filepath):
    """Convert xtc trajectories to dcd files using mdtraj.
    Returns path to output dcd file.
    """
    acc = basename(splitext(filepath)[0])
    try:
        import mdtraj
    except ImportError:
        raise ImportError('Please install mdtraj to convert to dcd.')
    else:
        top = mdtraj.load_psf(fetchBioexcelTopology(acc))
        traj = mdtraj.load_xtc(filepath, top=top)
        filepath = filepath.replace('xtc', 'dcd')
        traj.save_dcd(filepath)

    return filepath

def requestFromUrl(url, timeout):
    """Helper function to make a request from a url and return the response"""
    import requests

    response = None
    LOGGER.timeit('_bioexcel')
    response = None
    sleep = 2
    while LOGGER.timing('_bioexcel') < timeout:
        try:
            response = requests.get(url).content
        except Exception:
            pass
        else:
            break
        
        sleep = 20 if int(sleep * 1.5) >= 20 else int(sleep * 1.5)
        LOGGER.sleep(int(sleep), '. Trying to reconnect...')

    return response

def checkSelection(**kwargs):
    """Helper function to check selection"""
    selection = kwargs.get('selection', None)
    if selection is not None:
        if not isinstance(selection, str):
            raise TypeError('selection should be a string')

        if selection not in ['_C', 'backbone', 'backbone and _C']:
            raise ValueError("selection should be '_C', 'backbone' or 'backbone and _C'")

    return selection
