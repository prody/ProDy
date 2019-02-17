# -*- coding: utf-8 -*-
"""This module defines functions for interfacing with the EBI's
GO association (goa) databases for analysing gene/protein functions 
through the Gene Ontology (GO).

This module is based on the tutorial notebook at 
https://nbviewer.jupyter.org/urls/dessimozlab.github.io/go-handbook/GO%20Tutorial%20in%20Python%20-%20Solutions.ipynb"""

from io import BytesIO
from os.path import isdir, isfile, join, split, splitext, normpath
import os
import gzip
from ftplib import FTP

from collections import defaultdict
import re
import numpy as np
from prody.utilities import makePath, gunzip, relpath, copyFile, openURL
from prody.utilities import openFile, isListLike, sympath
from prody import LOGGER, PY3K
from prody.ensemble import Ensemble

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2


__all__ = ['GOADictList', 'parseOBO', 'parseGAF',
           'queryGOA', 'showGoLineage',
           'calcGoOverlap', 'calcDeepFunctionOverlaps',
           'calcEnsembleFunctionOverlaps']


class GOADictList:
    """A class for handling the list of GOA Dictionaries returned 
    by queryGOA
    """

    def __init__(self, parsingList, title='unnamed', **kwargs):
        go = kwargs.pop('go', None)
        if go is None:
            go = parseOBO(**kwargs)

        self._title = title
        self._list = parsingList
        go_ids = np.unique([entry['GO_ID'] for entry in parsingList])
        self._go_terms = [go[id] for id in go_ids]
        self.numTerms = len(self._go_terms)

        self._molecular = [
            term for term in self._go_terms if term.namespace == 'molecular_function']
        self._cellular = [
            term for term in self._go_terms if term.namespace == 'cellular_component']
        self._process = [
            term for term in self._go_terms if term.namespace == 'biological_process']

    def __getitem__(self, ind):
        return self._list[ind]

    def getTitle(self):
        return self._title

    def setTitle(self, value):
        self._title = value

    def getList(self):
        return self._list

    def __repr__(self):
        if self.numTerms == 1:
            return '<GOADictList: {0} (1 GOA dict)>'.format(self._title)
        return '<GOADictList: {0} ({1} GOA dicts)>'.format(self._title, self.numTerms)

    def __iter__(self):
        """Yield go_term instances."""
        for item in self._list:
            yield item

    def pop(self, index):
        """Pop dataBlock with the given index from the list of dataBlocks in GOADictList"""
        self._go_terms.pop(index)
        self.numTerms -= 1

    def getGoTerms(self):
        return self._go_terms

    def getMolecularFunctions(self):
        return self._molecular

    def getCellularComponents(self):
        return self._cellular

    def getBiologicalProcesses(self):
        return self._process


def parseOBO(**kwargs):
    """Parse a GO OBO file containing the GO itself.
    See `OBO`_ for more information on the file format.

    .. _OBO: http://owlcollab.github.io/oboformat/doc/obo-syntax.html
    """
    from goatools import obo_parser

    go_obo_url = kwargs.get('go_obo_url', None)
    if go_obo_url is None:
        go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'

    data_folder = kwargs.get('data_folder', None)
    if data_folder is None:
        data_folder = os.getcwd() + '/Data'

    # Check if we have the ./data directory already
    if(not os.path.isfile(data_folder)):
        # Emulate mkdir -p (no error if folder exists)
        try:
            os.mkdir(data_folder)
        except OSError as e:
            if(e.errno != 17):
                raise e
    else:
        raise Exception('Data path (' + data_folder + ') exists as a file. '
                        'Please rename, remove or change the desired location of the data path.')

    # Check if the file exists already
    if(not os.path.isfile(data_folder+'/go-basic.obo')):
        try:
            handle = openURL(go_obo_url)
        except Exception as err:
            LOGGER.warn('{0} download failed ({1}).'.format(
                go_obo_url, str(err)))
        else:
            data = handle.read()
            if len(data):
                filename = data_folder+'/go-basic.obo'

                with open(filename, 'w+b') as obofile:
                    obofile.write(data)

                LOGGER.debug('{0} downloaded ({1})'
                             .format(go_obo_url, sympath(filename)))
            else:
                LOGGER.warn('{0} download failed, reason unknown.'
                            .format(go_obo_url))

    else:
        go_obo = data_folder+'/go-basic.obo'

    return obo_parser.GODag(go_obo)


def parseGAF(database='PDB', **kwargs):
    """Parse a GO Association File (GAF) corresponding to
    a particular database collection into a dictionary 
    for ease of querying.

    See `GAF`_ for more information on the file format

    .. _GAF: http://geneontology.org/docs/go-annotation-file-gaf-format-20/

    :arg database: name of the database of interest
        default is PDB. Others include UNIPROT and 
        common names of many organisms.
    :type database: str

    :arg filename: filename for the gaf of interest
        default is goa_ and the database name in lower case
        and .gaf.gz
    :type filename: str
    """
    import Bio.UniProt.GOA as GOA

    if not isinstance(database, str):
        raise TypeError('database should be a string')

    database = database.upper()
    filename = kwargs.get('filename', None)
    if filename is None:
        if database == 'UNIPROT':
            filename = 'goa_' + database.lower() + '_all.gaf.gz'
        else:
            filename = 'goa_' + database.lower() + '.gaf'

    data_folder = kwargs.get('data_folder', os.getcwd())

    # If the file doesn't already exist, download it
    gaf = os.path.join(data_folder, filename)
    if not(os.path.exists(gaf) and os.path.getsize(gaf) > 0):
        LOGGER.info('Downloading file {0} to {1}'.format(filename, gaf))
        data_stream = BytesIO()
        ftp_host = 'ftp.ebi.ac.uk'
        ftp = FTP(ftp_host)
        ftp.login()

        try:
            ftp.cwd('pub/databases/GO/goa')
            ftp.cwd(database)
            ftp.retrbinary('RETR {}.gz'.format(filename), data_stream.write)
        except:
            raise ValueError('Cannot find the requested GO association file')

        # Logout from FTP server
        ftp.quit()

        zip_data = data_stream.getvalue()
        data_stream.close()

        rawdata = gunzip(zip_data)
        if PY3K:
            rawdata = rawdata.decode()

        with open(filename, 'w') as gaf_fp:
            gaf_fp.write(rawdata)

        LOGGER.info('Download completed for file {0}'.format(filename))

    with open(filename, 'rt') as gaf_fp:
        funcs = defaultdict(list)  # Initialise the dictionary of functions

        # Iterate on each function using Bio.UniProt.GOA library.
        LOGGER.info('Iterating through entries in {0}'.format(gaf))
        for entry in GOA.gafiterator(gaf_fp):
            id = entry.pop('DB_Object_ID')
            funcs[id].append(entry)

    return funcs


def queryGOA(*ids, **kwargs):
    """Query a GOA database by identifier.

    :arg ids: an identifier or a list-like of identifiers 
    :type ids: str, tuple, list, :class:`~numpy.ndarray`

    :arg database: name of the database of interest
        default is PDB. Others include UNIPROT and 
        common names of many organisms.
    :type database: str
    """
    database = kwargs.pop('database', 'PDB')

    gaf_dict = kwargs.pop('gaf_dict', None)
    if gaf_dict is None:
        gaf_dict = parseGAF(database=database, **kwargs)
        LOGGER.info('GAF parsing completed.')

    n_ids = len(ids)
    if n_ids == 1:
        if isListLike(ids[0]):
            ids = ids[0]
            n_ids = len(ids)

    if n_ids == 1:
        ids = list(ids)

    results = []
    unmapped = []
    LOGGER.progress('Querying GOA for {0} ids...'
                    .format(n_ids), n_ids, '_prody_queryGOA')
    for i, id in enumerate(ids):
        LOGGER.update(i, 'Querying GOA for id {0} of {1}...'
                      .format(i+1, n_ids), label='_prody_queryGOA')
        if not isinstance(id, str):
            raise TypeError('each ID should be a string')

        id = id.upper()

        if database == 'PDB':
            if not len(id) in [4, 5, 6]:
                raise ValueError('PDB IDs should be strings of length 4 to 6')

            if len(id) == 5 and str.isalpha(id[-1]):
                id = id[:4] + '_' + id[-1]

        if id in list(gaf_dict.keys()):
            results.append(gaf_dict[id])
        else:
            results.append([])
            unmapped.append(id)

    rets = []
    LOGGER.progress('Mapping GO terms back to GOA results for {0} ids...'
                    .format(n_ids), n_ids, '_prody_mapGO')
    for i, result in enumerate(results):
        LOGGER.update(i, 'Mapping GO terms back to GOA results id {0} of {1}...'
                      .format(i+1, n_ids), label='_prody_mapGO')
        rets.append(GOADictList(result, title=id))

    if n_ids == 1:
        rets = rets[0]

    return rets


def calcGoOverlap(*go_terms, **kwargs):
    """Calculate overlap between GO terms based on their distance
    in the graph. GO terms in different namespaces (molecular function,
    cellular component, and biological process) have undefined distances.

    :arg go_terms: a list of GO terms or GO IDs
    :type go_terms: list, tuple, `~numpy.ndarray`

    :arg pairwise: whether to calculate to a matrix of pairwise overlaps
        default is False
    :type pairwise: bool

    :arg distance: whether to return distances rather than calculating overlaps
        default is False
    :type distance: bool

    :arg go: GO graph. Default behaviour is to parse it with :func:`.parseOBO`.
    :type go: `~goatools.obo_parser.GODag`
    """
    pairwise = kwargs.pop('pairwise', False)
    distance = kwargs.pop('distance', False)

    go = kwargs.pop('go', None)
    if go is None:
        go = parseOBO(**kwargs)

    if not isListLike(go_terms):
        raise TypeError('please provide a list-like of go terms')

    try:
        go_terms = [go[term] for term in go_terms]
    except:
        try:
            go_terms = [term.id for term in go_terms]
        except:
            raise TypeError('go_terms should contain go terms or IDs')

    if isinstance(go_terms[0], str):
        go_terms = [go[term] for term in go_terms]

    if pairwise:
        distances = np.zeros((len(go_terms), len(go_terms)))

        for i in range(len(go_terms)):
            for j in range(i+1, len(go_terms)):
                dist = min_branch_length(go_terms[i], go_terms[j], go)
                distances[i, j] = distances[j, i] = dist

    else:
        distances = np.zeros((len(go_terms)-1))

        go_id1 = go_terms[0]
        for i, go_id2 in enumerate(go_terms[1:]):
            distances[i] = min_branch_length(go_id1, go_id2, go)

    if distance:
        return distances
    else:
        return 1. / distances


def min_branch_length(go_id1, go_id2, go):
    '''Find the minimum branch length between two terms in the GO DAG.

    :arg go_id1: the first GO ID
    :type go_id1: str

    :arg go_id2: the second GO ID
    :type go_id2:str

    :arg go: object containing a gene ontology (GO) directed acyclic graph (DAG)
    :type go: `~goatools.obo_parser.GODag`
    '''
    # First get the deepest common ancestor
    dca = deepest_common_ancestor([go_id1, go_id2], go)
    if dca is None:
        LOGGER.warn('There are no common ancestors between {0} and {1} so no meaningful distance can be calculated.'.format(
            go_id1, go_id2))
        return None

    # Then get the distance from the DCA to each term
    dca_depth = go[dca].depth
    d1 = go[go_id1].depth - dca_depth
    d2 = go[go_id2].depth - dca_depth

    # Return the total distance - i.e., to the deepest common ancestor and back.
    return d1 + d2


def deepest_common_ancestor(terms, go):
    '''Find the nearest common ancestor. 
    Only returns single most specific - assumes unique exists.

    :arg terms: a list of GO terms
    :type terms: tuple, list, :class:`~numpy.ndarray`

    :arg go: object containing a gene ontology (GO) directed acyclic graph (DAG)
    :type go: `~goatools.obo_parser.GODag`
    '''
    # Take the element at maximum depth.
    common_parent = common_parent_go_ids(terms, go)
    if common_parent == set():
        return None
    return max(common_parent, key=lambda t: go[t].depth)


def common_parent_go_ids(terms, go):
    '''This function finds the common ancestors in the GO 
    tree of the list of terms in the input.

    :arg terms: a list of GO terms
    :type terms: tuple, list, :class:`~numpy.ndarray`

    :arg go: object containing a gene ontology (GO) directed acyclic graph (DAG)
    :type go: `~goatools.obo_parser.GODag`
    '''
    # Find candidates from first
    rec = go[terms[0]]
    candidates = rec.get_all_parents()
    candidates.update({terms[0]})

    # Find intersection with second to nth term
    for term in terms[1:]:
        rec = go[term]
        parents = rec.get_all_parents()
        parents.update({term})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)

    return candidates


def showGoLineage(go_term, **kwargs):
    """Use pygraphviz and IPython notebook to show the lineage of a GO term

    :arg go: object containing a gene ontology (GO) directed acyclic graph (DAG) 
        default is to parse with :func:`.parseOBO`
    :type go: `~goatools.obo_parser.GODag`

    arg out_format: format for output. 
        Currently only output to file. This file will be displayed in Jupyter Notebook.
    type out_format: str

    arg filename: filename for output
        default behaviour is to use the GO term ID and append '_lineage.png'
    type filename: str
    """
    out_format = kwargs.pop('format', 'png')

    if out_format == 'png':
        filename = kwargs.pop('filename', '_'.join(
            go_term.id.split(':')) + '_lineage.png')

        go = kwargs.pop('go', None)
        if go is None:
            go = parseOBO(**kwargs)

        go.draw_lineage([go_term], lineage_img=filename)
        from IPython.display import Image
        Image(filename)


def calcDeepFunctionOverlaps(goa1, goa2, **kwargs):
    """Calculate function overlaps between the deep 
    (most detailed) molecular functions in particular 
    from two sets of GO terms.

    :arg goa1: the first set of GO terms
    :type goa1: tuple, list, :class:`~numpy.ndarray`

    :arg goa2: the second set of GO terms
    :type goa2: tuple, list, :class:`~numpy.ndarray`
    """
    return_funcs = kwargs.pop('return_funcs', False)

    xfuncs = findDeepestFunctions(goa1)
    yfuncs = findDeepestFunctions(goa2)
    overlaps = np.zeros((len(xfuncs), len(yfuncs)))
    for i, function_i in enumerate(xfuncs):
        for j, function_j in enumerate(yfuncs):
            overlaps[i, j] = calcGoOverlap(function_i, function_j, **kwargs)

    if return_funcs:
        return overlaps, xfuncs, yfuncs

    return overlaps


def findDeepestFunctions(go_terms, **kwargs):
    """Find the deepest (most detailed) molecular functions in 
    a list of GO terms.

    :arg go_terms: a list of GO terms
    :type go_terms: tuple, list, :class:`~numpy.ndarray`
    """
    if not isListLike(go_terms):
        raise TypeError('go_terms should be list-like')

    go = kwargs.pop('go', None)
    if go is None:
        go = parseOBO(**kwargs)

    try:
        go_terms = [go[term] for term in go_terms]
    except:
        try:
            go_terms = [term.id for term in go_terms]
        except:
            raise TypeError('go_terms should contain go terms or IDs')

    if isinstance(go_terms[0], str):
        go_terms = [go[term] for term in go_terms]

    deep_functions = []
    max_depth = 0
    for function in go_terms.getMolecularFunctions():
        if function.depth > max_depth:
            max_depth = function.depth

    for function in go_terms.getMolecularFunctions():
        if function.depth == max_depth:
            deep_functions.append(function)

    return deep_functions


def calcEnsembleFunctionOverlaps(ens, **kwargs):
    """Calculate function overlaps for an ensemble as the 
    mean of the value from :func:`calcDeepFunctionOverlaps`.

    :arg ens: an ensemble with labels
    :type ens: :class:`Ensemble`
    """
    if not isinstance(ens, Ensemble) and not isListLike(ens):
        raise TypeError('ens should be an ensemble or list-like')

    if isinstance(ens, Ensemble):
        ids = [label[:5] for label in ens.getLabels()]
    else:
        ids = ens

    if not isinstance(ids[0], str):
        raise TypeError('ens should have labels')

    distance = kwargs.pop('distance', False)
    pairwise = kwargs.pop('pairwise', False)

    goa_ens = queryGOA(ids, **kwargs)

    overlaps = np.zeros((len(goa_ens), len(goa_ens)))
    for i, funcs_i in enumerate(goa_ens):
        for j, funcs_j in enumerate(goa_ens[i:]):
            overlaps[i, j] = np.mean(calcDeepFunctionOverlaps(
                funcs_i, funcs_j, distance=distance, pairwise=pairwise, **kwargs))

    return overlaps
