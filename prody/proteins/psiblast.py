import os.path
import numpy as np
from prody import LOGGER, PY3K
from prody.utilities import dictElement, openURL, which

import platform, os, re, sys, time, urllib

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

import xml.etree.cElementTree as ET

from prody.sequence import Sequence
from prody.atomic import Atomic
from prody.proteins.pdbfile import parsePDB

__all__ = ['PsiBlastRecord', 'psiBlastCycle', 'psiBlastRun']

def psiBlastRun(sequence, cycles=2, filename=None, **kwargs):
    """Returns the results from a full PSI-BLAST run (multiple cycles).
    All arguments are the same as psiBlastCycle and are passed to it
    except for cycles.

    :arg cycles: the number of cycles to run
        default is 2
    :type cycles: int
    """
    psithr = kwargs.get('psithr', 1.0e-3)
    job_id = kwargs.get('previousjobid','') 
    selectedHits = kwargs.get('selectedHits','')

    cycles_done = 0
    results_list = []
    job_ids = []
    while cycles_done < cycles:
        if cycles_done > 0:
            selectedHits = 'http://www.ebi.ac.uk/Tools/services/rest/psiblast/result/' \
                      + job_id + '/preselected_seq'
            sequence = None
        job_id, results, sequence = psiBlastCycle(sequence, filename, \
                                                 previousjobid=job_id, \
                                                 selectedHits=selectedHits, \
                                                 **kwargs)
        results_list.append(results)
        job_ids.append(job_id)
        cycles_done += 1
        LOGGER.info('Finished cycle {0} with job ID {1}.'.format(cycles_done, job_id))

    return job_ids, results_list, sequence

def psiBlastCycle(sequence, filename=None, **kwargs):
    """Returns a :class:`PDBBlastRecord` instance that contains results from
    a single cycle of EBI psiblast.

    :arg sequence: an object with an associated sequence string 
         or a sequence string itself
    :type sequence: :class:`Atomic`, :class:`Sequence`, or str

    :arg filename: a *filename* to save the results in XML format
    :type filename: str

    The following search parameters can be adjusted by the user.
    We use the same default values as 
    http://www.ebi.ac.uk/Tools/services/rest/psiblast/parameterdetails/
    wherever applicable.

    :arg email: email address for reporting problems
        default is prody-devel@gmail.com
    :type email: str with an @ before a .

    :arg matrix: The comparison matrix to be used to score alignments when searching the database
        possible values are 'BLOSUM45', 'BLOSUM62', 'BLOSUM80', 'PAM30' and 'PAM70' 
        default is 'BLOSUM62'
    :type matrix: str

    :arg gapopen: Penalty taken away from the score when a gap is created in sequence alignments. 
        Increasing the gap opening penalty will decrease the number of gaps in the final alignment.
        Possible values range from 8 to 16 inclusive, default is 11
    :type gapopen: int

    :arg gapext: Penalty taken away from the score for each base or residue in the gap. 
        Increasing the gap extension penalty favors short gaps in the final alignment, 
        conversly decreasing the gap extension penalty favors long gaps in the final alignment. 
        Possible values range from 0 to 3, default is 1
    :type gapext: int

    :arg expthr: Expectation threshold that limits the number of scores and alignments reported. 
        This is the maximum number of times the match is expected to occur by chance.
        Possible values are 1.0e-200, 1.0e-100, 1.0e-50, 1.0e-10, 1.0e-5, 1.0e-4, 1.0e-3,
        1.0e-2, 0.1, 1.0, 10.0, 100, 1000
        default is 10.0
    :type expthr: float

    :arg psithr: Expectation value threshold for automatic selection of matched sequences for 
        inclusion in the PSSM at each iteration.
        Possible values are 1.0e-6, 1.0e-5, 1.0e-4, 2.0e-4, 5.0e-4, 1.0e-3, 2.0e-3, 5.0e-3,
        1.0e-2, 2.0e-2, 0.1, 0.3, 0.5, 1.0, 3.0, 10.0
        default is 1.0e-3
    :type psithr: float

    :arg scores: Maximum number of match score summaries reported in the result output.
        Possible values are 5, 10, 20, 50, 100, 200, 500, 750, 1000, or 5000
        Default is 500
    :type scores: int

    :arg alignments: Maximum number of match alignments reported in the result output.
        Possible values are 5, 10, 20, 50, 100, 200, 500, 750, 1000, or 5000
        Default is 500
    :type alignmets: int

    :arg dropoff: The amount a score can drop before extension of word hits is halted
        Possible values are 0, 2, 4, 6, 8, 10, 15, 20, 25, or 30
        Default is 15
    :type dropoff: int

    :arg finaldropoff: Dropoff value for final gapped alignment
        Possible values are 10, 12, 14, 16, 18, 20, 22, 24, 25, 26, 28, or 30
        Default is 25
    :type finaldropoff: int

    :arg filter: Filter regions of low sequence complexity. This can avoid issues with 
        low complexity sequences where matches are found due to composition rather than 
        meaningful sequence similarity. However, in some cases filtering also masks 
        regions of interest and so should be used with caution.
        Possible values are T and F, default is F
    :type filter: str

    :arg seqrange: Specify a range or section of the input sequence to use in the search.
        Example: Specifying '34-89' in an input sequence of total length 100, will tell BLAST 
        to only use residues 34 to 89, inclusive.
    :type seqrange: str of form START-END

    :arg database: a database name from those available. See
        http://www.ebi.ac.uk/Tools/services/rest/psiblast/parameterdetails/database
        default is pdb
    :type database: str

    :arg previousjobid: The job identifier for the previous PSI-BLAST iteration. 
        default is None
        You can change this if you want to continue from a previous run
    :type previousjobid: str

    :arg selectedHits: Name of a file containing a list of identifiers of the 
        hits from the previous iteration to use to construct the search PSSM 
        for this iteration.
        default is None
    :type selectedHits: str

    :arg cpfile: Name of a Checkpoint file from the previous iteration. 
        default is None
    :type cpfile: str

    :arg sleep: how long to wait to reconnect for status
         Sleep time is multiplied by 1.5 when results are not ready.
         default is 2 seconds
    :type sleep: float

    :arg timeout:  when to give up waiting for the results 
        default is 120 seconds
    :type timeout: float

    """
    if sequence == 'runexample':
        sequence = ('ASFPVEILPFLYLGCAKDSTNLDVLEEFGIKYILNVTPNLPNLFENAGEFKYKQIPI'
                    'SDHWSQNLSQFFPEAISFIDEARGKNCGVLVHSLAGISRSVTVTVAYLMQKLNLSMN'
                    'DAYDIVKMKKSNISPNFNFMGQLLDFERTL')

    elif isinstance(sequence, Atomic):
        sequence = sequence.calpha.getSequence()

    elif isinstance(sequence, Sequence):
        sequence = str(sequence)

    elif isinstance(sequence, str):
        if len(sequence) == 4 or len(sequence) == 5 or len(sequence) == 6:
            sequence = parsePDB(sequence)
            sequence = ag.calpha.getSequence()
        sequence = ''.join(sequence.split())

    elif sequence is None:
        pass

    else:
        raise TypeError('sequence must be Atomic, Sequence, or str not {0}'
                        .format(type(sequence)))

    query = [('sequence', sequence)]

    email = kwargs.get('email','prody-devel@gmail.com')
    if not isinstance(email, str):
        raise TypeError('email must be a string')
    elif email.find('@') == -1 or email.find('.') == -1 or len(email.split('@')) != 2:
        raise ValueError('email must be a valid email address with at least one . and exactly one @ sign')
    elif not email.find('@') < email.find(email.split('.')[-1]):
        raise ValueError('email must be a valid email address with a . after the @ sign')
    query.append(('email', email))
    query.append(('title', 'ProDy psiBlastPDB request'))
 
    matrix = kwargs.get('matrix', 'BLOSUM62')
    checkPsiBlastParameter('matrix', matrix)
    query.append(('matrix',matrix))

    gapopen = kwargs.get('gapopen',11)
    checkPsiBlastParameter('gapopen', gapopen)
    query.append(('gapopen',gapopen))

    gapext = kwargs.get('gapext',1)
    checkPsiBlastParameter('gapext', gapext)
    query.append(('gapext',gapext))

    expthr = kwargs.get('expthr', 10.)
    checkPsiBlastParameter('expthr', expthr)
    query.append(('expthr',expthr))
    
    psithr = kwargs.get('psithr',1.0e-3)
    checkPsiBlastParameter('psithr', psithr)
    query.append(('psithr',psithr))

    scores = kwargs.get('scores',500)
    checkPsiBlastParameter('scores', scores)
    query.append(('scores',scores))

    alignments = kwargs.get('alignments',500)
    checkPsiBlastParameter('alignments', alignments)
    query.append(('alignments',alignments))
    
    query.append(('alignView',0))
                    
    dropoff = kwargs.get('dropoff',15)
    checkPsiBlastParameter('dropoff', dropoff)
    query.append(('dropoff',dropoff))
        
    finaldropoff = kwargs.get('finaldropoff',25)
    checkPsiBlastParameter('finaldropoff', finaldropoff)
    query.append(('finaldropoff',finaldropoff))
        
    filter = kwargs.get('filter','F')
    checkPsiBlastParameter('filter', filter)
    query.append(('filter',filter))
            
    seqrange = kwargs.get('seqrange', None)
    if seqrange is None:
        seqrange = '0-' + str(len(sequence))
    elif not isinstance(seqrange, str):
        raise TypeError('seqrange should be a string')
    elif len(seqrange.split('-')) != 2:
        raise ValueError('seqrange should take the form START-END')
    try:
        start = int(seqrange.split('-')[0])
        end = int(seqrange.split('-')[1])
    except:
        raise ValueError('seqrange should be START-END with START and END being integers')
    query.append(('seqrange',seqrange))
    
    database = kwargs.get('database','pdb')
    checkPsiBlastParameter('database', database)
    query.append(('database',database))
 
    headers = { 'User-Agent' : 'ProDy' }
    
    try:
        import urllib.parse
        urlencode = lambda data: bytes(urllib.parse.urlencode(data), 'utf-8')
    except ImportError:
        from urllib import urlencode

    sleep = float(kwargs.pop('sleep', 2))
    timeout = float(kwargs.pop('timeout', 120))
    
    data = urlencode(query)

    # submit the job
    base_url = 'http://www.ebi.ac.uk/Tools/services/rest/psiblast/'
    url = base_url + 'run/'
    LOGGER.timeit('_prody_psi-blast')
    LOGGER.info('PSI-Blast searching NCBI PDB database for "{0}..."'
                .format(sequence[:5]))

    handle = openURL(url, data=data, headers=headers)
    job_id = handle.read()
    handle.close()

    # check the status
    url = base_url + 'status/' + job_id
    handle = openURL(url)
    status = handle.read()
    handle.close()
                    
    # keep checking the status until it's no longer running
    while status is 'RUNNING':
        LOGGER.sleep(int(sleep), 'to reconnect to EBI for status.')
        LOGGER.write('Connecting to EBI for status...')
        handle = openURL(url)
        status = handle.read()
        LOGGER.clear()
        sleep = int(sleep * 1.5)
        if LOGGER.timing('_prody_psi-blast') > timeout:
            LOGGER.warn('PSI-Blast search time out.')
            return None

    # check status once it's not running and tell the user
    LOGGER.sleep(int(sleep), 'to reconnect to EBI for status.')
    LOGGER.write('Connecting to EBI for status...')
    handle = openURL(url)
    status = handle.read()
    handle.close()
    LOGGER.info('The status is {0}'.format(status))
 
    # get the results
    url = base_url + 'result/' + job_id + '/xml'
    handle = openURL(url)
    results = handle.read()
    handle.close()

    LOGGER.clear()
    LOGGER.report('PSI-Blast search completed in %.1fs.', '_prody_psi-blast')
    
    try:
        ext_xml = filename.lower().endswith('.xml')
    except AttributeError:
        pass
    else:
        if not ext_xml:
            filename += '.xml'
        f_out = open(filename, 'w')
        f_out.write(results)
        f_out.close()
        LOGGER.info('Results are saved as {0}.'.format(repr(filename)))
    
    return job_id, PsiBlastRecord(results, sequence), sequence

def checkPsiBlastParameter(parameter, value):
    """Checks that the value provided for a parameter is in the xml page for that parameter
    and raises an error if it isn't.

    :arg parameter: parameter name
    :type parameter: str
    
    :arg value: value being checked
    :type value: any
    """
    info_file = urllib2.urlopen('http://www.ebi.ac.uk/Tools/services/rest/psiblast/parameterdetails/' + parameter)
    data = info_file.read()
    info_file.close()
    
    data = ET.XML(data)
    data = dictElement(data)

    if data['type'] == 'STRING':
        type = str
    elif data['type'] == 'INTEGER':
        type = int
    elif data['type'] == 'DOUBLE':
        type = float
    
    if not isinstance(value,type):
        raise TypeError(name + ' should be of type ' + \
                        str(type).split()[1].strip("'>"))
    
    values = []
    str_values = []
    for element in data['values']:
        values.append(type(dictElement(element)['value']))
        str_values.append(str(dictElement(element)['value']))

    if not value in values:
        raise ValueError(parameter + ' should be one of ' + \
                         ', '.join(str_values[:-1]) \
                         + ', or ' + str_values[-1])

    return


class PsiBlastRecord(object):

    """A class to store results from psi-blast searches."""


    def __init__(self, results, sequence=None):
        """Instantiate a PDBBlastRecord object instance.

        :arg result: psi-blast search results in XML format or an XML file
            that contains the results
        :type result: str

        :arg sequence: query sequence
        :type sequence: str
        """

        if sequence:
            try:
                sequence = ''.join(sequence.split())
                _ = sequence.isalpha()
            except AttributeError:
                raise TypeError('sequence must be a string')
            else:
                if not _:
                    raise ValueError('not a valid protein sequence')
        self._sequence = sequence

        tree = ET.fromstring(results)
        header = tree.getchildren()[0]
        parameters = dictElement(header.getchildren()[2], '{http://www.ebi.ac.uk/schema}')
        query_len = int(parameters.items()[0][1].getchildren()[0].items()[0][1])

        self._params = np.sum([parameters.items()[0][1].getchildren()[0].items(), \
                               parameters.items()[1:-2], \
                               parameters.items()[-2][1].getchildren()[0].items(), \
                               [(parameters.items()[-1])]])

        result = tree.getchildren()[1] 
        hits = []
        for hit in result.getchildren()[0]:
            alignments = dictElement(hit.getchildren()[0], '{http://www.ebi.ac.uk/schema}')
            data = dictElement(alignments['alignment'], '{http://www.ebi.ac.uk/schema}')
 
            for key in ['gaps', 'score']:
                if key in data.keys():
                    data[key] = int(data[key])
                else:
                    data[key] = 0

            data['query-len'] = query_len

            for key in ['bits', 'expectation', 'identity']:
                data[key] = float(data[key])

            p_identity = data['identity'] 
            data['percent_identity'] = p_identity
            p_overlap = (100.0 * (len(data['querySeq']) - data['gaps']) /
                         query_len)
            data['percent_coverage'] = p_overlap

            pdbch = dict(data)
            pdbch['pdb_id'] = hit.items()[0][1][:4]
            pdbch['chain_id'] = hit.items()[0][1][-1]
            pdbch['description'] = hit.items()[2][1]
            hits.append((p_identity, p_overlap, pdbch))

        hits.sort(key=lambda hit: hit[0], reverse=True)
        self._hits = hits
        
        if sequence and len(sequence) != query_len:
            raise ValueError('xml sequence length and the length of the provided '
                             'sequence do not match')
        
    def getSequence(self):
        """Returns the query sequence that was used in the search."""

        return self._sequence

    def getParameters(self):
        """Returns parameters used in blast search."""

        return dict(self._param)

    def getHits(self, percent_identity=0., percent_overlap=0., chain=False):
        """Returns a dictionary in which PDB identifiers are mapped to structure
        and alignment information.

        :arg percent_identity: PDB hits with percent sequence identity equal
            to or higher than this value will be returned, default is ``0.``
        :type percent_identity: float
        :arg percent_overlap: PDB hits with percent coverage of the query
          sequence equivalent or better will be returned, default is ``0.``
        :type percent_overlap: float
        :arg chain: if chain is **True**, individual chains in a PDB file
          will be considered as separate hits , default is **False**
        :type chain: bool"""

        assert isinstance(percent_identity, (float, int)), \
            'percent_identity must be a float or an integer'
        assert isinstance(percent_overlap, (float, int)), \
            'percent_overlap must be a float or an integer'
        assert isinstance(chain, bool), 'chain must be a boolean'

        hits = {}
        for p_identity, p_overlap, hit in self._hits:
            if p_identity < percent_identity:
                break
            if p_overlap < percent_overlap:
                continue
            if chain:
                key = (hit['pdb_id'], hit['chain_id'])
            else:
                key = hit['pdb_id']
            if not key in hits:
                hits[key] = hit
        return hits

    def getBest(self):
        """Returns a dictionary containing structure and alignment information
        for the hit with highest sequence identity."""

        return dict(self._hits[0][2])
 
