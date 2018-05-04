
from os import getcwd, remove
# from glob import glob
from os.path import sep as pathsep
from os.path import isdir, isfile, join, split, splitext, normpath

from prody import LOGGER, SETTINGS, PY3K
from prody.utilities import makePath, gunzip, relpath, copyFile, openURL
from prody.utilities import sympath

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

import xml.etree.ElementTree as ET

__all__ = ['CATHDB', 'fetchCATH', 'buildCATHNameDict', 'buildPDBChainCATHDict']

class CATHDB(object):
    def __init__(self, url=None, name_path=None, domain_path=None, boundary_path=None):
        self._url = url or 'http://download.cathdb.info'
        self._name_path = name_path or \
            '/cath/releases/latest-release/cath-classification-data/cath-names.txt'
        
        self._boundary_path = boundary_path or \
            '/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping.txt'

        self._domain_path = domain_path or \
            '/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt'

        self._raw_names = None
        self._raw_domains = None
        self._data = None

        self.update()
    
    def update(self):
        """Update data and files from CATH."""

        self._fetch()
        self._parse()

    def _fetch(self):
        """Download CATH files via HTTP."""

        try:
            response = urllib2.urlopen(self._url + self._name_path)
            name_data = response.read()
        except Exception as error:
            raise type(error)('HTTP connection problem when extracting CATH names: '\
                              + repr(error))

        try:
            response = urllib2.urlopen(self._url + self._domain_path)
            domain_data = response.read()
        except Exception as error:
            raise type(error)('HTTP connection problem when extracting domain information: '\
                              + repr(error))

        try:
            response = urllib2.urlopen(self._url + self._boundary_path)
            boundary_data = response.read()
        except Exception as error:
            raise type(error)('HTTP connection problem when extracting domain boundaries: '\
                              + repr(error))

        if PY3K:
            name_data = name_data.decode()
            domain_data = domain_data.decode()
            boundary_data = boundary_data.decode()

        self._raw_names = name_data
        self._raw_domains = domain_data
        self._raw_boundaries = boundary_data

        self._data = None
        
    def _parse(self):
        """Parse CATH files."""
        
        name_data = self._raw_names 
        domain_data = self._raw_domains
        boundary_data = self._raw_boundaries

        root = ET.Element('root')

        # parsing the name file
        lines = name_data.splitlines()
        for line in lines:
            text = line.strip()
            if not text:
                continue
            if text.startswith('#'):
                continue
            items = text.split('    ')

            cath_id, domain_id, cath_name = items
            pdb_id = domain_id[:4]
            chain = domain_id[4]

            level = getLevel(cath_id)
            attrib = {'rep': domain_id,
                      'name': cath_name,
                      'pdb': pdb_id,
                      'chain': chain}
            if level == 1:
                node = ET.SubElement(root, cath_id, attrib=attrib)
            else:
                parent_id = getParent(cath_id)
                parent = root.find(parent_id)
                if parent is None:
                    LOGGER.warn('error encountered when building the tree: {0}'.format(cath_id))
                    continue
                node = ET.SubElement(parent, cath_id, attrib=attrib)

        # parsing the domain file
        lines = domain_data.splitlines()
        for line in lines:
            # Column 1:  CATH domain name (seven characters)
            # Column 2:  Class number
            # Column 3:  Architecture number
            # Column 4:  Topology number
            # Column 5:  Homologous superfamily number
            # Column 6:  S35 sequence cluster number
            # Column 7:  S60 sequence cluster number
            # Column 8:  S95 sequence cluster number
            # Column 9:  S100 sequence cluster number
            # Column 10: S100 sequence count number
            # Column 11: Domain length
            # Column 12: Structure resolution (Angstroms)
            #            (999.000 for NMR structures and 1000.000 for obsolete PDB entries)
            text = line.strip()
            if not text:
                continue
            if text.startswith('#'):
                continue
            items = text.split()

            domain_id = items[0]
            pdb_id = domain_id[:4]
            chain = domain_id[4]
            
            cath_id = '.'.join(items[1:5])
            #S35, S60, S95, S100 = items[5:9]
            length = int(items[10])
            #resolution = float(items[11])

            attrib = {'cath': cath_id,
                      'pdb': pdb_id,
                      'chain': chain,
                      'length': length,
                      'range': ''}

            parent = root.find(cath_id)
            if parent is None:
                LOGGER.warn('error encountered when assigning domains: {0}'.format(domain_id))
                continue
            node = ET.SubElement(root, domain_id, attrib=attrib)

        # parsing the boundary file
        lines = domain_data.splitlines()
        for line in lines:
            # Example: 3hgnA02	13-111,228-240
            text = line.strip()
            if not text:
                continue
            if text.startswith('#'):
                continue
            items = text.split()

            domain_id, resrange = items

            node = root.find(cath_id)
            if node is None:
                LOGGER.warn('error encountered when assigning domain boundaries: {0}'.format(domain_id))
                continue
            node.attrib['range'] = resrange
        
        self._data = root

def getParent(cath_id):
    fields = cath_id.split('.')
    parent = '.'.join(f for f in fields[:-1])
    return parent

def getLevel(cath_id):
    return len(cath_id.split('.'))

def isLeaf(cath_id):
    level = getLevel(cath_id)

    if level > 4:
        raise ValueError('invalid cath_id: {0}'.format(cath_id))

    return level == 4

def fetchCATH(filename, ftp_host=None, ftp_path=None, **kwargs):
    """Downloads CATH file via FTP."""
    if ftp_host == None:
        ftp_host = 'orengoftp.biochem.ucl.ac.uk'
    if ftp_path == None:
        ftp_path = '/cath/releases/daily-release/newest/'
    from ftplib import FTP
    output_folder = kwargs.pop('folder', None)
    ftp_fn = filename
    try:
        ftp = FTP(ftp_host)
    except Exception as error:
        raise type(error)('FTP connection problem, potential reason: '
                          'no internet connectivity')
    else:
        success = 0
        failure = 0
        filenames = []
        ftp.login('')
        
        data = []
        try:
            ftp.cwd(ftp_path)
            ftp.retrbinary('RETR ' + ftp_fn, data.append)
        except Exception as error:
            if ftp_fn in ftp.nlst():
                LOGGER.warn('{0} download failed ({1}). It is '
                            'possible that you do not have rights to '
                            'download .gz files in the current network.'
                            .format(ftp_fn, str(error)))
            else:
                LOGGER.warn('{0} download failed. {1} does not exist '
                            'on {2}.'.format(ftp_fn, ftp_fn, ftp_host))
            failure += 1
            filenames.append(None)
        else:
            if len(data):
                if output_folder is None:
                    output_folder = getcwd()
                    filename_full = join(output_folder, ftp_fn)

                    with open(filename_full, 'w+b') as pdbfile:
                        write = pdbfile.write
                        [write(block) for block in data]

                    filename_full = normpath(relpath(filename_full))
                    LOGGER.debug('{0} downloaded ({1})'
                                    .format(ftp_fn, sympath(filename_full)))
                    success += 1
                    filenames.append(filename_full)
                else:
                    LOGGER.warn('{0} download failed, reason unknown.'
                                .format(ftp_fn))
                    failure += 1
                    filenames.append(None)
        ftp.quit()


def buildCATHNameDict(cath_file, iscommpressed=True):
    """Returns a dictionary for CATH names with key of CATH ID."""
    if iscommpressed:
        gunzip(cath_file, 'cath_b.names.temp')
        cath_file = 'cath_b.names.temp'
        
    cath_id2name = dict()
    with open(cath_file, 'r') as file_temp:
        for line in file_temp:
            ind_temp = line.find(' ')
            cath_id2name[line[:ind_temp]] = line[ind_temp:].strip()
    if iscommpressed:
        remove(cath_file) 
    return cath_id2name
    
def buildPDBChainCATHDict(cath_file, iscommpressed=True):
    """Returns a dictionary for CATH info (ID and version) with key of PDB chain."""
    if iscommpressed:
        gunzip(cath_file, 'cath_b.all.temp')
        cath_file = 'cath_b.all.temp'
    
    cath_dict_temp = dict()
    cath_i_dict = dict()
    with open(cath_file, 'r') as file_temp:
        for line in file_temp:
            line = line.strip()
            if line != '':
                line_list = line.split(' ')
                cath_dict_temp[line_list[0]] = line_list[1:]
                key, value = line[0:5], line[5:7]
                if key in cath_i_dict:
                    cath_i_dict[key].append(value)
                else:
                    cath_i_dict[key] = [value]
    pdbChain2CATH = dict()
    for key, values in cath_i_dict.items():
        pdbChain2CATH[key] = []
        for v in values:
            pdbChain2CATH[key].append(cath_dict_temp[key+v])
    if iscommpressed:
        remove(cath_file) 
    return pdbChain2CATH
    
    
# ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/daily-release/newest/

# fetchCATH('cath-b-newest-names.gz')
# cath_id2name = buildCATHNameDict('cath-b-newest-names.gz')

# fetchCATH('cath-b-newest-all.gz')
# pdbChain2CATH = buildPDBChainCATHDict('cath-b-newest-all.gz')
