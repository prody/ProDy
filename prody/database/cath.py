
from os import getcwd, remove
from os.path import sep as pathsep
from os.path import isdir, isfile, join, split, splitext, normpath

import re
from numbers import Integral

from prody import LOGGER, SETTINGS, PY3K
from prody.utilities import makePath, gunzip, relpath, copyFile, openURL
from prody.utilities import sympath, indentElement, isPDB

if PY3K:
    import urllib.parse as urllib
    import urllib.request as urllib2
else:
    import urllib
    import urllib2

import xml.etree.ElementTree as ET

__all__ = ['CATHDB']

CATH_URL = 'http://download.cathdb.info'
CATH_NAME_PATH = '/cath/releases/latest-release/cath-classification-data/cath-names.txt'
CATH_DOMAIN_PATH = '/cath/releases/latest-release/cath-classification-data/cath-domain-list.txt'
CATH_BOUNDARY_PATH = '/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping.txt'

isCATH = re.compile('^[0-9.]+$').match
# Note that if the chain id is blank, it will be represented by "0" in CATH domain ID.
isDomain = re.compile('^[A-Za-z0-9]{5}[0-9]{2}$').match

class CATHElement(ET.Element):
    def __init__(self, tag, attrib={}, **extra):
        super(CATHElement, self).__init__(tag, attrib, **extra)

    def __getitem__(self, index):
        if isinstance(index, Integral):
            return super(CATHElement, self).__getitem__(index)
        else:
            return self.attrib[index]

    def __repr__(self):
        if isDomain:
            return '<CATHElement: {0}>'.format(self.id)
        return '<CATHElement: {0} ({1} members)>'.format(self.id, len(self))

    def isDomain(self):
        return isDomain(self.id)

    def getID(self):
        return self.tag

    def setID(self, value):
        self.tag = value
    
    id = property(getID, setID)


class CATHDB(ET.ElementTree):
    def __init__(self, source=CATH_URL):
        self.root = CATHElement('root')
        super(CATHDB, self).__init__(self.root)
        self._source = source
        self._raw_names = None
        self._raw_domains = None
        self._raw_boundaries = None
        self._result = None

        if source is not None:
            self.update()
    
    def reset(self):
        # by removing the immediate children of root, 
        # all the descendants will be removed
        for child in self.root:
            self.root.remove(child)

    def update(self, source=None):
        """Update data and files from CATH."""

        self._source = source = self._source or source
        self.reset()
        if source is None:
            return

        LOGGER.timeit('_cath_update')
        
        type_ = 0
        tree = None
        if isinstance(source, str):
            if isfile(source):
                type_ = 1
            else:
                try:
                    tree = ET.fromstring(source)
                    type_ = 2
                except:
                    type_ = 0
        elif hasattr(source, 'read'):
                type_ = 1
        else:
            raise TypeError('source must be either an url, file name, file handle, '
                            'or text in xml format')

        if type_ == 0:
            LOGGER.info('Fetching data from CATH...')
            self._fetch()

            LOGGER.info('Parsing CATH files...')
            self._parse()
        elif type_ == 1:
            tree = ET.parse(source)

        # post-processing
        if type_ > 0:
            root = tree.getroot()
            nodes = root.findall('.//*')

            # remove prefix from node tags
            for node in nodes:
                node.tag = node.tag.lstrip('id.')

            # convert int to str
            length_nodes = root.findall('.//*[@length]')
            for node in length_nodes:
                node.attrib['length'] = int(node.attrib['length'])
            
            copy2(root, self.root)

        LOGGER.report('CATH local database built in %.2fs.', '_cath_update')

    def copy(self):
        root = CATHElement('root')
        db = CATHDB(source=None)
        copy2(self.root, root)
        return db

    def save(self, file='cath.xml'):
        """Write local CATH database to an XML file. *file* can either be a 
        file name or a handle."""

        LOGGER.timeit('_cath_write')
        LOGGER.info('Writing data to xml...')

        if not len(self.root):
            raise ValueError('local database has not been built, '
                             'please call update() first')

        tree = self.copy()
        root = tree.getroot()

        # convert int to str
        length_nodes = root.findall('.//*[@length]')
        for node in length_nodes:
            node.attrib['length'] = str(node.attrib['length'])
        
        # add prefix to node tags
        nodes = root.findall('.//*')
        for node in nodes:
            node.tag = 'id.' + node.tag

        # add indentation to nodes
        indentElement(root)

        if isinstance(file, str):
            file = open(file, 'wb')
        tree.write(file)
        file.close()

        LOGGER.report('CATH local database saved in %.2fs.', '_cath_write')

    def get(self, identifier='root'):
        """Get an :class:`xml.etree.Element` object given a tag name."""

        if len(self.root) is None:
            raise ValueError('local database has not been built, '
                             'please call update() first')
        if identifier != 'root':
            iter = self.iter(identifier)
            data = next(iter)
        else:
            data = self.root
        return data

    def _fetch(self):
        """Download CATH files via HTTP."""

        try:
            response = urllib2.urlopen(CATH_URL + CATH_NAME_PATH)
            name_data = response.read()
        except Exception as error:
            raise type(error)('HTTP connection problem when extracting CATH names: '\
                              + repr(error))

        try:
            response = urllib2.urlopen(CATH_URL + CATH_DOMAIN_PATH)
            domain_data = response.read()
        except Exception as error:
            raise type(error)('HTTP connection problem when extracting domain information: '\
                              + repr(error))

        try:
            response = urllib2.urlopen(CATH_URL + CATH_BOUNDARY_PATH)
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
        
    def _parse(self):
        """Parse CATH files into a :class:`~xml.etree.ElementTree` structure. 
        The leaf nodes are tagged by domain names, e.g. 3hgnA02, and other nodes 
        are tagged by their full CATH identifiers, e.g. 2.40.10.10. The parsing 
        process assumes the domain file is organized in a topdown fashion, i.e. 
        1 -> 1.20 -> 1.20.120 -> 1.20.120.1000. """
        
        name_data = self._raw_names 
        domain_data = self._raw_domains
        boundary_data = self._raw_boundaries

        root = self.root

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
            #pdb_id = domain_id[:4]
            #chain = domain_id[4]

            level = getLevel(cath_id)
            attrib = {'cath': cath_id,
                      'rep': domain_id,
                      'name': cath_name,}
                      #'pdb': pdb_id,
                      #'chain': chain}
            
            if level == 1:
                node = CATHElement(cath_id, attrib=attrib)
                root.append(node)
            else:
                parent_id = getParent(cath_id)
                parent = findNodeByCATH(root, parent_id)
                if parent is None:
                    LOGGER.warn('error encountered when building the tree: {0}'
                                .format(cath_id))
                    continue
                node = CATHElement(cath_id, attrib=attrib)
                parent.append(node)

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
            if chain == '0':
                chain = ''
            
            cath_id = '.'.join(items[1:5])
            #S35, S60, S95, S100 = items[5:9]
            length = int(items[10])
            #resolution = float(items[11])

            attrib = {'cath': cath_id,
                      'pdb': pdb_id,
                      'chain': chain,
                      'length': length,
                      'range': ''}

            parent = findNodeByCATH(root, cath_id)
            if parent is None:
                LOGGER.warn('error encountered when assigning domains: {0}'
                            .format(domain_id))
                continue
            node = CATHElement(domain_id, attrib=attrib)
            parent.append(node)

        # parsing the boundary file
        lines = boundary_data.splitlines()
        for line in lines:
            # Example: 3hgnA02	13-111,228-240
            text = line.strip()
            if not text:
                continue
            if text.startswith('#'):
                continue
            items = text.split()

            domain_id, resrange = items

            node = findNodeByCATH(root, cath_id)
            if node is None:
                LOGGER.warn('error encountered when assigning domain boundaries: {0}'
                            .format(domain_id))
                continue
            node.attrib['range'] = resrange
    
    def search(self, identifier):
        """Search the *identifier* against the local CATH database and return a list 
        of :class:`~xml.etree.ElementTree` instances. *identifier* can be a 4-digit 
        PDB ID plus an 1-digit chain ID (optional), 7-digit domain ID, or CATH ID."""

        ret = None
        if isCATH(identifier):
            ret = self.get(identifier)
        elif isDomain(identifier):
            identifier = identifier[:4].lower() + identifier[-3:]
            ret = self.get(identifier)
        elif isPDB(identifier):
            pdb_id = identifier[:4].lower()
            chain = identifier[4:]
            for gap in ' -_:':
                chain = chain.replace(gap, '')
            
            root = self.get()
            nodes = root.findall('.//*[@pdb="%s"]'%pdb_id)
            if chain.strip() != '':
                for i in range(len(nodes)-1, -1, -1):
                    node = nodes[i]
                    if node.attrib['chain'] != chain:
                        del nodes[i]
            ret = nodes
        else:
            raise ValueError('identifier can be a 4-digit PDB ID plus 1-digit chain ID '
                             '(optional), 7-digit domain ID, or CATH ID.')
        
        self._result = ret
        return ret

    def getCATH(self, node=None):
        if node is None:
            node = self._result
        return node.attrib['cath']

# Note: all the following functions should be kept hidden from users
def splitCATH(cath_id, cumul=False):
    fields = cath_id.split('.')

    if cumul:
        fields_ = list(fields)
        for i in range(len(fields)):
            fields_[i] = '.'.join(fields[:i+1])
        fields = fields_
    return fields

def findNodeByCATH(root, cath_id):
    fields = splitCATH(cath_id, cumul=True)

    if len(fields) == 3:
        pass
    node = root
    for field in fields:
        node = node.find(field)
        if node is None:
            return None
    return node

def getParent(cath_id):
    fields = splitCATH(cath_id)
    parent = '.'.join(f for f in fields[:-1])
    return parent

def getLevel(cath_id):
    return len(splitCATH(cath_id))

def isLeaf(cath_id):
    level = getLevel(cath_id)

    if level > 4:
        raise ValueError('invalid cath_id: {0}'.format(cath_id))

    return level == 4

def copy2(a, b):
    """Copy a tree structure from *a* to *b*"""
    for c in a:
        node = CATHElement(c.tag, c.attrib)
        b.append(node)
        copy2(c, node)
