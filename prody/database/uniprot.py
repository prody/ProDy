from prody import LOGGER
from prody.proteins import parsePDB
from prody.utilities import dictElement, dictElementLoop, openURL

import re
from xml.etree.cElementTree import XML

__all__ = ['UniprotRecord', 'searchUniprot', 'queryUniprot']

comma_splitter = re.compile(r'\s*,\s*').split

class UniprotRecord(object):
    """This class provides a wrapper for UniProt data including functions 
    for accessing particular fields and parsing associated PDB entries."""
    def __init__(self, data):
        self._rawdata = data
        self._pdbids = []
        # self._selstrs = []

        self._parse()

    def __repr__(self):
        return '<UniprotRecord: %s>'%self.getTitle()

    def __str__(self):
        return self.getTitle()

    def setData(self, value):
        self._rawdata = value
        self._parse()

    def getData(self):
        return self._rawdata

    def getPDBs(self):
        return self._pdbids

    # def getSelstrs(self):
    #     return self._selstrs

    def getSequence(self, index=0):
        return self.getEntry('sequence', index)
    
    def getAccession(self, index=0):
        return self.getEntry('accession', index)
    
    def getName(self, index=0):
        return self.getEntry('name', index)

    def getTitle(self):
        uid = self.getAccession()
        name = self.getName()
        return '%s (%s)'%(uid, name)

    def getEntry(self, item, index=0):
        key = '%s%4d'%(item, index)
        if key in self._rawdata:
            return self._rawdata[key]
        else:
            raise KeyError('%s does not exist in the Uniprot record'%key)

    def _parse(self):
        data = self._rawdata
        PDBIDs = []
        # SELSTRs = []
        for key, value in data.items():
            if not key.startswith('dbReference'):
                continue
            try:
                pdbid = value['PDB']
            except (KeyError, TypeError) as e:
                continue
            pdbchains = value['chains']
            
            # example chain strings: "A=27-139, B=140-150" or "A/B=27-150"
            chains = []
            ranges = []
            pdbchains = comma_splitter(pdbchains)
            for chain in pdbchains:
                chids, resrange = chain.split('=')
                chids = [chid.strip() for chid in chids.split('/')]
                resrange = resrange.split('-')

                for chid in chids:
                    chains.append(chid)
                    ranges.append(resrange)
            
            for chid, rng in zip(chains, ranges):
                pdbchid = pdbid + chid if chid != '@' else pdbid
                PDBIDs.append(pdbchid)
                # SELSTRs.append('resnum %s to %s'%tuple(rng))
        
        self._pdbids = PDBIDs
        # self._selstrs = SELSTRs

    def parsePDBs(self, **kwargs):
        """Load PDB into memory as :class:`.AtomGroup` instances using :func:`.parsePDB` and 
        perform selection based on residue ranges given by CATH."""
        
        pdbs = self.getPDBs()
        # selstrs = self.getSelstrs()
        header = kwargs.get('header', False)
        model = kwargs.get('model', None)

        LOGGER.timeit('_uniprot_parsePDB')
        LOGGER.info('Parsing {0} PDB files...'.format(len(pdbs)))
        ret = parsePDB(*pdbs, **kwargs)

        if model != 0:
            headers = None
            if header:
                prots, headers = ret
            else:
                prots = ret

            if not isinstance(prots, list):
                prots = [prots]

                if header:
                    headers = [headers]
                    ret = (prots, headers)
                else:
                    ret = prots
                    
            LOGGER.info('Extracting domains...')
            # for i in range(len(prots)):
            #     sel = prots[i].select(selstrs[i])
            #     prots[i] = sel
        LOGGER.report('Uniprot domains are parsed and extracted in %.2fs', '_uniprot_parsePDB')

        return ret

def queryUniprot(id, expand=[], regex=True):
    """Query Uniprot with *id* and return a `dict` containing the raw results. 
    Regular users should use :func:`searchUniprot` instead.
    
    :arg expand: entries through which you want to loop dictElements
        until there aren't any elements left
    :type expand: list
    """

    if not isinstance(id, str):
        raise TypeError('id should be a string')

    try:
        record_file = openURL('http://www.uniprot.org/uniprot/{0}.xml'.format(id))
    except:
        raise ValueError('No Uniprot record found with that id')
    
    data = record_file.read()
    record_file.close()
    data = XML(data)

    data = dictElement(data[0], '{http://uniprot.org/uniprot}', number_multiples=True)

    for key in data:
        value = data[key]
        if not key.startswith('dbReference'):
            continue
        
        try:
            if value.get('type') != 'PDB':
                continue
        except AttributeError:
            continue

        pdbid = value.get('id')
        refdata = {'PDB': pdbid}
        for prop in value:
            prop_key = prop.get('type')
            prop_val = prop.get('value')
            refdata[prop_key] = prop_val
        data[key] = refdata
            
    if expand:
        keys = []
        if regex:
            for lt in expand:
                lt_re = re.compile(lt)
                for key in data:
                    if lt_re.match(key):
                        keys.append(key)
        else:
            keys = expand
        data = dictElementLoop(data, keys, '{http://uniprot.org/uniprot}')
    
    return data

def searchUniprot(id):
    """Search Uniprot with *id* and return a :class:`UniprotRecord` containing the results. 
    """

    data = queryUniprot(id)
    
    return UniprotRecord(data)
