# -*- coding: utf-8 -*-
"""This module defines functions for parsing header data from PDB files."""

from collections import defaultdict
import os.path

import numpy as np

from prody import LOGGER
from prody.atomic import ATOMIC_FIELDS
from prody.atomic import Atomic, AtomGroup
from prody.atomic import getSequence
from prody.measure import Transformation
from prody.utilities import openFile, decToHybrid36

from .localpdb import fetchPDB

__all__ = ['Chemical', 'Polymer', 'DBRef', 'parsePDBHeader',
           'assignSecstr', 'buildBiomolecules']


class Chemical(object):

    """A data structure for storing information on chemical components
    (or heterogens) in PDB structures.

    A :class:`Chemical` instance has the following attributes:

    ===========  =====  =======================================================
    Attribute    Type   Description (RECORD TYPE)
    ===========  =====  =======================================================
    resname      str    residue name (or chemical component identifier) (HET)
    name         str    chemical name (HETNAM)
    chain        str    chain identifier (HET)
    resnum       int    residue (or sequence) number (HET)
    icode        str    insertion code (HET)
    natoms       int    number of atoms present in the structure (HET)
    description  str    description of the chemical component (HET)
    synonyms     list   synonyms (HETSYN)
    formula      str    chemical formula (FORMUL)
    pdbentry     str    PDB entry that chemical data is extracted from
    ===========  =====  =======================================================

    Chemical class instances can be obtained as follows:

    .. ipython:: python

       from prody import *
       chemical = parsePDBHeader('1zz2', 'chemicals')[0]
       chemical
       chemical.name
       chemical.natoms
       len(chemical)"""

    __slots__ = ['resname', 'name', 'chain', 'resnum', 'icode',
                 'natoms', 'description', 'synonyms', 'formula', 'pdbentry']

    def __init__(self, resname):

        #: residue name (or chemical component identifier)
        self.resname = resname
        #: chemical name
        self.name = None
        #: chain identifier
        self.chain = None
        #: residue (or sequence) number
        self.resnum = None
        #: insertion code
        self.icode = None
        #: number of atoms present in the structure
        self.natoms = None
        #: description of the chemical component
        self.description = None
        #: list of synonyms
        self.synonyms = None
        #: chemical formula
        self.formula = None
        #: PDB entry that chemical data is extracted from
        self.pdbentry = None

    def __str__(self):
        return self.resname

    def __repr__(self):
        return '<Chemical: {0} ({1}_{2}_{3})>'.format(self.resname,
                                                      self.pdbentry,
                                                      self.chain, self.resnum)

    def __len__(self):
        return self.natoms


class Polymer(object):

    """A data structure for storing information on polymer components
    (protein or nucleic) of PDB structures.

    A :class:`Polymer` instance has the following attributes:

    ==========  ======  ======================================================
    Attribute   Type    Description (RECORD TYPE)
    ==========  ======  ======================================================
    chid        str     chain identifier
    name        str     name of the polymer (macromolecule) (COMPND)
    fragment    str     specifies a domain or region of the molecule (COMPND)
    synonyms    list    synonyms for the polymer (COMPND)
    ec          list    associated Enzyme Commission numbers (COMPND)
    engineered  bool    indicates that the polymer was produced using
                        recombinant technology or by purely chemical synthesis
                        (COMPND)
    mutation    bool    indicates presence of a mutation (COMPND)
    comments    str     additional comments
    sequence    str     polymer chain sequence (SEQRES)
    dbrefs      list    sequence database records (DBREF[1|2] and SEQADV),
                        see :class:`DBRef`
    modified    list    | modified residues (MODRES)
                        | when modified residues are present, each will be
                          represented as: ``(resname, chid, resnum, icode, stdname,
                          comment)``
    pdbentry    str     PDB entry that polymer data is extracted from
    ==========  ======  ======================================================

    Polymer class instances can be obtained as follows:

    .. ipython:: python

       polymer = parsePDBHeader('2k39', 'polymers')[0]
       polymer
       polymer.pdbentry
       polymer.chid
       polymer.name
       polymer.sequence
       len(polymer.sequence)
       len(polymer)
       dbref = polymer.dbrefs[0]
       dbref.database
       dbref.accession
       dbref.idcode"""

    __slots__ = ['chid', 'name', 'fragment', 'synonyms', 'ec',
                 'engineered', 'mutation', 'comments', 'sequence', 'pdbentry',
                 'dbrefs', 'modified']

    def __init__(self, chid):

        #: chain identifier
        self.chid = chid
        #: name of the polymer (macromolecule)
        self.name = ''
        #: specifies a domain or region of the molecule
        self.fragment = None
        #: list of synonyms for the molecule
        self.synonyms = None
        #: list of associated Enzyme Commission numbers
        self.ec = None
        self.engineered = None
        """indicates that the molecule was produced using recombinant
        technology or by purely chemical synthesis"""
        #: sequence database reference records
        self.dbrefs = []
        #: indicates presence of a mutation
        self.mutation = None
        #: additional comments
        self.comments = None
        #: polymer chain sequence
        self.sequence = ''
        #: modified residues
        self.modified = None
        #: PDB entry that polymer data is extracted from
        self.pdbentry = None

    def __str__(self):
        return self.name

    def __repr__(self):
        return '<Polymer: {0} ({1}_{2})>'.format(self.name,
                                                 self.pdbentry, self.chid)

    def __len__(self):
        return len(self.sequence)

_PDB_DBREF = {
    'GB': 'GenBank',
    'PDB': 'PDB',
    'UNP': 'UniProt',
    'NORINE': 'Norine',
    'UNIMES': 'UNIMES',
    'EMDB': 'EMDB',
    'BMRB': 'BMRB'
}


class DBRef(object):

    """A data structure for storing reference to sequence databases for polymer
    components in  PDB structures.  Information if parsed from **DBREF[1|2]**
    and **SEQADV** records in PDB header."""

    __slots__ = ['database', 'dbabbr', 'idcode', 'accession',
                 'first', 'last', 'diff']

    def __init__(self):

        #: sequence database, one of UniProt, GenBank, Norine, UNIMES, or PDB
        self.database = None
        #: database abbreviation, one of UNP, GB, NORINE, UNIMES, or PDB
        self.dbabbr = None
        #: database identification code, i.e. entry name in UniProt
        self.idcode = None
        #: database accession code
        self.accession = None
        #: initial residue numbers, ``(resnum, icode, dbnum)``
        self.first = None
        #: ending residue numbers, ``(resnum, icode, dbnum)``
        self.last = None
        self.diff = []
        """list of differences between PDB and database sequences,
        ``(resname, resnum, icode, dbResname, dbResnum, comment)``"""

    def __str__(self):
        return self.accession

    def __repr__(self):
        return '<DBRef: {0} ({1})>'.format(self.accession, self.database)

_START_COORDINATE_SECTION = set(['ATOM  ', 'MODEL ', 'HETATM'])


def cleanString(string, nows=False):
    """*nows* is no white space."""

    if nows:
        return ''.join(string.strip().split())
    else:
        return ' '.join(string.strip().split())


def parsePDBHeader(pdb, *keys, **kwargs):
    """Returns header data dictionary for *pdb*.  This function is equivalent to
    ``parsePDB(pdb, header=True, model=0, meta=False)``, likewise *pdb* may be
    an identifier or a filename.

    List of header records that are parsed.

    ============ ================= ============================================
    Record type  Dictionary key(s)  Description
    ============ ================= ============================================
    HEADER       | classification  | molecule classification
                 | deposition_date | deposition date
                 | identifier      | PDB identifier
    TITLE        title             title for the experiment or analysis
    SPLIT        split             list of PDB entries that make up the whole
                                   structure when combined with this one
    COMPND       polymers          see :class:`Polymer`
    EXPDTA       experiment        information about the experiment
    NUMMDL       n_models          number of models
    MDLTYP       model_type        additional structural annotation
    AUTHOR       authors           list of contributors
    JRNL         reference         reference information dictionary:
                                     * *authors*: list of authors
                                     * *title*: title of the article
                                     * *editors*: list of editors
                                     * *issn*:
                                     * *reference*: journal, vol, issue, etc.
                                     * *publisher*: publisher information
                                     * *pmid*: pubmed identifier
                                     * *doi*: digital object identifier
    DBREF[1|2]   polymers          see :class:`Polymer` and :class:`DBRef`
    SEQADV       polymers          see :class:`Polymer`
    SEQRES       polymers          see :class:`Polymer`
    MODRES       polymers          see :class:`Polymer`
    HELIX        polymers          see :class:`Polymer`
    SHEET        polymers          see :class:`Polymer`
    HET          chemicals         see :class:`Chemical`
    HETNAM       chemicals         see :class:`Chemical`
    HETSYN       chemicals         see :class:`Chemical`
    FORMUL       chemicals         see :class:`Chemical`
    REMARK   2   resolution        resolution of structures, when applicable
    REMARK   4   version           PDB file version
    REMARK 350   biomoltrans       biomolecular transformation lines
                                   (unprocessed)
    REMARK 900	 related_entries   related entries in the PDB or EMDB
    ============ ================= ============================================

    Header records that are not parsed are: OBSLTE, CAVEAT, SOURCE, KEYWDS,
    REVDAT, SPRSDE, SSBOND, LINK, CISPEP, CRYST1, ORIGX1, ORIGX2, ORIGX3,
    MTRIX1, MTRIX2, MTRIX3, and REMARK X not mentioned above."""

    if not os.path.isfile(pdb):
        if len(pdb) == 4 and pdb.isalnum():
            filename = fetchPDB(pdb)
            if filename is None:
                raise IOError('PDB file for {0} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    pdb = openFile(pdb, 'rt')
    header, _ = getHeaderDict(pdb, *keys, **kwargs)
    pdb.close()
    return header


def getHeaderDict(stream, *keys, **kwargs):
    """Returns header data in a dictionary.  *stream* may be a list of PDB lines
    or a stream.
    
    Polymers have sequences that usually use one-letter residue name abbreviations by default. 
    To obtain long (usually three letter) abbrevations, set *longSeq* to **True**."""

    lines = defaultdict(list)
    loc = 0
    for loc, line in enumerate(stream):
        startswith = line[0:6]
        if startswith in _START_COORDINATE_SECTION:
            break
        lines[startswith].append((loc, line))
    if not loc:
        #raise ValueError('empty PDB file or stream')
        return None, loc
    for i, line in lines['REMARK']:
        lines[line[:10]].append((i, line))

    pdbid = _PDB_HEADER_MAP['identifier'](lines)
    lines['pdbid'] = pdbid
    if keys:
        keys = list(keys)
        for k, key in enumerate(keys):
            if key in _PDB_HEADER_MAP:
                if key == 'polymers':
                    value = _PDB_HEADER_MAP[key](lines, **kwargs)
                else:
                    value = _PDB_HEADER_MAP[key](lines)
                keys[k] = value
            else:
                raise KeyError('{0} is not a valid header data identifier'
                               .format(repr(key)))
            if key in ('chemicals', 'polymers'):
                for component in value:
                    component.pdbentry = pdbid
        if len(keys) == 1:
            return keys[0], loc
        else:
            return tuple(keys), loc
    else:
        header = {}
        for key, func in _PDB_HEADER_MAP.items():  # PY3K: OK
            value = func(lines)
            if value is not None:
                header[key] = value
        for chem in header.get('chemicals', []):
            chem.pdbentry = pdbid
            header[chem.resname] = chem
        for poly in header.get('polymers', []):
            poly.pdbentry = pdbid
            header[poly.chid] = poly
        return header, loc


def _getBiomoltrans(lines):

    
    biomolecule = defaultdict(list)
    for i, line in lines['REMARK 350']:
        if line[11:23] == 'BIOMOLECULE:':
            currentBiomolecule = line.split()[-1]
            applyToChains = []
        elif line[11:41] == 'APPLY THE FOLLOWING TO CHAINS:' \
        or line[30:41] == 'AND CHAINS:':
            applyToChains.extend(line[41:].replace(' ', '')
                                 .strip().strip(',').split(','))
        elif line[13:18] == 'BIOMT':
            biomt = biomolecule[currentBiomolecule]
            if line[13:19] == 'BIOMT1':
                if applyToChains == []:
                    applyToChains = biomt[0]
                biomt.append(applyToChains)
            elif line[13:19]:
                applyToChains = []
            biomt.append(line[23:])
    return dict(biomolecule)


def _getResolution(lines):

    for i, line in lines['REMARK   2']:
        if 'RESOLUTION' in line:
            try:
                return float(line[23:30])
            except:
                return None

def _getSCALE(lines):
    ctof = np.identity(4)
    if len(lines['SCALE1']) == 0: return {}
    ctof[0] = lines['SCALE1'][0][1].split()[1:5]
    ctof[1] = lines['SCALE2'][0][1].split()[1:5]
    ctof[2] = lines['SCALE3'][0][1].split()[1:5]
    ftoc = np.linalg.inv(ctof)
    return {'ctof':ctof, 'ftoc':ftoc}

def _getSpaceGroup(lines):
    rowct = 0
    mat = []
    xform = []
    sg = None
    for i, line in lines['REMARK 290']:
        if 'SYMMETRY OPERATORS FOR SPACE GROUP:' in line:
            try:
                sg = line.split('GROUP:')[1].strip()
            except:
                pass
        if line[13:18] == 'SMTRY':
            w = line.split()
            mat.append([float(x) for x in w[4:]])
            rowct += 1
            if rowct==3:
                mat.append( (0,0,0,1) )
                matnp = np.array(mat)
                #if np.sum(matnp-np.identity(4))!=0.0:
                xform.append(np.array(mat))
                rowct = 0
                mat = []
    return {'spaceGroup':sg, 'symMats': xform}

def _getRelatedEntries(lines):

    dbrefs = []
    for i, line in lines['REMARK 900']:
        if 'RELATED ID' in line:
            dbref = DBRef()
            end_of_id = line.find('RELATED DB')
            dbref.accession = line[23:end_of_id].strip()
            dbref.dbabbr = line[end_of_id+12:end_of_id+16].strip()
            dbref.database = _PDB_DBREF.get(dbref.dbabbr, 'Unknown')
            
            dbrefs.append(dbref)

    return dbrefs


def _getSpaceGroup(lines):

    for i, line in lines['REMARK 290']:
        if 'SYMMETRY OPERATORS FOR SPACE GROUP:' in line:
            try:
                return line.split('GROUP:')[1].strip()
            except:
                return None


def _getHelix(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = {}
    for i, line in lines['HELIX ']:
        try:
            chid = line[19]
            #        helix class,      serial number,   identifier
            value = (int(line[38:40]), int(line[7:10]), line[11:14].strip())
        except:
            continue

        initICode = line[25]
        initResnum = int(line[21:25])
        if initICode != ' ':
            for icode in alphas[alphas.index(initICode):]:
                helix[(chid, initResnum, icode)] = value
            initResnum += 1
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            for icode in alphas[:alphas.index(endICode)+1]:
                helix[(chid, endResnum, icode)] = value
            endResnum -= 1
        for resnum in range(initResnum, endResnum+1):
            helix[(chid, resnum, '')] = value
    return helix


def _getHelixRange(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = []
    for i, line in lines['HELIX ']:
        try:
            chid = line[19]
            Hclass=int(line[38:40])
            Hnr=int(line[7:10])
        except:
            continue

        initResnum = int(line[21:25])
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            endResnum -= 1
        helix.append(['H', chid, Hclass, Hnr, initResnum, endResnum]) 
     
    return helix


def _getSheet(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = {}
    for i, line in lines['SHEET ']:
        try:
            chid = line[21]
                     # sense           # strand num     # sheet id
            value = (int(line[38:40]), int(line[7:10]), line[11:14].strip())
        except:
            continue

        initICode = line[26]
        initResnum = int(line[22:26])
        if initICode != ' ':
            for icode in alphas[alphas.index(initICode):]:
                sheet[(chid, initResnum, icode)] = value
            initResnum += 1
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            for icode in alphas[:alphas.index(endICode)+1]:
                sheet[(chid, endResnum, icode)] = value
            endResnum -= 1
        for resnum in range(initResnum, endResnum+1):
            sheet[(chid, resnum, '')] = value
    return sheet


def _getSheetRange(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = []
    for i, line in lines['SHEET ']:
        try:
            chid = line[21]
            dir = int(line[38:40])
            Snr = int(line[7:10])
        except:
            continue

        initICode = line[26]
        initResnum = int(line[22:26])
        if initICode != ' ':
            initResnum += 1
        endICode = line[37]
        endResnum = int(line[33:37])
        if endICode != ' ':
            endResnum -= 1
        sheet.append(['E', chid, dir, Snr, initResnum, endResnum])
    return sheet


def _getReference(lines):
    """Returns a reference of the PDB entry."""

    ref = {}
    title = ''
    authors = []
    editors = []
    reference = ''
    publisher = ''
    for i, line in lines['JRNL  ']:
        try:
            what = line.split(None, 2)[1]
        except:
            continue
        if what == 'AUTH':
            authors.extend(line[19:].strip().split(','))
        elif what == 'TITL':
            title += line[19:]
        elif what == 'EDIT':
            editors.extend(line[19:].strip().split(','))
        elif what == 'REF':
            reference += line[19:]
        elif what == 'PUBL':
            publisher += line[19:]
        elif what == 'REFN':
            ref['issn'] = line[19:].strip()
        elif what == 'PMID':
            ref['pmid'] = line[19:].strip()
        elif what == 'DOI':
            ref['doi'] = line[19:].strip()
    ref['authors'] = authors
    ref['title'] = cleanString(title)
    ref['editors'] = editors
    ref['reference'] = cleanString(reference)
    ref['publisher'] = cleanString(publisher)

    return ref


def _getPolymers(lines, **kwargs):
    """Returns list of polymers (macromolecules).
    
    Polymers have sequences that usually use one-letter residue name abbreviations by default. 
    To obtain long (usually three letter) abbrevations, set *longSeq* to **True**."""

    pdbid = lines['pdbid']
    polymers = dict()
    for i, line in lines['SEQRES']:
        ch = line[11]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly

        longSeq = kwargs.get('longSeq', False)
        if longSeq:
            if poly.sequence != '':
                poly.sequence += ' '
            poly.sequence += getSequence(line[19:].split(), **kwargs)
        else:
            poly.sequence += ''.join(getSequence(line[19:].split(), **kwargs))

    for i, line in lines['DBREF ']:
        i += 1

        ch = line[12]
        if ch == ' ':
            if not len(polymers) == 1:
                LOGGER.warn('DBREF chain identifier is not specified '
                            '({0}:{1})'.format(pdbid, i))
                continue
            else:
                ch = list(polymers)[0]
        dbabbr = line[26:32].strip()
        dbref = DBRef()
        dbref.dbabbr = dbabbr
        dbref.database = _PDB_DBREF.get(dbabbr, 'Unknown')
        dbref.accession = line[33:41].strip()
        dbref.idcode = line[42:54].strip()

        try:
            first = int(line[14:18])
        except:
            LOGGER.warn('DBREF for chain {2}: failed to parse '
                        'initial sequence number of the PDB sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        try:
            last = int(line[20:24])
        except:
            LOGGER.warn('DBREF for chain {2}: failed to parse '
                        'ending sequence number of the PDB sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        try:
            dbref.first = (first, line[18], int(line[56:60]))
        except:
            LOGGER.warn('DBREF for chain {2}: failed to parse '
                        'initial sequence number of the database sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        try:
            dbref.last = (last, line[24].strip(), int(line[62:67]))
        except:
            LOGGER.warn('DBREF for chain {2}: failed to parse '
                        'ending sequence number of the database sequence '
                        '({0}:{1})'.format(pdbid, i, ch))

        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.dbrefs.append(dbref)

    dbref1 = lines['DBREF1']
    dbref2 = lines['DBREF2']
    if len(dbref1) != len(dbref2):
        LOGGER.warn('DBREF1 and DBREF1 records are not complete')
        dbref12 = []
    else:
        dbref12 = zip(dbref1, dbref2)  # PY3K: OK

    for dbref1, dbref2 in dbref12:
        i, line = dbref1
        i += 1
        ch = line[12]

        dbabbr = line[26:32].strip()
        dbref = DBRef()
        dbref.dbabbr = dbabbr
        dbref.database = _PDB_DBREF.get(dbabbr, 'Unknown')
        dbref.idcode = line[47:67].strip()

        try:
            first = int(line[14:18])
        except:
            LOGGER.warn('DBREF1 for chain {2}: failed to parse '
                        'initial sequence number of the PDB sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        try:
            last = int(line[20:24])
        except:
            LOGGER.warn('DBREF1 for chain {2}: failed to parse '
                        'ending sequence number of the PDB sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        i, line = dbref2
        i += 1
        if line[12] == ' ':
            LOGGER.warn('DBREF2 chain identifier is not specified '
                        '({0}:{1})'.format(pdbid, ch))
        elif line[12] != ch:
            LOGGER.warn('DBREF1 and DBREF2 chain id mismatch'
                        '({0}:{1})'.format(pdbid, ch))

        dbref.accession = line[18:40].strip()
        try:
            dbref.first = (first, line[18].strip(), int(line[45:55]))
        except:
            LOGGER.warn('DBREF2 for chain {2}: failed to parse '
                        'initial sequence number of the database sequence '
                        '({0}:{1})'.format(pdbid, i, ch))
        try:
            dbref.last = (last, line[24].strip(), int(line[57:67]))
        except:
            LOGGER.warn('DBREF2 for chain {2}: failed to parse '
                        'ending sequence number of the database sequence '
                        '({0}:{1})'.format(pdbid, i, ch))

        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.dbrefs.append(dbref)

    for poly in polymers.values():  # PY3K: OK
        resnum = []
        for dbref in poly.dbrefs:
            dbabbr = dbref.dbabbr
            if dbabbr == 'PDB':
                if not (pdbid == dbref.accession == dbref.idcode):
                    LOGGER.warn('DBREF for chain {2} refers to PDB '
                                'entry {3} ({0}:{1})'
                                .format(pdbid, i, ch, dbref.accession))
            else:
                if pdbid == dbref.accession or pdbid == dbref.idcode:
                    LOGGER.warn('DBREF for chain {2} is {3}, '
                                'expected PDB ({0}:{1})'
                                .format(pdbid, i, ch, dbabbr))
                    dbref.database = 'PDB'
            resnum.append((dbref.first[0], dbref.last[0]))
        resnum.sort()
        last = -10000
        for first, temp in resnum:
            if first <= last:
                LOGGER.warn('DBREF records overlap for chain {0} ({1})'
                            .format(poly.chid, pdbid))
            last = temp

    for i, line in lines['MODRES']:
        ch = line[16]
        if ch == ' ':
            if not len(polymers) == 1:
                LOGGER.warn('MODRES chain identifier is not specified '
                            '({0}:{1})'.format(pdbid, i))
                continue
            else:
                ch = list(polymers)[0]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        if poly.modified is None:
            poly.modified = []
        poly.modified.append((line[12:15].strip(), line[16],
                              line[18:22].strip() + line[22].strip(), 
                              line[24:27].strip(),
                              line[29:70].strip()))

    for i, line in lines['SEQADV']:
        i += 1
        ch = line[16]
        if ch == ' ':
            if not len(polymers) == 1:
                LOGGER.warn('SEQADV chain identifier is not specified '
                            '({0}:{1})'.format(pdbid, i))
                continue
            else:
                ch = list(polymers)[0]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        dbabbr = line[24:28].strip()
        resname = line[12:15].strip()
        try:
            resnum = int(line[18:22].strip())
        except:
            #LOGGER.warn('SEQADV for chain {2}: failed to parse PDB sequence '
            #            'number ({0}:{1})'.format(pdbid, i, ch))
            continue
        icode = line[22].strip()
        try:
            dbnum = int(line[43:48].strip())
        except:
            #LOGGER.warn('SEQADV for chain {2}: failed to parse database '
            #            'sequence number ({0}:{1})'.format(pdbid, i, ch))
            continue            

        comment = line[49:70].strip()
        match = False
        for dbref in poly.dbrefs:
            if not dbref.first[0] <= resnum <= dbref.last[0]:
                continue
            match = True
            if dbref.dbabbr != dbabbr:
                LOGGER.warn('SEQADV for chain {2}: reference database '
                            'mismatch, expected {3} parsed {4} '
                            '({0}:{1})'.format(pdbid, i, ch,
                            repr(dbref.dbabbr), repr(dbabbr)))
                continue
            dbacc = line[29:38].strip()
            if dbref.accession[:9] != dbacc[:9]:
                LOGGER.warn('SEQADV for chain {2}: accession code '
                            'mismatch, expected {3} parsed {4} '
                            '({0}:{1})'.format(pdbid, i, ch,
                            repr(dbref.accession), repr(dbacc)))
                continue
            dbref.diff.append((resname, resnum, icode, dbnum, dbnum, comment))
        if not match:
            LOGGER.warn('SEQADV for chain {2}: database sequence reference '
                        'not found ({0}:{1})'.format(pdbid, i, ch))
            continue

    string = ' '.join([line[10:].strip() for i, line in lines['COMPND']])
    if string.startswith('MOL_ID'):
        dict_ = {}
        for molecule in string[6:].split('MOL_ID'):
            dict_.clear()
            for token in molecule.split(';'):
                token = token.strip()
                if not token:
                    continue
                items = token.split(':', 1)
                if len(items) == 2:
                    key, value = items
                    dict_[key.strip()] = value.strip()

            chains = dict_.pop('CHAIN', '').strip()

            if not chains:
                continue
            for ch in chains.split(','):
                ch = ch.strip()
                poly = polymers.get(ch, Polymer(ch))
                polymers[ch] = poly
                poly.name = dict_.get('MOLECULE', '')

                poly.fragment = dict_.get('FRAGMENT', '')

                poly.comments = dict_.get('OTHER_DETAILS', '')

                val = dict_.get('SYNONYM', '')
                poly.synonyms = [s.strip() for s in val.split(',')
                                 ] if val else []

                val = dict_.get('EC', '')
                poly.ec = [s.strip() for s in val.split(',')] if val else []

                poly.engineered = dict_.get('ENGINEERED', '') == 'YES'
                poly.mutation = dict_.get('MUTATION', '') == 'YES'

    return list(polymers.values())


def _getChemicals(lines):
    """Returns list of chemical components (heterogens)."""

    chemicals = defaultdict(list)
    chem_names = defaultdict(str)
    chem_synonyms = defaultdict(str)
    chem_formulas = defaultdict(str)
    for i, line in lines['HET   ']:
        chem = Chemical(line[7:10].strip())
        chem.chain = line[12].strip()
        chem.resnum = int(line[13:17])
        chem.icode = line[17].strip()
        chem.natoms = int(line[20:25].strip() or '0')
        chem.description = line[30:70].strip()
        chemicals[chem.resname].append(chem)
    for i, line in lines['HETNAM']:
        chem = line[11:14].strip()
        chem_names[chem] += line[15:70].rstrip()
    for i, line in lines['HETSYN']:
        chem = line[11:14].strip()
        chem_synonyms[chem] += line[15:70].strip()
    for i, line in lines['FORMUL']:
        chem = line[12:15].strip()
        chem_formulas[chem] += line[18:70].rstrip()

    for chem, name in chem_names.items():  # PY3K: OK
        name = cleanString(name)
        for chem in chemicals[chem]:
            chem.name = cleanString(name, nows=True)
    for chem, formula in chem_formulas.items():  # PY3K: OK
        formula = cleanString(formula)
        for chem in chemicals[chem]:
            chem.formula = formula
    for chem, synonyms in chem_synonyms.items():  # PY3K: OK
        synonyms = cleanString(synonyms)
        synonyms = synonyms.split(';')
        for chem in chemicals[chem]:
            chem.synonyms = [syn.strip() for syn in synonyms]

    alist = []
    for chem in chemicals.values():  # PY3K: OK
        for chem in chem:
            alist.append(chem)
    return alist


def _getVersion(lines):

    for i, line in lines['REMARK   4']:
        if 'COMPLIES' in line:
            try:
                # Return a string, because floating makes 3.20, 3.2 or
                # may arise problems if wwPDB uses a version number like 3.30.1
                return line.split('V.')[1].split(',')[0].strip()
            except:
                return None


def _getNumModels(lines):

    # "NUMMDL", Integer, 11 - 14: Number of models.
    line = lines['NUMMDL']
    if line:
        i, line = line[0]
        try:
            return int(line[10:14])
        except:
            pass

def _getCRYST1(lines):
    line = lines['CRYST1']
    if line:
       i, line = line[0]
       try:
           return {'cellLength': (float(line[6:15]),
                                  float(line[15:24]),
                                  float(line[24:33])),
                   'cellAngles': (float(line[33:40]),
                                  float(line[40:47]),
                                  float(line[47:54])),
                   'spaceGroup': line[55:66].strip(),
                   'Z value' : int(line[66:70])
                   }
       except:
           pass

def _missingResidues(lines):
    """
    Parse REMARK 465, Missing residues records
    REMARK 465   M RES C SSSEQI
    REMARK 465     GLU B   448
    header['missing_residues'] = [(chid, resname, resnum, icode, modelNum)]
    """
    mr = []
    header = True
    for i, line in lines['REMARK 465']:
        #skip header records
        if line.startswith("REMARK 465   M RES C SSSEQI"):
            header=False
            continue
        if header: continue
        modelNumStr = line[10:14].strip()
        if modelNumStr: modelNum = int(modelNumStr)
        else: modelNum = None
        w = line[15:].split()
        icode = ''
        if w[2][-1].isalpha():
            icode = w[2][-1]
            resnum = int(w[2][:-1])
        else:
            resnum = int(w[2])
        mr.append((w[1], w[0], resnum, icode, modelNumStr))
    return mr
    
def _missingAtoms(lines):
    """
    Parse REMARK 470,  Missing Atom records
    REMARK 470   M RES CSSEQI  ATOMS
    REMARK 470     ARG A 412    CG   CD   NE   CZ   NH1  NH2  
    header['missing_atoms'] = [(chid, resname, resnum, icode, modelNum, [atom names])]
    """
    ma = []
    res_repr_to_ma_index = {} # {residue_string_repr: index of this residue in ma list}
    header = True
    for i, line in lines['REMARK 470']:
        #skip header records
        if line.startswith("REMARK 470   M RES CSSEQI  ATOMS"):
            header=False
            continue
        if header:
            continue
        modelNumStr = line[10:14].strip()
        if modelNumStr:
            modelNum = int(modelNumStr)
        else:
            modelNum = None
        resname = line[15:18]
        chid = line[19]
        resnum = int(line[20:24])
        icode = line[24]
        if icode == ' ':
            icode = ''
        key = '%s:%s%d%s'%(chid,resname,resnum,icode)
        if key in res_repr_to_ma_index: # missing atoms record for a residue can span multiple lines 1vzq.pdb
            ma[res_repr_to_ma_index[key]][-1].extend(line[25:].split())
        else:
            res_repr_to_ma_index[key] = len(ma)
            ma.append((chid, resname, resnum, icode, line[25:].split()))
    return ma
    
# Make sure that lambda functions defined below won't raise exceptions
_PDB_HEADER_MAP = {
    'helix': _getHelix,
    'helix_range': _getHelixRange,
    'sheet': _getSheet,
    'sheet_range': _getSheetRange,
    'chemicals': _getChemicals,
    'polymers': _getPolymers,
    'reference': _getReference,
    'resolution': _getResolution,
    'biomoltrans': _getBiomoltrans,
    'version': _getVersion,
    'deposition_date': lambda lines: lines['HEADER'][0][1][50:59].strip()
                                       if lines['HEADER'] else None,
    'classification': lambda lines: lines['HEADER'][0][1][10:50].strip()
                                        if lines['HEADER'] else None,
    'identifier': lambda lines: lines['HEADER'][0][1][62:66].strip()
                                    if lines['HEADER'] else None,
    'title': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['TITLE ']])
                ) if lines['TITLE '] else None,
    'experiment': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['EXPDTA']])
                ) if lines['EXPDTA'] else None,
    'authors': lambda lines: cleanString(
                ''.join([line[1][10:].rstrip() for line in lines['AUTHOR']]),
                True).split(',') if lines['AUTHOR'] else None,
    'split': lambda lines: (' '.join([line[1][11:].rstrip()
                                      for line in lines['SPLIT ']])).split()
                             if lines['SPLIT '] else None,
    'model_type': lambda lines: cleanString(
                   ''.join([line[1][10:].rstrip() for line in lines['MDLTYP']])
                   ) if lines['MDLTYP'] else None,
    'n_models': _getNumModels,
    'space_group': _getSpaceGroup,
    'related_entries': _getRelatedEntries,
    'CRYST1': _getCRYST1,
    'SCALE': _getSCALE,
    'missing_residues': _missingResidues,
    'missing_atoms': _missingAtoms,
}

mapHelix = {
    1: 'H',  # 4-turn helix (alpha helix)
    2: '',  # other helix, Right-handed omega
    3: 'I',  # 5-turn helix (pi helix)
    4: '',  # other helix, Right-handed gamma
    5: 'G',  # 3-turn helix (3-10 helix)
    6: '',  # Left-handed alpha
    7: '',  # Left-handed omega
    8: '',  # Left-handed gamma
    9: '',  # 2 - 7 ribbon/helix
    10: '',  # Polyproline
}

def isHelix(secstrs):
    torf = np.logical_and(secstrs=='', False)
    for h in mapHelix.values():
        if h != '':
            torf = np.logical_or(torf, secstrs==h)
    return torf

def isSheet(secstrs):
    torf = secstrs == 'E'
    return torf

def assignSecstr(header, atoms, coil=True):
    """Assign secondary structure from *header* dictionary to *atoms*.
    *header* must be a dictionary parsed using the :func:`.parsePDB`.
    *atoms* may be an instance of :class:`.AtomGroup`, :class:`.Selection`,
    :class:`.Chain` or :class:`.Residue`.  ProDy can be configured to
    automatically parse and assign secondary structure information using
    ``confProDy(auto_secondary=True)`` command.  See also :func:`.confProDy`
    function.

    The Dictionary of Protein Secondary Structure, in short DSSP, type
    single letter code assignments are used:

      * **G** = 3-turn helix (310 helix). Min length 3 residues.
      * **H** = 4-turn helix (alpha helix). Min length 4 residues.
      * **I** = 5-turn helix (pi helix). Min length 5 residues.
      * **T** = hydrogen bonded turn (3, 4 or 5 turn)
      * **E** = extended strand in parallel and/or anti-parallel
        beta-sheet conformation. Min length 2 residues.
      * **B** = residue in isolated beta-bridge (single pair beta-sheet
        hydrogen bond formation)
      * **S** = bend (the only non-hydrogen-bond based assignment).
      * **C** = residues not in one of above conformations.


    See http://en.wikipedia.org/wiki/Protein_secondary_structure#The_DSSP_code
    for more details.

    Following PDB helix classes are omitted:

      * Right-handed omega (2, class number)
      * Right-handed gamma (4)
      * Left-handed alpha (6)
      * Left-handed omega (7)
      * Left-handed gamma (8)
      * 2 - 7 ribbon/helix (9)
      * Polyproline (10)

    Secondary structures are assigned to all atoms in a residue.  Amino acid
    residues without any secondary structure assignments in the header
    section will be assigned coil (C) conformation.  This can be prevented
    by passing ``coil=False`` argument."""

    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    helix = header.get('helix', {})
    sheet = header.get('sheet', {})
    if len(helix) == 0 and len(sheet) == 0:
        #LOGGER.warn('header does not contain secondary structure data')
        return atoms

    ssa = atoms.getSecstrs()
    if ssa is None:
        if isinstance(atoms, AtomGroup):
            ag = atoms
        else:
            ag = atoms.getAtomGroup()
        ag.setSecstrs(np.zeros(ag.numAtoms(),
                      ATOMIC_FIELDS['secondary'].dtype))
        ag.setSecids(np.zeros(ag.numAtoms(),
                      ATOMIC_FIELDS['secid'].dtype))
        ag.setSecclasses(np.zeros(ag.numAtoms(),
                      ATOMIC_FIELDS['secclass'].dtype)) 
        ag.setSecindices(np.zeros(ag.numAtoms(),
                      ATOMIC_FIELDS['secindex'].dtype))  

    prot = atoms.select('protein')
    if prot is not None and coil:
        prot.setSecstrs('C')
    hierview = atoms.getHierView()
    count = 0
    getResidue = hierview.getResidue
    for key, value in helix.items():  # PY3K: OK
        res = getResidue(*key)
        if res is None:
            continue
        res.setSecids(value[2])
        res.setSecclasses(value[0])
        res.setSecindices(value[1])
        res.setSecstrs(mapHelix[value[0]])
        
        count += 1
    for key, value in sheet.items():  # PY3K: OK
        res = getResidue(*key)
        if res is None:
            continue
        res.setSecids(value[2])
        res.setSecclasses(value[0])
        res.setSecindices(value[1])
        res.setSecstrs('E')
        count += 1

    LOGGER.info('Secondary structures were assigned to {0} residues.'
                .format(count))

    return atoms


def buildBiomolecules(header, atoms, biomol=None):
    """Returns *atoms* after applying biomolecular transformations from *header*
    dictionary.  Biomolecular transformations are applied to all coordinate
    sets in the molecule.

    Some PDB files contain transformations for more than 1 biomolecules.  A
    specific set of transformations can be choosen using *biomol* argument.
    Transformation sets are identified by numbers, e.g. ``"1"``, ``"2"``, ...

    If multiple biomolecular transformations are provided in the *header*
    dictionary, biomolecules will be returned as
    :class:`.AtomGroup` instances in a :func:`list`.

    If the resulting biomolecule has more than 26 chains, the molecular
    assembly will be split into multiple :class:`.AtomGroup`
    instances each containing at most 26 chains.  These
    :class:`.AtomGroup` instances will be returned in a tuple.

    Note that atoms in biomolecules are ordered according to chain identifiers.
    When multiple chains in a biomolecule have the same chain identifier, they 
    are given different segment names to distinguish them.
    """

    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')

    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')

    biomt = header.get('biomoltrans')
    if not isinstance(biomt, dict) or len(biomt) == 0:
        LOGGER.warn("no biomolecular transformations found so original structure was used")
        return atoms

    if not isinstance(atoms, AtomGroup):
        atoms = atoms.copy()

    biomols = []
    if biomol is None:
        keys = list(biomt)
    else:
        biomol = str(biomol)
        if biomol in biomt:
            keys = [biomol]
        else:
            LOGGER.warn('Transformations for biomolecule {0} was not '
                        'found in the header dictionary.'.format(biomol))
            return None

    keys.sort()
    for i in keys:
        segnm = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'*20)
        ags = []
        mt = biomt[i]
        # mt is a list, first item is list of chain identifiers
        # following items are lines corresponding to transformation
        # mt must have 3n + 1 lines
        if (len(mt)) % 4 != 0:
            LOGGER.warn('Biomolecular transformations {0} were not '
                        'applied'.format(i))
            continue

        for times in range(int((len(mt)) / 4)):
            rotation = np.zeros((3, 3))
            translation = np.zeros(3)

            try:
                line0 = np.fromstring(mt[times*4+1], sep=' ')
            except:
                line0 = np.frombuffer(mt[times*4+1], sep=' ')
            rotation[0, :] = line0[:3]
            translation[0] = line0[3]

            try:
                line1 = np.fromstring(mt[times*4+2], sep=' ')
            except:
                line1 = np.frombuffer(mt[times*4+2], sep=' ')
            rotation[1, :] = line1[:3]
            translation[1] = line1[3]

            try:
                line2 = np.fromstring(mt[times*4+3], sep=' ')
            except:
                line2 = np.frombuffer(mt[times*4+3], sep=' ')
            rotation[2, :] = line2[:3]
            translation[2] = line2[3]
            
            t = Transformation(rotation, translation)

            newag = atoms.select('chain ' + ' or chain '.join(mt[times*4+0]))
            if newag is None:
                continue
            newag = newag.copy()
            segnames = newag.all.getSegnames()
            newag.all.setSegnames(np.array([segname + decToHybrid36(times+1, resnum=True) 
                                            for segname in segnames]))
            for acsi in range(newag.numCoordsets()):
                newag.setACSIndex(acsi)
                newag = t.apply(newag)
            newag.setACSIndex(0)
            ags.append(newag)

        if ags:
            newag = ags.pop(0)
            while ags:
                newag += ags.pop(0)
            newag.setTitle('{0} biomolecule {1}'
                           .format(atoms.getTitle(), i))
            biomols.append(newag)
            
    if biomols:
        if len(biomols) == 1:
            return biomols[0]
        else:
            return biomols
    else:
        return None
