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

"""This module defines functions for parsing header data from PDB files."""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

from collections import defaultdict
import os.path

import numpy as np

from prody import LOGGER
from prody.atomic import ATOMIC_FIELDS
from prody.atomic import Atomic, AtomGroup
from prody.atomic import getSequence
from prody.measure import Transformation
from prody.utilities import openFile

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
    modified    list    | modified residues (SEQMOD)
                        | when modified residues are present, each will be
                          represented as: ``(resname, resnum, icode, stdname,
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
    'UNIMES': 'UNIMES'
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


def parsePDBHeader(pdb, *keys):
    """Return header data dictionary for *pdb*.  This function is equivalent to
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
    REMARK 2     resolution        resolution of structures, when applicable
    REMARK 4     version           PDB file version
    REMARK 350   biomoltrans       biomolecular transformation lines
                                   (unprocessed)
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
    pdb = openFile(pdb)
    header, _ = getHeaderDict(pdb, *keys)
    pdb.close()
    return header


def getHeaderDict(stream, *keys):
    """Return header data in a dictionary.  *stream* may be a list of PDB lines
    or a stream."""

    lines = defaultdict(list)
    loc = 0
    for loc, line in enumerate(stream):
        startswith = line[0:6]
        if startswith in _START_COORDINATE_SECTION:
            break
        lines[startswith].append((loc, line))
    if not loc:
        raise ValueError('empty PDB file or stream')
    for i, line in lines['REMARK']:
        lines[line[:10]].append((i, line))

    pdbid = _PDB_HEADER_MAP['identifier'](lines)
    lines['pdbid'] = pdbid
    if keys:
        keys = list(keys)
        for k, key in enumerate(keys):
            if key in _PDB_HEADER_MAP:
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

    applyToChains = (' ')
    biomolecule = defaultdict(list)
    currentBiomolecule = '1'
    for i, line in lines['REMARK 350']:

        if line[13:18] == 'BIOMT':
            biomt = biomolecule[currentBiomolecule]
            if len(biomt) == 0:
                biomt.append(applyToChains)
            biomt.append(line[23:])
        elif line[11:41] == 'APPLY THE FOLLOWING TO CHAINS:':
            applyToChains = line[41:].replace(' ',
                                              '').strip().split(',')
        elif line[11:23] == 'BIOMOLECULE:':
            currentBiomolecule = line.split()[-1]
    return dict(biomolecule)


def _getResolution(lines):

    for i, line in lines['REMARK   2']:
        if 'RESOLUTION' in line:
            try:
                return float(line[23:30])
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


def _getSheet(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = {}
    for i, line in lines['SHEET ']:
        try:
            chid = line[21]
            value = (int(line[38:40]), int(line[7:10]),
                     line[11:14].strip())
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


def _getReference(lines):
    """Return a reference of the PDB entry."""

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


def _getPolymers(lines):
    """Return list of polymers (macromolecules)."""

    pdbid = lines['pdbid']
    polymers = dict()
    for i, line in lines['SEQRES']:
        ch = line[11]
        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        poly.sequence += ''.join(getSequence(line[19:].split()))

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
                        '({0}:{1})'.format(pdbid, i, ch))
        elif line[12] != ch:
            LOGGER.warn('DBREF1 and DBREF2 chain id mismatch'
                        '({0}:{1})'.format(pdbid, i, ch))

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
        poly.modified.append((line[12:15].strip(), line[18:22].strip() +
                              line[22].strip(), line[24:27].strip(),
                              line[29:70].strip()))

    for i, line in lines['SEQADV']:
        i += 1
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
        dbabbr = line[24:28].strip()
        resname = line[12:15].strip()
        try:
            resnum = int(line[18:22].strip())
        except:
            continue
            LOGGER.warn('SEQADV for chain {2}: failed to parse PDB sequence '
                        'number ({0}:{1})'.format(pdbid, i, ch))
        icode = line[22].strip()
        try:
            dbnum = int(line[43:48].strip())
        except:
            continue
            LOGGER.warn('SEQADV for chain {2}: failed to parse database '
                        'sequence number ({0}:{1})'.format(pdbid, i, ch))

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
            if dbref.accession != dbacc:
                LOGGER.warn('SEQADV for chain {2}: accession code '
                            'mismatch, expected {3} parsed {4} '
                            '({0}:{1})'.format(pdbid, i, ch,
                            repr(dbref.accession), repr(dbacc)))
                continue
            dbref.diff.append((resname, resnum, icode, dbnum, dbnum, comment))
        if not match:
            continue
            LOGGER.warn('SEQADV for chain {2}: database sequence reference '
                        'not found ({0}:{1})'.format(pdbid, i, ch))

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
    """Return list of chemical components (heterogens)."""

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
        chem_synonyms[chem] += line[15:70].rstrip()
    for i, line in lines['FORMUL']:
        chem = line[12:15].strip()
        chem_formulas[chem] += line[18:70].rstrip()

    for chem, name in chem_names.items():  # PY3K: OK
        name = cleanString(name)
        for chem in chemicals[chem]:
            chem.name = name
    for chem, formula in chem_formulas.items():  # PY3K: OK
        formula = cleanString(formula)
        for chem in chemicals[chem]:
            chem.formula = formula
    for chem, synonyms in chem_synonyms.items():  # PY3K: OK
        synonyms = cleanString(synonyms)
        synonyms = synonyms.split(';')
        for chem in chemicals[chem]:
            chem.synonyms = synonyms

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

# Make sure that lambda functions defined below won't raise exceptions
_PDB_HEADER_MAP = {
    'helix': _getHelix,
    'sheet': _getSheet,
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


def assignSecstr(header, atoms, coil=False):
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
        raise ValueError('header does not contain secondary structure data')

    ssa = atoms.getSecstrs()
    if ssa is None:
        if isinstance(atoms, AtomGroup):
            ag = atoms
        else:
            ag = atoms.getAtomGroup()
        ag.setSecstrs(np.zeros(ag.numAtoms(),
                      ATOMIC_FIELDS['secondary'].dtype))
    atoms.select('protein').setSecstrs('C')
    hierview = atoms.getHierView()
    count = 0
    getResidue = hierview.getResidue
    for key, value in helix.items():  # PY3K: OK
        res = getResidue(*key)
        if res is None:
            continue
        res.setSecstrs(mapHelix[value[0]])
        count += 1
    for key, res in sheet.items():  # PY3K: OK
        res = getResidue(*key)
        if res is None:
            continue
        res.setSecstrs('E')
        count += 1
    LOGGER.info('Secondary structures were assigned to {0} residues.'
                .format(count))
    return atoms


def buildBiomolecules(header, atoms, biomol=None):
    """Return *atoms* after applying biomolecular transformations from *header*
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
    """

    if not isinstance(header, dict):
        raise TypeError('header must be a dictionary')
    if not isinstance(atoms, Atomic):
        raise TypeError('atoms must be an Atomic instance')
    biomt = header.get('biomoltrans')
    if not isinstance(biomt, dict) or len(biomt) == 0:
        raise ValueError("header doesn't contain biomolecular transformations")

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
        if (len(mt) - 1) % 3 != 0:
            LOGGER.warn('Biomolecular transformations {0} were not '
                        'applied'.format(i))
            continue

        for times in range((len(mt) - 1) / 3):
            rotation = np.zeros((3, 3))
            translation = np.zeros(3)
            line = np.fromstring(mt[times*3+1], sep=' ')
            rotation[0, :] = line[:3]
            translation[0] = line[3]
            line = np.fromstring(mt[times*3+2], sep=' ')
            rotation[1, :] = line[:3]
            translation[1] = line[3]
            line = np.fromstring(mt[times*3+3], sep=' ')
            rotation[2, :] = line[:3]
            translation[2] = line[3]
            t = Transformation(rotation, translation)

            newag = atoms.select('chain ' + ' '.join(mt[0])).copy()
            if newag is None:
                continue
            newag.all.setSegnames(segnm.pop(0))
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
