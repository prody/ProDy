# -*- coding: utf-8 -*-
"""This module defines functions for parsing header data from PDB files."""

from collections import defaultdict, OrderedDict
import os.path

import numpy as np

from prody import LOGGER
from prody.atomic import ATOMIC_FIELDS
from prody.atomic import Atomic, AtomGroup
from prody.atomic import getSequence
from prody.atomic import flags
from prody.measure import Transformation
from prody.utilities import openFile, split

from .localpdb import fetchPDB
from .header import (Chemical, Polymer, DBRef, _PDB_DBREF,
                     cleanString)

from .starfile import parseSTARLines, parseSTARSection

__all__ = ['parseCIFHeader', 'getCIFHeaderDict']

_COMPND_KEY_MAPPINGS = {'_entity.id': 'MOL_ID',
                        '_entity.pdbx_description': 'MOLECULE',
                        '_entity.pdbx_fragment': 'FRAGMENT',
                        '_entity_name_com.name': 'SYNONYM',
                        '_entity.pdbx_ec': 'EC',
                        '_entity.pdbx_mutation': 'MUTATION',
                        '_entity.details': 'OTHER_DETAILS'}

_JRNL___KEY_MAPPINGS = {'_citation_author.name': 'AUTH',
                        '_citation.title': 'TITL',
                        #'_citation_editor.name': 'EDIT', # not normally found
                        '_citation.journal_abbrev': 'REF',
                        '_citation.journal_volume': 'REF',
                        '_citation.page_first': 'REF',
                        #'_citation.page_last': 'REF', # not in pdb
                        '_citation.year': 'REF',
                        #'_entity.pdbx_ec': 'PUBL', # publisher not supported
                        #'_citation.journal_id_ASTM': 'REFN',
                        '_citation.journal_id_ISSN': 'REFN',
                        #'_citation.journal_id_ISBN': 'REFN',
                        '_citation.pdbx_database_id_PubMed': 'PMID',
                        '_citation.pdbx_database_id_DOI': 'DOI'}

def _natomsFromFormula(formula, hydrogens=False):
    if ("(" in formula and ")" in formula):
        formula = formula.split("(")[1].split(")")[0]

    parts = formula.split()
    if not hydrogens:
        parts = [part for part in parts if part.find("H") == -1]

    return sum([_natomsFromFormulaPart(part) for part in parts])


def _natomsFromFormulaPart(part):
    digits = [s for s in part if s.isdigit()]
    if len(digits) == 0:
        return 1
    return int("".join(digits))

def parseCIFHeader(pdb, *keys):
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
            filename = fetchPDB(pdb, format='cif', compressed=False)
            if filename is None:
                raise IOError('PDB file for {0} could not be downloaded.'
                              .format(pdb))
            pdb = filename
        else:
            raise IOError('{0} is not a valid filename or a valid PDB '
                          'identifier.'.format(pdb))
    pdb = openFile(pdb, 'rt')
    header = getCIFHeaderDict(pdb, *keys)
    pdb.close()
    return header


def getCIFHeaderDict(stream, *keys):
    """Returns header data in a dictionary.  *stream* may be a list of PDB lines
    or a stream."""

    try:
        lines = stream.readlines()
        stream.close()
    except:
        lines = stream

    pdbid = _PDB_HEADER_MAP['identifier'](lines)
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
            return keys[0]
        else:
            return tuple(keys)
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
        return header


def _getBiomoltrans(lines):

    
    biomolecule = defaultdict(list)

    # 2 blocks are needed for this:
    # _pdbx_struct_assembly_gen: what to apply to which chains
    # _pdbx_struct_oper_list: everything else
    data1 = parseSTARSection(lines, '_pdbx_struct_assembly_gen')
    data2 = parseSTARSection(lines, '_pdbx_struct_oper_list')

    # extracting the data
    for n, item1 in enumerate(data1):
        currentBiomolecule = item1["_pdbx_struct_assembly_gen.assembly_id"]
        applyToChains = []

        chains = item1["_pdbx_struct_assembly_gen.asym_id_list"].split(',')
        applyToChains.extend(chains)

        biomt = biomolecule[currentBiomolecule]

        operators = item1["_pdbx_struct_assembly_gen.oper_expression"].split(',')
        for oper in operators:
            biomt.append(applyToChains)

            item2 = data2[int(oper)-1]

            for i in range(1,4):
                biomt.append(" ".join([
                    item2["_pdbx_struct_oper_list.matrix[{0}][1]".format(i)],
                    item2["_pdbx_struct_oper_list.matrix[{0}][2]".format(i)],
                    item2["_pdbx_struct_oper_list.matrix[{0}][3]".format(i)],
                    item2["_pdbx_struct_oper_list.vector[{0}]".format(i)],
                ]))

    return dict(biomolecule)


def _getResolution(lines):

    for line in lines:
        if line.startswith("_reflns.d_resolution_high"):
            try:
                return float(line.split()[1])
            except:
                pass
    return None


def _getRelatedEntries(lines):

    dbrefs = []

    try:
        key = "_pdbx_database_related"
        data = parseSTARSection(lines, key)
        for item in data:
            dbref = DBRef()
            dbref.accession = item[key + ".db_id"]
            dbref.dbabbr = item[key + ".db_name"]
            dbref.database = _PDB_DBREF.get(dbref.dbabbr, 'Unknown')
            
            dbrefs.append(dbref)
    except:
        LOGGER.warn("No related entries found")
        return None

    return dbrefs


def _getSpaceGroup(lines):

    for line in lines:
        if line.startswith("_symmetry.space_group_name_H-M"):
            try:
                return line.split()[1]
            except:
                pass
    return None


def _getHelix(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = {} 
    
    i = 0
    fields = OrderedDict()
    fieldCounter = -1
    foundHelixBlock = False
    doneHelixBlock = False
    start = 0
    stop = 0
    while not doneHelixBlock and i < len(lines):
        line = lines[i]
        if line[:13] == '_struct_conf.':
            fieldCounter += 1
            fields[line.split('.')[1].strip()] = fieldCounter

        if line.startswith('HELX_P'):
            if not foundHelixBlock:
                foundHelixBlock = True
                start = i
        else:
            if foundHelixBlock:
                doneHelixBlock = True
                stop = i
        i += 1

    if i < len(lines):
        for line in lines[start:stop]:
            data = line.split()

            try:
                chid = data[fields["beg_auth_asym_id"]]
                segn = data[fields["beg_label_asym_id"]]
                value = (int(data[fields["pdbx_PDB_helix_class"]]),
                        int(data[fields["id"]][6:]),
                        data[fields["pdbx_PDB_helix_id"]].strip())
            except:
                continue

            initICode = data[fields["pdbx_beg_PDB_ins_code"]]
            if initICode == '?':
                initICode = ' '

            initResnum = int(data[fields["beg_auth_seq_id"]])
            if initICode != ' ':
                for icode in alphas[alphas.index(initICode):]:
                    helix[(chid, initResnum, icode, segn)] = value
                initResnum += 1

            endICode = data[fields["pdbx_end_PDB_ins_code"]]
            if endICode == '?':
                endICode = ' '

            endResnum = int(data[fields["end_auth_seq_id"]])
            if endICode != ' ':
                for icode in alphas[:alphas.index(endICode)+1]:
                    helix[(chid, endResnum, icode, segn)] = value
                endResnum -= 1

            for resnum in range(initResnum, endResnum+1):
                helix[(chid, resnum, '', segn)] = value

    return helix


def _getHelixRange(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    helix = []

    i = 0
    fields = OrderedDict()
    fieldCounter = -1
    foundHelixBlock = False
    doneHelixBlock = False
    start = 0
    stop = 0
    while not doneHelixBlock and i < len(lines):
        line = lines[i]
        if line[:13] == '_struct_conf.':
            fieldCounter += 1
            fields[line.split('.')[1].strip()] = fieldCounter

        if line.startswith('HELX_P'):
            if not foundHelixBlock:
                foundHelixBlock = True
                start = i
        else:
            if foundHelixBlock:
                doneHelixBlock = True
                stop = i
        i += 1

    if i < len(lines):
        for line in lines[start:stop]:
            data = line.split()

            try:
                chid = data[fields["beg_auth_asym_id"]]
                Hclass=int(data[fields["pdbx_PDB_helix_class"]])
                Hnr=int(data[fields["id"]][6:])
            except:
                continue

            initResnum = int(data[fields["beg_auth_seq_id"]])
            endICode = data[fields["pdbx_end_PDB_ins_code"]]
            if endICode == '?':
                endICode = ' '
            endResnum = int(data[fields["end_auth_seq_id"]])
            if endICode != ' ':
                endResnum -= 1
            helix.append(['H', chid, Hclass, Hnr, initResnum, endResnum]) 
     
    return helix


def _getSheet(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = {}

    # mmCIF files have this data divided between 4 blocks

    # block 1 has how many strands are in each sheet - col 14:16 in PDB
    i = 0
    fields1 = OrderedDict()
    fieldCounter1 = -1
    foundSheetBlock1 = False
    foundSheetBlockData1 = False
    doneSheetBlock1 = False
    start1 = 0
    stop1 = 0
    while not doneSheetBlock1 and i < len(lines):
        line = lines[i]
        if line[:14] == '_struct_sheet.':
            fieldCounter1 += 1
            fields1[line.split('.')[1].strip()] = fieldCounter1
            if not foundSheetBlock1:
                foundSheetBlock1 = True

        if foundSheetBlock1:
            if not line.startswith('#'):
                if not foundSheetBlockData1:
                    start1 = i
                    foundSheetBlockData1 = True
            else:
                if foundSheetBlockData1:
                    doneSheetBlock1 = True
                    stop1 = i

        i += 1

    n_strands_dict = {}

    if i < len(lines):
        for line in lines[start1:stop1]:
            data = line.split()

            try:
                sheet_id = data[fields1["id"]]
                num_strands = int(data[fields1["number_strands"]])
            except:
                continue

            n_strands_dict[sheet_id] = num_strands

    # block 2 has the packing order and sense direction - col 38:40
    i = 0
    fields2 = OrderedDict()
    fieldCounter2 = -1
    foundSheetBlock2 = False
    foundSheetBlockData2 = False
    doneSheetBlock2 = False
    start2 = 0
    stop2 = 0
    while not doneSheetBlock2 and i < len(lines):
        line = lines[i]
        if line[:20] == '_struct_sheet_order.':
            fieldCounter2 += 1
            fields2[line.split('.')[1].strip()] = fieldCounter2
            if not foundSheetBlock2:
                foundSheetBlock2 = True

        if foundSheetBlock2:
            if not line.startswith('#'):
                if not foundSheetBlockData2:
                    start2 = i
                    foundSheetBlockData2 = True
            else:
                if foundSheetBlockData2:
                    doneSheetBlock2 = True
                    stop2 = i

        i += 1

    sense_dict = {}
    
    if i < len(lines):
        for sheet_id in n_strands_dict:
            sense_dict[sheet_id, 1] = 0

        for line in lines[start2:stop2]:
            data = line.split()

            try:
                sheet_id = data[fields2["sheet_id"]]
                strand2 = int(data[fields2["range_id_2"]])
                sense = data[fields2["sense"]]
            except:
                continue

            if sense == "parallel":
                sense_dict[sheet_id, strand2] = 1
            else:
                sense_dict[sheet_id, strand2] = -1

    # block 3 has the ranges - col 7:37 - NB: sheet and strand ID switched
    i = 0
    fields3 = OrderedDict()
    fieldCounter3 = -1
    foundSheetBlock3 = False
    foundSheetBlockData3 = False
    doneSheetBlock3 = False
    start3 = 0
    stop3 = 0
    while not doneSheetBlock3 and i < len(lines):
        line = lines[i]
        if line[:20] == '_struct_sheet_range.':
            fieldCounter3 += 1
            fields3[line.split('.')[1].strip()] = fieldCounter3
            if not foundSheetBlock3:
                foundSheetBlock3 = True

        if foundSheetBlock3:
            if not line.startswith('#'):
                if not foundSheetBlockData3:
                    start3 = i
                    foundSheetBlockData3 = True
            else:
                if foundSheetBlockData3:
                    doneSheetBlock3 = True
                    stop3 = i

        i += 1

    if i < len(lines):
        for line in lines[start3:stop3]:
            data = line.split()

            try:
                chid = data[fields3["beg_auth_asym_id"]]
                segn = data[fields3["beg_label_asym_id"]]

                strand_id = int(data[fields3["id"]])
                sheet_id = data[fields3["sheet_id"]]
                sense = sense_dict[sheet_id, strand_id]
                value = (sense, strand_id, sheet_id)
            except:
                continue

            initICode = data[fields3["pdbx_beg_PDB_ins_code"]]
            if initICode == '?':
                initICode = ' '

            initResnum = int(data[fields3["beg_auth_seq_id"]])
            if initICode != ' ':
                for icode in alphas[alphas.index(initICode):]:
                    sheet[(chid, initResnum, icode, segn)] = value
                initResnum += 1

            endICode = data[fields3["pdbx_end_PDB_ins_code"]]
            if endICode == '?':
                endICode = ' '

            endResnum = int(data[fields3["end_auth_seq_id"]])
            if endICode != ' ':
                for icode in alphas[:alphas.index(endICode)+1]:
                    sheet[(chid, endResnum, icode, segn)] = value
                endResnum -= 1

            for resnum in range(initResnum, endResnum+1):
                sheet[(chid, resnum, '', segn)] = value

    # block 4 has the hydrogen bonds - col 41:70 - NB: very different order
    # we don't actually use this info though

    return sheet


def _getSheetRange(lines):

    alphas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    sheet = []

    # mmCIF files have this data divided between 4 blocks

    # block 1 has how many strands are in each sheet - col 14:16 in PDB
    i = 0
    fields1 = OrderedDict()
    fieldCounter1 = -1
    foundSheetBlock1 = False
    foundSheetBlockData1 = False
    doneSheetBlock1 = False
    start1 = 0
    stop1 = 0
    while not doneSheetBlock1 and i < len(lines):
        line = lines[i]
        if line[:14] == '_struct_sheet.':
            fieldCounter1 += 1
            fields1[line.split('.')[1].strip()] = fieldCounter1
            if not foundSheetBlock1:
                foundSheetBlock1 = True

        if foundSheetBlock1:
            if not line.startswith('#'):
                if not foundSheetBlockData1:
                    start1 = i
                    foundSheetBlockData1 = True
            else:
                if foundSheetBlockData1:
                    doneSheetBlock1 = True
                    stop1 = i

        i += 1

    n_strands_dict = {}
    if i < len(lines):
        for line in lines[start1:stop1]:
            data = line.split()

            try:
                sheet_id = data[fields1["id"]]
                num_strands = int(data[fields1["number_strands"]])
            except:
                continue

            n_strands_dict[sheet_id] = num_strands

    # block 2 has the packing order and sense direction - col 38:40
    i = 0
    fields2 = OrderedDict()
    fieldCounter2 = -1
    foundSheetBlock2 = False
    foundSheetBlockData2 = False
    doneSheetBlock2 = False
    start2 = 0
    stop2 = 0
    while not doneSheetBlock2 and i < len(lines):
        line = lines[i]
        if line[:20] == '_struct_sheet_order.':
            fieldCounter2 += 1
            fields2[line.split('.')[1].strip()] = fieldCounter2
            if not foundSheetBlock2:
                foundSheetBlock2 = True

        if foundSheetBlock2:
            if not line.startswith('#'):
                if not foundSheetBlockData2:
                    start2 = i
                    foundSheetBlockData2 = True
            else:
                if foundSheetBlockData2:
                    doneSheetBlock2 = True
                    stop2 = i

        i += 1

    sense_dict = {}
    if i < len(lines):
        for sheet_id in n_strands_dict:
            sense_dict[sheet_id, 1] = 0

        for line in lines[start2:stop2]:
            data = line.split()

            try:
                sheet_id = data[fields2["sheet_id"]]
                strand2 = int(data[fields2["range_id_2"]])
                sense = data[fields2["sense"]]
            except:
                continue

            if sense == "parallel":
                sense_dict[sheet_id, strand2] = 1
            else:
                sense_dict[sheet_id, strand2] = -1

    # block 3 has the ranges - col 7:37 - NB: sheet and strand ID switched
    i = 0
    fields3 = OrderedDict()
    fieldCounter3 = -1
    foundSheetBlock3 = False
    foundSheetBlockData3 = False
    doneSheetBlock3 = False
    start3 = 0
    stop3 = 0
    while not doneSheetBlock3 and i < len(lines):
        line = lines[i]
        if line[:20] == '_struct_sheet_range.':
            fieldCounter3 += 1
            fields3[line.split('.')[1].strip()] = fieldCounter3
            if not foundSheetBlock3:
                foundSheetBlock3 = True

        if foundSheetBlock3:
            if not line.startswith('#'):
                if not foundSheetBlockData3:
                    start3 = i
                    foundSheetBlockData3 = True
            else:
                if foundSheetBlockData3:
                    doneSheetBlock3 = True
                    stop3 = i

        i += 1

    if i < len(lines):
        for line in lines[start3:stop3]:
            data = line.split()

            try:
                chid = data[fields3["beg_auth_asym_id"]]

                Snr = int(data[fields3["id"]])
                sheet_id = data[fields3["sheet_id"]]
                dir = sense_dict[sheet_id, Snr]
            except:
                continue

            initICode = data[fields3["pdbx_beg_PDB_ins_code"]]
            if initICode == '?':
                initICode = ' '

            initResnum = int(data[fields3["beg_auth_seq_id"]])
            if initICode != ' ':
                initResnum += 1

            endICode = data[fields3["pdbx_end_PDB_ins_code"]]
            if endICode == '?':
                endICode = ' '

            endResnum = int(data[fields3["end_auth_seq_id"]])
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

    # JRNL double block. Blocks 6 and 7 as copied from COMPND
    # Block 1 has most info. Block 2 has author info
    items1 = parseSTARSection(lines, "_citation")
    items2 = parseSTARSection(lines, "_citation_author")

    for row in items1:
        for k, value in row.items():
            try:
                what = _JRNL___KEY_MAPPINGS[k]
            except:
                continue
            if what == 'TITL':
                title += value.upper()
            elif what == 'REF':
                if k == "_citation.journal_volume":
                    reference += " V. "
                reference += ' ' + value.upper()
                reference = reference.strip()
            elif what == 'REFN':
                ref['issn'] = 'ISSN ' + value
            elif what == 'PMID':
                ref['pmid'] = value
            elif what == 'DOI':
                ref['doi'] = value.upper()

    for row in items2:
        for k, value in row.items(): 
            try:
                what = _JRNL___KEY_MAPPINGS[k]
            except:
                continue
            if what == 'AUTH':
                surname, initials = value.split(',')
                author = initials+surname
                authors.append(author.strip().upper())

        ref['authors'] = authors
        ref['title'] = cleanString(title)
        ref['editors'] = editors
        ref['reference'] = cleanString(reference)
        ref['publisher'] = cleanString(publisher)

    return ref


def _getPolymers(lines):
    """Returns list of polymers (macromolecules)."""

    pdbid = lines[0].split("_")[1].strip()
    polymers = dict()

    entities = defaultdict(list)

    # SEQRES block
    items1 = parseSTARSection(lines, '_entity_poly')

    for item in items1:
        chains = item['_entity_poly.pdbx_strand_id']
        entity = item['_entity_poly.entity_id']
        
        for ch in chains.split(","):
            entities[entity].append(ch)
            poly = polymers.get(ch, Polymer(ch))
            polymers[ch] = poly
            poly.sequence += ''.join(item[
                '_entity_poly.pdbx_seq_one_letter_code'][1:-2].split())

    # DBREF block 1
    items2 = parseSTARSection(lines, '_struct_ref')

    for item in items2:
        entity = item["_struct_ref.id"]
        chains = entities[entity]
        for ch in chains:
            dbabbr = item["_struct_ref.db_name"]
            dbref = DBRef()
            dbref.dbabbr = dbabbr
            dbref.database = _PDB_DBREF.get(dbabbr, 'Unknown')
            dbref.accession = item["_struct_ref.pdbx_db_accession"]
            dbref.idcode = item["_struct_ref.db_code"]

            poly = polymers[ch]
            poly.dbrefs.append(dbref)

    # DBREF block 2
    items3 = parseSTARSection(lines, "_struct_ref_seq")

    for i, item in enumerate(items3):
        i += 1

        ch = item["_struct_ref_seq.pdbx_strand_id"]
        poly = polymers[ch] 
        
        for dbref in poly.dbrefs:
            if item["_struct_ref_seq.pdbx_db_accession"] == dbref.accession:

                try:
                    first = int(item["_struct_ref_seq.pdbx_auth_seq_align_beg"])
                    initICode = item["_struct_ref_seq.pdbx_db_align_beg_ins_code"]
                    if initICode == '?':
                        initICode = ' '
                except:
                    try:
                        first = int(item["_struct_ref_seq.pdbx_auth_seq_align_beg"])
                        initICode = item["_struct_ref_seq.pdbx_seq_align_beg_ins_code"]
                        if initICode == '?':
                            initICode = ' '
                    except:
                        LOGGER.warn('DBREF for chain {2}: failed to parse '
                                    'initial sequence number of the PDB sequence '
                                    '({0}:{1})'.format(pdbid, i, ch))
                try:
                    last = int(item["_struct_ref_seq.pdbx_auth_seq_align_end"])
                    endICode = item["_struct_ref_seq.pdbx_db_align_end_ins_code"]
                    if endICode == '?':
                        endICode = ' '
                except:
                    try:
                        last = int(item["_struct_ref_seq.pdbx_auth_seq_align_end"])
                        endICode = item["_struct_ref_seq.pdbx_seq_align_end_ins_code"]
                        if endICode == '?':
                            endICode = ' '
                    except:
                        LOGGER.warn('DBREF for chain {2}: failed to parse '
                                    'ending sequence number of the PDB sequence '
                                    '({0}:{1})'.format(pdbid, i, ch))
                try:
                    first2 = int(item["_struct_ref_seq.db_align_beg"])
                    dbref.first = (first, initICode, first2)
                except:
                    LOGGER.warn('DBREF for chain {2}: failed to parse '
                                'initial sequence number of the database sequence '
                                '({0}:{1})'.format(pdbid, i, ch))
                try:
                    last2 = int(item["_struct_ref_seq.db_align_end"])
                    dbref.last = (last, endICode, last2)
                except:
                    LOGGER.warn('DBREF for chain {2}: failed to parse '
                                'ending sequence number of the database sequence '
                                '({0}:{1})'.format(pdbid, i, ch))

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
            
            try:
                resnum.append((dbref.first[0], dbref.last[0]))
            except:
                pass # we've already warned about this

        resnum.sort()
        last = -10000
        for first, temp in resnum:
            if first <= last:
                LOGGER.warn('DBREF records overlap for chain {0} ({1})'
                            .format(poly.chid, pdbid))
            last = temp

    # MODRES block
    data4 = parseSTARSection(lines, "_pdbx_struct_mod_residue")

    for data in data4:
        ch = data["_pdbx_struct_mod_residue.label_asym_id"]

        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        if poly.modified is None:
            poly.modified = []

        iCode = data["_pdbx_struct_mod_residue.PDB_ins_code"]
        if iCode == '?':
            iCode == '' # PDB one is stripped
        poly.modified.append((data["_pdbx_struct_mod_residue.auth_comp_id"],
                                data["_pdbx_struct_mod_residue.auth_asym_id"],
                                data["_pdbx_struct_mod_residue.auth_seq_id"] + iCode,
                                data["_pdbx_struct_mod_residue.parent_comp_id"],
                                data["_pdbx_struct_mod_residue.details"]))

    # SEQADV block
    data5 = parseSTARSection(lines, "_struct_ref_seq_dif")

    for i, data in enumerate(data5):
        ch = data["_struct_ref_seq_dif.pdbx_pdb_strand_id"]

        poly = polymers.get(ch, Polymer(ch))
        polymers[ch] = poly
        dbabbr = data["_struct_ref_seq_dif.pdbx_seq_db_name"]
        resname = data["_struct_ref_seq_dif.mon_id"]
        if resname == '?':
            resname = '' # strip for pdb

        try:
            resnum = int(data["_struct_ref_seq_dif.pdbx_auth_seq_num"])
        except:
            #LOGGER.warn('SEQADV for chain {2}: failed to parse PDB sequence '
            #            'number ({0}:{1})'.format(pdbid, i, ch))
            continue

        icode = data["_struct_ref_seq_dif.pdbx_pdb_ins_code"]
        if icode == '?':
            icode = '' # strip for pdb            
        
        try:
            dbnum = int(data["_struct_ref_seq_dif.pdbx_seq_db_seq_num"])
        except:
            #LOGGER.warn('SEQADV for chain {2}: failed to parse database '
            #            'sequence number ({0}:{1})'.format(pdbid, i, ch))
            continue            

        comment = data["_struct_ref_seq_dif.details"].upper()
        if comment == '?':
            comment = '' # strip for pdb 

        match = False
        for dbref in poly.dbrefs:
            if not dbref.first[0] <= resnum <= dbref.last[0]:
                continue
            match = True
            if dbref.dbabbr != dbabbr:
                LOGGER.warn('SEQADV for chain {2}: reference database '
                            'mismatch, expected {3} parsed {4} '
                            '({0}:{1})'.format(pdbid, i+1, ch,
                            repr(dbref.dbabbr), repr(dbabbr)))
                continue
            dbacc = data["_struct_ref_seq_dif.pdbx_seq_db_accession_code"]
            if dbref.accession[:9] != dbacc[:9]:
                LOGGER.warn('SEQADV for chain {2}: accession code '
                            'mismatch, expected {3} parsed {4} '
                            '({0}:{1})'.format(pdbid, i+1, ch,
                            repr(dbref.accession), repr(dbacc)))
                continue
            dbref.diff.append((resname, resnum, icode, dbnum, dbnum, comment))
        if not match:
            LOGGER.warn('SEQADV for chain {2}: database sequence reference '
                        'not found ({0}:{1})'.format(pdbid, i+1, ch))
            continue

    # COMPND double block. 
    # Block 6 has most info. Block 7 has synonyms
    data6 = parseSTARSection(lines, "_entity")
    data7 = parseSTARSection(lines, "_entity_name_com")

    dict_ = {}
    for molecule in data6:
        dict_.clear()
        for k, value in molecule.items():
            if k == '_entity.id':
                dict_['CHAIN'] = ', '.join(entities[value])

            try:
                key = _COMPND_KEY_MAPPINGS[k]
            except:
                continue
            val = value.strip()
            if val == '?':
                val = ''
            dict_[key.strip()] = val

        chains = dict_.pop('CHAIN', '').strip()

        if not chains:
            continue
        for ch in chains.split(','):
            ch = ch.strip()
            poly = polymers.get(ch, Polymer(ch))
            polymers[ch] = poly
            poly.name = dict_.get('MOLECULE', '').upper()

            poly.fragment = dict_.get('FRAGMENT', '').upper()

            poly.comments = dict_.get('OTHER_DETAILS', '').upper()

            val = dict_.get('EC', '')
            poly.ec = [s.strip() for s in val.split(',')] if val else []

            poly.mutation = dict_.get('MUTATION', '') != ''
            poly.engineered = dict_.get('ENGINEERED', poly.mutation)

    for molecule in data7:
        dict_.clear()
        for k, value in molecule.items():
            if k.find('entity_id') != -1:
                dict_['CHAIN'] = ', '.join(entities[value])

            try:
                key = _COMPND_KEY_MAPPINGS[k]
            except:
                continue
            dict_[key.strip()] = value.strip()

        chains = dict_.pop('CHAIN', '').strip()

        if not chains:
            continue
        for ch in chains.split(','):
            ch = ch.strip()
            poly = polymers.get(ch, Polymer(ch))
            polymers[ch] = poly

            val = dict_.get('SYNONYM', '')
            poly.synonyms = [s.strip().upper() for s in val.split(',')
                                ] if val else []

    return list(polymers.values())


def _getChemicals(lines):
    """Returns list of chemical components (heterogens)."""

    chemicals = defaultdict(list)
    chem_names = defaultdict(str)
    chem_synonyms = defaultdict(str)
    chem_formulas = defaultdict(str)
    chem_n_atoms = defaultdict(int)

    # Data is split across blocks again

    # 1st block we need is has info about location in structure

    # this instance only includes single sugars not branched structures
    items = parseSTARSection(lines, "_pdbx_nonpoly_scheme")

    for data in items:
        resname = data["_pdbx_nonpoly_scheme.mon_id"]
        if resname in flags.AMINOACIDS or resname == "HOH":
            continue

        chem = Chemical(resname)
        chem.chain = data["_pdbx_nonpoly_scheme.pdb_strand_id"]
        chem.resnum = int(data["_pdbx_nonpoly_scheme.pdb_seq_num"])

        icode = data["_pdbx_nonpoly_scheme.pdb_ins_code"]
        if icode == '.':
            icode = ''
        chem.icode = icode
        chem.description = '' # often empty in .pdb and not clearly here
        chemicals[chem.resname].append(chem)

    # next we get the equivalent one for branched sugars part
    items = parseSTARSection(lines, "_pdbx_branch_scheme")

    for data in items:
        resname = data["_pdbx_branch_scheme.mon_id"]
        if resname in flags.AMINOACIDS or resname == "HOH":
            continue

        chem = Chemical(resname)
        chem.chain = data["_pdbx_branch_scheme.pdb_asym_id"]
        chem.resnum = int(data["_pdbx_branch_scheme.pdb_seq_num"])

        chem.icode = '' # this part doesn't have this field
        chem.description = '' # often empty in .pdb and not clearly here
        chemicals[chem.resname].append(chem)

    # 2nd block to get has general info e.g. name and formula
    items = parseSTARSection(lines, "_chem_comp")

    for data in items:
        resname = data["_chem_comp.id"]
        if resname in flags.AMINOACIDS or resname == "HOH":
            continue

        chem_names[resname] += data["_chem_comp.name"].upper()

        synonym = data["_chem_comp.pdbx_synonyms"]
        if synonym == '?':
            synonym = ' '
        synonym = synonym.rstrip()
        if synonym.startswith(';') and synonym.endswith(';'):
            synonym = synonym[1:-1]
        chem_synonyms[resname] += synonym
        
        chem_formulas[resname] += data["_chem_comp.formula"]


    for key, name in chem_names.items():  # PY3K: OK
        name = cleanString(name)
        for chem in chemicals[key]:
            chem.name = name

    for key, formula in chem_formulas.items():  # PY3K: OK
        formula = cleanString(formula)
        repeats = len(chemicals[key])
        for chem in chemicals[key]:
            chem.formula = '{0}({1})'.format(repeats, formula)
            chem.natoms = _natomsFromFormula(formula)

    for key, synonyms in chem_synonyms.items():  # PY3K: OK
        synonyms = cleanString(synonyms)
        synonyms = synonyms.split(';')
        for chem in chemicals[key]:
            if synonyms != ['']:
                chem.synonyms = [syn.strip() for syn in synonyms]

    alist = []
    for chem in chemicals.values():  # PY3K: OK
        for chem in chem:
            alist.append(chem)
    return alist


def _getVersion(lines):

    for line in lines:
        if line.startswith("_audit_conform.dict_version"):
            return line.split()[1]
    return None


def _getNumModels(lines):

    for line in lines:
        if line.startswith("_pdbx_nmr_ensemble.conformers_submitted_total_number"):
            data = line.split()[1]
            try:
                return int(data)
            except:
                pass
    return None


def _getTitle(lines):

    title = ''

    try:
        data = parseSTARSection(lines, "_struct")
        for item in data:
            title += item['_struct.title'].upper()
    except:
        pass

    if len(title) == 0:
        title = None
        
    return title


def _getAuthors(lines):

    authors = []

    try:
        data = parseSTARSection(lines, "_audit_author")
        for item in data:
            author = ''.join(item['_audit_author.name'].split(', ')[::-1])
            authors.append(author.upper())
    except:
        pass

    if len(authors) == 0:
        authors = None
        
    return authors


def _getSplit(lines):

    split = []

    key = "_pdbx_database_related"

    try:
        data, _ = parseSTARSection(lines, key)
        for item in data:
            if item[key + '.content_type'] == 'split':
                split.append(item[key + '.db_id'])
    except:
        pass

    if len(split) == 0:
        split = None
        
    return split


def _getModelType(lines):

    model_type = ''

    model_type += [line.split()[1] for line in lines
                  if line.find("_struct.pdbx_model_type_details") != -1][0]

    if model_type == '?':
        model_type = None

    return model_type


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
    'deposition_date': lambda lines: [line.split()[1]
                                      if line.find("initial_deposition_date") != -1 else None
                                      for line in lines][0],
    'classification': lambda lines: [line.split()[1]
                                     if line.find("_struct_keywords.pdbx_keywords") != -1 else None
                                     for line in lines][0],
    'identifier': lambda lines: lines[0].split("_")[1].strip() if len(lines[0].split("_")) else '',
    'title': _getTitle,
    'experiment': lambda lines: [line.split()[1]
                                 if line.find("_exptl.method") != -1 else None
                                 for line in lines][0],
    'authors': _getAuthors,
    'split': _getSplit,
    'model_type': _getModelType,
    'n_models': _getNumModels,
    'space_group': _getSpaceGroup,
    'related_entries': _getRelatedEntries,
}
