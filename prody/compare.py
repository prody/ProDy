# ProDy: A Python Package for Protein Dynamics Analysis
# 
# Copyright (C) 2010-2011 Ahmet Bakan
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
"""This module defines functions for comparing and mapping polypeptide chains.

Functions
---------

**Compare chains**:

  * :func:`matchChains`
  * :func:`mapOntoChain`
        
**Adjust settings**:
        
  * :func:`getPairwiseAlignmentMethod`
  * :func:`setPairwiseAlignmentMethod`
  * :func:`getPairwiseMatchScore`
  * :func:`setPairwiseMatchScore`
  * :func:`getPairwiseMismatchScore`
  * :func:`setPairwiseMismatchScore`
  * :func:`getPairwiseGapOpeningPenalty`
  * :func:`setPairwiseGapOpeningPenalty`
  * :func:`getPairwiseGapExtensionPenalty`
  * :func:`setPairwiseGapExtensionPenalty`
    
**Miscellaneous**:
    
  * :func:`getIntAsStr`
  
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2011 Ahmet Bakan'

import numpy as np
pairwise2 = None

import prody
from prody import ProDyLogger as LOGGER

from . import AtomMap, select

__all__ = ['matchChains',
           'mapOntoChain',
           'getPairwiseMatchScore', 'setPairwiseMatchScore',
           'getPairwiseMismatchScore', 'setPairwiseMismatchScore',
           'getPairwiseGapOpeningPenalty', 'setPairwiseGapOpeningPenalty',
           'getPairwiseGapExtensionPenalty', 'setPairwiseGapExtensionPenalty',
           'getPairwiseAlignmentMethod', 'setPairwiseAlignmentMethod',
           'getIntAsStr',
           ]

PAIRWISE_MATCH_SCORE = 1.0
PAIRWISE_MISMATCH_SCORE = 0.0
PAIRWISE_GAP_OPENING_PENALTY = -1.
PAIRWISE_GAP_EXTENSION_PENALTY = -0.1
PAIRWISE_ALIGNMENT_METHOD = 'global'
PAIRWISE_ALIGNMENT_GAP = '-'


GAPCHARS = ['-', '.']
NONE_A = '_'

_aaa2a = {
'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
'HSD': 'H', 'HSP': 'H', 'HSE': 'H',
'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O', 'XLE': 'J'
}

_a2aaa = {
'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 
'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

def _getSequence(resnames):
    """Return sequence of 1-letter codes for a given list of 3-letter amino acid 
    codes."""
    
    sequence = ''
    for rn in resnames:
        sequence += _aaa2a.get(rn, 'X')
    return sequence


def getPairwiseMatchScore():
    """Return match score used to align sequences."""
    
    return PAIRWISE_MATCH_SCORE


def setPairwiseMatchScore(pairwise_match_score):
    """Set match score used to align sequences."""
    
    if isinstance(pairwise_match_score, (float, int)) and pairwise_match_score >= 0:
        global PAIRWISE_MATCH_SCORE 
        PAIRWISE_MATCH_SCORE = pairwise_match_score
    else:
        raise TypeError('pairwise_match_score must be a positive number or zero')


def getPairwiseMismatchScore():
    """Return mismatch score used to align sequences."""
    
    return PAIRWISE_MISMATCH_SCORE


def setPairwiseMismatchScore(pairwise_mismatch_score):
    """Set mismatch score used to align sequences."""
    
    if isinstance(pairwise_mismatch_score, (float, int)) and pairwise_mismatch_score >= 0:
        global PAIRWISE_MISMATCH_SCORE
        PAIRWISE_MISMATCH_SCORE = pairwise_mismatch_score
    else:
        raise TypeError('pairwise_mismatch_score must be a positive number or zero')


def getPairwiseGapOpeningPenalty():
    """Return gap opening penalty used for pairwise alignment."""
    
    return PAIRWISE_GAP_OPENING_PENALTY


def setPairwiseGapOpeningPenalty(pairwise_gap_opening_penalty):
    """Set gap opening penalty used for pairwise alignment."""
    
    if isinstance(pairwise_gap_opening_penalty, (float, int)) and pairwise_gap_opening_penalty <= 0:
        global PAIRWISE_GAP_OPENING_PENALTY
        PAIRWISE_GAP_OPENING_PENALTY = pairwise_gap_opening_penalty
    else:
        raise TypeError('pairwise_gap_opening_penalty must be a negative number or zero')


def getPairwiseGapExtensionPenalty():
    """Return gap extension penalty used for pairwise alignment"""
    
    return PAIRWISE_GAP_EXTENSION_PENALTY


def setPairwiseGapExtensionPenalty(pairwise_gap_extension_penalty):
    """Set gap extension penalty used for pairwise alignment"""
    
    if isinstance(pairwise_gap_extension_penalty, (float, int)) and pairwise_gap_extension_penalty <= 0:
        global PAIRWISE_GAP_EXTENSION_PENALTY
        PAIRWISE_GAP_EXTENSION_PENALTY = pairwise_gap_extension_penalty
    else:
        raise TypeError('pairwise_gap_extension_penalty must be a negative number or zero')


def getPairwiseAlignmentMethod():
    """Return pairwise alignment method."""
    
    return PAIRWISE_ALIGNMENT_METHOD


def setPairwiseAlignmentMethod(method):
    """Set pairwise alignment method ("global" or "local")."""
    
    if method in ('local', 'global'):
        global PAIRWISE_ALIGNMENT_METHOD
        PAIRWISE_ALIGNMENT_METHOD = method
    else:
        raise ValueError('method must be "local" or "global"')


def getIntAsStr(lint, sep=' ', rng=' to '):
    """Return a structured string for a given list of ordered integers.
    
    :arg lint: list of integers
    :arg sep: range or number separator         
    :arg rng: inclusive range symbol

    E.g. for sep=' ' and rng = ' to ' 
        [1, 2, 3, 4, 10, 15, 16, 17] -> "1 to 4 10 15 to 17"
    for sep=',' and rng = '-'
        [1, 2, 3, 4, 10, 15, 16, 17] -> "1-4,10,15-17"
    
    """
    
    strint = ''
    i = -1
    for j in lint:
        if j < 0:
            continue
        if i < 0:
            i = j
        diff = j - i
        if diff == 0:
            strint += str(j)
        elif diff > 1: 
            strint += rng + str(i) + sep + str(j)
            k = j 
        i = j
    if diff == 1: 
        strint += rng + str(i)
    elif diff > 1 and k != j: 
        strint += rng + str(i) + sep + str(j)
    return strint


def matchChains(atoms1, atoms2, **kwargs):
    """Returns pairs of polypeptide chains sharing sequence identity.
    
    Makes an all-to-all comparison of chains in *atoms1* and *atoms2*.
    
    This function returns a list of matches. Each match is a tuple
    that contains 4 items:

      * Matched chain from *atoms1* as a :class:`~prody.atomic.AtomMap` instance, 
      * Matched chain from *atoms2* as a :class:`~prody.atomic.AtomMap` instance,
      * Percent sequence identitity of the match,
      * Percent sequence coverage of the match.
    
    List of matches are sorted in decreasing percent sequence identity order. 
    AtomMap can be used to calculate RMSD values and superpose atom groups.
    
    If *subset* is set to *calpha* or *backbone*, only alpha carbon
    atoms or backbone atoms will be paired. If set to *all*, all atoms
    common to matched residues will be returned.
    
    This function tries to match chains based on residue numbers and types. 
    All chains in *atoms1* is compared to all *chains* in *atoms2*. 
    This works well for different structures of the same
    protein. When it fails, :mod:`Bio.pairwise2` is used for sequence
    alignment, and matching is based on the sequence alignment.
    User can control, whether sequence alignment is performed or not with
    *pwalign* keyword. 
    
    :arg atoms1: atoms that will be mapped to the target *chain*
    :type atoms1: :class:`~prody.atomic.Chain`, :class:`~prody.atomic.AtomGroup`, 
                 or :class:`~prody.atomic.Selection`
    
    :arg atoms2: chain to which atoms will be mapped
    :type atoms2: :class:`~prody.atomic.Chain`, :class:`~prody.atomic.AtomGroup`, 
                 or :class:`~prody.atomic.Selection`
    
    :keyword subset: "calpha" (or "ca"), "backbone" (or "bb"), or "all", 
        default is "calpha"
    :type subset: string
    
    :keyword seqid: percent sequence identity, default is 90.
    :type seqid: float

    :keyword coverage: percent coverage, default is 90.
    :type coverage: float

    :keyword pwalign: perform pairwise sequence alignment 
    :type pwalign: bool
    
    .. versionchanged:: 0.6.1
       ``"ca"`` and and ``*bb*`` are accepted as *subset* argument. 
    
    """
    
    if not isinstance(atoms1, (prody.AtomGroup, prody.Chain, prody.Selection)):
        raise TypeError('atoms1 must be an AtomGroup, Chain, or Selection')
    if not isinstance(atoms2, (prody.AtomGroup, prody.Chain, prody.Selection)):
        raise TypeError('atoms2 must be an AtomGroup, Chain, or Selection')
    
    subset = kwargs.get('subset', 'calpha')
    if subset not in _SUBSETS:
        raise ValueError('{0:s} is not a valid subset argument'
            .format(str(subset)))
    seqid = kwargs.get('seqid', 90.)
    coverage  = kwargs.get('coverage', 90.)
    pwalign = kwargs.get('pwalign', None)
    
    if isinstance(atoms1, prody.Chain):
        simpch = SimpleChain()
        simpch.buildFromChain(atoms1)
        chains1 = [simpch]
        atoms1 = atoms1.getAtomGroup()
    else:
        chains1 = list()
        for ch in atoms1.getHierView().iterChains():
            simpch = SimpleChain()
            simpch.buildFromChain(ch)
            if len(simpch) > 0:
                chains1.append(simpch)
        if not isinstance(atoms1, prody.AtomGroup):
            atoms1 = atoms1.getAtomGroup()
        LOGGER.debug('Checking {0:s}: {1:d} chains are identified'.format(str(atoms1), len(chains1)))
    if isinstance(atoms2, prody.Chain):
        simpch = SimpleChain()
        simpch.buildFromChain(atoms2)
        chains2 = [simpch]
        atoms2 = atoms2.getAtomGroup()
    else:
        chains2 = list()
        for ch in atoms2.getHierView().iterChains():
            simpch = SimpleChain()
            simpch.buildFromChain(ch)
            if len(simpch) > 0:
                chains2.append(simpch)
        if not isinstance(atoms2, prody.AtomGroup):
            atoms2 = atoms2.getAtomGroup()
        LOGGER.debug('Checking {0:s}: {1:d} chains are identified'.format(str(atoms2), len(chains2)))
    
    matches = []
    unmatched = []
    LOGGER.debug('Trying to match chains based on residue numbers and identities:')
    for simpch1 in chains1:
        for simpch2 in chains2:
            LOGGER.debug('  Comparing {0:s} (len={1:d}) and {2:s} (len={3:d}):'
                        .format(simpch1.getName(), len(simpch1), simpch2.getName(), len(simpch2)))
            
            match1, match2, nmatches = _getTrivialMatch(simpch1, simpch2)
            _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
            _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))

            if _seqid >= seqid and _cover >= coverage:
                LOGGER.debug('\tMatch: {0:d} residues match with {1:.0f}% sequence identity and {2:.0f}% coverage.'
                            .format(len(match1), _seqid, _cover))
                matches.append((match1, match2, _seqid, _cover))
            else:
                LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, cover={1:.0f}%).'.format(_seqid, _cover))
                unmatched.append((simpch1, simpch2))


    if (not matches and (pwalign is None or pwalign)) or pwalign: 
        if pairwise2 is None: prody.importBioPairwise2()
        if pairwise2:
            for simpch1, simpch2 in unmatched:
                LOGGER.debug('Trying to match chains based on {0:s} sequence alignment:'.format(PAIRWISE_ALIGNMENT_METHOD))
                match1, match2, nmatches = _getAlignedMatch(simpch1, simpch2)
                _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
                _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))
                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMatch: {0:d} residues match with {1:.0f}% sequence identity and {2:.0f}% coverage.'
                                .format(len(match1), _seqid, _cover))
                    matches.append((match1, match2, _seqid, _cover))
                else:
                    LOGGER.debug('\tThese chains do not match.')
        else:
            LOGGER.warning('Pairwise alignment skipped (Biopython is not found).')
    
    if subset == 'calpha':
        subset = 'ca' 
    elif subset == 'backbone':
        subset = 'bb' 
        
    for mi, result in enumerate(matches):
        match1, match2, _seqid, _cover = result
        
        indices1 = []
        indices2 = []
        
        for i in xrange(len(match1)):
            ares = match1[i]
            bres = match2[i]

            if subset == 'ca':
                try:
                    aid = ares.getAtomNames().tolist().index('CA')
                except ValueError:
                    aid = None
                try:
                    bid = bres.getAtomNames().tolist().index('CA')
                    if aid is not None:
                        indices1.append(ares._indices[aid])
                        indices2.append(bres._indices[bid])
                except ValueError:
                    pass
            elif subset == 'bb':
                for bban in select.BACKBONE_ATOM_NAMES:
                    try:
                        aid = ares.getAtomNames().tolist().index(bban)
                    except ValueError:
                        continue
                    try:
                        bid = bres.getAtomNames().tolist().index(bban)
                        indices1.append(ares._indices[aid])
                        indices2.append(bres._indices[bid])
                    except ValueError:
                        continue
            elif subset is None or subset is 'all':
                aans = ares.getAtomNames()
                bans = bres.getAtomNames().tolist()

                aids = ares.getIndices()
                #bids = bres.getIndices()
                
                for j in xrange(len(aans)):
                    try:
                        bid = bres._indices[ bans.index( aans[j] ) ]
                        indices1.append(aids[j])
                        indices2.append(bid)
                    except ValueError:
                        pass
                
                
        indices1 = np.array(indices1)
        indices2 = np.array(indices2)
        lengh = len(indices1)
        match1 = AtomMap(atoms1, indices1, np.arange(lengh), np.array([]),
                               str(atoms1) + ' -> ' + str(atoms2),
                               atoms1._acsi)#'index ' + getIntAsStr(indices1), 
                                 
        match2 = AtomMap(atoms2, indices2, np.arange(lengh), np.array([]),
                               str(atoms2) + ' -> ' + str(atoms1),
                               atoms2._acsi)#'index ' + getIntAsStr(indices2), 
                                 
        matches[mi] = (match1, match2, _seqid, _cover)
    if len(matches) > 0:
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        matches.sort(compare, reverse=True)
    return matches

def _getTrivialMatch(ach, bch):
    """Return lists of matching residues (match is based on residue number).
    
    """
    #if not isinstance(ach, SimpleChain):
    #    raise TypeError('ach must be a SimpleChain instance')
    #if not isinstance(bch, SimpleChain):
    #    raise TypeError('bch must be a SimpleChain instance')
    amatch = []
    bmatch = []
    match = 0.0
    for ares in ach:
        bres = bch[ares.getNumber()]
        if bres is not None:
            if ares.getName() == bres.getName():
                match += 1
            amatch.append(ares.getResidue())
            bmatch.append(bres.getResidue())
    
    return amatch, bmatch, match
    
def _getAlignedMatch(ach, bch):
    """Return list of matching residues (match is based on sequence alignment).
    
    """
    if pairwise2 is None: prody.importBioPairwise2()
    if not pairwise2: return None
    if PAIRWISE_ALIGNMENT_METHOD == 'local':
        alignment = pairwise2.align.localms(ach.getSequence(), bch.getSequence(), 
                                                PAIRWISE_MATCH_SCORE, 
                                                PAIRWISE_MISMATCH_SCORE,
                                                PAIRWISE_GAP_OPENING_PENALTY, 
                                                PAIRWISE_GAP_EXTENSION_PENALTY,
                                                one_alignment_only=1)
    else:
        alignment = pairwise2.align.globalms(ach.getSequence(), bch.getSequence(), 
                                                 PAIRWISE_MATCH_SCORE, 
                                                 PAIRWISE_MISMATCH_SCORE,
                                                 PAIRWISE_GAP_OPENING_PENALTY, 
                                                 PAIRWISE_GAP_EXTENSION_PENALTY,
                                                 one_alignment_only=1)
               
    this = alignment[0][0]
    that = alignment[0][1]
    amatch = []
    bmatch = []
    aiter = ach.__iter__()
    biter = bch.__iter__()
    match = 0.0
    for i in xrange(len(this)):
        a = this[i]
        b = that[i]
        if a != PAIRWISE_ALIGNMENT_GAP:
            ares = next(aiter)
        if b != PAIRWISE_ALIGNMENT_GAP:
            bres = next(biter)
            if a != PAIRWISE_ALIGNMENT_GAP:
                amatch.append(ares.getResidue())
                bmatch.append(bres.getResidue())
                if a == b:
                    match += 1
    return amatch, bmatch, match

_SUBSETS = set(['ca', 'calpha', 'bb', 'backbone', 'all'])

def mapOntoChain(atoms, chain, **kwargs):
    """Map *atoms* onto *chain*. 
    
    This function returns a list of mappings. Each mapping is a tuple
    that contains 4 items:

      * Mapped chain as a :class:`~prody.atomic.AtomMap` instance, 
      * *chain* as it is
      * Percent sequence identitity,
      * Percent sequence coverage)
         
    Mappings are returned in decreasing percent sequence identity order.
    AtomMap that keeps mapped atom indices contains dummy atoms in place of 
    unmapped atoms.
    
    This function tries to map *atoms* to *chain* based on residue
    numbers and types. Each individual chain in *atoms* is compared to
    target *chain*. This works well for different structures of the same
    protein. When it fails, :mod:`Bio.pairwise2` is used for sequence
    alignment, and mapping is performed based on the sequence alignment.
    User can control, whether sequence alignment is performed or not with
    *pwalign* keyword. 
    
    :arg atoms: atoms that will be mapped to the target *chain*
    :type atoms: :class:`~prody.atomic.Chain`, :class:`~prody.atomic.AtomGroup`, 
                 or :class:`~prody.atomic.Selection`
    
    :arg chain: chain to which atoms will be mapped
    :type chain: :class:`~prody.atomic.Chain`
    
    :keyword subset: "calpha" (or "ca"), "backbone" (or "bb"), or "all", 
        default is "calpha"  
    :type subset: string
    
    :keyword seqid: percent sequence identity, default is 90.
    :type seqid: float

    :keyword coverage: percent coverage, default is 90.
    :type coverage: float

    :keyword pwalign: perform pairwise sequence alignment 
    :type pwalign: bool
    
    .. versionchanged:: 0.6.1
       Mapping can be performed for backbone or all atoms.
       ``"ca"`` and and ``*bb*`` are accepted as *subset* argument. 
    
    """
    target_chain = chain
    if not isinstance(atoms, (prody.AtomGroup, prody.Chain, prody.Selection)):
        raise TypeError('atoms must be an AtomGroup, a Chain, or a Selection instance')
    if not isinstance(target_chain, prody.Chain):
        raise TypeError('target_chain must be Chain instance')
        
    subset = str(kwargs.get('subset', 'calpha')).lower()
    if subset not in _SUBSETS:
        raise ValueError('{0:s} is not a valid subset argument'
            .format(str(subset)))
    seqid = kwargs.get('seqid', 90.)
    coverage  = kwargs.get('coverage', 90.) 
    pwalign = kwargs.get('pwalign', None)
    
    if isinstance(atoms, prody.Chain):
        chains = [atoms]
        map_ag = atoms.getAtomGroup()
    else:
        if isinstance(atoms, prody.AtomGroup):
            map_ag = atoms
        else:
            map_ag = atoms.getAtomGroup()
        chains = list(atoms.getHierView().iterChains())
        LOGGER.debug('Evaluating "{0:s}": {1:d} chains are identified'
                     .format(str(atoms), len(chains)))
    
    if subset != 'all':
        target_chain = target_chain.select(subset
                                ).getHierView()[target_chain.getIdentifier()]
    
    mappings = []
    unmapped = []
    target_ag = target_chain._ag
    simple_target = SimpleChain(True)
    simple_target.buildFromChain(target_chain)
    LOGGER.debug('Trying to map atoms based on residue numbers and identities:')
    for chain in chains:
        simple_chain = SimpleChain(True)
        simple_chain.buildFromChain(chain)
        if len(simple_chain) == 0:
            LOGGER.debug('  Skipping {0:s}, which does not contain any amino acid residues.'
                         .format(simple_chain))
            continue
        
        LOGGER.debug('  Comparing {0:s} (len={2:d}) with {1:s}:'.format(
                     simple_chain.getName(), simple_target.getName(), len(simple_chain)))
        
        target_list, chain_list, n_match, n_mapped = _getTrivialMapping(simple_target, simple_chain)
        if n_mapped > 0:
            _seqid = n_match * 100 / n_mapped
            #_cover = n_mapped * 100 / len(simple_target)
            _cover = n_mapped * 100 / max(len(simple_target), len(simple_chain))
        else:
            _seqid = 0
            _cover = 0
        
        if _seqid >= seqid and _cover >= coverage:
            LOGGER.debug('\tMapped: {0:d} residues match with {1:.0f}% sequence identity '
                         'and {2:.0f}% coverage.'.format(n_mapped, _seqid, _cover))
            mappings.append((target_list, chain_list, _seqid, _cover))
        else:
            LOGGER.debug('\tFailed to match chains based on residue numbers (seqid={0:.0f}%, cover={1:.0f}%).'.format(_seqid, _cover))
            unmapped.append(simple_chain)


    if (not mappings and (pwalign is None or pwalign)) or pwalign: 
        LOGGER.debug('Trying to map atoms based on {0:s} sequence alignment:'.format(PAIRWISE_ALIGNMENT_METHOD))
        for simple_chain in unmapped:
            LOGGER.debug('  Comparing {0:s} (len={2:d}) with {1:s}:'.format(
                         simple_chain.getName(), simple_target.getName(), len(simple_chain)))
            result = _getAlignedMapping(simple_target, simple_chain)
            if result is not None:
                target_list, chain_list, n_match, n_mapped = result
                if n_mapped > 0:
                    _seqid = n_match * 100 / n_mapped
                    #_cover = n_mapped * 100 / len(simple_target)
                    _cover = n_mapped * 100 / max(len(simple_target), len(simple_chain))
                else:
                    _seqid = 0
                    _cover = 0
                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMapped: {0:d} residues match with {1:.0f}% sequence identity and {2:.0f}% coverage.'
                                .format(n_mapped, _seqid, _cover))
                    mappings.append((target_list, chain_list, _seqid, _cover))
                #else:
                #    LOGGER.debug('\tThese two chains do not match.')
                #    return None
    
    for mi, result in enumerate(mappings):
        residues_target, residues_chain, _seqid, _cover = result
        indices_target = []
        indices_chain = []
        indices_mapping = []
        indices_dummies = []
        counter = 0
        for i in xrange(len(residues_target)):
            res_tar = residues_target[i]
            res_chn = residues_chain[i]
            
            for atom_tar in res_tar:
                indices_target.append(atom_tar.getIndex())
                if res_chn is not None:
                    try:
                        atom_chn = res_chn.getAtom(atom_tar.getName())
                        indices_chain.append(atom_chn.getIndex())
                        indices_mapping.append(counter)
                    except KeyError:
                        indices_dummies.append(counter)
                        pass
                else:
                    indices_dummies.append(counter)
                counter += 1
        #n_atoms = len(indices_target)   
        atommap = AtomMap(map_ag, 
                          indices_chain,
                          indices_mapping,
                          indices_dummies,
                          simple_chain.getName() + ' -> ' + simple_target.getName(),
                          chain._acsi)

        selection = AtomMap(target_ag, 
                                    indices_target,
                                    np.arange(len(indices_target)),
                                    np.array([]),
                                    simple_target.getName() + ' -> ' + simple_chain.getName(),
                                    target_chain._acsi)
                                    
        mappings[mi] = (atommap, selection, _seqid, _cover)
    if len(mappings) > 0:
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        mappings.sort(compare, reverse=True)
    return mappings


def _getMapping(target, chain, seqid, coverage):
    

    target_list, chain_list, n_match, n_mapped = _getTrivialMapping(target, chain)
    if n_mapped > 0:
        _seqid = n_match * 100 / n_mapped
        _cover = n_mapped * 100 / len(target)
    else:
        _seqid = 0
        _cover = 0
    
    if _seqid >= seqid and _cover >= coverage:
        LOGGER.debug('\t{0:d} residues match with {1:.0f}% sequence identity '
                     'and {2:.0f}% coverage based on residue numbers.'
                    .format(n_mapped, _seqid, _cover))
    else:
        LOGGER.debug('\tFailed to match chains based on residue numbers (seqid={0:.0f}%, cover={1:.0f}%).'
                     .format(_seqid, _cover))
        result = _getAlignedMapping(target, chain)
        if result is not None:
            target_list, chain_list, n_match, n_mapped = result
            if n_mapped > 0:
                _seqid = n_match * 100 / n_mapped
                _cover = n_mapped * 100 / len(target)
            else:
                _seqid = 0
                _cover = 0
            if _seqid >= seqid and _cover >= coverage:
                LOGGER.debug('\tMatch: {0:d} residues match with {1:.0f}% sequence identity and {2:.0f}% coverage based on {3:s} alignment using Bio.pairwise2.'
                            .format(n_mapped, _seqid, _cover, PAIRWISE_ALIGNMENT_METHOD))
        
    
    
    if _seqid < seqid or _cover < coverage:
        LOGGER.debug('\tThese two chains do not match.')
        return None
    else:
        return target_list, chain_list, _seqid, _cover

def _getTrivialMapping(target, chain):
    """Return lists of matching residues (match is based on residue number).
    
    """
    target_list = []
    chain_list = []
    n_match = 0
    n_mapped = 0
    chain_dict_get = chain._dict.get
    append = target_list.append 
    for target_residue in target:
        append(target_residue.getResidue())

        chain_residue = chain_dict_get(target_residue.getNumber(), None)
        if chain_residue is None:
            chain_list.append(chain_residue)
        else:
            if target_residue.getName() == chain_residue.getName():
                n_match += 1
            chain_list.append(chain_residue.getResidue())
            n_mapped += 1
            
    return target_list, chain_list, n_match, n_mapped

def _getAlignedMapping(target, chain):
    if pairwise2 is None: prody.importBioPairwise2()
    if not pairwise2: return None
    if PAIRWISE_ALIGNMENT_METHOD == 'local':
        alignment = pairwise2.align.localms(target.getSequence(), chain.getSequence(), 
                                                PAIRWISE_MATCH_SCORE, 
                                                PAIRWISE_MISMATCH_SCORE,
                                                PAIRWISE_GAP_OPENING_PENALTY, 
                                                PAIRWISE_GAP_EXTENSION_PENALTY,
                                                one_alignment_only=1)
    else:
        alignment = pairwise2.align.globalms(target.getSequence(), chain.getSequence(), 
                                                 PAIRWISE_MATCH_SCORE,
                                                 PAIRWISE_MISMATCH_SCORE,
                                                 PAIRWISE_GAP_OPENING_PENALTY, 
                                                 PAIRWISE_GAP_EXTENSION_PENALTY,
                                                 one_alignment_only=1)
               
    this = alignment[0][0]
    that = alignment[0][1]
    
    amatch = []
    bmatch = []
    aiter = target.__iter__()
    biter = chain.__iter__()
    n_match = 0
    n_mapped = 0
    for i in xrange(len(this)):
        a = this[i]
        b = that[i]
        if a not in (PAIRWISE_ALIGNMENT_GAP, NONE_A):
            ares = next(aiter)
            amatch.append(ares.getResidue())
            if b not in (PAIRWISE_ALIGNMENT_GAP, NONE_A):
                bres = next(biter)
                bmatch.append(bres.getResidue())
                if a == b:
                    n_match += 1
                n_mapped += 1
            else:
                bmatch.append(None)
        elif b not in (PAIRWISE_ALIGNMENT_GAP, NONE_A):
                bres = next(biter)
    return amatch, bmatch, n_match, n_mapped

class SimpleResidue(object):
    
    __slots__ = ['_res', '_name', '_num']
    
    def __init__(self, number, name, residue=None):
        self._num = number
        self._name = name
        self._res = residue
        
    def __repr__(self):
        return '<SimpleResidue: {0:s}{1:d}>'.format(self._name, self._num)
        
    def getResidue(self):
        return self._res
    
    def getNumber(self):
        return self._num
    
    def getName(self):
        return self._name


class SimpleChain(object):
    
    """An internal class used to compare two polypeptide chains.


    SimpleChain instances can be indexed using residue numbers. If a residue
    with given number is not found in the chain, ``None`` is returned.    
    
    """
    
    __slots__ = ['_list', '_seq', '_name', '_dict', '_gaps']
    
    def __init__(self, allow_gaps=False):
        """Initialize SimpleChain with a chain id and a sequence (available).
        
        :arg allow_gaps: allow gaps in the sequence of simple chain instance  
        :type allow_gaps: bool, default is False
        
        """
        self._dict = dict()
        self._list = list()
        self._seq = ''
        self._name = None
        self._gaps = allow_gaps

    def __len__(self):
        return len(self._list)
    
    def __iter__(self):
        return self._list.__iter__()
    
    def __repr__(self):
        return '<SimpleChain: {0:s} with {1:d} residues>'.format(
                    self._name, len(self._list))

    def __str__(self):
        return '{0:s} with {1:d} residues'.format( self._name, len(self._list))

    def __getitem__(self, index):
        try:
            return self._dict[index]
        except KeyError:
            return None
    
    def getSequence(self):
        return self._seq
    
    def getName(self):
        return self._name
    
    def buildFromSequence(self, sequence, resnums=None):        
        """Build from amino acid sequence.
        
        "-" or "." are acceptable amino acid types and are treated as gaps.

        :arg sequence: sequence of single letter amino acid codes 
        :type sequence: str
        :arg resnums: residue numbers corresponding the sequence
        :type resnums: a list of numbers, or a string representation of numbers
        
        Examples of *resnums* are:
            
            * 1:200 250:300
            
        """
        if not isinstance(sequence, str):
            raise TypeError('sequence must be string')
        self._dict = dict()
        self._list = list()
        self._seq = ''

        sqlen = len(sequence) - sum([sequence.count(gap) for gap in GAPCHARS])
        
        if resnums is None:
            resnums = range(1, sqlen+1)
        resid = 0
        for i, aa in enumerate(sequence):
            if self._gaps and aa in GAPCHARS:
                self._sequence += NONE_A
            else:
                resid = resnums[i]
                simpres = SimpleResidue(resid, aa)
                self._list.append(simpres)
                self._dict[resid] = simpres 
                self._sequence += aa
    
    def buildFromChain(self, chain):
        """Build from a :class:`prody.atomic.Chain`."""
        if not isinstance(chain, prody.Chain):
            raise TypeError('chain must be a Chain instance')
        self._dict = dict()
        self._list = list()
        self._seq = ''

        temp = 0
        for res in chain.iterResidues():
            if res.getAtom('CA') is None:
                continue
            resid = res.getNumber()
            aa = _aaa2a.get(res.getName(), 'X')
            simpres = SimpleResidue(resid, aa, res)
            if self._gaps:
                diff = resid - temp - 1
                if diff > 0:
                    self._seq += NONE_A * diff
                temp = resid
            self._seq += aa
            self._list.append(simpres)
            self._dict[resid] = simpres
        self._name = 'Chain {0:s} from {1:s}'.format(chain.getIdentifier(),
                                                     chain._ag.getName())

