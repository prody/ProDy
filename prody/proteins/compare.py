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

"""This module defines functions for comparing and mapping polypeptide chains.
"""

__author__ = 'Ahmet Bakan'
__copyright__ = 'Copyright (C) 2010-2012 Ahmet Bakan'

import numpy as np
PW2 = None

from prody.atomic import AtomMap as AM
from prody.atomic import Chain, AtomGroup, Selection
from prody.atomic import AAMAP
from prody.atomic import flags
from prody.measure import calcTransformation, calcRMSD, printRMSD, calcDistance
from prody import LOGGER, SELECT

__all__ = ['matchChains', 'matchAlign', 'mapOntoChain',
           'getMatchScore', 'setMatchScore',
           'getMismatchScore', 'setMismatchScore',
           'getGapPenalty', 'setGapPenalty',
           'getGapExtPenalty', 'setGapExtPenalty',
           'getAlignmentMethod', 'setAlignmentMethod']

MATCH_SCORE = 1.0
MISMATCH_SCORE = 0.0
GAP_PENALTY = -1.
GAP_EXT_PENALTY = -0.1
ALIGNMENT_METHOD = 'global'
GAP = '-'


GAPCHARS = ['-', '.']
NONE_A = '_'

_a2aaa = {
'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 
'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 
'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 
'Y': 'TYR', 'V': 'VAL'
}

def importBioPairwise2():

    global PW2
    if PW2 is None:    
        try:
            import pairwise2
        except ImportError:
            try:
                from Bio import pairwise2
            except ImportError:
                raise ImportError('pairwise2 module could not be imported. '
                                  'Reinstall ProDy or install BioPython '
                                  'to solve the problem.')
        PW2 = pairwise2
    return PW2

def getMatchScore():
    """Return match score used to align sequences."""
    
    return MATCH_SCORE

def setMatchScore(match_score):
    """Set match score used to align sequences."""
    
    if isinstance(match_score, (float, int)) and match_score >= 0:
        global MATCH_SCORE 
        MATCH_SCORE = match_score
    else:
        raise TypeError('match_score must be a positive number or zero')

def getMismatchScore():
    """Return mismatch score used to align sequences."""
    
    return MISMATCH_SCORE

def setMismatchScore(mismatch_score):
    """Set mismatch score used to align sequences."""
    
    if isinstance(mismatch_score, (float, int)) and mismatch_score >= 0:
        global MISMATCH_SCORE
        MISMATCH_SCORE = mismatch_score
    else:
        raise TypeError('mismatch_score must be a positive number or zero')

def getGapPenalty():
    """Return gap opening penalty used for pairwise alignment."""
    
    return GAP_PENALTY

def setGapPenalty(gap_penalty):
    """Set gap opening penalty used for pairwise alignment."""
    
    if isinstance(gap_penalty, (float, int)) and gap_penalty <= 0:
        global GAP_PENALTY
        GAP_PENALTY = gap_penalty
    else:
        raise TypeError('gap_penalty must be a negative number')

def getGapExtPenalty():
    """Return gap extension penalty used for pairwise alignment."""
    
    return GAP_EXT_PENALTY

def setGapExtPenalty(gap_ext_penalty):
    """Set gap extension penalty used for pairwise alignment."""
    
    if isinstance(gap_ext_penalty, (float, int)) and gap_ext_penalty <= 0:
        global GAP_EXT_PENALTY
        GAP_EXT_PENALTY = gap_ext_penalty
    else:
        raise TypeError('gap_ext_penalty must be a negative number or zero')

def getAlignmentMethod():
    """Return pairwise alignment method."""
    
    return ALIGNMENT_METHOD

def setAlignmentMethod(method):
    """Set pairwise alignment method (global or local)."""
    
    if method in ('local', 'global'):
        global ALIGNMENT_METHOD
        ALIGNMENT_METHOD = method
    else:
        raise ValueError('method must be "local" or "global"')

class SimpleResidue(object):
    
    __slots__ = ['_res', '_name', '_num', '_inc']
    
    def __init__(self, number, name, insertioncode='', residue=None):
        self._num = number
        self._name = name
        self._inc = insertioncode
        self._res = residue
        
    def __repr__(self):
        return '<SimpleResidue: {0:s}{1:d}>'.format(self._name, self._num)
        
    def getResidue(self):
        return self._res
    
    def getResnum(self):
        return self._num
    
    def getIcode(self):
        return self._inc

    def getResname(self):
        return self._name


class SimpleChain(object):
    
    """An internal class used to compare two polypeptide chains.


    SimpleChain instances can be indexed using residue numbers. If a residue
    with given number is not found in the chain, **None** is returned."""
    
    __slots__ = ['_list', '_seq', '_title', '_dict', '_gaps']
    
    def __init__(self, chain=None, allow_gaps=False):
        """Initialize SimpleChain with a chain id and a sequence (available).
        
        :arg chain: chain instance or single-letter amino acid sequence  
        :type chain: str, :class:`.Chain` 
        
        :arg allow_gaps: allow gaps in the sequence of simple chain instance, 
            default is False  
        :type allow_gaps: bool"""
        
        self._dict = dict()
        self._list = list()
        self._seq = ''
        self._title = None
        self._gaps = allow_gaps
        if isinstance(chain, Chain): 
            self.buildFromChain(chain)
        elif isinstance(chain, str):
            self.buildFromSequence(chain)
    def __len__(self):
        return len(self._list)
    
    def __iter__(self):
        return self._list.__iter__()
    
    def __repr__(self):
        return '<SimpleChain: {0:s} with {1:d} residues>'.format(
                    self._title, len(self._list))

    def __str__(self):
        return '{0:s} with {1:d} residues'.format(self._title, len(self._list))

    def __getitem__(self, index):
        if isinstance(index, int):
            self._dict.get((index, ''))
        return self._dict.get(index)
    
    def getSequence(self):
        return self._seq
    
    def getTitle(self):
        return self._title
    
    def buildFromSequence(self, sequence, resnums=None):        
        """Build from amino acid sequence.
        
        "-" or "." are acceptable amino acid types and are treated as gaps.

        :arg sequence: sequence of single letter amino acid codes 
        :type sequence: str
        :arg resnums: residue numbers corresponding the sequence
        :type resnums: a list of numbers, or a string representation of numbers
        
        Examples of *resnums* are:
            
            * 1:200 250:300"""
            
        assert isinstance(sequence, str), 'sequence must be string'
        assert sequence.isalpha(), 'sequence must be all alpha'

        if resnums is None:
            resnums = range(1, len(sequence)+1)
        resid = 0
        gaps = self._gaps
        for i, aa in enumerate(sequence):
            if gaps and aa in GAPCHARS:
                self._seq += NONE_A
            else:
                resid = resnums[i]
                simpres = SimpleResidue(resid, aa)
                self._list.append(simpres)
                self._dict[resid] = simpres 
                self._seq += aa
    
    def buildFromChain(self, chain):
        """Build from a :class:`.Chain`."""
        
        assert isinstance(chain, Chain), 'chain must be a Chain instance'
        gaps = self._gaps
        residues = list(chain.iterResidues())
        temp = residues[0].getResnum()-1
        protein_resnames = flags.AMINOACIDS
        for res in chain:
            if not res.getResname() in protein_resnames:
                continue
            resid = res.getResnum()
            incod = res.getIcode()
            aa = AAMAP.get(res.getResname(), 'X')
            simpres = SimpleResidue(resid, aa, incod, res)
            if gaps:
                diff = resid - temp - 1
                if diff > 0:
                    self._seq += NONE_A * diff
                temp = resid
            self._seq += aa
            self._list.append(simpres)
            self._dict[(resid, incod)] = simpres
        self._title = 'Chain {0:s} from {1:s}'.format(chain.getChid(),
                                             chain.getAtomGroup().getTitle())


def countUnpairedBreaks(chone, chtwo, resnum=True):
    """This function is under development.
    Return number of unpaired breaks in aligned chains *chone* and *chtwo*, 
    which are expected to be :class:`.AtomMap` instances obtained from one of 
    :func:`.matchChains` or :func:`.mapOntoChain` functions. 
    
    Pairwise global or local alignment of chains with missing residues may be
    problematic, as in the following illustrations for example.  This function
    helps identifying some of these problems.
     
    Breaks in a chain are determined using Cα-Cα distances between consecutive
    residues, i.e. Cα to Cα distance larger than 4 Å corresponds to a break or
    gap of 1 or more residues.  This function counts such breaks in *chone* or 
    *chtwo* that is not paired with a break in the other.
    
    The following example illustrates a problem that may occur when aligning
    two structures of the same protein chain where one of the chains have
    a few missing residues::
        
      Correct alignment: A.L.R.S - - V.W.Y.K.L  -> no unpaired breaks
      Target chain     : A.L.R.S.V.T.V.W.Y.K.L
      Wrong alignment  : A.L.R.S_V - - W.Y.K.L  
                                |
                                --> 1 unpaired break, counted
      
      Key:
          - (dash) is a gap in the alignment
          . (dot) is a peptide bond
          _ (underscore) is a break
      
    In this case, one unpaired break is an indicator of the problem in the 
    alignment.

    The following example illustrates a case where an unpaired break is due to
    an insertion in the homologous chain::

      Target chain     : 1A.2L.3R.4S.5V.6T.7V
      Homologous chain : 1A.2L.3K.4S.6V_9S.10L  
                                       |
                                       --> 1 unpaired break, not counted
      
    In this case, residue numbers are used to determine whether the unpaired
    break is due to an insertion/deletion in the chain or misalignment."""

    try:    
        if chone.numAtoms() != chtwo.numAtoms():
            raise ValueError('number of atoms do not match')
    except AttributeError:
        raise TypeError('one and two must be Atomic instances')
    
    mapped = chone.getMappedFlags() * chtwo.getMappedFlags()
    if mapped.sum() == 0:
        raise ValueError('chains do not have common mapped atoms')
    chone = chone[mapped]
    chtwo = chtwo[mapped]
    
    rnone = chone.getResnums()
    rntwo = chtwo.getResnums()
    
    brone = calcDistance(chone[1:], chone[:-1]) > 4.
    brtwo = calcDistance(chtwo[1:], chtwo[:-1]) > 4.
    from code import interact
    #interact(local=locals())
    
    brone[(rnone[1:] - rnone[:-1]) > 1] = False
    brtwo[(rntwo[1:] - rntwo[:-1]) > 1] = False
    
    brone = set(brone.nonzero()[0])
    brtwo = set(brtwo.nonzero()[0])
    #interact(local=locals())
    return len(brone.union(brtwo)) - len(brone.intersection(brtwo))
    

_SUBSETS = set(['ca', 'calpha', 'bb', 'backbone', 'all'])

def matchAlign(mobile, target, **kwargs):
    """Superpose *mobile* onto *target* based on best matching pair of chains.
    This function uses :func:`matchChains` for matching chains and returns a 
    tuple that contains the following items:
      
      * *mobile* after it is superposed,
      * matching chain from *mobile* as a :class:`.AtomMap` instance, 
      * matching chain from *target* as a :class:`.AtomMap` instance,
      * percent sequence identity of the match,
      * percent sequence overlap of the match.
      
    :arg mobile: atoms that contain a protein chain
    :type mobile: :class:`.Chain`, :class:`.AtomGroup`, :class:`.Selection`
    
    :arg target: atoms that contain a protein chain
    :type target: :class:`.Chain`, :class:`.AtomGroup`, :class:`.Selection`
     
    :arg tarsel: *target* atoms that will be used for alignment,
        default is ``'calpha'``
    :type tarsel: str
     
    :arg allcsets: align all coordinate sets of *mobile*, default is **True**
    :type allcsets: bool

    :keyword seqid: percent sequence identity, default is 90
    :type seqid: float

    :keyword overlap: percent overlap, default is 90
    :type overlap: float

    :keyword pwalign: perform pairwise sequence alignment 
    :type pwalign: bool"""
    
    selstr = kwargs.pop('tarsel', 'calpha')
    if selstr == 'calpha':
        selstr = None
    subset = 'calpha'
    if selstr:
        if selstr in _SUBSETS:
            subset = selstr
        else: 
            subset = 'all'
        sel = target.select(selstr)
        if sel is None:
            raise ValueError('selection {0:s} did not match any atoms'
                             .format(repr(selstr)))
        chid = set(sel.getChids())
        if len(chid) == 1:
            chid = chid.pop()
            target = target.select('chain ' + chid)
    
    match = matchChains(mobile, target, subset=subset, **kwargs)
    if not match:
        return
    match = match[0]
    mob = match[0]
    tar = match[1]
    if selstr:
        which = SELECT.getIndices(tar, selstr)
        n_atoms = len(which)
    else:
        which = slice(None)
        n_atoms = len(tar)
        selstr = 'calpha'
    
    if kwargs.get('allcets', True):
        cslabel = kwargs.get('cslabel', 'Coordinate set')
        csincr = kwargs.get('csincr', 0)
        indent = (lambda acsi, cslabel=cslabel, csincr=csincr:
                  '  {0:s} {1:d}: '.format(cslabel, acsi + csincr))
        csets = range(mobile.numCoordsets())
    else:
        indent = lambda acsi: ''
        csets = [mobile.getACSIndex()]
    
    LOGGER.info('Alignment is based on {0:d} atoms matching {1:s}.'
                .format(n_atoms, repr(selstr)))
    printRMSD(tar._getCoords()[which], mob._getCoordsets()[:, which], 
              msg='Before alignment ')
    for acsi in csets:
        mob.setACSIndex(acsi)
        mobile.setACSIndex(acsi)
        calcTransformation(mob._getCoords()[which], 
                           tar._getCoords()[which]).apply(mobile)
    printRMSD(tar._getCoords()[which], mob._getCoordsets()[:, which], 
              msg='After alignment  ')
    return (mobile,) + match

def matchChains(atoms1, atoms2, **kwargs):
    """Return pairs of chains matched based on sequence similarity.  Makes an 
    all-to-all comparison of chains in *atoms1* and *atoms2*.  Chains are 
    obtained from hierarchical views (:class:`.HierView`) of atom groups.  
    This function returns a list of matching chains in a tuples that contain
    4 items:

      * matching chain from *atoms1* as a :class:`.AtomMap` 
        instance, 
      * matching chain from *atoms2* as a :class:`.AtomMap` 
        instance,
      * percent sequence identity of the match,
      * percent sequence overlap of the match.
    
    List of matches are sorted in decreasing percent sequence identity order. 
    :class:`.AtomMap` instances can be used to calculate RMSD values and 
    superpose atom groups.
    
    :arg atoms1: atoms that contain a chain
    :type atoms1: :class:`.Chain`, :class:`.AtomGroup`, :class:`.Selection`
    
    :arg atoms2: atoms that contain a chain
    :type atoms2: :class:`.Chain`, :class:`.AtomGroup`, :class:`.Selection`
    
    :keyword subset: ``"calpha"`` (or ``"ca"``), ``"backbone"`` (or ``"bb"``), 
        or ``"all"``, default is ``"calpha"``
    :type subset: string
    
    :keyword seqid: percent sequence identity, default is 90
    :type seqid: float

    :keyword overlap: percent overlap, default is 90
    :type overlap: float

    :keyword pwalign: perform pairwise sequence alignment 
    :type pwalign: bool

    If *subset* is set to *calpha* or *backbone*, only alpha carbon
    atoms or backbone atoms will be paired. If set to *all*, all atoms
    common to matched residues will be returned.
    
    This function tries to match chains based on residue numbers and names.  
    All chains in *atoms1* is compared to all chains in *atoms2*.  This works 
    well for different structures of the same protein.  When it fails, 
    :mod:`Bio.pairwise2` is used for pairwise sequence alignment, and matching 
    is performed based on the sequence alignment.  User can control, whether 
    sequence alignment is performed or not with *pwalign* keyword.  If 
    ``pwalign=True`` is passed, pairwise alignment is enforced."""
    
    if not isinstance(atoms1, (AtomGroup, Chain, Selection)):
        raise TypeError('atoms1 must be an AtomGroup, Chain, or Selection')
    if not isinstance(atoms2, (AtomGroup, Chain, Selection)):
        raise TypeError('atoms2 must be an AtomGroup, Chain, or Selection')
    
    subset = kwargs.get('subset', 'calpha')
    if subset not in _SUBSETS:
        raise ValueError('{0:s} is not a valid subset argument'
                         .format(str(subset)))
    seqid = kwargs.get('seqid', 90.)
    assert isinstance(seqid, (float, int)), 'seqid must be float'
    assert 0 < seqid <= 100, 'seqid must be in the range from 0 to 100'
    coverage = kwargs.get('overlap')
    if coverage is None:
        coverage = kwargs.get('coverage', 90.)
    assert isinstance(coverage, (float, int)), 'overlap must be float'
    assert 0 < coverage <= 100, 'overlap must be in the range from 0 to 100'
    pwalign = kwargs.get('pwalign', None)
    
    if isinstance(atoms1, Chain):
        chains1 = [atoms1]
        atoms1 = atoms1.getAtomGroup()
    else:
        chains1 = list(atoms1.getHierView().iterChains())
        if not isinstance(atoms1, AtomGroup):
            atoms1 = atoms1.getAtomGroup()
    chains = list()
    for ch in chains1:
        simpch = SimpleChain(ch)
        if len(simpch) > 0:
            chains.append(simpch)
    chains1 = chains
    if not isinstance(atoms1, Chain):
        LOGGER.debug('Checking {0:s}: {1:d} chains are identified'
                     .format(str(atoms1), len(chains1)))
        
    if isinstance(atoms2, Chain):
        chains2 = [atoms2]
        atoms2 = atoms2.getAtomGroup()
    else:
        chains2 = list(atoms2.getHierView().iterChains())
        if not isinstance(atoms2, AtomGroup):
            atoms2 = atoms2.getAtomGroup()
    chains = list()
    for ch in chains2:
        simpch = SimpleChain(ch)
        if len(simpch) > 0:
            chains.append(simpch)
    chains2 = chains
    if not isinstance(atoms2, Chain):
        LOGGER.debug('Checking {0:s}: {1:d} chains are identified'
                     .format(str(atoms2), len(chains2)))

    matches = []
    unmatched = []
    LOGGER.debug('Trying to match chains based on residue numbers and names:')
    for simpch1 in chains1:
        for simpch2 in chains2:
            LOGGER.debug('  Comparing {0:s} (len={1:d}) and {2:s} (len={3:d}):'
                         .format(simpch1.getTitle(), len(simpch1), 
                                 simpch2.getTitle(), len(simpch2)))
            
            match1, match2, nmatches = getTrivialMatch(simpch1, simpch2)
            _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
            _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))

            if _seqid >= seqid and _cover >= coverage:
                LOGGER.debug('\tMatch: {0:d} residues match with {1:.0f}% '
                             'sequence identity and {2:.0f}% overlap.'
                            .format(len(match1), _seqid, _cover))
                matches.append((match1, match2, _seqid, _cover, simpch1, simpch2))
            else:
                LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, '
                             'overlap={1:.0f}%).'.format(_seqid, _cover))
                unmatched.append((simpch1, simpch2))
            

    if pwalign or (not matches and (pwalign is None or pwalign)): 
        pairwise2 = importBioPairwise2()
        if pairwise2:
            LOGGER.debug('Trying to match chains based on {0:s} sequence '
                         'alignment:'.format(ALIGNMENT_METHOD))
            for simpch1, simpch2 in unmatched:
                LOGGER.debug('  Comparing {0:s} (len={1:d}) and {2:s} '
                             '(len={3:d}):'
                             .format(simpch1.getTitle(), len(simpch1), 
                                     simpch2.getTitle(), len(simpch2)))
                match1, match2, nmatches = getAlignedMatch(simpch1, simpch2)
                _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
                _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))
                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMatch: {0:d} residues match with {1:.0f}% '
                                 'sequence identity and {2:.0f}% overlap.'
                                 .format(len(match1), _seqid, _cover))
                    matches.append((match1, match2, _seqid, _cover, 
                                    simpch1, simpch2))
                else:
                    LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, '
                                 'overlap={1:.0f}%).'
                                 .format(_seqid, _cover))
        else:
            LOGGER.warning('Pairwise alignment could not be performed.')
    if not matches:
        return None
    if subset == 'calpha':
        subset = 'ca' 
    elif subset == 'backbone':
        subset = 'bb'
    for mi, result in enumerate(matches):
        match1, match2, _seqid, _cover, simpch1, simpch2 = result
        
        indices1 = []
        indices2 = []
        
        for i in xrange(len(match1)):
            ares = match1[i]
            bres = match2[i]

            if subset == 'ca':
                try:
                    aid = ares.getNames().tolist().index('CA')
                except ValueError:
                    aid = None
                try:
                    bid = bres.getNames().tolist().index('CA')
                    if aid is not None:
                        indices1.append(ares._indices[aid])
                        indices2.append(bres._indices[bid])
                except ValueError:
                    pass
            elif subset == 'bb':
                for bban in ('N', 'CA', 'C', 'O'):
                    try:
                        aid = ares.getNames().tolist().index(bban)
                    except ValueError:
                        continue
                    try:
                        bid = bres.getNames().tolist().index(bban)
                        indices1.append(ares._indices[aid])
                        indices2.append(bres._indices[bid])
                    except ValueError:
                        continue
            elif subset is None or subset is 'all':
                aans = ares.getNames()
                bans = bres.getNames().tolist()

                aids = ares.getIndices()
                #bids = bres.getIndices()
                
                for j in xrange(len(aans)):
                    try:
                        bid = bres._indices[ bans.index( aans[j] ) ]
                        indices1.append(aids[j])
                        indices2.append(bid)
                    except ValueError:
                        pass

        indices1 = np.array(indices1, int)
        indices2 = np.array(indices2, int)
        lengh = len(indices1)
        
        match1 = AM(atoms1, indices1, np.arange(lengh), None,
                    simpch1.getTitle() + ' -> ' + simpch2.getTitle(),
                    acsi=atoms1.getACSIndex(), allintarrays=True) 
        match2 = AM(atoms2, indices2, np.arange(lengh), None,
                    simpch2.getTitle() + ' -> ' + simpch1.getTitle(),
                    acsi=atoms2.getACSIndex(), allintarrays=True) 
                                 
        matches[mi] = (match1, match2, _seqid, _cover)
    if len(matches) > 1:
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        matches.sort(compare, reverse=True)
    return matches

def getTrivialMatch(ach, bch):
    """Return lists of matching residues (match based on residue number).
    
    """
    #if not isinstance(ach, SimpleChain):
    #    raise TypeError('ach must be a SimpleChain instance')
    #if not isinstance(bch, SimpleChain):
    #    raise TypeError('bch must be a SimpleChain instance')
    amatch = []
    bmatch = []
    match = 0.0
    for ares in ach:
        bres = bch[(ares.getResnum(), ares.getIcode())]
        if bres is not None:
            if ares.getResname() == bres.getResname():
                match += 1
            amatch.append(ares.getResidue())
            bmatch.append(bres.getResidue())
    
    return amatch, bmatch, match
    
def getAlignedMatch(ach, bch):
    """Return list of matching residues (match is based on sequence alignment).
    """
    
    pairwise2 = importBioPairwise2()
    if ALIGNMENT_METHOD == 'local':
        alignment = pairwise2.align.localms(ach.getSequence(), 
                                            bch.getSequence(), 
                                            MATCH_SCORE, MISMATCH_SCORE,
                                            GAP_PENALTY, GAP_EXT_PENALTY,
                                            one_alignment_only=1)
    else:
        alignment = pairwise2.align.globalms(ach.getSequence(), 
                                             bch.getSequence(), 
                                             MATCH_SCORE, MISMATCH_SCORE,
                                             GAP_PENALTY, GAP_EXT_PENALTY,
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
        if a != GAP:
            ares = next(aiter)
        if b != GAP:
            bres = next(biter)
            if a != GAP:
                amatch.append(ares.getResidue())
                bmatch.append(bres.getResidue())
                if a == b:
                    match += 1
    return amatch, bmatch, match

def mapOntoChain(atoms, chain, **kwargs):
    """Map *atoms* onto *chain*.  This function returns a list of mappings. 
    Each mapping is a tuple that contains 4 items:

      * Mapped chain as an :class:`.AtomMap` instance, 
      * *chain* as an :class:`.AtomMap` instance,
      * Percent sequence identitity,
      * Percent sequence overlap
         
    Mappings are returned in decreasing percent sequence identity order.
    :class:`.AtomMap` that keeps mapped atom indices contains dummy atoms 
    in place of unmapped atoms.
    
    :arg atoms: atoms that will be mapped to the target *chain*
    :type atoms: :class:`.Chain`, :class:`.AtomGroup`, :class:`.Selection`
    
    :arg chain: chain to which atoms will be mapped
    :type chain: :class:`.Chain`
    
    :keyword seqid: percent sequence identity, default is 90
    :type seqid: float

    :keyword overlap: percent overlap, default is 90
    :type overlap: float

    :keyword pwalign: perform pairwise sequence alignment 
    :type pwalign: bool
    
    This function tries to map *atoms* to *chain* based on residue
    numbers and types. Each individual chain in *atoms* is compared to
    target *chain*. This works well for different structures of the same
    protein. When it fails, :mod:`Bio.pairwise2` is used for sequence
    alignment, and mapping is performed based on the sequence alignment.
    User can control, whether sequence alignment is performed or not with
    *pwalign* keyword. If ``pwalign=True`` is passed, pairwise alignment is 
    enforced."""
    
    """
    :keyword subset: "calpha" (or "ca"), "backbone" (or "bb"), or "all", 
        default is "calpha"  
    :type subset: string
    """
    
    target_chain = chain
    if not isinstance(atoms, (AtomGroup, Chain, Selection)):
        raise TypeError('atoms must be an AtomGroup, a Chain, or a '
                        'Selection instance')
    if not isinstance(target_chain, Chain):
        raise TypeError('chain must be Chain instance')
        
    subset = str(kwargs.get('subset', 'calpha')).lower()
    if subset not in _SUBSETS:
        raise ValueError('{0:s} is not a valid subset argument'
                         .format(str(subset)))
    seqid = kwargs.get('seqid', 90.)
    coverage = kwargs.get('overlap')
    if coverage is None:
        coverage = kwargs.get('coverage', 90.)
    pwalign = kwargs.get('pwalign', None)
    
    if isinstance(atoms, Chain):
        chains = [atoms]
        map_ag = atoms.getAtomGroup()
    else:
        if isinstance(atoms, AtomGroup):
            map_ag = atoms
        else:
            map_ag = atoms.getAtomGroup()
        chains = list(atoms.getHierView().iterChains())
        LOGGER.debug('Evaluating {0:s}: {1:d} chains are identified'
                     .format(str(atoms), len(chains)))
    
    if subset != 'all':
        target_chain = target_chain.select(subset
                                ).getHierView()[target_chain.getChid()]
    
    mappings = []
    unmapped = []
    target_ag = target_chain.getAtomGroup()
    simple_target = SimpleChain(target_chain, True)
    LOGGER.debug('Trying to map atoms based on residue numbers and '
                 'identities:')
    for chain in chains:
        simple_chain = SimpleChain(True)
        simple_chain.buildFromChain(chain)
        if len(simple_chain) == 0:
            LOGGER.debug('  Skipping {0:s}, which does not contain any amino '
                         'acid residues.'.format(simple_chain))
            continue
        LOGGER.debug('  Comparing {0:s} (len={1:d}) with {2:s}:'
                     .format(simple_chain.getTitle(), len(simple_chain), 
                             simple_target.getTitle()))
        
        target_list, chain_list, n_match, n_mapped = getTrivialMapping(
                                                simple_target, simple_chain)
        if n_mapped > 0:
            _seqid = n_match * 100 / n_mapped
            _cover = n_mapped * 100 / max(len(simple_target), len(simple_chain))
        else:
            _seqid = 0
            _cover = 0
        
        if _seqid >= seqid and _cover >= coverage:
            LOGGER.debug('\tMapped: {0:d} residues match with {1:.0f}% '
                         'sequence identity and {2:.0f}% overlap.'
                         .format(n_mapped, _seqid, _cover))
            mappings.append((target_list, chain_list, _seqid, _cover))
        else:
            LOGGER.debug('\tFailed to match chains based on residue numbers '
                         '(seqid={0:.0f}%, overlap={1:.0f}%).'
                        .format(_seqid, _cover))
            unmapped.append(simple_chain)


    if pwalign or (not mappings and (pwalign is None or pwalign)): 
        LOGGER.debug('Trying to map atoms based on {0:s} sequence alignment:'
                     .format(ALIGNMENT_METHOD))
        for simple_chain in unmapped:
            LOGGER.debug('  Comparing {0:s} (len={1:d}) with {2:s}:'
                         .format(simple_chain.getTitle(), len(simple_chain), 
                                 simple_target.getTitle()))
            result = getAlignedMapping(simple_target, simple_chain)
            if result is not None:
                target_list, chain_list, n_match, n_mapped = result
                if n_mapped > 0:
                    _seqid = n_match * 100 / n_mapped
                    _cover = n_mapped * 100 / max(len(simple_target), 
                                                  len(simple_chain))
                else:
                    _seqid = 0
                    _cover = 0
                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMapped: {0:d} residues match with {1:.0f}%'
                                 ' sequence identity and {2:.0f}% overlap.'
                                 .format(n_mapped, _seqid, _cover))
                    mappings.append((target_list, chain_list, _seqid, _cover))
                else:
                    LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, '
                                 'overlap={1:.0f}%).'
                                 .format(_seqid, _cover))
    
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
                    atom_chn = res_chn.getAtom(atom_tar.getName())
                    if atom_chn is not None:
                        indices_chain.append(atom_chn.getIndex())
                        indices_mapping.append(counter)
                    else:
                        indices_dummies.append(counter)
                else:
                    indices_dummies.append(counter)
                counter += 1
        #n_atoms = len(indices_target)   
        atommap = AM(map_ag, 
                     indices_chain,
                     indices_mapping,
                     indices_dummies,
                     simple_chain.getTitle() + ' -> ' + 
                                                    simple_target.getTitle(),
                     acsi=chain.getACSIndex())
        selection = AM(target_ag,
                       indices_target,
                       np.arange(len(indices_target)),
                       None,
                       simple_target.getTitle() + ' -> ' + 
                                       simple_chain.getTitle(),
                       acsi=target_chain.getACSIndex())
                                    
        mappings[mi] = (atommap, selection, _seqid, _cover)
    if len(mappings) > 1:
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        mappings.sort(compare, reverse=True)
    return mappings

def getTrivialMapping(target, chain):
    """Return lists of matching residues (map based on residue number)."""
    
    target_list = []
    chain_list = []
    n_match = 0
    n_mapped = 0
    chain_dict_get = chain._dict.get
    append = target_list.append 
    for target_residue in target:
        append(target_residue.getResidue())

        chain_residue = chain_dict_get((target_residue.getResnum(), 
                                        target_residue.getIcode()))
        if chain_residue is None:
            chain_list.append(chain_residue)
        else:
            if target_residue.getResname() == chain_residue.getResname():
                n_match += 1
            chain_list.append(chain_residue.getResidue())
            n_mapped += 1
            
    return target_list, chain_list, n_match, n_mapped

def getAlignedMapping(target, chain):
    pairwise2 = importBioPairwise2()
    if ALIGNMENT_METHOD == 'local':
        alignment = pairwise2.align.localms(target.getSequence(), 
                                            chain.getSequence(), 
                                            MATCH_SCORE, MISMATCH_SCORE,
                                            GAP_PENALTY,  GAP_EXT_PENALTY,
                                            one_alignment_only=1)
    else:
        alignment = pairwise2.align.globalms(target.getSequence(), 
                                             chain.getSequence(), 
                                             MATCH_SCORE, MISMATCH_SCORE,
                                             GAP_PENALTY, GAP_EXT_PENALTY,
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
        if a not in (GAP, NONE_A):
            ares = next(aiter)
            amatch.append(ares.getResidue())
            if b not in (GAP, NONE_A):
                bres = next(biter)
                bmatch.append(bres.getResidue())
                if a == b:
                    n_match += 1
                n_mapped += 1
            else:
                bmatch.append(None)
        elif b not in (GAP, NONE_A):
                bres = next(biter)
    return amatch, bmatch, n_match, n_mapped


if __name__ == '__main__':
    
    from prody import *
    p = parsePDB('1p38')
    r = parsePDB('1r39')
    chtwo, chone = mapOntoChain(r, p['A'])[0][:2]
    print countUnpairedBreaks(chone, chtwo)
