# -*- coding: utf-8 -*-
"""This module defines functions for comparing and mapping polypeptide chains.
"""

from numbers import Integral

import numpy as np
from numpy import arange

from prody.atomic import AtomMap as AM
from prody.atomic import AtomGroup, Chain, AtomSubset, Selection
from prody.atomic import AAMAP
from prody.atomic import flags
from prody.measure import calcTransformation, printRMSD, calcDistance, calcRMSD, superpose
from prody import LOGGER, SELECT, PY2K, PY3K
from prody.sequence import MSA
from prody.utilities import cmp, pystr, isListLike, multilap, SolutionDepletionException, index
from prody.utilities import MATCH_SCORE, MISMATCH_SCORE, GAP_PENALTY, GAP_EXT_PENALTY, ALIGNMENT_METHOD

from Bio import pairwise2

if PY2K:
    range = xrange

if PY3K:
    basestring = str

__all__ = ['matchChains', 'matchAlign', 'mapChainOntoChain', 'mapOntoChain', 'alignChains',
           'mapOntoChains', 'bestMatch', 'sameChid', 'userDefined', 'sameChainPos',
           'mapOntoChainByAlignment', 'getMatchScore', 'setMatchScore',
           'getMismatchScore', 'setMismatchScore', 'getGapPenalty', 
           'setGapPenalty', 'getGapExtPenalty', 'setGapExtPenalty',
           'getGoodSeqId', 'setGoodSeqId', 'getGoodCoverage', 'combineAtomMaps',
           'setGoodCoverage', 'getAlignmentMethod', 'setAlignmentMethod']

GOOD_SEQID = 90.
GOOD_COVERAGE = 90.

GAPCHARS = ['-', '.']
NONE_A = '_'

_a2aaa = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
    'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
    'Y': 'TYR', 'V': 'VAL'
}

def calcScores(n_match, n_mapped, n_total):
    if n_mapped > 0:
        seqid = n_match * 100 / n_mapped
        cover = n_mapped * 100 / n_total
    else:
        seqid = cover = 0
    return seqid, cover


def getGoodSeqId():
    """Returns good sequence identity."""

    return GOOD_SEQID


def setGoodSeqId(seqid):
    """Set good sequence identity."""

    if isinstance(seqid, (float, int)) and seqid >= 0:
        global GOOD_SEQID
        GOOD_SEQID = seqid
    else:
        raise TypeError('seqid must be a positive number or zero')

def getGoodCoverage():
    """Returns good sequence coverage."""

    return GOOD_COVERAGE


def setGoodCoverage(coverage):
    """Set good sequence coverage."""

    if isinstance(coverage, (float, int)) and coverage >= 0:
        global GOOD_COVERAGE
        GOOD_COVERAGE = coverage
    else:
        raise TypeError('coverage must be a positive number or zero')

def getMatchScore():
    """Returns match score used to align sequences."""

    return MATCH_SCORE


def setMatchScore(match_score):
    """Set match score used to align sequences."""

    if isinstance(match_score, (float, int)) and match_score >= 0:
        global MATCH_SCORE
        MATCH_SCORE = match_score
    else:
        raise TypeError('match_score must be a positive number or zero')


def getMismatchScore():
    """Returns mismatch score used to align sequences."""

    return MISMATCH_SCORE


def setMismatchScore(mismatch_score):
    """Set mismatch score used to align sequences."""

    if isinstance(mismatch_score, (float, int)) and mismatch_score >= 0:
        global MISMATCH_SCORE
        MISMATCH_SCORE = mismatch_score
    else:
        raise TypeError('mismatch_score must be a positive number or zero')


def getGapPenalty():
    """Returns gap opening penalty used for pairwise alignment."""

    return GAP_PENALTY


def setGapPenalty(gap_penalty):
    """Set gap opening penalty used for pairwise alignment."""

    if isinstance(gap_penalty, (float, int)) and gap_penalty <= 0:
        global GAP_PENALTY
        GAP_PENALTY = gap_penalty
    else:
        raise TypeError('gap_penalty must be a negative number')


def getGapExtPenalty():
    """Returns gap extension penalty used for pairwise alignment."""

    return GAP_EXT_PENALTY


def setGapExtPenalty(gap_ext_penalty):
    """Set gap extension penalty used for pairwise alignment."""

    if isinstance(gap_ext_penalty, (float, int)) and gap_ext_penalty <= 0:
        global GAP_EXT_PENALTY
        GAP_EXT_PENALTY = gap_ext_penalty
    else:
        raise TypeError('gap_ext_penalty must be a negative number or zero')


def getAlignmentMethod():
    """Returns pairwise alignment method."""

    return ALIGNMENT_METHOD


def setAlignmentMethod(method):
    """Set pairwise alignment method (global or local)."""

    if method in ('local', 'global'):
        global ALIGNMENT_METHOD
        ALIGNMENT_METHOD = method
    else:
        raise ValueError('method must be "local" or "global"')


class SimpleResidue(object):

    __slots__ = ['_chain', '_index', '_res', '_name', '_num', '_inc']

    def __init__(self, chain, index, number, name, insertioncode='', residue=None):
        self._num = number
        self._name = name
        self._inc = insertioncode
        self._res = residue
        self._chain = chain
        self._index = index

    def __repr__(self):
        return '<SimpleResidue: {0}{1}>'.format(self._name, self._num)

    def __iter__(self):
        res = self.getResidue()

        if res is self:
            # iterate the residue itself as atoms in it, assuming the 
            # residue contains a single alpha carbon
            for atom in [self]:
                yield atom
        else:
            for atom in res:
                yield atom

    def getResidue(self):
        if self._res is not None:
            return self._res
        else:
            return self

    def getChain(self):
        return self._chain

    def getResnum(self):
        return self._num

    def getIndex(self):
        return self._index

    def getIcode(self):
        return self._inc

    def getResname(self):
        return self._name

    def getCoords(self):
        if self._res:
            return self._res._getCoords()
        return None

    def getName(self):
        "A dummy method that returns the atom name of alpha carbon of this residue."

        return 'CA'

class SimpleChain(object):

    """An internal class used to compare two polypeptide chains.


    SimpleChain instances can be indexed using residue numbers. If a residue
    with given number is not found in the chain, **None** is returned."""

    __slots__ = ['_list', '_seq', '_title', '_dict', '_gaps', '_chain']

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
        self._title = ''
        self._gaps = allow_gaps
        self._chain = None
        if isinstance(chain, Chain):
            self.buildFromChain(chain)
        elif isinstance(chain, str):
            self.buildFromSequence(chain)

    def __len__(self):
        return len(self._list)

    def __iter__(self):
        return self._list.__iter__()

    def __repr__(self):
        return '<SimpleChain: {0} with {1} residues>'.format(
            self._title, len(self._list))

    def __str__(self):
        return '{0} with {1} residues'.format(self._title, len(self._list))

    def __getitem__(self, index):
        if isinstance(index, Integral):
            return self._dict.get((index, ''))
        return self._dict.get(index)

    def getSequence(self):
        return self._seq

    def getTitle(self):
        return self._title

    def getCoords(self, calpha=False):
        if self._chain is None:
            return None
        if calpha:
            return self._chain.ca._getCoords()
        else:
            return self._chain._getCoords()

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
            resnums = arange(1, len(sequence)+1)
        resid = 0
        gaps = self._gaps
        for i, aa in enumerate(sequence):
            resid = resnums[i]
            if gaps and aa in GAPCHARS:
                aa = NONE_A
                simpres = None 
            else:
                simpres = SimpleResidue(self, i, resid, aa)

            self._list.append(simpres)
            self._dict[resid] = simpres
            self._seq += aa
        
        self._title = 'built from sequence %s...'%sequence[:5]

    def buildFromChain(self, chain):
        """Build from a :class:`.Chain`."""

        assert isinstance(chain, Chain), 'chain must be a Chain instance'
        self._chain = chain
        gaps = self._gaps
        residues = list(chain.iterResidues())
        temp = residues[0].getResnum()-1
        protein_resnames = flags.AMINOACIDS
        for i, res in enumerate(chain):
            resid = res.getResnum()
            incod = res.getIcode()
            aa = AAMAP.get(res.getResname(), 'X')
            if aa == '-':
                aa = 'X'
            if aa == 'X':
                LOGGER.warn("no one-letter mapping found for %s " % repr(res))
            simpres = SimpleResidue(self, i, resid, aa, incod, res)
            if gaps:
                diff = resid - temp - 1
                if diff > 0:
                    self._seq += NONE_A * diff
                temp = resid
            self._seq += aa
            self._list.append(simpres)
            self._dict[(resid, incod)] = simpres
        self._title = 'Chain {0} from {1}'.format(chain.getChid(),
                                                  chain.getAtomGroup()
                                                  .getTitle())


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

    mapped = chone.getFlags('mapped') * chtwo.getFlags('mapped')
    if mapped.sum() == 0:
        raise ValueError('chains do not have common mapped atoms')
    chone = chone[mapped]
    chtwo = chtwo[mapped]

    rnone = chone.getResnums()
    rntwo = chtwo.getResnums()

    brone = calcDistance(chone[1:], chone[:-1]) > 4.
    brtwo = calcDistance(chtwo[1:], chtwo[:-1]) > 4.

    brone[(rnone[1:] - rnone[:-1]) > 1] = False
    brtwo[(rntwo[1:] - rntwo[:-1]) > 1] = False

    brone = set(brone.nonzero()[0])
    brtwo = set(brtwo.nonzero()[0])
    return len(brone.union(brtwo)) - len(brone.intersection(brtwo))


_SUBSETS = {
    'ca': 'ca',
    'calpha': 'ca',
    'bb': 'bb',
    'backbone': 'bb',
    'heavy': 'noh',
    'noh': 'noh',
    'all': 'all'
}


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
            raise ValueError('selection {0} did not match any atoms'
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
        csets = range(mobile.numCoordsets())  # PY3K: OK
    else:
        csets = [mobile.getACSIndex()]

    LOGGER.info('Alignment is based on {0} atoms matching {1}.'
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
    """Returns pairs of chains matched based on sequence similarity.  Makes an
    all-to-all comparison of chains in *atoms1* and *atoms2*.  Chains are
    obtained from hierarchical views (:class:`.HierView`) of atom groups.
    This function returns a list of matching chains in a tuple that contain
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

    :keyword subset: one of the following well-defined subsets of atoms:
        ``"calpha"`` (or ``"ca"``), ``"backbone"`` (or ``"bb"``),
        ``"heavy"`` (or ``"noh"``), or ``"all"``, default is ``"calpha"``
    :type subset: str

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
        raise ValueError('{0} is not a valid subset argument'
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
        LOGGER.debug('Checking {0}: {1} chains are identified'
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
        LOGGER.debug('Checking {0}: {1} chains are identified'
                     .format(str(atoms2), len(chains2)))

    matches = []
    unmatched = []
    if pwalign:
        for simpch1 in chains1:
            for simpch2 in chains2:
                unmatched.append((simpch1, simpch2))
    else:
        LOGGER.debug('Trying to match chains based on residue numbers and names:')
        for simpch1 in chains1:
            for simpch2 in chains2:
                LOGGER.debug('  Comparing {0} (len={1}) and {2} (len={3}):'
                            .format(simpch1.getTitle(), len(simpch1),
                                    simpch2.getTitle(), len(simpch2)))

                match1, match2, nmatches = getTrivialMatch(simpch1, simpch2)
                _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
                _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))

                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMatch: {0} residues match with {1:.0f}% '
                                'sequence identity and {2:.0f}% overlap.'
                                .format(len(match1), _seqid, _cover))
                    matches.append((match1, match2, _seqid, _cover, simpch1, simpch2))
                else:
                    LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, '
                                'overlap={1:.0f}%).'.format(_seqid, _cover))
                    unmatched.append((simpch1, simpch2))

    if pwalign or (not matches and (pwalign is None or pwalign)):
        if pairwise2:
            if unmatched:
                LOGGER.debug('Trying to match chains based on {0} sequence '
                            'alignment:'.format(ALIGNMENT_METHOD))
                for simpch1, simpch2 in unmatched:
                    LOGGER.debug(' Comparing {0} (len={1}) and {2} '
                                '(len={3}):'
                                .format(simpch1.getTitle(), len(simpch1),
                                        simpch2.getTitle(), len(simpch2)))
                    match1, match2, nmatches = getAlignedMatch(simpch1, simpch2)
                    _seqid = nmatches * 100 / min(len(simpch1), len(simpch2))
                    _cover = len(match2) * 100 / max(len(simpch1), len(simpch2))
                    if _seqid >= seqid and _cover >= coverage:
                        LOGGER.debug('\tMatch: {0} residues match with {1:.0f}% '
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
    subset = _SUBSETS[subset]
    for mi, result in enumerate(matches):
        match1, match2, _seqid, _cover, simpch1, simpch2 = result

        indices1 = []
        indices2 = []

        for i in range(len(match1)):
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
                    except ValueError:
                        continue
                    else:
                        indices1.append(ares._indices[aid])
                        indices2.append(bres._indices[bid])
            elif subset == 'noh':
                for han, aid, noh in zip(ares.getNames(), ares._indices,
                                         ares.getFlags('noh')):
                    if not noh:
                        continue
                    try:
                        bid = bres.getNames().tolist().index(han)
                    except ValueError:
                        continue
                    else:
                        indices1.append(aid)
                        indices2.append(bres._indices[bid])
            elif subset is None or subset == 'all':
                aans = ares.getNames()
                bans = bres.getNames().tolist()

                aids = ares.getIndices()
                #bids = bres.getIndices()

                for j in range(len(aans)):
                    try:
                        bid = bres._indices[bans.index(aans[j])]
                        indices1.append(aids[j])
                        indices2.append(bid)
                    except ValueError:
                        pass

        indices1 = np.array(indices1, int)
        indices2 = np.array(indices2, int)

        match1 = AM(atoms1, indices1, atoms1.getACSIndex(),
                    title=simpch1.getTitle() + ' -> ' + simpch2.getTitle(),
                    intarrays=True)
        match2 = AM(atoms2, indices2, atoms2.getACSIndex(),
                    title=simpch2.getTitle() + ' -> ' + simpch1.getTitle(),
                    intarrays=True)

        matches[mi] = (match1, match2, _seqid, _cover)
    if len(matches) > 1:
        matches.sort(key=lambda m: m[-2:], reverse=True)
    return matches


def getTrivialMatch(ach, bch):
    """Returns lists of matching residues (match based on residue number).

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
    """Returns list of matching residues (match is based on sequence alignment).
    """

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

    amatch = []
    bmatch = []
    match = 0.0
    try:
        this = alignment[0][0]
        that = alignment[0][1]
    except IndexError:
        LOGGER.warning('Matching chains resulted in empty alignment.')
        return amatch, bmatch, match
    aiter = ach.__iter__()
    biter = bch.__iter__()
    for i in range(len(this)):
        a = this[i]
        b = that[i]
        if a not in GAPCHARS:
            ares = next(aiter)
        if b not in GAPCHARS:
            bres = next(biter)
            if a not in GAPCHARS:
                amatch.append(ares.getResidue())
                bmatch.append(bres.getResidue())
                if a == b:
                    match += 1
    return amatch, bmatch, match


def mapOntoChain(atoms, chain, **kwargs):
    """Map *atoms* onto *chain*. This function is a wrapper of 
    :func:`.mapChainOntoChain` that manages to map chains onto target *chain*. 
    The function returns a list of mappings. Each mapping is a tuple that 
    contains 4 items:

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

    :keyword subset: one of the following well-defined subsets of atoms:
        ``"calpha"`` (or ``"ca"``), ``"backbone"`` (or ``"bb"``),
        ``"heavy"`` (or ``"noh"``), or ``"all"``, default is ``"calpha"``
    :type subset: str

    See :func:`.mapChainOntoChain` for other keyword arguments. 
    This function tries to map *atoms* to *chain* based on residue
    numbers and types. Each individual chain in *atoms* is compared to
    target *chain*.
    
    .. [IS98] Shindyalov IN, Bourne PE. Protein structure alignment by 
       incremental combinatorial extension (CE) of the optimal path. 
       *Protein engineering* **1998** 11(9):739-47.
    """

    if not isinstance(atoms, (AtomGroup, AtomSubset)):
        raise TypeError('atoms must be an AtomGroup or a AtomSubset (Chain, '
                        'Segment, etc.) instance')
    if not isinstance(chain, Chain):
        raise TypeError('chain must be Chain instance')

    subset = str(kwargs.get('subset', 'all')).lower()
    if subset not in _SUBSETS:
        raise ValueError('{0} is not a valid subset argument'
                         .format(str(subset)))

    if subset != 'all':
        chid = chain.getChid()
        segname = chain.getSegname()
        chain_subset = chain.select(subset)
        target_chain = chain_subset.getHierView()[segname, chid]
        
        mobile = atoms.select(subset)
    else:
        target_chain = chain
        mobile = atoms

    if isinstance(mobile, Chain):
        chains = [mobile]
    else:
        chains = list(mobile.getHierView().iterChains())
        LOGGER.debug('Evaluating {0}: {1} chains are identified'
                     .format(str(atoms), len(chains))) 

    mappings = []
    simple_target = SimpleChain(target_chain, False)
    LOGGER.debug('Trying to map atoms based on residue numbers and '
                'identities:')
    for chain in chains:
        simple_chain = SimpleChain(chain, False)
        mapping = mapChainOntoChain(simple_chain, simple_target, **kwargs)
        if mapping is not None:
            mappings.append(mapping)

    if len(mappings) > 1:
        mappings.sort(key=lambda m: m[-2]*m[-1], reverse=True)
    return mappings

def mapChainOntoChain(mobile, target, **kwargs):
    """Map *mobile* chain onto *target* chain.  This function returns a mapping that 
    contains 4 items:

      * Mapped chain as an :class:`.AtomMap` instance,
      * *chain* as an :class:`.AtomMap` instance,
      * Percent sequence identitity,
      * Percent sequence overlap

    Mappings are returned in decreasing percent sequence identity order.
    :class:`.AtomMap` that keeps mapped atom indices contains dummy atoms
    in place of unmapped atoms.

    :arg mobile: mobile that will be mapped to the *target* chain
    :type mobile: :class:`.Chain`

    :arg target: chain to which atoms will be mapped
    :type target: :class:`.Chain`

    :keyword seqid: percent sequence identity, default is **90**. Note that this parameter is 
        only effective for sequence alignment
    :type seqid: float

    :keyword overlap: percent overlap with *target*, default is **70**
    :type overlap: float

    :keyword mapping: what method will be used if the trivial mapping based on residue numbers 
        fails. If ``"ce"`` or ``"cealign"``, then the CE algorithm [IS98]_ will be 
        performed. It can also be a list of prealigned sequences, a :class:`.MSA` instance,
        or a dict of indices such as that derived from a :class:`.DaliRecord`.
        If set to **True** then the sequence alignment from :mod:`~Bio.pairwise2` 
        will be used. If set to **False**, only the trivial mapping will be performed. 
        Default is **"auto"**
    :type mapping: list, str, bool

    :keyword pwalign: if **True**, then pairwise sequence alignment will 
        be performed. If **False** then a simple mapping will be performed 
        based on residue numbers (as well as insertion codes). This will be 
        overridden by the *mapping* keyword's value.
    :type pwalign: bool
    
    .. [IS98] Shindyalov IN, Bourne PE. Protein structure alignment by 
       incremental combinatorial extension (CE) of the optimal path. 
       *Protein engineering* **1998** 11(9):739-47.
    """

    if isinstance(mobile, Chain):
        simple_mobile = SimpleChain(mobile, False)
    elif isinstance(mobile, SimpleChain):
        simple_mobile = mobile
        mobile = mobile._chain
    else:
        raise TypeError('mobile must be a Chain instance') 

    if len(simple_mobile) == 0:
        LOGGER.debug('\tCannot process {0}, which does not contain any amino '
                    'acid residues.'.format(simple_mobile))
        return None

    if isinstance(target, Chain):
        simple_target = SimpleChain(target, False)
    elif isinstance(target, SimpleChain):
        simple_target = target
        target = target._chain
    else:
        raise TypeError('target must be a Chain instance') 

    seqid = kwargs.get('seqid', 90.) 
    coverage = kwargs.get('overlap', 70.)
    coverage = kwargs.get('coverage', coverage) 
    pwalign = kwargs.get('pwalign', 'auto')
    pwalign = kwargs.get('mapping', pwalign)
    alignment = None

    if isinstance(pwalign, basestring):
        pwalign = pystr(pwalign).strip().lower()
    elif not isinstance(pwalign, bool):
        alignment = pwalign
        pwalign = True
    
    map_ag = mobile.getAtomGroup() if mobile is not None else None
    target_ag = target.getAtomGroup() if target is not None else None

    if map_ag is None and target_ag is None:
        raise ValueError('At least one of mobile and target should be a Chain object '
                         'or a SimpleChain object associated with a Chain object.')

    mapping = None
    LOGGER.debug('Trying to map atoms based on residue numbers and '
            'identities:')
    LOGGER.debug('  Comparing {0} (len={1}) with {2}:'
                .format(simple_mobile.getTitle(), len(simple_mobile),
                        simple_target.getTitle()))

    # trivial mapping serves as a first simple trial of alignment the two 
    # sequences based on residue number, therefore the sequence identity 
    # (GOOD_SEQID) criterion is strict.
    target_list, chain_list, n_match, n_mapped = getTrivialMapping(
        simple_target, simple_mobile)
    _seqid, _cover = calcScores(n_match, n_mapped, len(simple_target))

    trivial_seqid = GOOD_SEQID if pwalign else seqid
    trivial_cover = GOOD_COVERAGE if pwalign else coverage
    if _seqid >= trivial_seqid and _cover >= trivial_cover:
        LOGGER.debug('\tMapped: {0} residues match with {1:.0f}% '
                'sequence identity and {2:.0f}% overlap.'
                .format(n_mapped, _seqid, _cover))
        mapping = (target_list, chain_list, _seqid, _cover)
    else:
        if not pwalign:
            LOGGER.debug('\tFailed to match chains based on residue numbers '
                    '(seqid={0:.0f}%, overlap={1:.0f}%).'
                    .format(_seqid, _cover))

    if pwalign and mapping is None:
        SEQ_ALIGNMENT = ('seq', ALIGNMENT_METHOD + ' sequence alignment', seqid, coverage)
        CE_ALIGNMENT = ('ce', 'CEalign', 0., coverage)

        if not 'seqid' in kwargs:
            tar_seqid = 0.
        else:
            tar_seqid = seqid
        PREDEF_ALIGNMENT = ('predef', 'predefined alignment', tar_seqid, coverage)

        if alignment is None:
            if pwalign in ['ce', 'cealign']:
                methods = [CE_ALIGNMENT]
            elif pwalign == 'auto': 
                methods = [SEQ_ALIGNMENT, 
                           CE_ALIGNMENT]
            else:
                methods = [SEQ_ALIGNMENT]
        else:
            methods = [PREDEF_ALIGNMENT]

        for method, desc, seqid, coverage in methods:
            LOGGER.debug('Trying to map atoms based on {0}:'
                        .format(desc))

            LOGGER.debug('  Comparing {0} (len={1}) with {2}:'
                        .format(simple_mobile.getTitle(), len(simple_mobile),
                                simple_target.getTitle()))
            if method == 'ce':
                result = getCEAlignMapping(simple_target, simple_mobile)
            elif method == 'seq':
                result = getAlignedMapping(simple_target, simple_mobile)
            else:
                if isinstance(alignment, dict):
                    result = getDictMapping(simple_target, simple_mobile, alignment)
                else:
                    result = getAlignedMapping(simple_target, simple_mobile, alignment)

            if result is not None:
                target_list, chain_list, n_match, n_mapped = result
                _seqid, _cover = calcScores(n_match, n_mapped, max(len(simple_target),
                                                                   len(simple_mobile)))

                if _seqid >= seqid and _cover >= coverage:
                    LOGGER.debug('\tMapped: {0} residues match with {1:.0f}%'
                                    ' sequence identity and {2:.0f}% overlap.'
                                    .format(n_mapped, _seqid, _cover))
                    mapping = (target_list, chain_list, _seqid, _cover)
                    break
                else:
                    LOGGER.debug('\tFailed to match chains (seqid={0:.0f}%, '
                                    'overlap={1:.0f}%).'
                                    .format(_seqid, _cover))

    if mapping is not None:
        residues_target, residues_chain, _seqid, _cover = mapping
        indices_target = []
        indices_chain = []
        indices_mapping = []
        indices_dummies = []
        counter = 0
        for i in range(len(residues_target)):
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

        ch_tar = next((r for r in residues_target if r is not None)).getChain()
        ch_chn = next((r for r in residues_chain if r is not None)).getChain()

        if isinstance(ch_tar, Chain):
            title_tar = 'Chain {0} from {1}'.format(ch_tar.getChid(), ch_tar.getAtomGroup().getTitle())
        else:
            title_tar = 'SimpleChain {0}'.format(ch_tar.getTitle())

        if isinstance(ch_chn, Chain):
            title_chn = 'Chain {0} from {1}'.format(ch_chn.getChid(), ch_chn.getAtomGroup().getTitle())
        else:
            title_chn = 'SimpleChain {0}'.format(ch_tar.getTitle())

        # note that chain here is from atoms
        if map_ag is not None:
            atommap = AM(map_ag, indices_chain, mobile.getACSIndex(),
                         mapping=indices_mapping, dummies=indices_dummies,
                         title=title_chn + ' -> ' + title_tar)
        else:
            atommap = None

        if target_ag is not None:
            selection = AM(target_ag, indices_target, target.getACSIndex(),
                           title=title_tar + ' -> ' + title_chn, intarrays=True)
        else:
            selection = None

        mapping = (atommap, selection, _seqid, _cover)
    return mapping

def userDefined(chain1, chain2, correspondence):
    id1 = chain1.getTitle()
    id2 = chain2.getTitle()

    if not isinstance(correspondence, dict):
        chmap = {}
        try:
            chmap[id1] = correspondence[0]
            chmap[id2] = correspondence[1]
        except (IndexError, TypeError):
            raise TypeError('correspondence should be a dict with keys being titles of atoms and ref, '
                            'and values are str indicating chID correspondences')

    corr1 = correspondence[id1]
    corr2 = correspondence[id2]

    if len(corr1) != len(corr2):
        raise ValueError('%s and %s have different number of chain identifiers '
                         'in the correspondence'%(id1, id2))

    try:
        i = corr1.index(chain1.getChid())
        chid = corr2[i]
    except:
        return False

    return chain2.getChid() == chid

def sameChid(chain1, chain2):
    return chain1.getChid() == chain2.getChid()

def sameChainPos(chain1, chain2):
    chids_arr1 = np.array([chain.getChid() 
                           for chain in list(chain1.getAtomGroup().getHierView())])
    position1 = np.where(chids_arr1 == chain1.getChid())[0][0]
    
    chids_arr2 = np.array([chain.getChid() 
                           for chain in list(chain2.getAtomGroup().getHierView())])
    position2 = np.where(chids_arr2 == chain2.getChid())[0][0]
                         
    return position1 == position2

def bestMatch(chain1, chain2):
    return True

def mapOntoChains(atoms, ref, match_func=bestMatch, **kwargs):
    """This function is a generalization and wrapper of :func:`.mapOntoChain` that 
    manages to map chains onto chains (instead of a single chain). 
    
    :arg atoms: atoms to map onto the reference
    :type atoms: :class:`.Atomic`
    
    :arg ref: reference structure for mapping
    :type ref: :class:`.Atomic`

    :arg match_func: function determines which chains from ``ref`` and ``atoms`` are matched.
        Default is to use the best match.
    :type match_func: func
    """
    
    if not isinstance(atoms, (SimpleChain, AtomGroup, AtomSubset)):
        raise TypeError('atoms must be an AtomGroup or a AtomSubset (Chain, '
                        'Segment, etc.) instance')
    if not isinstance(ref, (SimpleChain, AtomGroup, AtomSubset)):
        raise TypeError('ref must be an AtomGroup or a AtomSubset (Chain, '
                        'Segment, etc.) instance')

    subset = str(kwargs.get('subset', 'all')).lower()
    if subset not in _SUBSETS:
        raise ValueError('{0} is not a valid subset argument'
                         .format(str(subset)))

    if subset != 'all':
        if not isinstance(ref, SimpleChain):
            target = ref.select(subset)
        else:
            target = ref
        if not isinstance(atoms, SimpleChain):
            mobile = atoms.select(subset)
        else:
            mobile = atoms
    else:
        target = ref
        mobile = atoms

    if isinstance(mobile, (SimpleChain, Chain)):
        chs_atm = [mobile]
    else:
        chs_atm = [chain for chain in mobile.getHierView().iterChains()]
    
    if isinstance(target, (SimpleChain, Chain)):
        chs_ref = [target]
    else:
        chs_ref = [chain for chain in target.getHierView().iterChains()]

    # iterate through chains of both target and mobile
    mappings = np.empty((len(chs_ref), len(chs_atm)), dtype='O')
    for i, chain in enumerate(chs_ref):
        simple_chain = chain if isinstance(chain, SimpleChain) else SimpleChain(chain, False) 
        for j, target_chain in enumerate(chs_atm):
            if not match_func(chain, target_chain):
                continue

            simple_target = target_chain if isinstance(target_chain, SimpleChain) else SimpleChain(target_chain, False)
            mappings[i, j] = mapChainOntoChain(simple_target, simple_chain, **kwargs)

    return mappings

def mapOntoChainByAlignment(atoms, chain, **kwargs):
    """This function is similar to :func:`.mapOntoChain` but correspondence 
    of chains is found by alignment provided. 
    
    :arg alignments: A list of predefined alignments. It can be also a 
        dictionary or :class:`MSA` instance where the keys or 
        labels are the title of *atoms* or *chains*. 
    :type alignments: list, dict, :class:`MSA`
    """

    alignments = kwargs.pop('alignments', None)
    if alignments is None:
        return mapOntoChain(atoms, chain, **kwargs)
    else:
        if isinstance(alignments, (MSA, dict)):
            refseq = str(alignments[chain.getTitle()])
            tarseq = str(alignments[atoms.getTitle()])
            alignment = [refseq, tarseq]
        else:
            index = kwargs.pop('index', 0)
            alignment = alignments[index]

        tar_aligned_seq = alignment[-1]
        for char in GAPCHARS:
            tar_aligned_seq = tar_aligned_seq.replace(char, '').upper()
        hv = atoms.getHierView()
        for target_chain in hv.iterChains():
            tar_seq = target_chain.getSequence().upper()
            if tar_seq == tar_aligned_seq:
                mappings = mapOntoChain(target_chain, chain, pwalign=alignment, **kwargs)
                return mappings
        LOGGER.warn('The sequence of chain does not match that in alignment (%s).'%atoms.getTitle())
    return []

def getTrivialMapping(target, chain):
    """Returns lists of matching residues (map based on residue number)."""

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

def getDictMapping(target, chain, map_dict):
    """Returns lists of matching residues (based on *map_dict*)."""

    pdbid = chain._chain.getTitle()[:4].lower()
    chid = chain._chain.getChid().upper()
    key = pdbid + chid

    mapping = map_dict.get(key)
    if mapping is None:
        LOGGER.warn('map_dict does not have the mapping for {0}'.format(key))
        return None

    tar_indices = mapping[0]
    chn_indices = mapping[1]

    chain_res_list = [res for res in chain]

    amatch = []
    bmatch = []
    n_match = 0
    n_mapped = 0
    for i, a in enumerate(target):  
        ares = a.getResidue()
        amatch.append(ares)
        if i in tar_indices:
            try:
                n = index(tar_indices, i)
            except IndexError:
                LOGGER.warn('\nthe number of residues in the map_dict ({0} residues) is inconsistent with {2} ({1} residues)'
                            .format(max(tar_indices)+1, len(chain_res_list), target.getTitle()))
                return None
            try:
                b = chain_res_list[chn_indices[n]]
            except IndexError:
                LOGGER.warn('\nthe number of residues in the map_dict ({0} residues) is inconsistent with {2} ({1} residues)'
                            .format(max(chn_indices)+1, len(chain_res_list), chain.getTitle()))
                return None
            bres = b.getResidue()
            bmatch.append(bres)
            if a.getResname() == b.getResname():
                n_match += 1
            n_mapped += 1
        else:
            bmatch.append(None)

    return amatch, bmatch, n_match, n_mapped

def getAlignedMapping(target, chain, alignment=None):
    """Returns lists of matching residues (map based on pairwise 
    alignment or predefined alignment)."""

    if alignment is None:
        if ALIGNMENT_METHOD == 'local':
            alignments = pairwise2.align.localms(target.getSequence(),
                                                chain.getSequence(),
                                                MATCH_SCORE, MISMATCH_SCORE,
                                                GAP_PENALTY,  GAP_EXT_PENALTY,
                                                one_alignment_only=1)
        else:
            alignments = pairwise2.align.globalms(target.getSequence(),
                                                chain.getSequence(),
                                                MATCH_SCORE, MISMATCH_SCORE,
                                                GAP_PENALTY, GAP_EXT_PENALTY,
                                                one_alignment_only=1)
        alignment = alignments[0]
        this, that = alignment[:2]
    else:
        def _findAlignment(sequence, alignment):
            for seq in alignment:
                strseq = str(seq).upper()
                for gap in GAPCHARS:
                    strseq = strseq.replace(gap, '')

                if sequence.upper() == strseq:
                    return str(seq)
            return None

        this = _findAlignment(target.getSequence(), alignment)
        if this is None:
            LOGGER.warn('alignment does not contain the target ({0}) sequence'
                        .format(target.getTitle()))
            return None

        that = _findAlignment(chain.getSequence(), alignment)
        if that is None:
            LOGGER.warn('alignment does not contain the chain ({0}) sequence'
                        .format(chain.getTitle()))
            return None

    amatch = []
    bmatch = []
    aiter = target.__iter__()
    biter = chain.__iter__()
    n_match = 0
    n_mapped = 0
    gap_chars = list(GAPCHARS)
    gap_chars.append(NONE_A)
    for i in range(len(this)):
        a = this[i]
        b = that[i]
        if a not in gap_chars:
            ares = next(aiter)
            amatch.append(ares.getResidue())
            if b not in gap_chars:
                bres = next(biter)
                bmatch.append(bres.getResidue())
                if a == b:
                    n_match += 1
                n_mapped += 1
            else:
                bmatch.append(None)
        elif b not in gap_chars:
            bres = next(biter)
    return amatch, bmatch, n_match, n_mapped

def getCEAlignMapping(target, chain):
    try:
        from .ccealign import ccealign
    except ImportError:
        LOGGER.warn('Could not import ccealign C/C++ extension.'
                    'It may not be installed properly.')
        return None

    if not ("X" in target.getSequence() or "X" in chain.getSequence()):
        calpha=True
    else:
        calpha=False

    tar_coords = target.getCoords(calpha=calpha).tolist()
    mob_coords = chain.getCoords(calpha=calpha).tolist()

    if len(tar_coords) < 8:
        LOGGER.warn('target ({1}) is too small to be aligned '
                    'by CE algorithm (at least {0} residues)'
                    .format(8, repr(target)))
        return None

    if len(mob_coords) < 8:
        LOGGER.warn('chain ({1}) is too small to be aligned '
                    'by CE algorithm (at least {0} residues)'
                    .format(8, repr(chain)))
        return None

    try:
        aln_info = ccealign((tar_coords, mob_coords))
    except:
        LOGGER.warn('cealign could not align {0} and {1}'.format(repr(target), repr(chain)))
        return None

    paths, bestIdx, nres, rmsd = aln_info[:4]
    path = paths[bestIdx]

    tar_indices = []
    chn_indices = []
    for i, j in path:
        # if i not in tar_dummies:
        #     if j not in mob_dummies:
        tar_indices.append(i)
        chn_indices.append(j)

    chain_res_list = [res for res in chain]

    amatch = []
    bmatch = []
    n_match = 0
    n_mapped = 0
    for i, a in enumerate(target):
        ares = a.getResidue()
        amatch.append(ares)
        if i in tar_indices:
            n = tar_indices.index(i)
            try:
                b = chain_res_list[chn_indices[n]]
            except IndexError:
                bmatch.append(None)
                continue
            bres = b.getResidue()
            bmatch.append(bres)
            if a.getResname() == b.getResname():
                n_match += 1
            n_mapped += 1
        else:
            bmatch.append(None)

    return amatch, bmatch, n_match, n_mapped

def combineAtomMaps(mappings, target=None, **kwargs):
    """Builds a grand :class:`.AtomMap` instance based on *mappings* obtained from 
    :func:`.mapOntoChains`. The function also accepts the output :func:`.mapOntoChain` 
    but will trivially return all the :class:`.AtomMap` in *mappings*. 
    *mappings* should be a list or an array of matching chains in a tuple that contain
    4 items:

      * matching chain from *atoms1* as a :class:`.AtomMap` instance,
      * matching chain from *atoms2* as a :class:`.AtomMap` instance,
      * percent sequence identity of the match,
      * percent sequence overlap of the match.

    :arg mappings: a list or an array of matching chains in a tuple, or just the tuple
    :type mappings: tuple, list, :class:`~numpy.ndarray`

    :arg target: reference structure for superposition and checking RMSD
    :type target: :class:`.Atomic`

    :arg drmsd: amount deviation of the RMSD with respect to the top ranking atommap. 
        This is to allow multiple matches when *mobile* has more chains than *target*. 
        Default is 3.0
    :type drmsd: float

    :arg rmsd_reject: upper RMSD cutoff that rejects an atommap. Default is 15.0
    :type rmsd_reject: float

    :arg least: the least number of atommaps requested. If **None**, it will be automatically 
        determined by the number of chains present in *target* and *mobile*. 
        Default is **None**
    :type least: int

    :arg debug: a container (dict) that saves the following information for debugging purposes:
        * coverage: original coverage matrix, rows and columns correspond to the reference and the 
        mobile, respectively,
        * solutions: matched index groups that obtained by modeling the coverage matrix as a linear 
        assignment problem,
        * rmsd: a list of ranked RMSDs of identified atommaps.
    :type debug: dict

    """

    BIG_NUMBER = 1e6
    drmsd = kwargs.pop('drmsd', 3.)
    debug = kwargs.pop('debug', {})
    reject_rmsd = kwargs.pop('rmsd_reject', 15.)
    least_n_atommaps = kwargs.pop('least', None)
    
    def _build(mappings, nodes=[]):
        m, n = mappings.shape
        cov_matrix = np.zeros((m, n), dtype=float)
        cost_matrix = np.zeros((m, n), dtype=float) 
        for i in range(m):
            for j in range(n):
                mapping = mappings[i, j]
                if mapping is None:
                    cov_matrix[i, j] = 0  
                    cost_matrix[i, j] = BIG_NUMBER
                else:
                    cov_matrix[i, j] = mapping[3] / 100.
                    cost_matrix[i, j] = 1 - cov_matrix[i, j]

        # uses LAP to find the optimal mappings of chains
        atommaps = []
        (R, C), crrpds = multilap(cost_matrix, nodes, BIG_NUMBER)

        for row_ind, col_ind in crrpds:
            if len(row_ind) != m:
                continue
            atommap = None
            title = ''
            for r, c in zip(row_ind, col_ind):
                # if one of the chains failed to match then discard the entire atommap
                if mappings[r, c] is None: 
                    atommap = None
                    break
                atommap_ = mappings[r, c][0]
                title_ = '(' + atommap_.getTitle() + ')' 
                if atommap is None:
                    atommap = atommap_
                    title = title_
                else:
                    atommap += atommap_
                    title = title_ + ' + ' + title
            if atommap is not None:
                atommap.setTitle(title)
                atommaps.append(atommap)

        return atommaps, cov_matrix, (R, C)

    def _optimize(atommaps):
        # extract nonoverlaping mappings
        if len(atommaps):
            atommaps, rmsds = rankAtomMaps(atommaps, target)

            if rmsds is not None:
                debug['rmsd'] = list(rmsds)

                # if rmsd_cutoff is not None:
                #     for i in reversed(range(len(atommaps))):
                #         if rmsds[i] > rmsd_cutoff:
                #             atommaps.pop(i)
                #             rmsds.pop(i)

                # pre-store chain IDs of atommaps
                atommap_segchids = []
                for atommap in atommaps:
                    nodummies = atommap.select('not dummy')
                    chids = nodummies.getChids()
                    segids = nodummies.getSegnames()
                    segchids = []
                    for segid, chid in zip(segids, chids):
                        if (segid, chid) not in segchids:
                            segchids.append((segid, chid))
                    atommap_segchids.append(segchids)
                
                atommaps_ = []
                rmsd_standard = rmsds[0]
                while len(atommaps):
                    atommap = atommaps.pop(0)
                    rmsd = rmsds.pop(0)
                    segchids = atommap_segchids.pop(0)

                    if reject_rmsd is not None:
                        if rmsd > reject_rmsd:
                            break

                    if rmsd > rmsd_standard + drmsd:
                        break

                    atommaps_.append(atommap)

                    # remove atommaps that share chains with the popped atommap
                    for i in reversed(range(len(atommap_segchids))):
                        amsegchids = atommap_segchids[i]

                        for segchid in amsegchids:
                            if segchid in segchids:
                                atommaps.pop(i)
                                atommap_segchids.pop(i)
                                break

                atommaps = atommaps_
            else:
                debug['rmsd'] = None

        return atommaps

    # checkers
    if not isListLike(mappings):
        raise TypeError('mappings should be list-like')
    
    if len(mappings) == 0:
        raise ValueError('mappings cannot be empty')

    if isinstance(mappings, tuple):
        am, am_r, s, c = mappings
        return am
    
    mappings = np.atleast_2d(mappings)
    
    if mappings.ndim != 2:
        raise ValueError('mappings can only be either an 1-D or 2-D array')
    
    # build atommaps
    LOGGER.debug('Finding the atommaps based on their coverages...')
    nodes = []
    atommaps, cov_matrix, (R, C) = _build(mappings, nodes)
    if least_n_atommaps is None:
        n_mapped = 0
        for r, c in zip(R, C):
            if cov_matrix[r, c] > 0:
                n_mapped += 1
        least_n_atommaps = int(np.floor(float(n_mapped) / mappings.shape[0]))
        LOGGER.debug('Identified that there exists %d atommap(s) potentially.'%least_n_atommaps)

    debug['coverage'] = cov_matrix
    debug['solution'] = [1]

    # optimize atommaps based on superposition if target is given
    if target is not None and len(nodes):
        atommaps = _optimize(atommaps)
        i = 2

        if len(atommaps) < least_n_atommaps:
            LOGGER.debug('At least %d atommaps requested. '
                         'Finding alternative solutions.'%least_n_atommaps)

            LOGGER.progress('Solving for %d-best solution...', None, label='_atommap_lap')
            while len(atommaps) < least_n_atommaps:
                LOGGER.update(i, label='_atommap_lap')
                try:
                    more_atommaps, _, _ = _build(mappings, nodes)
                except SolutionDepletionException:
                    break
                more_atommaps = _optimize(more_atommaps)
                for j in reversed(range(len(more_atommaps))):
                    if more_atommaps[j] in atommaps:
                        more_atommaps.pop(j)
                if len(more_atommaps):
                    debug['solution'].append(i)
                atommaps.extend(more_atommaps)

                i += 1
            LOGGER.finish()
            LOGGER.report('%d atommaps were found in %%.2fs. %d requested'%(len(atommaps), least_n_atommaps), 
                          label='_atommap_lap')
    
    if len(atommaps) == 0:
        if np.count_nonzero(cov_matrix) == 0:
            LOGGER.warn('no atommaps were available. Consider adjusting accepting criteria')
        else:
            LOGGER.warn('no atommaps were found. Consider inceasing rmsd_reject or drmsd')
    return atommaps

def rankAtomMaps(atommaps, target):
    """Ranks :class:`.AtomMap` instances from *atommaps* based on its RMSD 
    with *target*.
    """
    
    rmsds = []
    coords0 = target.getCoords()
    if coords0 is None:
        rmsds = None
    else:
        for atommap in atommaps:
            weights = atommap.getFlags('mapped')
            coords = atommap.getCoords()
            rcoords, t = superpose(coords, coords0, weights)
            rmsd = calcRMSD(rcoords, coords0, weights)

            rmsds.append(rmsd)
        
        I = np.argsort(rmsds)

        atommaps = [atommaps[i] for i in I]
        rmsds = [rmsds[i] for i in I]
        
    return atommaps, rmsds


def alignChains(atoms, target, match_func=bestMatch, **kwargs):
    """Aligns chains of *atoms* to those of *target* using :func:`.mapOntoChains` 
    and :func:`.combineAtomMaps`. Please check out those two functions for details 
    about the parameters.
    """

    mappings = mapOntoChains(atoms, target, match_func, **kwargs)
    m, n = mappings.shape
    if m > n:
        LOGGER.warn('%s has fewer chains than %s'%(atoms.getTitle(), target.getTitle()))
        return []

    atommaps = combineAtomMaps(mappings, target, **kwargs)

    return atommaps


if __name__ == '__main__':

    from prody import *
    p = parsePDB('1p38')
    r = parsePDB('1r39')
    chtwo, chone = mapOntoChain(r, p['A'])[0][:2]
