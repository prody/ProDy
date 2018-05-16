# -*- coding: utf-8 -*-
"""This module defines functions for comparing and mapping polypeptide chains.
"""

from numbers import Integral

import numpy as np
from numpy import arange
PW2 = None

from prody.atomic import AtomMap as AM
from prody.atomic import Chain, AtomGroup, Selection
from prody.atomic import AAMAP
from prody.atomic import flags
from prody.measure import calcTransformation, printRMSD, calcDistance
from prody import LOGGER, SELECT, PY2K, PY3K
from prody.sequence import MSA
from prody.utilities import cmp

if PY2K:
    range = xrange

if PY3K:
    basestring = str

__all__ = ['matchChains', 'matchAlign', 'mapOntoChain', 'mapChainByChain', 
           'mapOntoChainByAlignment', 'getMatchScore', 'setMatchScore',
           'getMismatchScore', 'setMismatchScore', 'getGapPenalty', 
           'setGapPenalty', 'getGapExtPenalty', 'setGapExtPenalty',
           'getTrivialSeqId', 'setTrivialSeqId', 'getTrivialCoverage', 
           'setTrivialCoverage', 'getAlignmentMethod', 'setAlignmentMethod']

TRIVIAL_SEQID = 90.
TRIVIAL_COVERAGE = 50.
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
            from . import pairwise2
        except ImportError:
            try:
                from Bio import pairwise2
            except ImportError:
                raise ImportError('pairwise2 module could not be imported. '
                                  'Reinstall ProDy or install Biopython '
                                  'to solve the problem.')
        PW2 = pairwise2
    return PW2


def getTrivialSeqId():
    """Returns sequence identity used in the trivial mapping."""

    return TRIVIAL_SEQID


def setTrivialSeqId(seqid):
    """Set sequence identity used in the trivial mapping."""

    if isinstance(seqid, (float, int)) and seqid >= 0:
        global TRIVIAL_SEQID
        TRIVIAL_SEQID = seqid
    else:
        raise TypeError('seqid must be a positive number or zero')

def getTrivialCoverage():
    """Returns sequence coverage used in the trivial mapping."""

    return TRIVIAL_COVERAGE


def setTrivialCoverage(coverage):
    """Set sequence coverage used in the trivial mapping."""

    if isinstance(coverage, (float, int)) and coverage >= 0:
        global TRIVIAL_COVERAGE
        TRIVIAL_COVERAGE = coverage
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

    __slots__ = ['_res', '_name', '_num', '_inc']

    def __init__(self, number, name, insertioncode='', residue=None):
        self._num = number
        self._name = name
        self._inc = insertioncode
        self._res = residue

    def __repr__(self):
        return '<SimpleResidue: {0}{1}>'.format(self._name, self._num)

    def getResidue(self):
        return self._res

    def getResnum(self):
        return self._num

    def getIcode(self):
        return self._inc

    def getResname(self):
        return self._name

    def getCoords(self):
        if self._res:
            return self._res._getCoords()
        return None

class SimpleChain(object):

    """An internal class used to compare two polypeptide chains.


    SimpleChain instances can be indexed using residue numbers. If a residue
    with given number is not found in the chain, **None** is returned."""

    __slots__ = ['_list', '_seq', '_title', '_dict', '_gaps', '_coords', '_chain']

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
        self._coords = None
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

    def getCoords(self):
        return self._coords

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
        self._chain = chain
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
        self._coords = chain._getCoords()
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

    :keyword subset: one of the following well-defined subsets of atoms:
        ``"calpha"`` (or ``"ca"``), ``"backbone"`` (or ``"bb"``),
        ``"heavy"`` (or ``"noh"``), or ``"all"``, default is ``"calpha"``
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
        pairwise2 = importBioPairwise2()
        if pairwise2:
            LOGGER.debug('Trying to match chains based on {0} sequence '
                         'alignment:'.format(ALIGNMENT_METHOD))
            for simpch1, simpch2 in unmatched:
                LOGGER.debug('  Comparing {0} (len={1}) and {2} '
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
            elif subset is None or subset is 'all':
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
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        matches.sort(compare, reverse=True)
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
    for i in range(len(this)):
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

    :keyword subset: one of the following well-defined subsets of atoms:
        ``"calpha"`` (or ``"ca"``), ``"backbone"`` (or ``"bb"``),
        ``"heavy"`` (or ``"noh"``), or ``"all"``, default is ``"calpha"``
    :type subset: string

    :keyword seqid: percent sequence identity, default is 90
    :type seqid: float

    :keyword overlap: percent overlap, default is 90
    :type overlap: float

    :keyword mapping: if ``"ce"`` or ``"cealign"``, then the CE algorithm [IS98]_ will be 
        performed. It can also be a list of prealigned sequences, a :class:`.MSA` instance,
        or a dict of indices such as that derived from a :class:`.DaliRecord`.
        If set to anything other than the options listed above, including the default value 
        (**None**), a simple mapping will be first attempted and if that failed 
        then sequence alignment with a function from :mod:`~Bio.pairwise2` will be used 
        unless *pwalign* is set to **False**, in which case the mapping will fail.
    :type mapping: list, str

    :keyword pwalign: if **True**, then pairwise sequence alignment will 
        be performed. If **False** then a simple mapping will be performed 
        based on residue numbers (as well as insertion codes). This will be 
        overridden by the *mapping* keyword's value. 
    :type pwalign: bool

    This function tries to map *atoms* to *chain* based on residue
    numbers and types. Each individual chain in *atoms* is compared to
    target *chain*.
    
    .. [IS98] Shindyalov IN, Bourne PE. Protein structure alignment by 
       incremental combinatorial extension (CE) of the optimal path. 
       *Protein engineering* **1998** 11(9):739-47.
    """

    if not isinstance(atoms, (AtomGroup, Chain, Selection)):
        raise TypeError('atoms must be an AtomGroup, a Chain, or a '
                        'Selection instance')
    if not isinstance(chain, Chain):
        raise TypeError('chain must be Chain instance')

    subset = str(kwargs.get('subset', 'calpha')).lower()
    if subset not in _SUBSETS:
        raise ValueError('{0} is not a valid subset argument'
                         .format(str(subset)))
    seqid = kwargs.get('seqid', 90.)
    coverage = kwargs.get('overlap')
    if coverage is None:
        coverage = kwargs.get('coverage', 70.)
    pwalign = kwargs.get('pwalign', None)
    pwalign = kwargs.get('mapping', pwalign)
    alignment = None
    if pwalign is not None:
        if isinstance(pwalign, basestring):
            pwalign = str(pwalign).strip().lower()
        elif not isinstance(pwalign, bool):
            alignment = pwalign
            pwalign = True

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
        map_ag = mobile.getAtomGroup()
    else:
        if isinstance(mobile, AtomGroup):
            map_ag = mobile
        else:
            map_ag = mobile.getAtomGroup()
        chains = list(mobile.getHierView().iterChains())
        LOGGER.debug('Evaluating {0}: {1} chains are identified'
                     .format(str(atoms), len(chains))) 

    mappings = []
    unmapped = []
    unmapped_chids = []
    target_ag = target_chain.getAtomGroup()
    simple_target = SimpleChain(target_chain, False)
    LOGGER.debug('Trying to map atoms based on residue numbers and '
                'identities:')
    for chain in chains:
        simple_chain = SimpleChain(chain, False)
        if len(simple_chain) == 0:
            LOGGER.debug('  Skipping {0}, which does not contain any amino '
                        'acid residues.'.format(simple_chain))
            continue
        LOGGER.debug('  Comparing {0} (len={1}) with {2}:'
                    .format(simple_chain.getTitle(), len(simple_chain),
                            simple_target.getTitle()))

        # trivial mapping serves as a first simple trial of alignment the two 
        # sequences based on residue number, therefore the sequence identity 
        # (TRIVIAL_SEQID) criterion is strict.
        _seqid = _cover = -1
        target_list, chain_list, n_match, n_mapped = getTrivialMapping(
            simple_target, simple_chain)
        if n_mapped > 0:
            _seqid = n_match * 100 / n_mapped
            _cover = n_mapped * 100 / max(len(simple_target), len(simple_chain))

        trivial_seqid = TRIVIAL_SEQID if pwalign else seqid
        trivial_cover = TRIVIAL_COVERAGE if pwalign else coverage
        if _seqid >= trivial_seqid and _cover >= trivial_cover:
            LOGGER.debug('\tMapped: {0} residues match with {1:.0f}% '
                    'sequence identity and {2:.0f}% overlap.'
                    .format(n_mapped, _seqid, _cover))
            mappings.append((target_list, chain_list, _seqid, _cover))
        else:
            if not pwalign:
                LOGGER.debug('\tFailed to match chains based on residue numbers '
                        '(seqid={0:.0f}%, overlap={1:.0f}%).'
                        .format(_seqid, _cover))
            unmapped.append(simple_chain)
            unmapped_chids.append(chain.getChid())

    if not mappings and pwalign is None:
        pwalign = True

    if pwalign and unmapped:
        if alignment is None:
            if pwalign in ['ce', 'cealign']:
                aln_type = 'structure alignment'
                method = 'CE'
            else:
                aln_type = 'sequence alignment'
                method = ALIGNMENT_METHOD
        else:
            aln_type = 'alignment'
            method = 'predefined'
        LOGGER.debug('Trying to map atoms based on {0} {1}:'
                     .format(method, aln_type))

        for chid, simple_chain in zip(unmapped_chids, unmapped):
            LOGGER.debug('  Comparing {0} (len={1}) with {2}:'
                        .format(simple_chain.getTitle(), len(simple_chain),
                                simple_target.getTitle()))
            if method == 'CE':
                result = getCEAlignMapping(simple_target, simple_chain)
            else:
                if isinstance(alignment, dict):
                    result = getDictMapping(simple_target, simple_chain, map_dict=alignment)
                else:
                    result = getAlignedMapping(simple_target, simple_chain, alignment=alignment)

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
                    LOGGER.debug('\tMapped: {0} residues match with {1:.0f}%'
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
        #n_atoms = len(indices_target)

        ch_tar = next((r for r in residues_target if r is not None)).getChain()
        ch_chn = next((r for r in residues_chain if r is not None)).getChain()
        title_tar = 'Chain {0} from {1}'.format(ch_tar.getChid(), ch_tar.getAtomGroup().getTitle())
        title_chn = 'Chain {0} from {1}'.format(ch_chn.getChid(), ch_chn.getAtomGroup().getTitle())

        # note that chain here is from atoms
        atommap = AM(map_ag, indices_chain, chain.getACSIndex(),
                     mapping=indices_mapping, dummies=indices_dummies,
                     title=title_chn + ' -> ' + title_tar )
        selection = AM(target_ag, indices_target, target_chain.getACSIndex(),
                       title=title_tar + ' -> ' + title_chn, intarrays=True)

        mappings[mi] = (atommap, selection, _seqid, _cover)
    if len(mappings) > 1:
        def compare(m1, m2):
            return cmp(m1[2], m2[2])
        mappings.sort(compare, reverse=True)
    return mappings

def mapChainByChain(atoms, ref, **kwargs):
    """This function is similar to :func:`.mapOntoChain` but correspondence 
    of chains is found by their chain identifiers. """
    hv = atoms.getHierView()
    for chain in ref.getHierView().iterChains():
        for target_chain in hv.iterChains():
            if target_chain.getChid() == chain.getChid():
                mappings = mapOntoChainByAlignment(target_chain, chain, **kwargs)
                return mappings
    return []

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
                n = tar_indices.index(i)
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
        pairwise2 = importBioPairwise2()
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

    def _findAlignment(sequence, alignment):
        for seq in alignment:
            strseq = str(seq).upper()
            for gap in GAPCHARS:
                strseq = strseq.replace(gap, '')

            if sequence.upper() == strseq:
                return seq
        return None

    this = _findAlignment(target.getSequence(), alignment)
    if this is None:
        LOGGER.warn('alignment does not contain the target ({0}) sequence'
                    .format(this.getTitle()))
        return None

    that = _findAlignment(chain.getSequence(), alignment)
    if that is None:
        LOGGER.warn('alignment does not contain the chain ({0}) sequence'
                    .format(that.getTitle()))
        return None

    amatch = []
    bmatch = []
    aiter = target.__iter__()
    biter = chain.__iter__()
    n_match = 0
    n_mapped = 0
    for i in range(len(this)):
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

def getCEAlignMapping(target, chain):
    from .ccealign import ccealign

    tar_coords = target.getCoords().tolist()
    mob_coords = chain.getCoords().tolist()
    
    def add_tail_dummies(coords, window=8):
        natoms = len(coords)
        if natoms < window:
            return None
        rest = natoms % window

        tail_indices = []
        for i in range(rest):
            tail_indices.append(len(coords))
            coords.append(coords[-(i+1)])
        
        return tail_indices

    window = 8
    tar_dummies = add_tail_dummies(tar_coords, window)
    if tar_dummies is None:
        LOGGER.warn('target ({1}) is too small to be aligned '
                    'by CE algorithm (at least {0} residues)'
                    .format(window, repr(target)))
        return None

    mob_dummies = add_tail_dummies(mob_coords, window)
    if mob_dummies is None:
        LOGGER.warn('chain ({1}) is too small to be aligned '
                    'by CE algorithm (at least {0} residues)'
                    .format(window, repr(chain)))
        return None

    aln_info = ccealign((tar_coords, mob_coords))

    paths, bestIdx, nres, rmsd = aln_info[:4]
    path = paths[bestIdx]

    tar_indices = []
    chn_indices = []
    for i, j in path:
        if i not in tar_dummies:
            if j not in mob_dummies:
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
            b = chain_res_list[chn_indices[n]]
            bres = b.getResidue()
            bmatch.append(bres)
            if a.getResname() == b.getResname():
                n_match += 1
            n_mapped += 1
        else:
            bmatch.append(None)

    return amatch, bmatch, n_match, n_mapped

if __name__ == '__main__':

    from prody import *
    p = parsePDB('1p38')
    r = parsePDB('1r39')
    chtwo, chone = mapOntoChain(r, p['A'])[0][:2]
