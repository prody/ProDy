# -*- coding: utf-8 -*-

"""This module is called CaviTracer and defines functions for calculating 
channels, tunnels, pores, and surface cavities within protein structure.
"""

__author__ = 'Karolina Mikulska-Ruminska', 'Jan Brezovsky', 'Eryk Trzcinski'
__credits__ = ['Karolina Mikulska-Ruminska', 'Jan Brezovsky', 'Eryk Trzcinski']
__email__ = ['karolamik@fizyka.umk.pl']

import os
import logging
from collections import namedtuple
from contextlib import contextmanager
from collections import defaultdict
import numpy as np
from prody import LOGGER, PY3K
from prody.atomic import Atomic
from prody.utilities import getCoords, isListLike
from prody.proteins import writePDB, parsePDB, parsePQR
from prody.ensemble import Ensemble
from prody.trajectory import Trajectory
from prody.measure import calcCenter

if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

__all__ = ['getVmdModel', 'calcChannels', 'calcChannelsMultipleFrames', 
           'getChannelParameters', 'getChannelAtoms', 'showChannels', 
           'showCavities', 'showSurfaceCavities', 'selectChannelBySelection', 
           'getChannelResidueNames',
           'calcChannelSurfaceOverlaps', 'calcSurfaceCavities', 
           'calcSurfaceCavitiesMultipleFrames', 'getSurfaceCavityParameters',
           'getSurfaceCavityResidueNames', 'selectSurfaceCavityBySelection',
           'calcSurfaceCavityOverlaps', 'getLiningResindexMask',
           'getLiningResidues','getSiteScores','getSiteScoreStats', 'getLiningResIndices',
           'getSurfaceCavityResidueNamesMultipleFrames',
           'getSurfaceCavityParametersMultipleFrames', 
           'getChannelParametersMultipleFrames', '_reportAtomsInputComposition',
           'getObjectLiningInfo', 'getObjectLiningInfoMultipleFrames',
           'getSurfaceCavityLiningInfo','getSurfaceCavityLiningInfoMultipleFrames',
           'getConnectedCavityChannelLiningInfoMultipleFrames',
           'getConnectedCavityChannelLiningInfo',
           'getChannelResidueNamesMultipleFrames', 'calcPoresFromChannels',
           'showPores', 'getPoreParameters', 'getPoreResidueNames',
           'calcPoresFromChannelsMultipleFrames', 'getPoreParametersMultipleFrames',
           'getPoreResidueNamesMultipleFrames', 'scanChannelParameters',
           'getLinkParameters', 'getLinkResidueNames',
           'getLinkParametersMultipleFrames', 'getLinkResidueNamesMultipleFrames',
           'scanSurfaceCavityParameters',
           'calcFrequentObjectResidues', 'showFrequentObjectResidues',
           'connectChannelsToSurfaceCavities',
           'calcFrequentObjectResidues', 'showFrequentObjectResidues',
           'writeChannelsCIF',
           'buildLiningSequenceDistanceMatrix', 'buildCenterlineDistanceMatrix',
           'buildVoxelOverlapDistanceMatrix','prepareChannelRecords',
           'compareClusterings', 'computeBootstrapStability', 'coassociationMatrix',
           'clusterConsensusStrength',
            'renameChannelFilesWithClusterLabels', 'clusterChannels',
           'writeChannelsCIF', 'writePyMolCaviTracerScript',
           'mergeFramesPQR',
           'writeChimeraXMultiModelScript', 'writePyMolMultiModelScript',
           'writeVmdMultiModelScript',
           'ClusterCalculator','ClusteringRecords','Cluster',
           'writeClusterLabels']

# Van der Waals radii in Angstrom, by element symbol (upper case). The radii the
# tessellation is built on, and the ones the lining report measures a Voronoi
# vertex against, so both read them from here.
# Sources:
# BD  - Bondi family - J. Phys. Chem. 1964, 68, 441-451.,
#                      J. Phys. Chem. 1966, 70, 3006.,
#                      J. Phys. Chem. A 2009, 113, 5806-5812.
# AZ  - Alvarez - Dalton Trans. 2013, 42, 8617-8636.
# C&T - Charry and Tkatchenko - J. Chem. Theory Comput. 2024, 20, 7469-7478.

VDW_RADII = {
    # BD + AZ agree within 0.1 A, or the AZ value is uncertain: keeping BD.
    # Note that the values for Si, Br and Te are corrections from the later
    # papers of the BD family.
    'H': 1.20, 'HE': 1.40, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'F': 1.47, 'NE': 1.54,'NA': 2.27, 'SI': 2.22, 'P': 1.80, 'S': 1.80, 
    'CL': 1.75, 'AR': 1.88,'K': 2.75, 'GE': 2.11, 'AS': 1.85, 'SE': 1.90, 
    'BR': 1.83, 'KR': 2.02,'RB': 3.03, 'SB': 2.06, 'TE': 2.00, 'I': 1.98, 
    'XE': 2.16, 'CS': 3.43,

    # AZ values here, given the known issues of BD: some of its values are
    # effectively metallic radii, and a bare-atom probe artifact. In ambiguous
    # cases AZ shows the best agreement with the C&T radii.
    'LI': 2.12, 'BE': 1.98, 'MG': 2.51, 'AL': 2.25, 'CA': 2.62, 'SC': 2.58, 
    'TI': 2.46, 'V': 2.42, 'CR': 2.45, 'MN': 2.45, 'FE': 2.44, 'CO': 2.40,
    'NI': 2.40, 'CU': 2.38, 'ZN': 2.39, 'GA': 2.32, 'SR': 2.84, 'Y': 2.75, 
    'ZR': 2.52, 'NB': 2.56, 'MO': 2.45, 'TC': 2.44, 'RU': 2.46, 'RH': 2.44,
    'PD': 2.15, 'AG': 2.53, 'CD': 2.49, 'IN': 2.43, 'SN': 2.42, 'BA': 3.03, 
    'LA': 2.98, 'CE': 2.88, 'PR': 2.92, 'ND': 2.95, 'SM': 2.90, 'EU': 2.87, 
    'GD': 2.83, 'TB': 2.79, 'DY': 2.87, 'HO': 2.81, 'ER': 2.83, 'TM': 2.79,
    'YB': 2.80, 'LU': 2.74, 'HF': 2.63, 'TA': 2.53, 'W': 2.57, 'RE': 2.49, 
    'OS': 2.48, 'IR': 2.41, 'PT': 2.29, 'AU': 2.32, 'HG': 2.45, 'TL': 2.47, 
    'PB': 2.60, 'BI': 2.54, 'AC': 2.80, 'TH': 2.93, 'PA': 2.88, 'U': 2.71, 
    'NP': 2.82, 'PU': 2.81, 'AM': 2.83, 'CM': 3.05, 'BK': 3.40, 'CF': 3.05, 
    'ES': 2.70,

    # BD only: no contesting reference for these, but the values are likely
    # underestimated.
    'PO': 1.97, 'AT': 2.02, 'RN': 2.20, 'FR': 3.48, 'RA': 2.83,

    # Not an element: the radius the lining queries fall back on where the table
    # covers no entry, named so that the fallback is visible here rather than
    # buried as a literal at the point of use. 2.0.
    'UNKNOWN': 2.00
}

# Per-residue property scales, for the physicochemical categories of the sb-ncbr
# tunnels mmCIF schema (see writeChannelsCIF). The schema names no scale, but the
# MOLE method page its items were written from does, and it pins each one by its
# extremes; those extremes are quoted below and agree with the tables here, which
# is what makes the choice checkable rather than asserted.
# Sources:
# KD  - Kyte and Doolittle - J. Mol. Biol. 1982, 157, 105-132.
# CID - Cid, Bunster, Canales and Gazitua - Protein Eng. 1992, 5, 373-375.
#       (AAindex CIDH920105; note the AAindex record is laid out in A/L R/K ...
#        column pairs, so a naive row-wise read of it transposes the scale.)
# ZIM - Zimmerman, Eliezer and Simha - J. Theor. Biol. 1968, 21, 170-201.

# Hydropathy. The schema pins ARG -4.5 as most hydrophilic and ILE 4.5 as most
# hydrophobic; both hold here.
_KYTE_DOOLITTLE = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5,
    'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4, 'HIS': -3.2, 'ILE': 4.5,
    'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6,
    'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}

# Hydrophobicity, as an average of normalized scales. The schema pins GLU -1.140
# as most hydrophilic and ILE 1.810 as most hydrophobic; both hold here, and they
# are what identifies this scale as Cid's rather than one of the several other
# normalized consensus scales the description would otherwise fit.
_CID_HYDROPHOBICITY = {
    'ALA': 0.02, 'ARG': -0.42, 'ASN': -0.77, 'ASP': -1.04, 'CYS': 0.77,
    'GLN': -1.10, 'GLU': -1.14, 'GLY': -0.80, 'HIS': 0.26, 'ILE': 1.81,
    'LEU': 1.14, 'LYS': -0.41, 'MET': 1.00, 'PHE': 1.35, 'PRO': -0.09,
    'SER': -0.97, 'THR': -0.77, 'TRP': 1.71, 'TYR': 1.11, 'VAL': 1.13
}

# Polarity. The schema pins ALA and GLY at 0.00, SER at 1.67, GLU at 49.90 and
# ARG at 52.00; all four hold here.
_ZIMMERMAN_POLARITY = {
    'ALA': 0.00, 'ARG': 52.00, 'ASN': 3.38, 'ASP': 49.70, 'CYS': 1.48,
    'GLN': 3.53, 'GLU': 49.90, 'GLY': 0.00, 'HIS': 51.60, 'ILE': 0.13,
    'LEU': 0.13, 'LYS': 49.50, 'MET': 1.43, 'PHE': 0.35, 'PRO': 1.58,
    'SER': 1.67, 'THR': 1.66, 'TRP': 2.10, 'TYR': 1.61, 'VAL': 0.13
}

# Relative mutability, from empirical substitution matrices; the scale's reference
# is ALA at 100. Beware that Dayhoff's relative mutability is a different table
# under the same name - it is the one ProtScale serves - and differs by up to 40
# points, so a lookup by name alone lands on the wrong one. MOLE's own method page
# quotes thirteen of these values, and all thirteen agree with the table here.
# JTT - Jones, Taylor and Thornton - Bioinformatics 1992, 8, 275-282.
_JTT_MUTABILITY = {
    'ALA': 100, 'ARG': 83, 'ASN': 104, 'ASP': 86, 'CYS': 44,
    'GLN': 84, 'GLU': 77, 'GLY': 50, 'HIS': 91, 'ILE': 103,
    'LEU': 54, 'LYS': 72, 'MET': 93, 'PHE': 51, 'PRO': 58,
    'SER': 117, 'THR': 107, 'TRP': 25, 'TYR': 50, 'VAL': 98
}

# Lipophilicity and solubility of the Cbeta side-chain fragments, predicted rather
# than measured: MOLE obtained them from chemicalize.org and publishes them only as
# prose in its method page, which is where these were read from. Glycine has no
# side-chain fragment, so all three give it 0 rather than leaving it out.
#
# logP is the octanol/water partition coefficient of the fragment and logD the
# distribution coefficient at pH 7.4, which differ only where the fragment ionises
# - the two tables are equal for the fifteen neutral residues and lower in logD for
# ASP, GLU, ARG, LYS and HIS. That relationship is worth keeping in mind if these
# are ever re-typed, since it is the one internal check on them.
_MOLE_LOGP = {
    'ALA': 1.08, 'ARG': -0.08, 'ASN': -1.03, 'ASP': -0.22, 'CYS': 0.84,
    'GLN': -0.33, 'GLU': 0.48, 'GLY': 0.00, 'HIS': -0.01, 'ILE': 2.24,
    'LEU': 2.08, 'LYS': 0.70, 'MET': 1.48, 'PHE': 2.49, 'PRO': 1.80,
    'SER': -0.52, 'THR': -0.16, 'TRP': 2.59, 'TYR': 2.18, 'VAL': 1.80
}

_MOLE_LOGD = {
    'ALA': 1.08, 'ARG': -2.49, 'ASN': -1.03, 'ASP': -3.00, 'CYS': 0.84,
    'GLN': -0.33, 'GLU': -2.12, 'GLY': 0.00, 'HIS': -0.11, 'ILE': 2.24,
    'LEU': 2.08, 'LYS': -1.91, 'MET': 1.48, 'PHE': 2.49, 'PRO': 1.80,
    'SER': -0.52, 'THR': -0.16, 'TRP': 2.59, 'TYR': 2.18, 'VAL': 1.80
}

# Water solubility at pH 7.4, as a unit-stripped base-10 logarithm of mol/litre.
_MOLE_LOGS = {
    'ALA': 0.59, 'ARG': 1.63, 'ASN': 0.54, 'ASP': 2.63, 'CYS': 0.16,
    'GLN': 0.13, 'GLU': 2.23, 'GLY': 0.00, 'HIS': -0.20, 'ILE': -1.85,
    'LEU': -1.79, 'LYS': 1.46, 'MET': -0.72, 'PHE': -1.81, 'PRO': -1.30,
    'SER': 1.11, 'THR': 0.77, 'TRP': -2.48, 'TYR': -1.44, 'VAL': -1.30
}

# The residues every scale above covers - they share one set of twenty keys. What
# lies outside it is what a lining report has to be honest about: a nucleotide, a
# cofactor, a modified or differently protonated residue, an ion.
_SCALE_RESIDUES = frozenset(_KYTE_DOOLITTLE)

# Charge, as the schema itself defines it rather than by way of ProDy's acidic
# and basic flags: the item description spells the formula out, so taking it
# literally is what makes the number reproducible from the schema alone.
_CHARGED_RESIDUES = {'ARG': 1, 'LYS': 1, 'HIS': 1, 'ASP': -1, 'GLU': -1}

# Residues that can carry a charge at some pH, which is a wider set than the one
# charged at pH 7 - the schema's wording is "can go into an ionic state".
_IONIZABLE_RESIDUES = frozenset(['ASP', 'GLU', 'HIS', 'CYS', 'TYR', 'LYS',
                                 'ARG'])

_OVERLAP_OFFSET_CACHE = {}

# The tag a per-object file carries in its name, against the word the same
# object is written under inside the residue and parameter text files. The two
# differ for every object but the cavities, which would broke the reading side
# of selectChannelBySelection when the per-object files labels change e.g., from
# _channel0.pqr to _chl0.pqr
_OBJECT_TAGS = {'chl': 'channel', 'lnk': 'link', 'pore': 'pore',
                'cavity': 'cavity'}

@contextmanager
def _warningsDelivered():
    """Let WARNING records through for the duration of the block.

    Importing ProDy installs a logging filter that drops every WARNING record from
    the package logger (prody.dynamics.adaptive2, at module scope), so every
    ``LOGGER.warn`` in ProDy is silently discarded. Whether that filter should exist
    at all is a question for the package as a whole; until it is settled, this
    module at least delivers its own warnings.

    Offending filters are found by asking them, rather than by importing the class
    and matching on type: any filter that rejects a synthetic WARNING record is
    detached for the block and reinstated afterwards. That keeps this working if the
    filter is renamed or moved, and it leaves the filter in force everywhere else,
    so nothing outside this module changes behaviour."""
    logger = LOGGER._logger
    probe = logging.LogRecord(logger.name, logging.WARNING, __file__, 0,
                              '', (), None)
    muting = [f for f in list(logger.filters)
              if not (f.filter(probe) if hasattr(f, 'filter') else f(probe))]
    for f in muting:
        logger.removeFilter(f)
    try:
        yield
    finally:
        for f in muting:
            logger.addFilter(f)


def _warn(message):
    """``LOGGER.warn``, but actually emitted. See :func:`_warningsDelivered`."""
    with _warningsDelivered():
        LOGGER.warn(message)


def checkAndImport(package_name):
    """Check for package and import it if possible and return **True**.
    Otherwise, return **False
        
    :arg package_name: name of package
    :type package_name: str

    :arg import_command: optional command to import submodules or with an alias
        default **None** means use "import {0}".format(package_name)
    :type import_command: None, str """
    
    if not isinstance(package_name, str):
        raise TypeError('package_name should be a string')

    if PY3K:
        import importlib.util
        if importlib.util.find_spec(package_name) is None:
            _warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    else:
        try:
            __import__(package_name)
        except ImportError:
            _warn("Package " + str(package_name) + " is not installed. "
            "Please install it to use this function.")
            return False
    
    return True


def _requireCoords(atoms):
    """Raise :exc:`TypeError` unless *atoms* can supply coordinates.

    :func:`~prody.utilities.getCoords` is the whole check: it reads the
    coordinates through ``_getCoords``/``getCoords``, accepts a plain array, and
    raises :exc:`TypeError` for anything else. Each entry point used to wrap that
    call in a ``try``/``except AttributeError`` handler that re-implemented its
    body -- and could not run, since ``getCoords`` converts the
    :exc:`AttributeError` into the :exc:`TypeError` it is meant to raise."""

    getCoords(atoms)


def _numberedPath(filename, tag, index, suffix='', stem=None):
    """``channels.pqr`` with ``('chl', 0)`` becomes ``channels_chl0.pqr``.

    *suffix* is appended after the number, for the objects that have a far end
    to name as well as a near one: a link out of start point 13 into start point
    9 comes out as ``channels_sp13_lnk0_sp9.pqr``.

    *stem* overrides the part before the tag. ``None`` takes it from the file
    name, as above; ``''`` drops it and its separator, giving ``chl0.pqr``,
    which is what a run told only a directory writes - there the stem carries no
    information, being the same placeholder for every run.

    Only the file's own name is numbered. The per-object files used to be named
    by ``filename.replace('.pqr', ...)``, which rewrites every occurrence of the
    extension anywhere in the path: an output directory called ``run.pqr/`` was
    renamed along with the file and the write then failed, a name carrying no
    extension was not numbered at all -- so each object overwrote the previous
    one and the combined file with it -- and a name holding both extensions was
    numbered twice."""

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    path = Path(filename)
    if stem is None:
        stem = path.stem
    return path.with_name("{0}{1}{2}{3}{4}".format(
        stem + '_' if stem else '', tag, index, suffix, path.suffix))


def _frameBounds(n_frames, start_frame=0, stop_frame=-1):
    """Half-open ``[first, last)`` over the 0-based frames a run covers.

    The one definition of what ``start_frame`` and ``stop_frame`` mean, because
    having several was how they came to disagree. ``stop_frame`` is **inclusive** -
    the last frame analysed, as every docstring in this module says - and ``-1``
    means "through the last one". Frames are numbered from 0, matching
    :meth:`~.AtomGroup.setACSIndex` and :meth:`~.AtomGroup.getCoordsets`; note that
    a multi-model PDB numbers its ``MODEL`` records from 1, so ``stop_frame=0`` is
    ``MODEL 1``.

    *n_frames* is how many there are, or ``None`` where that is not known yet - a
    trajectory read frame by frame - in which case ``last`` comes back ``None``,
    which slices to the end just as an omitted bound does.

    The reason this exists rather than a slice at each call site: ``-1`` is a
    sentinel here, but it is also a perfectly good Python index, and
    ``coordsets[start:-1]`` therefore silently drops the last frame instead of
    keeping all of them. A slice is also exclusive where this API is inclusive, so
    the same expression is off by one again for every other value. Both mistakes
    read as correct code, which is why they lasted.

    An out-of-range *stop_frame* is clamped rather than raising: callers wrote
    ``stop_frame=atoms.numCoordsets()`` to work around the exclusive behaviour, and
    that should keep meaning "all of them" rather than becoming an IndexError."""

    first = max(0, int(start_frame))

    if stop_frame is None or int(stop_frame) < 0:
        last = n_frames
    else:
        last = int(stop_frame) + 1
        if n_frames is not None:
            last = min(last, n_frames)

    if last is not None and last < first:
        last = first

    return first, last


def _frameOutputPath(output_path, index, name, suffix='.pqr'):
    """Where one frame of a multi-frame run writes.

    A directory takes the files inside it, named after what they hold and the
    frame they came from: ``2kid`` gives ``2kid/cavities3.pqr``, and the
    per-object files under it carry that stem, so the frame stays in every name
    and frames cannot overwrite one another. Anything else is used as a prefix,
    as it always has been: ``2kid`` with no such directory gives ``2kid3.pqr``
    beside it.

    *suffix* is what the frame's file is named with, so that an mmCIF run numbers
    its frames exactly as a PQR one does. The schema has no frame of its own - it
    describes one structure - so a frame per file is what keeps each written file
    something the schema can actually describe."""

    import os

    path = str(output_path)
    if os.path.isdir(path):
        return os.path.join(path, "{0}{1}{2}".format(name, index, suffix))
    return path + "{0}{1}".format(index, suffix)


def _splitObjectFileName(filename):
    """``run1_sp0_chl3.pqr`` becomes ``('run1', 'channel', 3)``.

    The inverse of :func:`_numberedPath`: it reads back the prefix, the kind of
    object and its number from a per-object file name. Start points are not part
    of the identity - an object is numbered over the whole run, not within its
    start point - so both the ``_sp<N>`` before the tag and the ``_sp<M>`` a
    link carries after its number are dropped. The prefix is empty for the names
    a run told only a directory writes (``sp0_chl3.pqr``). Returns ``None`` when
    the name is not one of ours."""

    import re

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    match = re.match(r'^(?:(?P<head>.*)_)?(?P<tag>{0})(?P<index>\d+)'
                     r'(?:_sp\d+)?$'.format('|'.join(_OBJECT_TAGS)),
                     Path(filename).stem)
    if match is None:
        return None
    # A start point sits between the prefix and the tag, so it lands in the head
    # and is taken off here rather than in the pattern, where an optional group
    # before a greedy one would just as happily read "sp0" as the prefix itself.
    head = re.sub(r'(?:^|_)sp\d+$', '', match.group('head') or '')
    return head, _OBJECT_TAGS[match.group('tag')], int(match.group('index'))


def _selectObjectRows(files, suffix, object_name, file_prefix=None):
    """The rows of ``<prefix><suffix>`` belonging to *files*, in their order.

    One text file holds every object of a run, one row each, keyed by the prefix
    the run was written under and the object's number. Each text file is read
    once however many objects are looked up in it."""

    import os

    rows = []
    cache = {}
    unnamed = []
    for path in files:
        name = os.path.basename(str(path))
        parsed = _splitObjectFileName(name)
        if parsed is None:
            # The combined file is one of these every time, so they are counted
            # and reported once rather than warned about one by one.
            unnamed.append(name)
            continue
        prefix, found_name, index = parsed
        if found_name != object_name:
            continue                   # a file of another kind, not an error
        if file_prefix is not None:
            prefix = file_prefix

        text_file = prefix + suffix
        if text_file not in cache:
            try:
                with open(text_file, 'r') as handle:
                    cache[text_file] = handle.readlines()
            except (IOError, OSError):
                cache[text_file] = None
                _warn("{0} not found, so nothing was extracted for {1}. It is "
                      "written by the get*Parameters and get*ResidueNames "
                      "functions under the name passed to them; pass "
                      "file_prefix if that name is not {2!r}.".format(
                          text_file, name, prefix))
        lines = cache[text_file]
        if lines is None:
            continue

        key = '{0}_{1}{2}:'.format(prefix, object_name, index)
        matched = [line for line in lines if line.startswith(key)]
        if not matched:
            _warn("no row starting {0!r} in {1}, so {2} contributed "
                  "nothing.".format(key, text_file, name))
        rows.extend(matched)

    if unnamed:
        LOGGER.info("{0} file(s) hold no single numbered object, so no {1} row "
                    "was looked up for them: {2}.".format(
                        len(unnamed), object_name, ', '.join(unnamed)))

    return rows


def _kdTree(coords):
    """A :class:`~scipy.spatial.cKDTree` over *coords*.

    Accepts either an :class:`.Atomic` object or an ``(n, 3)`` array, so the one
    helper serves the lining queries, the enclosure test and the route
    comparison. scipy stays a function-level import, as everywhere else in this
    module, so that importing ProDy does not require it."""

    from scipy.spatial import cKDTree

    if hasattr(coords, 'getCoords'):
        coords = coords.getCoords()
    return cKDTree(np.asarray(coords, dtype=float))


def _getOverlapSphereOffsets(radius, resolution):
    """Return integer voxel offsets inside a sphere of a given radius."""

    key = (round(float(radius), 3), float(resolution))

    if key in _OVERLAP_OFFSET_CACHE:
        return _OVERLAP_OFFSET_CACHE[key]

    n = int(np.ceil(radius / resolution))
    grid = np.arange(-n, n + 1, dtype=int)

    dx, dy, dz = np.meshgrid(grid, grid, grid, indexing='ij')
    offsets = np.vstack((dx.ravel(), dy.ravel(), dz.ravel())).T
    xyz = offsets.astype(float) * resolution
    mask = np.sum(xyz * xyz, axis=1) <= radius * radius
    offsets = offsets[mask]
    _OVERLAP_OFFSET_CACHE[key] = offsets

    return offsets


def _surfaceFromPqrWorker(args):
    """Create voxelized FIL surface for one PQR file."""

    pqr_file, resolution = args
    atoms = parsePQR(pqr_file)
    fil = atoms.select('resname FIL')

    if fil is None:
        return set()

    coords = fil.getCoords()
    radii = fil.getRadii()
    surface = set()

    for center, radius in zip(coords, radii):
        center_idx = np.rint(center / resolution).astype(int)
        offsets = _getOverlapSphereOffsets(float(radius), resolution)
        voxels = offsets + center_idx
        surface.update(map(tuple, voxels))

    return surface


def _calcChannelsMultipleFramesWorker(args):
    """Compute channels. Supporting function for muliprocessing in :func:`calcChannelsMultipleFrames`."""
    frame_nr, atoms, frame_coords, frame_output_path, separate, start_point, return_details, kwargs = args

    # LOGGER.info("Frame/model: {0}".format(frame_nr))
    atoms_copy = atoms.copy()
    atoms_copy.setCoords(frame_coords)
    try:
        if return_details:
            channels, surfaces, details =  calcChannels(atoms_copy,
                                                output_path=frame_output_path,
                                                separate=separate,
                                                    start_point=start_point,
                                                    return_details=return_details,
                                                    **kwargs)
        else:
            channels, surfaces =  calcChannels(atoms_copy,
                                    output_path=frame_output_path,
                                    separate=separate,
                                        start_point=start_point,
                                        return_details=return_details,
                                        **kwargs)
        if not channels:
            if return_details:
                return frame_nr, [], [], {}
            return frame_nr,[],[]
    except Exception as e:
        print("Frame {0} failed: {1}".format(frame_nr, e))
        if return_details:
            return frame_nr, [], [], {}
        return frame_nr, [], []
    if return_details:
        return frame_nr, channels, surfaces, details
    return frame_nr, channels, surfaces

def _calcSurfaceCavitiesMultipleFramesWorker(args):
    """Compute surface cavities. Supporting function for multiprocessing in :func:`calcSurfaceCavitiesMultipleFrames`."""
    frame_nr, atoms, frame_coords, frame_output_path, separate, kwargs = args

    LOGGER.info("Frame/model: {0}".format(frame_nr))
    atoms_copy = atoms.copy()
    atoms_copy.setCoords(frame_coords)

    return calcSurfaceCavities(atoms_copy, output_path=frame_output_path, separate=separate, **kwargs)
    
def _calcPoresFromChannelsWorker(args):
    """Reconstruct pores from channels. Supporting function for multiprocessing
    in :func:`calcPoresFromChannelsMultipleFrames`."""
    frame_nr, channels, details, output_path, separate, kwargs = args
    LOGGER.info("Frame/model: {0}".format(frame_nr))
    return calcPoresFromChannels(channels, details, output_path=output_path,
                                 separate=separate, **kwargs)

def _findContactRun(mask, min_contact_points):
    """Return the first index of the last valid contiguous True run in *mask*.

    The path is assumed to be oriented from the protein interior towards the
    surface cavity. The run nearest the surface end is preferred. """

    mask = np.asarray(mask, dtype=bool)

    if len(mask) == 0:
        return None

    end = len(mask) - 1

    while end >= 0:
        if not mask[end]:
            end -= 1
            continue

        start = end

        while start > 0 and mask[start - 1]:
            start -= 1

        if end - start + 1 >= min_contact_points:
            return start

        end = start - 1

    return None


def _selectLocalSurfaceCavity(cavity, cavity_surface, channel,
    tolerance=1.0, cavity_margin=2.0, num_samples=5):
    """Select the connected part of a surface cavity associated with a channel.

    Surface-cavity tetrahedra are first filtered according to their proximity
    to the channel surface. Only the connected component containing the actual
    channel-cavity contact region is retained. """

    cavity_tetrahedra = np.asarray(cavity.tetrahedra, dtype=np.intp)

    if len(cavity_tetrahedra) == 0:
        return cavity_tetrahedra
    
    if cavity_margin is None:
        return cavity_tetrahedra.copy()
    
    cavity_simplices = np.asarray(cavity_surface[3])
    cavity_vertices = np.asarray(cavity_surface[4])

    cavity_xyz = cavity_vertices[cavity_tetrahedra]

    # Sample the reconstructed channel including its radius.
    channel_centers, channel_radii = _sampleObjectSpheres(channel, num_samples)

    # Distance of every cavity Voronoi vertex from its nearest sampled channel center.
    channel_tree = _kdTree(channel_centers)
    distances, nearest = channel_tree.query(cavity_xyz)
    nearest_radii = channel_radii[nearest]

    # Seed tetrahedra are in direct contact with the channel centerline.
    contact_mask = distances <= tolerance

    if not np.any(contact_mask):
        return np.array([], dtype=np.intp)

    # Candidate cavity tetrahedra must remain close to the actual channel
    # surface, not merely to one connection point.
    candidate_mask = distances <= nearest_radii + cavity_margin

    # Contact tetrahedra must always remain candidates.
    candidate_mask |= contact_mask
    candidate_local = np.where(candidate_mask)[0]
    contact_local = np.where(contact_mask)[0]

    if len(candidate_local) == 0:
        return np.array([], dtype=np.intp)

    # Build topology only inside this cavity. Two Delaunay tetrahedra are
    # neighbours if they share one triangular face, i.e. three atom indices.
    face_to_tetrahedra = {}

    for local_index in candidate_local:
        tetra = int(cavity_tetrahedra[local_index])
        simplex = cavity_simplices[tetra]

        faces = (tuple(sorted((simplex[0], simplex[1], simplex[2]))),
            tuple(sorted((simplex[0], simplex[1], simplex[3]))),
            tuple(sorted((simplex[0], simplex[2], simplex[3]))),
            tuple(sorted((simplex[1], simplex[2], simplex[3]))))

        for face in faces:
            face_to_tetrahedra.setdefault(face, []).append(local_index)

    adjacency = {int(i): set() for i in candidate_local}

    for tetrahedra in face_to_tetrahedra.values():
        if len(tetrahedra) < 2:
            continue

        for i in tetrahedra:
            for j in tetrahedra:
                if i != j:
                    adjacency[int(i)].add(int(j))

    # Flood-fill from every genuine channel-cavity contact tetrahedron.
    seeds = [int(i) for i in contact_local if candidate_mask[i]]

    if not seeds:
        return np.array([], dtype=np.intp)

    visited = set(seeds)
    stack = list(seeds)

    while stack:
        current = stack.pop()

        for neighbor in adjacency.get(current, ()):
            if neighbor in visited:
                continue

            visited.add(neighbor)
            stack.append(neighbor)

    selected_local = np.array(sorted(visited), dtype=np.intp)

    return cavity_tetrahedra[selected_local]


def _saveConnectedCavityChannels(connected, cavity_surface, filename,
    separate=False, num_samples=5):
    """Save local surface-cavity regions together with connected channels."""

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    filename = Path(filename)

    if filename.is_dir():
        filename = filename / 'connected_cavities_channels.pqr'
        separate_stem = ''
    else:
        separate_stem = None

        if filename.suffix not in ('.pdb', '.pqr'):
            filename = filename.with_suffix('.pqr')

    cavity_vertices = cavity_surface[4]


    def records(result, atom_index=1):
        cavity_index = result['cavity_index']
        channel_index = result['channel_index']
        cavity = result['cavity']
        channel = result['trimmed_channel']
        cavity_tetrahedra = result['cavity_tetrahedra']
        connection_point = result['connection_point']

        lines = []
        lines.append("REMARK   connected cavity %d channel %d  "
            "connection=%.3f %.3f %.3f A\n" %
            (cavity_index, channel_index, connection_point[0],
             connection_point[1], connection_point[2]))

        cavity_xyz = cavity_vertices[cavity_tetrahedra]

        lines.append("REMARK   local cavity %d  tetrahedra=%d  "
            "original_tetrahedra=%d\n" %
            (cavity_index, len(cavity_tetrahedra), len(cavity.tetrahedra)))

        cavity_radius = ChannelCalculator.CAVITY_MARKER_RADIUS

        for x, y, z in cavity_xyz:
            lines.append("ATOM  %5d  H   FIL C%4d    "
                "%8.3f%8.3f%8.3f%6.2f%6.2f\n" %
                (atom_index, cavity_index + 1, x, y, z, 1.00, cavity_radius))
            atom_index += 1

        channel_lines, samples = calculator_records(channel_index, channel, atom_index, num_samples)
        lines.extend(channel_lines)
        atom_index += samples
        lines.append("\n")

        return lines, atom_index

    def calculator_records(channel_index, channel, atom_index, num_samples):
        centers, radii = _sampleObjectSpheres(channel, num_samples)

        lines = [ChannelCalculator._channelRemark(channel_index, channel, label='channel')]

        for i, (x, y, z, radius) in enumerate(
                zip(centers[:, 0], centers[:, 1],
                    centers[:, 2], radii),
                start=atom_index):

            lines.append("ATOM  %5d  H   FIL H%4d    "
                "%8.3f%8.3f%8.3f%6.2f%6.2f\n" %
                (i, channel_index + 1, x, y, z, 1.00, radius))

        for i in range(atom_index, atom_index + len(centers) - 1):
            lines.append("CONECT%5d%5d\n" % (i, i + 1))

        return lines, len(centers)

    with open(str(filename), 'w') as handle:
        atom_index = 1

        for result in connected:
            lines, atom_index = records(result, atom_index)
            handle.writelines(lines)

    LOGGER.info("Connected surface cavities and channels saved to {0}.".format(filename))

    if separate:
        for pair_index, result in enumerate(connected):

            pair_filename = _numberedPath(filename, 'cavchl', pair_index, stem=separate_stem)

            with open(str(pair_filename), 'w') as handle:
                lines, _ = records(result, 1)
                handle.writelines(lines)

        LOGGER.info("Saved {0} individual connected cavity-channel file(s).".format(len(connected)))


def _reportAtomsInputComposition(atoms, inner_radius=None, diagram=None):
    """Report the composition of atoms supplied for channel analysis.

    This function checks whether the input atomic structure contains only
    protein atoms or also includes water, non-water HETATM records, or other
    non-protein components. If non-protein atoms are present, a warning is
    issued indicating that all supplied atoms will be included in the channel
    calculation.

    It also reads how far the structure is protonated, and reports where that
    sits badly with the run about to be made: a probe smaller than water on a
    structure missing its hydrogens, hydrogens present but unused because the
    probe was left at the size an unprotonated structure needs, or a
    radius-blind diagram on a structure whose radii vary most.

    The function does not modify or filter the input structure. To analyze only
    the protein, the user should provide an appropriate ProDy selection, for
    example ``atoms.select('protein')``.

    :arg inner_radius: probe radius the run will use, if known. Judged against
        the hydrogen content; omit it to skip that pair of checks.
    :type inner_radius: float

    :arg diagram: Voronoi diagram the run will use, if known. Omit it to skip
        the check on radius-blind tessellation.
    :type diagram: str """

    if not isinstance(atoms, Atomic):
        raise TypeError(
            "atoms must be a ProDy Atomic object, such as an AtomGroup "
            "or Selection")

    protein = atoms.select('protein')
    nucleic = atoms.select('nucleic')
    water = atoms.select('water')
    hetero = atoms.select('hetero and not water')
    other = atoms.select('not protein and not nucleic and not hetero')
    # Nucleic acid is part of the biomolecule the channels run through, not something
    # that found its way into the selection, so it is named alongside the protein
    # rather than counted among the components worth warning about.
    foreign = atoms.select('not protein and not nucleic')

    if foreign is None:
        LOGGER.info("The atoms supplied to calcChannels contain {0} atoms only.".format(
            " and ".join(name for name, selection in (('protein', protein),
                                                      ('nucleic acid', nucleic))
                         if selection is not None)))
    else:
        components = []

        if water is not None:
            components.append(
                "water: {0} atoms in {1} residues".format(
                    water.numAtoms(),
                    len(np.unique(water.getResindices()))))

        if hetero is not None:
            components.append(
                "non-water hetero components: {0} atoms "
                "(resnames: {1})".format(
                    hetero.numAtoms(),
                    ", ".join(sorted(np.unique(hetero.getResnames())))))

        if other is not None:
            components.append(
                "other components: {0} atoms "
                "(resnames: {1})".format(
                    other.numAtoms(),
                    ", ".join(sorted(np.unique(other.getResnames())))))

        # The advice names only what the structure actually holds, so a protein with
        # a ligand is still pointed at 'protein' and only a complex is told about
        # the wider selection.
        present = [(name, keyword) for name, keyword, selection in
                   (('protein', 'protein', protein),
                    ('nucleic acid', 'nucleic', nucleic)) if selection is not None]
        _warn("The atoms supplied to calcChannels() contain components other than "
            "{1}: {0}. All supplied atoms except waters will be used for channel "
            "analysis. To analyze only the {1}, provide an appropriate selection, "
            "for example atoms.select('{2}').".format(
                "; ".join(components),
                " and ".join(name for name, _ in present),
                " or ".join(keyword for _, keyword in present)))

    # How far the structure is protonated, measured on the biomolecule alone: waters
    # are dropped before the diagram is built, and a solvated but otherwise bare
    # structure would look protonated through its water hydrogens. Each kind of chain
    # is scored against the hydrogens per heavy atom it should carry, since those
    # differ: the standard amino acids hold about as many hydrogens as heavy atoms,
    # so a complete protein sits near 1.0, while a nucleotide is much richer in heavy
    # atoms -- phosphate oxygens and ring nitrogens bear none -- so a complete nucleic
    # acid sits near 0.55. A mixed structure counts as protonated only when every kind
    # of chain in it is, because hydrogens missing anywhere open interstices wherever
    # the channel happens to run.
    chains = [(selection, name, expected) for selection, name, expected in
              ((protein, 'protein', 1.0), (nucleic, 'nucleic acid', 0.55))
              if selection is not None]
    if not chains:
        rest = atoms.select('not water')
        if rest is None:
            return
        chains = [(rest, 'structure', 1.0)]

    n_hydrogen = 0
    scarcest = None
    for selection, name, expected in chains:
        elements = np.char.upper(np.asarray(selection.getElements(), dtype=str))
        count = int(np.count_nonzero(elements == 'H'))
        heavy = int(np.count_nonzero(elements != 'H'))
        n_hydrogen += count
        # Fraction of the hydrogens a complete chain of this kind would carry, which
        # puts protein and nucleic acid on one scale.
        fraction = count / float(max(heavy, 1)) / expected
        if scarcest is None or fraction < scarcest[0]:
            scarcest = (fraction, name, count)

    # A complete file scores within a few percent of 1.0, one that kept only its polar
    # hydrogens about a fifth of that, and one straight from a refinement nothing, so
    # the threshold has a wide margin either side and is not delicate.
    fully_protonated = scarcest[0] >= 0.7

    # An experimental structure generally carries no hydrogens -- X-ray and cryo-EM
    # alike, since neither resolves them except at the very highest resolutions -- and
    # its carbons keep their full vdW radius, so the ~0.6 A the missing H occupied is
    # left as void, around every heavy atom at once, including buried contacts that
    # never come apart. That is usually harmless, and is often defended as standing in
    # for thermal motion: a probe of water size cannot enter those interstices anyway,
    # and protonated and unprotonated runs agree from about 1.2 A upwards. Below that
    # the probe is small enough to thread them and the interior percolates into a
    # sponge rather than merely widening. Those routes might be fictitious, not the real
    # ones made wider. So a sub-water probe needs real hydrogens -- all of them, since
    # it is the apolar C-H that fill those interstices, and a structure holding only
    # its polar hydrogens leaves them just as open as one holding none.
    if inner_radius is not None and inner_radius < 1.2 and not fully_protonated:
        _warn("inner_radius={0:.2f} is below 1.2 Å but the {1} {2}: the space "
              "left by the missing H is then wide enough for the probe to pass, and "
              "channels will be found through interstices that do not exist in the "
              "real structure (their number can rise several-fold). Either add "
              "hydrogens, or raise inner_radius to 1.2 Å or more, where protonated and "
              "unprotonated structures give similar channels.".format(
                  inner_radius, scarcest[1],
                  "carries no hydrogens" if not scarcest[2] else
                  "carries only {0:.0f}% of the hydrogens a complete one would".format(
                      100.0 * scarcest[0])))

    # The reverse mismatch. The 1.2 A floor above is a workaround for absent
    # hydrogens, not a property of the probe, so a structure that carries them is
    # being measured with a probe coarser than its own detail: narrow connections
    # are reported closed rather than measured. Only a note -- nothing is wrong with
    # the result, it is simply more conservative than the input requires.
    if inner_radius is not None and inner_radius >= 1.2 and fully_protonated:
        LOGGER.info("The structure carries its hydrogens ({0:.0f}% of what a complete "
            "{1} would hold), so inner_radius={2:.2f} is more conservative than it needs "
            "to be: the 1.2 Å floor exists only to keep a sub-water probe out of the "
            "space that missing hydrogens leave open, and here that space is filled. A "
            "smaller probe, down to about 0.9 Å, measures the narrow connections instead "
            "of reporting them closed.".format(
                100.0 * scarcest[0], scarcest[1], inner_radius))

    # 'simple' builds an *unweighted* Delaunay of the atom centres, i.e. it
    # ignores the differences between atomic radii. That approximation is worst
    # when the radius spread is largest -- which is exactly when hydrogens (small
    # vdW) are present -- so warn there and steer the user to a radius-aware mode.
    # With H absent the heavy-atom radii are much closer, so 'simple' is more
    # defensible and matches the heavy-atom-only input most tools accept (at the
    # cost of over-large empty space where the missing H would sit).
    if diagram == "simple" and n_hydrogen:
        _warn("diagram='simple' with hydrogens present: the unweighted "
            "Voronoi diagram ignores radius differences, which are largest when H "
            "are present, so its topology and clearances are significantly less "
            "accurate. Consider diagram='homogenized' (or 'weighted'), which "
            "account for per-atom radii.")


def getVmdModel(vmd_path, atoms, representation='NewCartoon'):
    """Generates a 3D model of molecular structures using VMD and returns 
    it as an Open3D TriangleMesh.

    This function creates a temporary PDB file from the provided atomic data
    and uses VMD (Visual Molecular Dynamics) to render this data into an STL
    file, which is then loaded into Open3D as a TriangleMesh. The function
    handles the creation and cleanup of temporary files and manages the
    subprocess call to VMD.
    
    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d

    If problem with `ipykernel.comm.Comm` class appeared while using getVmdModel
    please update dash: python -m pip install -U dash


    :arg vmd_path: Path to the VMD executable. This is required to run VMD and 
        execute the TCL script.
    :type vmd_path: str

    :arg atoms: Atomic data to be written to a PDB file. This should be an 
        object or data structure that is compatible with the `writePDB` function.
    :type atoms: object

    :raises ImportError: If required libraries ('subprocess', 'pathlib', 
        'tempfile', 'open3d') are not installed, an ImportError is raised, 
        specifying which libraries are missing.

    :raises ValueError: If the STL file is not created or is empty, or if the 
        STL file cannot be read as a TriangleMesh,
        a ValueError is raised.

    :returns: An Open3D TriangleMesh object representing the 3D model generated
        from the PDB data.
    :rtype: open3d.geometry.TriangleMesh

    Example usage:
    model = getVmdModel('/path/to/vmd', atoms) """

    required = ['subprocess', 'pathlib', 'tempfile', 'open3d']
    missing = []
    errorMsg = None
    for name in required:
        if not checkAndImport(name):
            missing.append(name)
            if errorMsg is None:
                errorMsg = 'To run getVmdModel, ' \
                'please install {0}'.format(missing[0])
            else:
                errorMsg += ', ' + name

    if len(missing) > 0:
        if len(missing) > 1:
            errorMsg = ', '.join(errorMsg.split(', ')[:-1]) + ' and ' + errorMsg.split(', ')[-1]
        raise ImportError(errorMsg)

    import subprocess
    import tempfile
    import open3d as o3d
    import os

    representation_map = {
    'newcartoon': 'NewCartoon',
    'cartoon': 'NewCartoon',
    'vdw': 'VDW 1.0 20.0',
    'surf': 'Surf',
    'quicksurf': 'QuickSurf 1.0 0.5 0.25 2.0',
    'cpk': 'CPK 1.0 0.3 20.0 20.0'}

    rep_key = representation.lower()
    if rep_key not in representation_map:
        raise ValueError(
            "representation must be one of: 'NewCartoon', 'VDW', 'Surf', " \
            "'QuickSurf', or 'CPK'")
    representation_style = representation_map[rep_key]

    if PY3K:
        from pathlib import Path
    else:
        Path = lambda x: x 

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as temp_pdb:
        temp_pdb_path = Path(temp_pdb.name)
        writePDB(temp_pdb.name, atoms)

    with tempfile.NamedTemporaryFile(suffix=".tcl", delete=False) as temp_script:
        temp_script_path = Path(temp_script.name)

        if PY3K:
            output_path = temp_script_path.parent / "output.stl"
        else:
            output_path = os.path.join(os.path.dirname(temp_script.name), 
                                       "output.stl")

        vmd_script = """
        set file_path [lindex $argv 0]
        set output_path [lindex $argv 1]

        mol new $file_path
        mol modstyle 0 0 %s

        set id_matrix {{1 0 0 0} {0 1 0 0} {0 0 1 0} {0 0 0 1}}
        molinfo top set center_matrix [list $id_matrix]
        molinfo top set rotate_matrix [list $id_matrix]
        molinfo top set scale_matrix [list $id_matrix]

        rendering_method stl
        render STL $output_path

        exit
        """ % representation_style
        
        temp_script.write(vmd_script.encode('utf-8'))

    command = [vmd_path, '-e', str(temp_script_path), '-args', 
               str(temp_pdb_path), str(output_path)]

    try:
        if PY3K:
            subprocess.run(command, check=True)
        else:
            returncode = subprocess.call(command)
            if returncode != 0:
                LOGGER.info("VMD exited with status " + str(returncode) + ".")
    except Exception as e:
        _warn("An unexpected error occurred: " + str(e))
    finally:
        if os.path.exists(temp_script_path):
            os.unlink(temp_script_path)
        if os.path.exists(temp_pdb_path):
            os.unlink(temp_pdb_path)

        if not os.path.exists(output_path) or os.stat(output_path).st_size == 0:
            raise ValueError("STL file was not created or is empty.")

        stl_mesh = o3d.io.read_triangle_mesh(str(output_path))

        if stl_mesh.is_empty():
            raise ValueError("Failed to read the STL file as a TriangleMesh.")

        if os.path.exists(output_path):
            os.unlink(output_path)

        LOGGER.info("Model created successfully.")
    return stl_mesh


def showChannels(channels, model=None, surface=None):
    """Visualizes the channels or pores, and optionally, the molecular model and 
    surface, using Open3D.
    
    This function renders a 3D visualization of molecular channels based on 
    their spline representations. It can also display a molecular model (e.g., 
    the protein structure) and a surface (e.g., cavity surface) in the same 
    visualization. The function utilizes the Open3D library to create and 
    render the 3D meshes.

    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d
    
:arg channels: A list of channel objects or a single channel object. Each
        channel should have a `getSplines()` method that returns two
        interpolators over one parameter domain: one for the centerline and one
        for the radii.
    :type channels: list or single channel object
    
    :arg model: An optional Open3D TriangleMesh object representing the
        molecular model, such as a protein. If provided, this model will be 
        rendered in the visualization.
        Model can be generated using getVmdModel() function.
    :type model: open3d.geometry.TriangleMesh, optional
    
    :arg surface: An optional list containing the surface data. The list should
         have two elements:
        - `points`: The coordinates of the vertices on the surface.
        - `simp`: The simplices that define the surface (e.g., triangles or 
           tetrahedra).
        If provided, the surface will be rendered as a wireframe overlay in the
        visualization.
    :type surface: list (with two numpy arrays), optional
    
    :raises ImportError: If the Open3D library is not installed, an ImportError
        is raised, prompting the user to install Open3D.
    
    :returns: None. This function only renders the visualization.
    
    Example usage:
    showChannels(channels, model=protein_mesh, surface=surface_data) """
    
    if not checkAndImport('open3d'):
        errorMsg = 'To run showChannels, please install open3d (version 0.19.0).'
        raise ImportError(errorMsg)
            
    import open3d as o3d
    
    def create_mesh_from_spline(centerline_spline, radius_spline, n=5):
        N = n * len(centerline_spline.x)
        t = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], N)
        centers = centerline_spline(t)
        radii = radius_spline(t)

        spheres = [
            o3d.geometry.TriangleMesh.create_sphere(radius=r,
                                                    resolution=20).translate(c)
            for r, c in zip(radii, centers)
        ]
        mesh = spheres[0]
        for sphere in spheres[1:]:
            mesh += sphere

        return mesh
    
    if not isinstance(channels, list):
        channels = [channels]
    
    channel_meshes = [create_mesh_from_spline(*channel.getSplines()) for channel in channels]
    meshes_to_visualize = [o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1, origin=[0, 0, 0])]

    if model is not None:
        model.compute_vertex_normals()
        model.paint_uniform_color([0.1, 0.7, 0.3])
        meshes_to_visualize.append(model)
            
    if channel_meshes is not None:
        if not isinstance(channel_meshes, list):
            channel_meshes = [channel_meshes]
        for channel_mesh in channel_meshes:
            channel_mesh.compute_vertex_normals()
            channel_mesh.paint_uniform_color([0.5, 0.0, 0.5])
        meshes_to_visualize.extend(channel_meshes)
        
    if surface is not None:
        points = surface[0]
        simp = surface[1]
        
        triangles = []
        for tetra in simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
            
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        
        meshes_to_visualize.append(line_set)
        
    if len(meshes_to_visualize) > 1:
        o3d.visualization.draw_geometries(meshes_to_visualize)
    else:
        LOGGER.info("Nothing to visualize.")

showPores = showChannels


def showCavities(surface, show_surface=False):
    """Visualizes the cavities within a molecular surface using Open3D.

    This function displays a 3D visualization of cavities detected in a 
    molecular structure.
    It uses the Open3D library to render the cavities as a triangle mesh. 
    Optionally, it can also display the molecular surface as a wireframe 
    overlay.

    To install Open3D use: 
    conda install open3d (for Anaconda users; version open3d-0.19.0 was used 
    during the development) or pip install open3d

    :arg surface: A list containing three elements:
        - `points`: The coordinates of the vertices (atoms) in the molecular 
        structure.
        - `surf_simp`: The simplices that define the molecular surface.
        - `simp_cavities`: The simplices corresponding to the detected cavities.
    :type surface: list (with three numpy arrays)

    :arg show_surface: A boolean flag indicating whether to display the 
        molecular surface
        as a wireframe overlay in the visualization. If True, the surface will 
        be displayed in addition to the cavities. Default is False.
    :type show_surface: bool

    :raises ImportError: If the Open3D library is not installed, an ImportError
        is raised, prompting the user to install Open3D.

    :returns: None

    Example usage:
    showCavities(surface_data, show_surface=True) """
    
    if not checkAndImport('open3d'):
        errorMsg = 'To run showChannels, please install open3d.'
        raise ImportError(errorMsg)
            
    import open3d as o3d
    
    points = surface[0]
    surf_simp = surface[1]
    simp_cavities = surface[2]
    
    triangles = []
    for tetra in simp_cavities:
        triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                          sorted([tetra[0], tetra[1], tetra[3]]),
                          sorted([tetra[0], tetra[2], tetra[3]]),
                          sorted([tetra[1], tetra[2], tetra[3]])])
        
    surface_triangles = np.unique(np.array(triangles), axis=0, 
                                  return_counts=True)[0]
        
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(points)
    mesh.triangles = o3d.utility.Vector3iVector(surface_triangles)
        
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color([0.1, 0.7, 0.3])
        
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.add_geometry(mesh)
        
    if show_surface == True:
        triangles = []
        for tetra in surf_simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
            
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        
        vis.add_geometry(line_set)
            
    vis.get_render_option().mesh_show_back_face = True
    vis.get_render_option().background_color = np.array([1, 1, 1])
    vis.update_renderer()
    vis.run()
    vis.destroy_window()


def showSurfaceCavities(surface, cavities=None, model=None, show_surface=False, 
    cavity_atoms=None, mode='tetra', alpha=4.0, smoothing=0):
    """Visualize surface cavities together with an optional protein model.
    
    This function displays surface cavities calculated by :func:`calcSurfaceCavities`.
    Cavities can be visualized either directly from the tetrahedral/Voronoi
    representation returned by the calculation, or from pseudoatoms loaded from
    a PDB/PQR file through `cavity_atoms`.

    - `mode='tetra'` displays cavities as a tetrahedron-derived mesh.
      This mode follows the original geometric representation most closely,
      but the resulting surface may appear faceted.
    - `mode='smooth'` builds a smoother surface from Voronoi vertices
      assigned to each cavity using an alpha-shape reconstruction. This mode
      requires `cavities` to be provided.

    If `cavity_atoms` is provided, cavities are visualized from the
    coordinates of the supplied pseudoatoms instead of from `surface[2]` or
    `cavities`. 
    
    :arg surface: Surface data returned by :func:`calcSurfaceCavities`.
        Required for `mode='tetra'`, for `mode='smooth'` when
        `cavity_atoms` is not provided, and for displaying the molecular
        surface when `show_surface=True`. The expected list contains:

        - `surface[0]`: atomic coordinates used for the calculation,
        - `surface[1]`: simplices defining the molecular surface,
        - `surface[2]`: merged cavity simplices,
        - `surface[3]`: simplices after surface-layer removal,
        - `surface[4]`: Voronoi vertices used to represent surface cavities.

    :type surface: list or None

    :arg cavities: List of :class:`Cavity` objects returned by
        :func:`calcSurfaceCavities`. Required when `mode='smooth'` and
        `cavity_atoms` is not provided, because the function uses
        `cavity.tetrahedra` to select the corresponding Voronoi vertices.
    :type cavities: list or None

    :arg model: Optional Open3D `TriangleMesh` representing the protein or
        another molecular model. The model can be generated with
        :func:`getVmdModel`.
    :type model: open3d.geometry.TriangleMesh, or None

    :arg show_surface: If `True`, display the molecular surface wireframe
        derived from `surface[1]` in addition to the cavity representation.
        This requires `surface` to be provided. Default is `False`.
    :type show_surface: bool

    :arg mode: Visualization mode used when `cavity_atoms` is not provided.
        Accepted values are `'tetra'` and `'smooth'`. Default is `'tetra'`.
    :type mode: str

    :arg alpha: Alpha value used for alpha-shape surface reconstruction in
        `mode='smooth'` and when visualizing `cavity_atoms`. Smaller values
        produce tighter surfaces, while larger values may connect more distant
        points and generate broader surfaces. Default is 4.0.
    :type alpha: float

    :arg smoothing: Number of Taubin smoothing iterations applied to the
        reconstructed cavity mesh. If `0` or `None`, no smoothing is
        applied. Default is 0.
    :type smoothing: int or None

    :arg cavity_atoms: Optional pseudoatom representation of surface
        cavities. This can be either a path to a PDB/PQR file or a parsed ProDy
        `AtomGroup`, or an Open3D `TriangleMesh` generated, for example, with 
        :func:`getVmdModel`.
    :type cavity_atoms: str, :class:`.AtomGroup`, open3d.geometry.TriangleMesh,
         or None
    
    Examples:
    p = parsePDB('1tqn')
    protein = p.select('protein')
    cavities, surface = calcSurfaceCavities(protein, output_path='cavities.pqr')
    model_protein = getVmdModel(vmd_path, protein)
    cav_model = getVmdModel(vmd_path, parsePQR('cavities.pqr'), representation='QuickSurf')
    
    showSurfaceCavities(surface, model=model_protein, show_surface=True)
    or
    showSurfaceCavities(surface, model=model_protein, cavity_atoms='cavities.pqr', show_surface=True)
    or
    showSurfaceCavities(surface, model=model_protein, cavity_atoms=cav_model, show_surface=True)
    """

    if not checkAndImport('open3d'):
        raise ImportError('To run showSurfaceCavities, please install open3d.')

    import open3d as o3d
    import os

    points = surface[0]
    surf_simp = surface[1]
    simp_cavities = surface[2]
    
    if mode not in ['tetra', 'smooth']:
        raise ValueError("mode must be 'tetra' or 'smooth'")

    meshes_to_visualize = []

    if model is not None:
        model.compute_vertex_normals()
        model.paint_uniform_color([0.8, 0.8, 0.8])
        meshes_to_visualize.append(model)
    
    if cavity_atoms is not None:
        cavity_atoms_given = True
    
        if isinstance(cavity_atoms, o3d.geometry.TriangleMesh):
            cavity_atoms.compute_vertex_normals()
            cavity_atoms.paint_uniform_color([0.1, 0.7, 0.3])
            meshes_to_visualize.append(cavity_atoms)
        
        else:  
            if isinstance(cavity_atoms, str):
                ext = os.path.splitext(cavity_atoms)[1].lower()
                if ext == '.pqr':
                    cavity_atoms = parsePQR(cavity_atoms)
                else:
                    cavity_atoms = parsePDB(cavity_atoms)

            if not hasattr(cavity_atoms, 'getCoords'):
                raise TypeError("cavity_atoms must be a PDB/PQR filename, a " 
                                "ProDy AtomGroup,or an Open3D TriangleMesh.")

            resnums = np.unique(cavity_atoms.getResnums())

            for resnum in resnums:
                sele = cavity_atoms.select('resnum {0}'.format(resnum))
                if sele is None:
                    continue

                pts = sele.getCoords()
                if pts is None or len(pts) < 4:
                    continue

                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(pts)
                cavity_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)

                if smoothing is not None and smoothing > 0:
                    cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)

                cavity_mesh.compute_vertex_normals()
                cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
                meshes_to_visualize.append(cavity_mesh)

    else:
        if surface is None:
            raise ValueError("surface must be provided when cavity_atoms is not given")

        points = surface[0]
        surf_simp = surface[1]
        simp_cavities = surface[2]
    
        if mode == 'tetra':
            triangles = []
            for tetra in simp_cavities:
                triangles.extend([
                    sorted([tetra[0], tetra[1], tetra[2]]),
                    sorted([tetra[0], tetra[1], tetra[3]]),
                    sorted([tetra[0], tetra[2], tetra[3]]),
                    sorted([tetra[1], tetra[2], tetra[3]])])

            surface_triangles = np.unique(np.array(triangles), axis=0, 
                                          return_counts=True)[0]
            cavity_mesh = o3d.geometry.TriangleMesh()
            cavity_mesh.vertices = o3d.utility.Vector3dVector(points)
            cavity_mesh.triangles = o3d.utility.Vector3iVector(surface_triangles)
            cavity_mesh.compute_vertex_normals()
            
            if smoothing is not None and smoothing > 0:
                cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)
                cavity_mesh.compute_vertex_normals()
            
            cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
            meshes_to_visualize.append(cavity_mesh)
            
        elif mode == 'smooth':
            if cavities is None:
                raise ValueError("cavities must be provided for mode='smooth'")
            
            verti = surface[4]
            for cavity in cavities:
                if cavity.tetrahedra is None or len(cavity.tetrahedra) < 4:
                    continue

                pts = verti[cavity.tetrahedra]
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(pts)
                cavity_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)

                if smoothing is not None and smoothing > 0:
                    cavity_mesh = cavity_mesh.filter_smooth_taubin(number_of_iterations=smoothing)

                cavity_mesh.compute_vertex_normals()
                cavity_mesh.paint_uniform_color([0.1, 0.7, 0.3])
                meshes_to_visualize.append(cavity_mesh)

    if show_surface == True:
        triangles = []
        for tetra in surf_simp:
            triangles.extend([sorted([tetra[0], tetra[1], tetra[2]]),
                                  sorted([tetra[0], tetra[1], tetra[3]]),
                                  sorted([tetra[0], tetra[2], tetra[3]]),
                                  sorted([tetra[1], tetra[2], tetra[3]])])
            
        triangles = np.array(triangles)
        triangles.sort(axis=1)
            
        triangles_tuple = [tuple(tri) for tri in triangles]
        unique_triangles, counts = np.unique(triangles_tuple, 
                                             return_counts=True, axis=0)
        surface_triangles = unique_triangles[counts == 1]
            
        lines = []
        for simplex in surface_triangles:
            for i in range(3):
                for j in range(i + 1, 3):
                    lines.append([simplex[i], simplex[j]])
            
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        line_set.lines = o3d.utility.Vector2iVector(lines)
        meshes_to_visualize.append(line_set)

    o3d.visualization.draw_geometries(meshes_to_visualize)

def calcChannels(atoms, output_path=None, separate=False, start_point=None,
    start_point_search=3.0, surf_radius=15, inner_radius=1.2, min_depth=5,
    min_volume=None, max_volume=None, max_depth=None, sparsity=6,
    cavities_only=False, diagram="homogenized", max_deviation=0.1,
    route_divergence=0.2, return_details=False, output_format='pqr', **kwargs):
    """Computes and identifies channels within a molecular structure using 
    Voronoi and Delaunay tessellations.

    This function analyzes the provided atomic structure to detect channels, 
    which are voids or pathways within the molecular structure. It employs 
    Voronoi and Delaunay tessellations to identify these regions, then filters
    and refines the detected channels based on various parameters such as the 
    minimum depth and bottleneck size. The results can be saved to a PQR file 
    (PDB is optional) if an output path is provided. The `separate` parameter 
    controls whether each detected channel is saved to a separate file or if all 
    channels are saved in a single file.

    The implementation is inspired by the methods described in the following 
    publications:
    "MOLE 2.0: advanced approach for analysis of biomacromolecular channels" by
     D. Sehnal, et al., published in J Chemoinform, 5 (39) 2013.
     
     "CAVER: Algorithms for Analyzing Dynamics of Tunnels in Macromolecules". by
     A. Pavelka, et al., published in IEEE ACM T COMPUT BI, (13) 2016.

    "Software Tools for Identification, Visualization and Analysis of Protein 
    Tunnels and Channels". by J. Brezovsky, et al., Biotechnol Adv (31) 2013.


    :arg atoms: An object representing the molecular structure, typically 
        containing atomic coordinates and element types.
    :type atoms: `Atoms` object

    :arg output_path: Optional path to save the resulting channels and
        associated data in PQR (or PDB) format. If None, results are not saved.
        Default is None.

        Naming a directory names nothing after the run: the files are named
        after what they hold - ``channels.pqr`` beside ``links.pqr``, and with
        ``separate=True`` per-object files carrying no stem at all
        (``sp0_chl3.pqr``). Naming a file puts its stem on all of them
        (``run1.pqr`` gives ``run1_links.pqr`` and ``run1_sp0_chl3.pqr``), which
        is what tells two runs apart when they are written side by side.
    :type output_path: str or None

    :arg separate: If True, each detected channel is saved to a separate PDB
        file. If False, all channels are saved in a single PDB file. Default is
        False.

        The per-object files carry the start point they were traced from, so
        that everything belonging to one void globs together: ``out.pqr``
        yields ``out_sp0_chl3.pqr`` for the fourth channel overall, traced from
        start point 0. A chamber link carries the void it arrives at as well,
        after the number - ``out_sp13_lnk0_sp9.pqr`` runs from start point 13
        into start point 9 - so either end can be searched for. Start points are
        numbered largest void first and listed in the log; see ``seed_radius``.
        A search that ran from a single start point tags nothing with it, every
        object having the same one, and writes ``out_chl3.pqr``.

        Ignored for ``output_format="mmcif"``, which always writes one file.
    :type separate: bool

    :arg output_format: What ``output_path`` is written as; ``"pqr"`` (the
        default) or ``"mmcif"``. The format is a property of the one output path
        rather than a second path of its own, so a run has one place its results
        go. Anything else raises, rather than falling back on the default: a
        misspelled format that quietly wrote PQR would be found only by opening
        the file.

        There is no ``"pdb"``. The flag picks the format family and the path
        picks the extension, which is how PDB has always been chosen here -
        ``output_path`` ending in ``.pdb``.

        ``"mmcif"`` writes the sb-ncbr tunnels schema
        (https://github.com/sb-ncbr/tunnels-schema) through
        :func:`writeChannelsCIF`: the profile, the layers and the lining
        residues of every channel in one file, with keys joining them, which is
        what the PQR and the residue text files cannot express between them.
        Chamber links go into the same file under
        ``_sb_ncbr_channel.type`` ``Path``, and a directory takes
        ``channels.cif``.
    :type output_format: str

    :arg start_point: Optional starting point for channel search. This can be
        either a 3D coordinate point or an atomic selection/AtomGroup. If the
        3D coordinate point will be provided, the algorithm seeds the channel
        search near this point (overriding the default automatic seed selection
        based on the deepest tetrahedron); see ``start_point_search`` for how
        the seed tetrahedron itself is picked. Coordinates must be given in Å.
        If an atomic selection is provided, its geometric center is used as the
         starting point.

        A start point names a site, so it does more than move a seed: the search
        is restricted to the single cavity holding the tetrahedron nearest the
        point, and channels are reported for that site alone rather than one
        bundle per cavity in the structure. The automatic passes that decide
        where to start are skipped with it -- the chamber seeding and the
        ``seed_volume`` floor both defer to the point.
    :type start_point: list, tuple, or ndarray (length 3), :class:`.Atomic`, or None

    :arg start_point_search: Only used when ``start_point`` is provided. Radius,
        in Angstrom, of the neighbourhood of ``start_point`` searched for the seed
        tetrahedron. The tetrahedron nearest ``start_point`` is often a tight one,
        and since every channel of the cavity starts there, its inscribed radius
        caps all of their bottlenecks and appears as one shared bottleneck at the
        joint beginning of the bundle. Seeded instead is the widest tetrahedron within
        ``start_point_search`` of ``start_point`` that belongs to the same cavity, is
        at least ``min_depth`` below the surface (so the seed cannot drift out towards
        the mouth) and is reachable from the nearest one through that neighbourhood (so
        it stays in the void the point sits in rather than crossing a wall). The seed
        may therefore sit a little shallower than ``start_point`` itself, which is
        usually placed on a ligand or a catalytic residue and often lies deeper than
        the widest part of the pocket around it. The nearest tetrahedron is kept,
        whatever its own depth, unless a wider one qualifies, so the seed is never
        narrower than the tetrahedron at the point.

        It is a requirement and not only a search budget: the search must begin
        within this distance of the point, and if no cavity has a tetrahedron that
        close, no channels are computed and a warning reports how far the nearest one
        is and where it lies, so that ``start_point`` can be corrected or
        ``start_point_search`` raised. Seeding the nearest tetrahedron however far
        away it sits would answer a point that misses the void -- the centroid of a
        residue selection often lands inside an atom -- with channels through whatever
        cavity happens to lie nearest, and nothing in the result would say so.

        Default is 3.0; use 0 to seed the nearest tetrahedron as-is, which asks for no
        neighbourhood and so imposes no distance requirement either.
    :type start_point_search: float

    :arg surf_radius: Radius, in Angstrom, of the probe that says what counts as
        the outside. The tessellation is eroded from the boundary inward wherever
        this probe fits, so everything it can reach is exterior and what it cannot
        enter is the void the channels are traced through. Default is 15.

        It is deliberately large. The erosion is followed by the local peel of
        ``min_enclosure``, which strips the shell of exterior the probe bridged
        over rather than entered, so where the erosion *stops* is decided by
        burial and not by this radius. Above about 3.5 Angstrom the reported
        channels stop moving: across the test structures the result is identical,
        channel for channel, from 3.5 up to 20, and the run is no slower at the
        top of that range than at the bottom. The default sits high in the
        plateau rather than at its edge.

        What the value has to be large enough for is a wide pore. A probe smaller
        than the pore passes through it, so the lumen is classified as outside and
        carved away, and only the pockets around it are reported; the structure is
        turned inside out with no sign that anything went wrong. The probe must
        therefore exceed the radius of the widest opening that should count as
        interior - for a channel-forming protein or a large assembly that is well
        above the 3 Angstrom values common for compact enzymes, hence the default.

        For a globular protein any value in the plateau gives the same answer, so
        lowering it is safe but buys nothing. Note that ``calcSurfaceCavities``
        does not use this default: a surface pocket is shallow and open by
        definition, and it sets its own, smaller radius.

        One caveat applies to ``diagram="weighted"`` only. That path truncates the
        Apollonius diagram at ``max(2 * surf_radius, 8)`` Angstrom of clearance, so
        a large radius here makes an already expensive tessellation much more so.
        The weighted diagram is experimental; with the default radius, expect it
        to be impractical on anything but small structures.
    :type surf_radius: float

    :arg inner_radius: The second radius threshold used to define the inner surface of
        the channels. Default is 1.2, which is the smallest value safe on a
        structure without hydrogens.

        Below about 1.2 Angstrom the probe is smaller than a water molecule, and
        then the structure must carry explicit hydrogens. Without them, every
        carbon keeps its full vdW radius while the space its hydrogens occupied is
        left empty, and a sub-water probe is small enough to thread those
        interstices: the interior percolates into a sponge and the channel count
        can rise several-fold. At 1.2 Angstrom and above, protonated and
        unprotonated structures give similar channels, and an X-ray file may
        be used as it comes. A warning is issued for the unsafe combination, so
        go below the default only on a fully protonated structure -- partial
        protonation does not count, since it is the apolar hydrogens that fill
        those interstices.

        The floor is a workaround for absent hydrogens rather than a property of
        the probe, so on a structure that carries them the default is coarser
        than the input deserves and narrow connections are reported closed
        instead of measured. That case is noted too, and about 0.9 Angstrom is
        then a reasonable probe.

        Note that this sets where channels are traced, not how wide the reported
        ones end up being: a channel can be narrower than ``inner_radius`` at its
        tightest point. Use ``bottleneck`` to put a floor on that.
    :type inner_radius: float

    :arg min_depth: The minimum depth, in Angstrom, a cavity must reach to be
        considered as a channel. Depth is the geodesic distance from the cavity's
        surface opening to its farthest point along the Voronoi network, so it is a
        physical length independent of the tessellation density. Default is 5.
    :type min_depth: float

    :arg max_depth: Maximum cavity depth, in Angstrom. Portions of a cavity deeper
        than this value are trimmed away. Default is None (no trimming).
    :type max_depth: float

    :arg min_volume: Minimum volume, in cubic Angstrom, required for a channel or
        cavity to be retained. Default is None.

        Like ``bottleneck``, it drops entries from the finished list and does not
        change the search: where a search may *start* is set by ``seed_volume``,
        the floor on the void - cavity or chamber - a search runs from, so raising
        this to report only the large channels never costs you the seed that finds
        them. The two are also measured on different scales - a channel volume is
        the probe sphere swept along the centerline, hence free space, while a
        cavity or chamber volume sums the Delaunay tetrahedra it holds - so the
        same number does not mean the same thing on both.
    :type min_volume: float

    :arg max_volume: Maximum volume allowed for a channel/cavity to be 
        retained. Default is None.
    :type max_volume: float

    :arg sparsity: Smallest centre-to-centre distance, in Angstrom, at which two
        channel surface openings (mouths) still count as separate. Two channels
        leaving closer than this are treated as sharing one opening and are merged
        if they also take the same corridor (see ``route_divergence``). It acts as a
        floor: where the tessellation already measures the two mouths as wider
        than ``sparsity``, their own radii decide instead, so the value only sets
        the minimum separation a reported pair of openings can have. Being applied
        *after* the search, it can only merge channels, never hide one, and is a
        reporting preference rather than part of the geometry. A higher value
        reports fewer channels. It has no effect on the cavities, which are found
        from the exit tetrahedra before any merging. Default is 6.
    :type sparsity: float

    :arg diagram: 
        "homogenized" (default) - every atom is substituted by a set of 
        homogeneous balls whose common radius equals the smallest van der Waals
        radius present in the structure, before building the Voronoi and 
        Delaunay tessellations. This yields an accurate estimate of the 
        additively weighted Voronoi diagram from an ordinary one, as done in 
        MolAxis and CAVER 3
        "simple" - the original atoms are used with their individual van der 
        Waals radii directly. This is very inaccurate and should be avoided in 
        almost all cases.
        "weighted" - the true additively-weighted (Apollonius) Voronoi diagram is
        built directly from the atoms with their individual van der Waals radii,
        using the third-party ``vorpy`` package (with a compiled kernel when
        ``numba`` is available). This is the exact diagram the "homogenized" mode
        approximates, at a higher cost (~ 100x slower, for 4000 atoms). To use this
        approach install ``vorpy`` library using ``pip install vorpy3``.

        Treat "weighted" as experimental. Besides the cost of the tessellation
        itself, it is the one mode whose diagram depends on ``surf_radius``: it is
        truncated at ``max(2 * surf_radius, 8)`` Angstrom of clearance, so the
        default radius makes it far more expensive again, and lowering
        ``surf_radius`` to keep it affordable is the trade-off to be aware of.
    :type diagram: str

    :arg max_deviation: Maximum tolerated deviation, in Angstrom, between the 
        union surface of the substitute balls and the original van der Waals 
        surface when ``diagram = homogenized`` . It controls the trade-off
        between surface accuracy and the number of balls generated: an atom 
        whose radius exceeds the smallest radius (``rho``) by more than 
        ``max_deviation`` is filled with several balls, otherwise it is kept as
         a single ``rho`` ball. Default is 0.1. Guideline values:

        * ``0.1`` fine accurate surface with minimal errors, but on 13-15x 
            more balls than original
        * ``0.15`` in heavy-atom-only structures it starts filling carbon, which
             is otherwise left as a single ball with a uniform ~0.18 A inset).
        * ``0.2`` speed optimized ; e.g. carbon fills to ~15 balls when 
            hydrogens are present (``rho``=1.2), resulting roughly to ~8 times 
            more balls. Without hydrogens, carbons are single balls.

        Only used when ``diagram = homogenized``.
    :type max_deviation: float

    :arg route_divergence: How far, in Angstrom, two channels may part company
        per Angstrom of channel, and still count as one. Each channel is a tube -
        the centerline plus the clearance around it - and at every point the
        measurement is the gap between the two tube *surfaces*: negative where
        they overlap, since two clearance balls that share a point have no atom
        between them and the routes are locally the same corridor. The gaps that
        are positive are integrated along both routes and divided by the average
        of the two channel lengths, so how far apart they get and how much of the
        route is involved enter one number: a two Angstrom arm off a ninety
        Angstrom channel reads differently from the same arm off a ten Angstrom
        one. ``0.0`` merges only pairs whose tubes touch along their whole
        length; larger values tolerate more divergence.

        Measuring surfaces rather than centerlines is what lets a single value
        serve everywhere, the question having no absolute scale: two paths a
        couple of Angstrom apart in a wide chamber have nothing between them,
        while the same distance in a narrow throat spans a wall.

        Two channels are merged (cheapest kept) only when they take the same
        corridor **and** leave through the same opening (see ``sparsity``); a
        corridor that forks near the surface and exits twice is one tunnel, but
        two different corridors to one opening, or one corridor reaching two
        openings, are two tunnels. The comparison is geometric rather than a
        shared prefix of tetrahedra, so it is unaffected by *where* two routes
        diverge - variants that split and rejoin still count as one - and by
        ``max_deviation``, a tetrahedron count not being mesh-invariant where
        Angstrom are.

        Default is 0.2, inside the band over which the answer does not change:
        over every comparison the dedup makes on the test corpus, the pairs it
        merges score at most 0.161 and the pairs it keeps at least 0.208. That
        band is narrow, and the two ways out of it are not equally bad - set too
        high it merges a corridor away, and a channel that is not reported
        cannot be noticed, while set too low it reports a duplicate, which can
        be seen and judged.
    :type route_divergence: float

    :arg return_details: If True return an additional dictionary containing
        internal calculation data, including the channel calculator, simplices,
        neighboring tetrahedra, Voronoi vertices, atomic coordinates, and van der
        Waals radii. Default is False.
    :type return_details: bool

    The remaining options below are not part of the signature and are accepted
    **as keyword arguments only**; a normal run never touches them. Any other
    keyword raises :exc:`TypeError` rather than being ignored, so a misspelled
    option is reported instead of silently falling back to its default.

    :arg bottleneck: Minimum bottleneck radius, in Angstrom, a channel must have
        to be reported. Default is ``inner_radius``, so a channel the traversal
        probe itself would not fit through is not reported; pass ``bottleneck=0``
        to switch the filter off and see every channel the search traced.

        Set it whenever the question is what can actually pass through, because
        ``inner_radius`` alone does not guarantee it: channels regularly come out
        somewhat narrower than the probe that found them. The difference is minor
        with ``diagram="homogenized"`` and a small ``max_deviation``, but with
        ``diagram="simple"`` a channel can be several times narrower, and there
        this is the only real width control. Useful values are 0,
        ``inner_radius`` itself, or the radius of the ligand or ion of interest
        if that is what you are screening for.

        Unlike ``inner_radius``, it does not change the search: it drops entries
        from the finished list, before they are numbered and written to file,
        rather than rerouting them. Filtering is therefore cheap, but it cannot
        recover a wide route that the search did not take - if raising it leaves
        you with too few channels, raise ``inner_radius`` instead and let the
        channels be traced afresh, or set ``prune_narrow_edges``, which makes
        this floor bound the search as well.
    :type bottleneck: float

    :arg prune_narrow_edges: Whether to leave Voronoi edges narrower than
        ``bottleneck`` out of the search graph, instead of routing through them
        and discarding the resulting channel afterwards. Default ``False``.

        The two behaviours are the two established approaches: MOLE filters the
        finished channels, CAVER prunes the graph. Filtering alone can lose a
        corridor. Where a cheap route to a mouth pinches below the floor and a
        wider route to the same mouth costs more, the search returns the cheap
        one, the filter drops it, and the wider route - open, and the one asked
        for - is never reported, because only one path per exit is returned.
        Pruning takes those edges out of the search, so the wider way round is
        what the search finds.

        Pruning cannot cost a channel. A path over a dropped edge pinches at or
        below that edge's gate, so the floor rejects it either way, while a path
        that survives the floor crosses no dropped edge and keeps its cost, so it
        remains the cheapest route to its exit. The difference is one-sided: the
        same channels, plus any the filter was masking.

        Off by default because it changes what is reported, and the channels it
        adds are ones no earlier run listed. ``bottleneck=0`` prunes nothing
        whatever this is set to, there being no floor to prune at.
    :type prune_narrow_edges: bool

    :arg seed_radius: Probe radius, in Angstrom, that decides where a channel may
        *start*, as opposed to ``inner_radius``, which decides what it may pass
        *through*. Default is ``max(1.4, inner_radius)``.

        The cleared void is carved a second time with this larger probe, and each
        connected piece that survives - a *chamber* - seeds one channel search, so
        a large or branched cavity is reported from each of its lobes. Without
        chambers the whole cavity gets a single seed at its deepest tetrahedron,
        which is right for a compact pocket but arbitrary for a branched one,
        since the necks between the real pockets are themselves cavity under a
        small probe and the deepest point can land in any one of them. Every
        channel is then at least as long as the cavity is deep, and the openings
        of the other lobes go unreported, because a route from that far seed is
        absorbed at a nearer opening before it reaches them.

        A larger value separates the lobes more strictly, at the price of
        discarding the narrow ones entirely; a smaller one merges them back into a
        single chamber. It has to exceed ``inner_radius`` to prune anything, and
        at or below it the chamber pass is skipped and the single deepest seed is
        used instead. Passing ``start_point`` overrides it completely, since an
        explicit starting point is where the search begins.
    :type seed_radius: float

    :arg seed_volume: Smallest void, in cubic Angstrom, that may seed a channel
        search, measured on the cavity scale (Delaunay tetrahedra summed).
        Default is 50.

        It applies to both kinds of search site, on that one scale: to every
        chamber (see ``seed_radius``), and to every cavity, since a cavity holding
        no chamber is itself searched whole from a seed of its own. A void under it
        is taken for tessellation debris rather than a site.

        The two cases differ in what that costs. A *chamber* under the floor, in a
        cavity that clears it, costs nothing: the cavity is searched from its other
        seeds either way, and a chamber that seeds nothing still conducts, so the
        routes through it are still found. A *cavity* under the floor is not
        searched at all, and neither is anything inside it - which is the point:
        the tessellation leaves single-tetrahedron slivers lying wholly in the
        surface layer, and such a sliver is its own mouth. It reports either
        nothing (it has no target to path to, and comes out as sealed) or a
        one-step channel that is a facet of the surface rather than a tunnel. A low
        ``min_depth`` is what exposes them in numbers, but nothing except their size
        tells them from a real void at any depth.

        Raise it to seed only the roomy lobes of a branched cavity; ``None``
        applies no floor and seeds every cavity and chamber found. Passing
        ``start_point`` overrides it completely, since an explicit starting point
        is where the search begins.

        It is deliberately not ``min_volume``: that one is a floor on what gets
        *reported*, and tying the two would mean that asking for only the large
        channels also removed the seeds that trace them.
    :type seed_volume: float or None

    :arg max_seeds: Largest number of chambers one cavity may be seeded at, widest
        chamber first; ``None`` for no cap. Default is 20. Every seed adds
        candidate routes that have to be compared against all the others, and a
        percolated interior (a small ``inner_radius`` on a structure without
        hydrogens) can offer dozens of chambers, so the cap keeps the run time
        bounded on exactly those structures whose chambers are least trustworthy.
    :type max_seeds: int or None

    :arg chamber_links: Whether to report how the chambers of one cavity connect
        to each other. Default is ``True``; has no effect unless the cavity was
        seeded per chamber (see ``seed_radius``).

        A deep chamber often has no opening of its own and reaches the solvent
        only through a shallower chamber. Its own route out then duplicates the
        shallower chamber's shorter channel and is dropped as such, which would
        leave the deep chamber's access unreported. A *link* is that access: the
        route from the deep chamber, cut where it first joins a shallower one, so
        that the whole way out reads as ``link(deep -> shallow) +
        channel(shallow -> surface)``. The link's bottleneck is the neck between
        the two chambers, which is what governs whether anything can pass. Links
        run from deep to shallow only, so a pair of chambers yields one link,
        and the two are ranked by how close each comes to the surface anywhere
        along itself rather than by the depth of its seed - a wide chamber that
        opens straight to the solvent can still hold the deepest seed of the
        cavity, and it is the pockets hanging off it that need the link.

        Links are returned in ``details['links']``, which requires
        ``return_details=True``, and not in the returned channel list. They are
        written beside the channels as ``<stem>_links.pqr``, or as
        ``<stem>_sp<from>_lnk<i>_sp<to>.pqr`` per link when ``separate=True``,
        both ends being named because a link is the one object that has two of
        them; the same pair is written into the record as ``from spN to spM``.
        :func:`getLinkParameters` and :func:`getLinkResidueNames` measure and
        line them. They are filtered by ``bottleneck`` like the
        channels - a neck the probe cannot pass is not a way through - but not by
        ``min_volume``, which is a floor meant for objects that reach the surface,
        whereas a link is by construction only the neck.

        With ``chamber_links=False`` a route that leaves its own chamber is no
        longer cut where it joins another one, and is reported whole as a channel
        instead. Switching links off therefore gives more channels rather than the
        same ones minus the links.

        A chamber whose channels were all dropped as duplicates and whose link is
        then dropped by ``bottleneck`` reports nothing at all. That is a finding,
        not an omission: the chamber is sealed at this width. Because an empty
        report is easy to miss, the number of such chambers is logged, and
        lowering ``bottleneck`` shows how they connect.
    :type chamber_links: bool

    :arg weighted_cache: Cache the raw additively-weighted Voronoi diagram to disk so
        that re-running ``diagram="weighted"`` on the same structure skips the
        expensive vorpy tessellation (which dominates the ~10 min run time). The
        diagram only depends on the atoms and ``surf_radius``, so re-runs that change only
        ``inner_radius``, ``bottleneck``, ``sparsity``, ``start_point`` etc. reuse it. ``True``
        (default) caches next to ``output_path`` (or, absent one, under the structure
        title in the current directory); pass a path to place it explicitly, or
        ``False`` to disable. The cache is keyed by content, so editing the structure
        or ``surf_radius`` transparently forces a recompute. Only used for ``diagram="weighted"``.
    :type weighted_cache: bool or str

    :arg weighted_mouth_depth: Only used for ``diagram="weighted"``. The additively-
        weighted (Apollonius) tessellation is not a clean simplicial complex, so it
        leaves false interior boundary faces that the pipeline would misread as
        surface openings, truncating channels to stubs. To repair this a *homogenized*
        diagram of the same atoms is built as an interior/exterior oracle, and only
        exit tetrahedra whose Voronoi vertex lies within ``weighted_mouth_depth``
        Angstrom (geodesic distance below the molecular surface) are treated as mouths.
        Default 2.5 (the value at which the recovered channels match the
        ``"homogenized"`` result); ``None`` disables the relabeling.
    :type weighted_mouth_depth: float or None

    :arg min_enclosure: Fraction of directions that must be blocked by protein for
        a tetrahedron to count as interior, in ``[0, 1]``. Default is 0.70.

        Once the surf_radius surface is built it is eroded inward with the inner_radius probe, to
        strip the shell of true exterior that an surf_radius probe bridges over rather than
        enters (the "moat"); that shell would otherwise join the cavity and offer
        wide, low-cost routes along the outside of the protein. Erosion continues
        while the tetrahedra at the front are *open*, meaning that fewer than
        ``min_enclosure`` of the directions leaving them run into protein within
        the reach of :meth:`.ChannelCalculator.calcEnclosure`, and halts at the
        first buried layer.

        Bounding the erosion by size instead does not work. A count of tetrahedron
        layers is not mesh-invariant, as a layer is one tetrahedron thick and
        tetrahedra shrink as ``max_deviation`` is lowered. A depth in Angstrom is
        not ``surf_radius``-invariant, as the moat is as deep as ``surf_radius - inner_radius`` in a concavity
        but vanishes on a flat face, so a depth that clears it where it is thick
        also marches down the channel mouths and erodes the channels themselves.
        Testing burial locally instead leaves ``surf_radius`` to decide only where the
        erosion starts, not where it stops, so results are independent of it, and
        ``surf_radius`` is left doing the one job it should: capping the mouths.

        It is bounded from both sides, and the window is narrow.

        Too low and the moat is not fully stripped. Too high and the erosion 
        never meets a layer buried enough to stop it, so it percolates down the
          channels and eats the cavity. The default sits at the floor, which is
        the safe end: over-peeling deletes real channels, whereas
        under-peeling shows up as surface-riding routes that can be recognised.
    :type min_enclosure: float

    :arg max_peel_depth: Optional hard cap, **in Angstrom**, on how far the peel
        above may advance from the surf_radius surface. ``None`` (default) is uncapped, and
        the enclosure test alone decides where erosion stops. Set it only as a
        backstop on a structure where the peel misbehaves; it is deliberately not
        tied to ``surf_radius``, since a cap that scales with ``surf_radius`` reintroduces exactly
        the ``surf_radius`` dependence that ``min_enclosure`` exists to remove.
    :type max_peel_depth: float or None
        
    :arg edge_cost: How each Voronoi edge is priced in the Dijkstra tunnel search.
        ``"integral"`` prices each edge by the integral of its clearance profile
        along the edge, which is mesh-invariant. ``"bottleneck"`` uses the legacy
        ``length / (gate**2 + eps)``, charging the whole edge at its single
        narrowest point, whose value (and the routing it produces) drifts as
        ``max_deviation`` coarsens. ``None`` (default) selects ``"integral"`` for
        ``diagram="homogenized"``/``"simple"`` (straight edges, where the integral
        is exact) and ``"bottleneck"`` for ``diagram="weighted"``. ``"integral"``
        is rejected for ``diagram="weighted"`` (Apollonius edges are arcs the
        straight-chord integral cannot price). The reported bottleneck radius is
        unaffected by this choice.
    :type edge_cost: str or None

    :returns: A tuple containing two elements:
        - `channels`: A list of detected channels, where each channel is an 
          object containing information about its path and geometry.
        - `surface`: A list containing additional information for further 
          visualization, including the atomic coordinates, simplices defining 
          the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein 
        and calculates van der Waals radii. When ``homogenize`` is True, each 
        atom is replaced by homogeneous balls of the smallest radius present 
        so that an ordinary tessellation approximates the additively weighted 
        Voronoi diagram. It then performs 3D Delaunay triangulation and Voronoi
         tessellation on the resulting coordinates.
    2. **State Management:** Creates and updates different stages of channel 
        detection of the protein structure to filter out simplices based on the
         given radii.
    3. **Surface Layer Calculation:** Determines the surface and second-layer 
        simplices from the filtered results.
    4. **Cavity and Channel Detection:** Finds and filters cavities based on 
        their depth and calculates channels using Dijkstra's algorithm.
    5. **Visualization and Saving:** Generates meshes for  detected channels, 
        filters them by bottleneck size, and either saves the results to a PDB 
        file or visualizes them based on the specified parameters.
       
    Example usage:
    channels, surface = calcChannels(atoms, output_path="channels", separate=True)
    
    channels, surface = calcChannels(atoms, output_path="all_channels.pdb", 
                                     start_point=[-22.312, -20.065, -11.144])
    
    start_sel = protein.select('resid 212 309 483')
    channels, surface = calcChannels(atoms, output_path="all_channels.pdb", 
                                     start_point=start_sel)
    
    To save the results as PDB file:
    channels, surface = calcChannels(atoms, output_path="channels.pdb",
                                     separate=False, surf_radius=15, inner_radius=1.2, min_depth=5,
                                     bottleneck=1, sparsity=6) """

    # Advanced options, accepted as keyword arguments only and kept out of the
    # signature above, which is long enough already. These are settings a normal
    # run never touches.
    CHANNELS_ADVANCED_OPTIONS = {
        'bottleneck': inner_radius,
        'prune_narrow_edges': False,
        'seed_radius': max(1.4, inner_radius),
        'seed_volume': 50.0,
        'max_seeds': 20,
        'chamber_links': True,
        'min_enclosure': 0.70,
        'max_peel_depth': None,
        'edge_cost': None,
        'weighted_cache': True,
        'weighted_mouth_depth': 2.5,
        'min_tetrahedra': None, 
        'max_tetrahedra': None, 
    }

    # Unknown keywords are an error rather than silently ignored: a misspelled
    # option would otherwise be dropped without a trace and the run would quietly
    # proceed on the default, which is the failure mode an explicit signature
    # prevents. Checked first, so a typo fails before any of the work.
    unexpected = sorted(set(kwargs) - set(CHANNELS_ADVANCED_OPTIONS))
    if unexpected:
        raise TypeError('calcChannels() got an unexpected keyword argument {0}. '
                        'Keyword-only options are: {1}.'.format(
                            ', '.join(repr(key) for key in unexpected),
                            ', '.join(sorted(CHANNELS_ADVANCED_OPTIONS))))

    options = dict(CHANNELS_ADVANCED_OPTIONS, **kwargs)
    bottleneck = options['bottleneck']
    prune_narrow_edges = options['prune_narrow_edges']
    seed_radius = options['seed_radius']
    seed_volume = options['seed_volume']
    max_seeds = options['max_seeds']
    chamber_links = options['chamber_links']
    min_enclosure = options['min_enclosure']
    max_peel_depth = options['max_peel_depth']
    edge_cost = options['edge_cost']
    weighted_cache = options['weighted_cache']
    weighted_mouth_depth = options['weighted_mouth_depth']
    min_tetrahedra = options['min_tetrahedra']
    max_tetrahedra = options['max_tetrahedra'] 
    

    required = ['heapq', 'collections', 'scipy', 'pathlib', 'warnings']
    missing = []
    errorMsg = None
    for name in required:
        if not checkAndImport(name):
            missing.append(name)
            if errorMsg is None:
                errorMsg = 'To run calcChannels, please install {0}'.format(missing[0])
            else:
                errorMsg += ', ' + name

    if len(missing) > 0:
        if len(missing) > 1:
            errorMsg = ', '.join(errorMsg.split(', ')[:-1]) + ' and ' + errorMsg.split(', ')[-1]
        raise ImportError(errorMsg)

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path
    
    if start_point is not None:
    
        if hasattr(start_point, 'getCoords'):
            if start_point.numAtoms() == 0:
                raise ValueError("start_point selection contains no atoms")
            start_point = calcCenter(start_point)

        elif not isListLike(start_point):
            raise TypeError("start_point must be a selection/AtomGroup or a list/tuple/ndarray "
                "with three numeric values")

        if len(start_point) != 3:
            raise ValueError(
                "start_point must be a selection/AtomGroup or a list of three numbers, e.g. "
                "start_point=[-12.312, 5.065, -1.144]")
        
        start_point = np.array(start_point, dtype=float)        
        
        LOGGER.info("Using user-provided start_point for channel seed: [{:.3f}, {:.3f}, {:.3f}] Å"
            .format(start_point[0], start_point[1], start_point[2]))

    LOGGER.timeit('_prody_calcChannels')

    # Edge-cost mode for the Dijkstra routing (buildSparseGraph). Default: the
    # mesh-invariant profile integral for the straight-edged homogenized/simple
    # diagrams, the legacy l/(d^2+b) for weighted (Apollonius edges are arcs,
    # where the straight-chord integral is only approximate; weighted is
    # experimental). An explicit value overrides the per-diagram default.
    if edge_cost is None:
        edge_cost = 'bottleneck' if diagram == 'weighted' else 'integral'
    elif edge_cost not in ('integral', 'bottleneck'):
        raise ValueError("edge_cost must be 'integral', 'bottleneck' or None, "
                         "got {0!r}".format(edge_cost))
    elif edge_cost == 'integral' and diagram == 'weighted':
        raise ValueError("edge_cost='integral' is only valid for the straight-edge "
                         "diagrams ('homogenized'/'simple'); the weighted "
                         "(Apollonius) diagram has arc edges the straight-chord "
                         "integral cannot price. Use edge_cost='bottleneck' (the "
                         "default for diagram='weighted') or None.")
    
    _reportAtomsInputComposition(atoms, inner_radius, diagram)
    atoms = atoms.select('not water') # water is excluded from the selection
    calculator = ChannelCalculator(atoms, inner_radius=inner_radius, sparsity=sparsity,
                                   edge_cost=edge_cost)

    coords = atoms.getCoords()
    vdw_radii = calculator.getVdwRadii(atoms.getElements())
    # Burial is a property of the protein, not of the tessellation, so the enclosure
    # test that strips the moat runs against the real atoms rather than the balls the
    # diagram happens to be built on. Homogenization would otherwise make it a
    # function of max_deviation, which min_enclosure must not be.
    atom_coords = coords
    # For diagram="weighted" only: a homogenized-surface depth oracle used to relabel
    # the additively-weighted diagram's leaky surface mouths (see getSurfaceCavities).
    mouth_oracle = None

    if diagram == "homogenized":
        LOGGER.timeit('_prody_channels_homogenize')
        coords, vdw_radii = calculator.homogenizeAtoms(coords, vdw_radii, max_deviation)
        LOGGER.report("Substituted {0} atoms with {1} homogeneous balls of radius {2:.2f} Å in %.2fs.".format(
            atoms.numAtoms(), len(coords), float(vdw_radii[0])), '_prody_channels_homogenize')

    LOGGER.timeit('_prody_channels_tessellation')
    if diagram == "weighted":
        # True additively-weighted (Apollonius) Voronoi network via the third-party
        # vorpy package: van der Waals radii are baked into the diagram exactly,
        # instead of being approximated by homogenising atoms into uniform balls.
        # buildAwTessellation returns the same (simplices, neighbors, vertices) triple
        # a scipy Delaunay would, so the downstream erosion/cavity pipeline is
        # untouched. Because every AW vertex is equidistant (additively) to its 4
        # tangent atoms, the sum-based clearance test in deleteSimplices3d reduces
        # exactly to the per-atom clearance, so the returned clearances are not needed.
        if not checkAndImport('vorpy'):
            raise ImportError('diagram="weighted" requires the vorpy package for the '
                'additively-weighted (Apollonius) Voronoi diagram. Install vorpy, or '
                'use diagram="homogenized"/"simple".')
        # The compiled calc_vert kernel needs numba; without it the weighted path
        # still works but is ~5x slower, so fall back with a warning rather than fail.
        accelerate = checkAndImport('numba')
        if not accelerate:
            _warn('numba is not installed; the additively-weighted tessellation '
                'will run without the compiled kernel and may be very slow.')
        from ._vorpy_aw import buildAwTessellation, resolveCachePath
        try:
            title = atoms.getTitle()
        except Exception:
            title = None
        cache_path = resolveCachePath(weighted_cache, output_path, title)
        simplices, neighbors, verts, _ = buildAwTessellation(
            coords, vdw_radii, max_vert=max(2.0 * surf_radius, 8), accelerate=accelerate,
            cache=cache_path)
        LOGGER.report('Additively-weighted (Apollonius) tessellation of {0} atoms '
            'constructed in %.2fs.'.format(len(coords)),
            '_prody_channels_tessellation')
        # The AW->simplicial mapping leaves false interior boundary faces that the
        # pipeline would misread as surface mouths (collapsing channels to stubs).
        # Build a homogenized diagram of the same atoms as an interior/exterior depth
        # oracle; getSurfaceCavities then keeps only exit tetrahedra within
        # weighted_mouth_depth Angstrom (geodesic) of the true molecular surface.
        if weighted_mouth_depth is not None:
            LOGGER.timeit('_prody_channels_mouth_oracle')
            mouth_oracle = calculator.buildSurfaceDepthOracle(
                coords, vdw_radii, surf_radius, max_deviation, weighted_mouth_depth)
            LOGGER.report('Homogenized surface oracle (weighted mouth relabeling) '
                'built in %.2fs.', '_prody_channels_mouth_oracle')
    else:
        from scipy.spatial import Delaunay
        # We deliberately do NOT joggle/jitter the input (no QJ), unlike CAVER,
        # which perturbs by ~0.001 A to dodge the cospherical "T5" degeneracy it
        # reports in nearly every structure. scipy's default Qhull options
        # (Qbb Qc Qz) merge cospherical facets instead of joggling, so the
        # circumcenters stay finite even at degeneracies (verified: 0 NaN/inf
        # across millions of tetrahedra, even under heavy homogenized refinement).
        # Not joggling keeps the pipeline exactly reproducible (homogenizeAtoms is
        # a fixed Fibonacci lattice, and nothing here uses an RNG). The only
        # residual is coincident circumcenters at true degeneracies, handled where
        # it matters by the twin-tetrahedron guard in _edgeBottleneck.
        dela = Delaunay(coords)
        # circumcenters straight from the Delaunay paraboloid lifting, so we
        # skip the redundant second Qhull pass (scipy Voronoi). Numerically identical
        # to voro.vertices for points in general position.
        simplices = dela.simplices
        neighbors = dela.neighbors
        verts = calculator.calcCircumcenters(dela)
        LOGGER.report('Delaunay tessellation of {0} points constructed in %.2fs.'.format(
            len(coords)), '_prody_channels_tessellation')

    LOGGER.timeit('_prody_channels_surface')
    s_prt = State(simplices, neighbors, verts)
    
    if PY3K:
        s_tmp = State(*s_prt.getState())
        s_prv = State(None, None, None)
    else:
        s_tmp = apply(State, s_prt.getState())
        s_prv = State(None, None, None) 
        
    while True:
        s_prv.setState(*s_tmp.getState())
        
        if PY3K:
            s_tmp.setState(*calculator.deleteSimplices3d(coords, *(s_tmp.getState() + tuple([vdw_radii, surf_radius, True]))))
        else:
            tmp_state = calculator.deleteSimplices3d(coords, *(s_tmp.getState() + [vdw_radii, surf_radius, True]))
            s_tmp.setState(*tmp_state)

        if s_tmp == s_prv:
            break
        
    s_srf = State(*s_tmp.getState())

    # Moat removal: erode the surf_radius surface inward with the inner_radius probe, stripping the shell
    # of true exterior that a large surf_radius probe bridges over instead of entering (it would
    # otherwise join the cavity and offer wide, low-cost routes along the outside).
    # Erosion stops where the tetrahedra stop being open to the solvent, which is a
    # local criterion, so neither the mesh nor surf_radius sets how deep the peel goes.
    s_srf = State(*calculator.peelSurfaceByEnclosure(
        coords, *(s_srf.getState() + tuple([vdw_radii, inner_radius, atom_coords,
                                            min_enclosure, max_peel_depth]))))

    s_inr = State(*calculator.deleteSimplices3d(coords, *(s_srf.getState() + tuple([vdw_radii, inner_radius, False]))))

    l_first_layer_simp, l_second_layer_simp = calculator.surfaceLayer(s_srf.simp, s_inr.simp, s_srf.neigh)
    s_clr = State(*calculator.deleteSection(l_first_layer_simp, *s_inr.getState()))
    LOGGER.report('Surface and inner simplices filtered in %.2fs.', '_prody_channels_surface')

    LOGGER.timeit('_prody_channels_cavities')
    c_cavities = calculator.findGroups(s_clr.neigh)
    c_surface_cavities = calculator.getSurfaceCavities(c_cavities, s_clr.simp,
                                                       l_second_layer_simp,
                                                       s_clr, mouth_oracle)

    calculator.findDeepestTetrahedra(c_surface_cavities, s_clr.neigh, s_clr.verti,
                                     coords, s_clr.simp)
    if start_point is not None:
        c_surface_cavities = calculator.setStartingTetrahedraFromPoint(
            c_surface_cavities, s_clr.verti, start_point, coords, vdw_radii,
            s_clr.simp, s_clr.neigh, start_point_search, min_depth)

    c_filtered_cavities = calculator.filterCavities(c_surface_cavities, min_depth)

    # seed_volume is a floor on where a search may start, and a cavity is a search
    # site exactly as a chamber is: a cavity holding no chamber is searched whole,
    # from a seed of its own, and one whose chambers all fall under the floor falls
    # back to the same. So the floor is applied to the cavities here as well as to
    # the chambers in setStartingTetrahedraFromChambers, and means one thing
    # throughout - no void below it seeds a search, whether it is a lobe of a large
    # cavity or a cavity entire. Both are measured on the Delaunay scale the site
    # table reports.
    #
    # Without this a single-tetrahedron sliver of tessellation debris became a
    # search site bottleneck = min(gatesof its own whenever min_depth let it through: such a sliver lies
    # wholly in the surface layer, so it is its own mouth, and it reports either
    # nothing (sealed, having no target to path to) or a one-step channel that is
    # a facet of the surface rather than a tunnel. A low min_depth is what exposes
    # them - measured at 12 of 15 sites on 1grm and 300+ of 336 on LinB with
    # min_depth=0 - but nothing except their size distinguishes them at any depth.
    #
    # Skipped when start_point is given: the user has said where the search
    # begins, and that wins over any floor.
    debris = 0
    if not cavities_only:
        calculator.calculate_cavity_volumes(c_filtered_cavities, s_clr.simp,
                                            coords)

        if start_point is None and seed_volume is not None:
            searched = calculator.filterCavitiesByVolume(
                c_filtered_cavities, min_volume=seed_volume)
            debris = len(c_filtered_cavities) - len(searched)
            if debris and not searched:
                _warn('every cavity is smaller than seed_volume={0:g} Å³, so no '
                      'search site is left and no channel can be found. Lower '
                      'seed_volume, or set it to None, to search the small '
                      'cavities too.'.format(float(seed_volume)))
            c_filtered_cavities = searched

    if debris:
        report = ('Cavities: {0} found, {1} deeper than min_depth={2:.1f} Å, {3} '
                  'of them at least seed_volume={4:g} Å³ and searched for '
                  'channels; the {5} smaller ones are tessellation debris and are '
                  'left unsearched, in %.2fs.').format(
                      len(c_surface_cavities), len(c_filtered_cavities) + debris,
                      float(min_depth), len(c_filtered_cavities),
                      float(seed_volume), debris)
    else:
        report = ('Cavities: {0} found, {1} deeper than min_depth={2:.1f} Å and '
                  '{3}, in %.2fs.').format(
                      len(c_surface_cavities), len(c_filtered_cavities),
                      float(min_depth),
                      'kept' if cavities_only else 'searched for channels')
    LOGGER.report(report, '_prody_channels_cavities')

    if cavities_only:
        if max_depth is not None:
            calculator.trimCavitiesByDepth(c_filtered_cavities, max_depth)

        if min_tetrahedra is not None or max_tetrahedra is not None:
            c_filtered_cavities = calculator.filterCavitiesByTetrahedra(
                c_filtered_cavities, min_tetrahedra, max_tetrahedra)

        calculator.calculate_cavity_volumes(c_filtered_cavities, s_clr.simp, coords)

        if min_volume is not None or max_volume is not None:
            c_filtered_cavities = calculator.filterCavitiesByVolume(
                c_filtered_cavities, min_volume, max_volume)
    # The channel path has its volumes already: they are what the seed_volume
    # floor above was applied to.

    # Largest first. Cavities come out of findGroups in connected-component
    # order, which follows tetrahedron indices, so cavity 0 was whichever void
    # happened to hold the lowest-numbered tetrahedron - an identity tag with
    # nothing behind it, and one that contradicted the start points, which are
    # numbered by volume. Sorted once here, where the volumes are final in both
    # paths, so that every consumer sees one order: the cavity files written by
    # calcSurfaceCavities, the cavity a start point names in the report, and the
    # order the cavities are searched in.
    c_filtered_cavities = calculator.orderCavitiesByVolume(c_filtered_cavities)

    merged_cavities = calculator.mergeCavities(c_filtered_cavities, s_clr.simp)

    # Early-return for the calcSurfaceCavities function:
    if cavities_only:
        LOGGER.info("Returning surface cavities")
        
        if output_path:
            output_path = Path(output_path)
            
            # A directory names no run, so nothing here is named after one: the
            # file is named after what it holds, cavities.pqr beside the
            # channels.pqr, links.pqr and pores.pqr the other entry points
            # write, and the stem is not carried into the per-cavity files,
            # which come out as cavity0.pqr.
            separate_stem = None
            if output_path.is_dir():
                output_path = output_path / "cavities.pqr"
                separate_stem = ''
            elif not (output_path.suffix == ".pdb" or output_path.suffix == ".pqr"):
                output_path = output_path.with_suffix(".pqr")

            if not separate:
                LOGGER.info("Saving surface cavities to " + str(output_path) + ".")
            else:
                LOGGER.info("Saving multiple surface cavities to directory " + str(output_path.parent) + ".")

            calculator.saveCavitiesToPdb(c_filtered_cavities, s_clr.verti,
                                         output_path, separate, separate_stem)

        LOGGER.report('Surface cavity calculation completed in %.2fs.', '_prody_calcChannels')
        return c_filtered_cavities, [coords, s_srf.simp, merged_cavities, s_clr.simp, s_clr.verti]

    LOGGER.timeit('_prody_channels_pathfinding')
    # build the weighted adjacency matrix once for the whole cleared
    # state, then run a single multi-target Dijkstra per cavity (scipy csgraph),
    # instead of one heap Dijkstra per (seed, exit) pair.
    simplices, neighbors, vertices = s_clr.getState()
    graph = calculator.buildSparseGraph(
        simplices, neighbors, vertices, coords, vdw_radii,
        gate_floor=bottleneck if prune_narrow_edges else 0.0)

    # Seed each cavity at its chambers rather than at its single deepest
    # tetrahedron. Placed after buildSparseGraph because the seed of a chamber is
    # its widest qualifying tetrahedron and the per-tetrahedron clearance that
    # measures "widest" is cached there; computing it again here would repeat the
    # same min over the whole cleared state. Skipped when start_point is given -
    # the user has said where to start, and that wins over any chamber - and when
    # seed_radius does not exceed inner_radius, where the carve cannot prune
    # anything and the cavity would come back as one chamber (which would not be
    # a no-op: it would silently move the seed from the deepest tetrahedron to
    # the widest one).
    chamber_labels = None
    if start_point is None and seed_radius > inner_radius:
        LOGGER.timeit('_prody_channels_chambers')
        # Depth of every tetrahedron lying in a cavity, gathered into one array so
        # that the carve can stop where the cavities and the channels stop.
        # -inf elsewhere: outside a cavity there is no depth to speak of, and
        # nothing out there should be chamber in any case.
        chamber_depths = np.full(len(simplices), -np.inf)
        for cavity in c_filtered_cavities:
            for tetra, depth in cavity.tetrahedra_depths.items():
                chamber_depths[int(tetra)] = depth
        chamber_labels, chamber_volumes = calculator.findChambers(
            simplices, neighbors, vertices, coords, vdw_radii, seed_radius,
            chamber_depths, min_depth)
        reseeded, chamber_notes = calculator.setStartingTetrahedraFromChambers(
            c_filtered_cavities, chamber_labels, chamber_volumes, min_depth,
            seed_volume, max_seeds)
        # Heading, then the per-cavity detail under it, then the count of what
        # the search will actually run over. The detail is logged here rather
        # than where it is produced so that it never precedes its own heading.
        with_chambers = sum(1 for cavity in c_filtered_cavities
                            if cavity.chambers_found)
        whole = len(c_filtered_cavities) - reseeded
        LOGGER.info('Chambers (probe {0:.2f} Å): {1} of the {2} searched '
            'cavities have them{3}.'.format(
                seed_radius, with_chambers, len(c_filtered_cavities),
                '' if not whole else
                '; the other {0} {1} searched whole'.format(
                    whole, 'is' if whole == 1 else 'are')))
        for note in chamber_notes:
            LOGGER.info('    ' + note)
        LOGGER.report('{0} search sites (sp) in %.2fs: one per seeded chamber, '
            'one per cavity searched whole.'.format(
                sum(len(cavity.starting_tetrahedron)
                    for cavity in c_filtered_cavities)),
            '_prody_channels_chambers')

    for cavity in c_filtered_cavities:
        calculator.dijkstra(cavity, graph, simplices, neighbors, vertices,
                            coords, vdw_radii, route_divergence,
                            chamber_labels if chamber_links else None)
    # Sites, not cavities: one Dijkstra runs per site, so the number of sites is
    # what the time divides into.
    LOGGER.report('Channel search (Dijkstra) over {0} search sites in {1} '
        'cavities completed in %.2fs.'.format(
            sum(len(cavity.starting_tetrahedron)
                for cavity in c_filtered_cavities),
            len(c_filtered_cavities)), '_prody_channels_pathfinding')

    calculator.filterChannelsByBottleneck(c_filtered_cavities, bottleneck)
    
    if min_volume is not None or max_volume is not None:
        calculator.filterChannelsByVolume(c_filtered_cavities, min_volume, 
                                          max_volume)
    
    channels = [channel for cavity in c_filtered_cavities for channel in cavity.channels]
    # Order channels by ascending Dijkstra cost so that channel 0 is the best
    # tunnel (a short path through wide tetrahedra). This ordering drives both
    # the returned list and the channel numbering in the saved PQR/PDB files.
    channels.sort(key=lambda ch: ch.cost if ch.cost is not None else float('inf'))

    # Links are filtered by the same width rule as the channels - a neck the
    # traversal probe would not fit through is not a way from one chamber to the
    # next - but not by min_volume, which is a size floor written for objects that
    # reach the surface; a link is by construction only the neck.
    links = [link for cavity in c_filtered_cavities for link in cavity.links
             if link.bottleneck >= bottleneck]
    links.sort(key=lambda ln: ln.cost if ln.cost is not None else float('inf'))

    # Start points, numbered over the cavities in the order they were searched:
    # one per chamber where the cavity was seeded per chamber, one for the whole
    # cavity where no chamber qualified. Every object carries the index of the one
    # it was traced from, which is all the labelling needed - the voids nest
    # strictly, a chamber never spanning two cavities, so a start point names one
    # void unambiguously and the objects sharing it are the ways out of it.
    # Volumes for the report: a chamber's own where the cavity was seeded per
    # chamber, the whole cavity's where it was not. Both on the Delaunay scale,
    # which is not the swept-sphere scale a channel volume is on.
    #
    # Numbered largest first, over all the start points at once rather than
    # within each cavity, so that sp0 is the biggest void the search ran from
    # and the numbering is the order the sites are worth looking at. The two
    # kinds are ranked against each other: each row reports the size of the void
    # its start point actually names, the whole cavity where the cavity is the
    # site and the one lobe where a lobe is. Ties are broken by cavity and then
    # by seed, so the numbering is reproducible.
    seed_rows = []
    for cavity_index, cavity in enumerate(c_filtered_cavities):
        seeds = [int(seed) for seed in np.asarray(cavity.starting_tetrahedron)]
        for seed in seeds:
            seed_rows.append((cavity_index, cavity, seed, len(seeds),
                              cavity.tetrahedra_depths.get(seed, 0.0),
                              cavity.seed_volumes.get(seed, cavity.volume)))
    seed_rows.sort(key=lambda row: (-row[5], row[0], row[2]))

    origins = {}
    origin_rows = []
    ordinals = {}
    for cavity_index, cavity, seed, count, depth, volume in seed_rows:
        origins[seed] = len(origins)
        # The chamber's rank inside its own cavity, taken from this same order,
        # so that "chamber 1 of 20 in cavity 0" is that cavity's largest and
        # never contradicts the start point numbers.
        ordinal = ordinals.get(cavity_index, 0)
        ordinals[cavity_index] = ordinal + 1
        origin_rows.append((cavity_index,
                            ordinal if cavity.seed_chambers else None,
                            count, depth, volume, cavity.chambers_found))
    for entry in channels + links:
        entry.origin = origins.get(int(np.asarray(entry.tetrahedra)[0]))

    # An sp<n> tag tells one search site's objects from another's, so a run that
    # searched from a single site has nothing to tell: every object would carry
    # the same sp0 in its file name and its REMARK. The origins themselves are
    # kept - they are what the table below counts by - and only the naming is
    # dropped, so the files come out as chl0.pqr with a REMARK reading
    # "channel 0" alone, as they did before multi-seed searching.
    name_sites = len(origins) > 1

    # The far end of a link, named the way its near end is. Chamber labels are
    # internal to the carve, so they are translated here, where the start points
    # they stand for have just been numbered.
    for cavity in c_filtered_cavities:
        reached = {label: origins[int(seed)]
                   for seed, label in cavity.seed_chambers.items()}
        for link in cavity.links:
            link.destination = reached.get(link.joined_chamber)

    LOGGER.info("Found {0} channel{1}{2}.".format(
        len(channels), '' if len(channels) == 1 else 's',
        '' if not links else
        " and {0} link{1} (a link joins a deep chamber to a shallower one and "
        "never reaches the surface)".format(
            len(links), '' if len(links) == 1 else 's')))

    sealed = 0
    if origin_rows:
        # One row per site, in the order the sites are numbered. Laid out as a
        # table because every row says the same six things: written as sentences
        # they came out in three different grammars ("chamber 1 of 20 seeded in
        # cavity 0 (of 41 found)", "the only chamber of cavity 1", "cavity 2
        # whole (no chamber found)"), which read as three unrelated remarks
        # rather than as one list. How many chambers a cavity has and how many
        # of them were seeded is said once, in the chamber lines above, instead
        # of being repeated on every row that belongs to that cavity.
        rows = []
        for index, (cavity_index, chamber, seeded, depth,
                    volume, found) in enumerate(origin_rows):
            if chamber is None:
                where = 'cavity {0}, whole'.format(cavity_index)
            else:
                where = 'cavity {0}, chamber {1}/{2}'.format(
                    cavity_index, chamber + 1, seeded)
            traced = [entry for entry in channels if entry.origin == index]
            linked = [entry for entry in links if entry.origin == index]
            # Where the links go, so that the connections can be followed from
            # the report without opening the files. Named only on the deep end,
            # which is the end a link is reported from.
            note = '  -> {0}'.format(', '.join(
                'sp{0}'.format(entry.destination) for entry in linked
                if entry.destination is not None)) if linked else ''
            if not (traced or linked):
                # Why the site reports nothing; the bottleneck it failed and
                # what to do about it are in the summary line below the table.
                # Counted here rather than recomputed there, so that the number
                # in that line is the number of rows marked here by definition.
                note = '  sealed'
                sealed += 1
            # The seed's own Voronoi vertex, written so that it can be pasted back
            # as start_point to search this one site again - the automatic passes
            # otherwise report which void they started from but never where, and
            # the placement cannot be reproduced or adjusted without it.
            seed_vertex = vertices[seed_rows[index][2]]
            rows.append(['sp{0}'.format(index),
                         '[{0:.3f}, {1:.3f}, {2:.3f}]'.format(*seed_vertex),
                         where, '{0:.0f}'.format(volume),
                         '{0:.1f}'.format(depth),
                         str(len(traced)) if traced else '-',
                         str(len(linked)) if linked else '-', note])

        header = ['site', 'start_point [Å]', 'void', 'volume [Å³]', 'depth [Å]',
                  'channels', 'links']

        # A single search site names nothing with an sp<n> - every row would be
        # sp0, and nothing written carries the tag either - and it can hold no
        # link, which needs a second site to arrive at. Both columns are dropped
        # there, leaving the one row to say where the search began, what the void
        # is and what came out of it. The notes are kept out of the widths, being
        # ragged text at the end of the row rather than a column of it, so the
        # coordinates sit ahead of them rather than being run into by an arrow.
        keep = range(len(header)) if name_sites else [1, 2, 3, 4, 5]
        notes = [row[len(header)] for row in rows]
        rows = [[row[column] for column in keep] for row in rows]
        header = [header[column] for column in keep]
        # site, start_point and void are text and read left-aligned; the counts
        # and measurements after them are compared down the page and go right.
        labels = 3 if name_sites else 2
        widths = [max(len(row[column]) for row in [header] + rows)
                  for column in range(len(header))]

        def formatRow(row):
            # The labels left, the numbers right, so that a column of volumes or
            # counts can be compared by eye down the page.
            return '    ' + '  '.join(
                text.ljust(width) if column < labels else text.rjust(width)
                for column, (text, width) in enumerate(zip(row, widths)))

        LOGGER.info("Search sites (sp), the void each search ran from, largest "
                    "first; sp<n> tags every channel, link and output file:"
                    if name_sites else
                    "The void the search ran from:")
        LOGGER.info(formatRow(header))
        for row, note in zip(rows, notes):
            LOGGER.info(formatRow(row) + note)
        LOGGER.info("    (site volumes measure the void itself and are not on "
                    "the swept-sphere scale of the channel volumes)")

    # A site can end up reporting nothing at all: its own channels dropped by the
    # dedup as duplicates of a shallower site's shorter ones, and its link then
    # dropped by the bottleneck filter because the neck it would leave through is
    # too tight for the probe. That is a real finding - the void is sealed at this
    # width - but an invisible one, because nothing is written about a site that
    # reports no object, and the void simply goes missing from the output.
    #
    # A cavity searched whole is a site on the same terms as a chamber, so it is
    # counted here on the same terms: the count is the number of rows the table
    # marked sealed, whichever kind of void they name. Counting only the chamber
    # seeds (cavity.seed_chambers) made the summary contradict the table it
    # summarises - "The 1 site marked sealed above" under a table marking
    # hundreds.
    if sealed:
        LOGGER.info("The {0} site{1} marked sealed above report neither a "
                    "channel nor a link: no route out of them survived - either "
                    "narrower than bottleneck={2:.2f} Å, or dropped as a "
                    "duplicate of a shallower site's, or the void is its own "
                    "mouth and has nowhere to path to. Lower bottleneck to see "
                    "how the first kind connect.".format(
                        sealed, '' if sealed == 1 else 's', bottleneck))

    if output_path and not (channels or links):
        # Nothing found, so nothing is written. An empty file would say only
        # that a run happened, which the count reported above already says, and
        # it cannot be told apart from a run that failed while writing. What an
        # earlier run left behind is worth a word, though: it survives now, and
        # goes on looking like this run's output.
        _warnStaleOutputs(output_path, output_format, separate)

    elif output_path and _isMmcifFormat(output_format, separate):
        writeChannelsCIF(output_path, channels, atoms, links=links,
                         auto=start_point is None)

    elif output_path:
        output_path, links_path, _, separate_stem = \
            _pqrOutputPaths(output_path)

        # One line for the whole of what was written, so that the reader sees
        # where the channels went and where the links went in one place.
        if not separate:
            LOGGER.info("Saving {0} channels to {1}{2}.".format(
                len(channels), output_path,
                '' if not links else
                " and {0} links to {1}".format(len(links), links_path)))
        else:
            LOGGER.info("Saving {0} channels{1} to directory {2}, one file per "
                        "object named {3}chl<n>{4}.".format(
                            len(channels),
                            '' if not links else
                            " and {0} links".format(len(links)),
                            output_path.parent,
                            'sp<site>_' if name_sites else '',
                            '' if not links else
                            " and sp<site>_lnk<n>" if name_sites
                            else " and lnk<n>"))
        calculator.saveChannelsToPdb(channels, output_path, separate,
                                     separate_stem=separate_stem,
                                     name_sites=name_sites)
        if links:
            calculator.saveChannelsToPdb(links, links_path, separate,
                                         tag='lnk', label='link',
                                         separate_path=output_path,
                                         separate_stem=separate_stem,
                                         name_sites=name_sites)
    else:
        LOGGER.info("No output path given.")

    LOGGER.report('Channel calculation completed in %.2fs.', '_prody_calcChannels')

    # Additional information can be obtained
    if return_details:
        details = {'calculator': calculator,
                    'simplices': s_clr.simp,
                    'neighbors': s_clr.neigh,
                    'vertices': s_clr.verti,
                    'coords': coords,
                    'vdw_radii': vdw_radii,
                    # Chamber links ride in details rather than in the returned
                    # tuple, so that `channels, surface = calcChannels(...)` keeps
                    # working. Empty unless the cavities were seeded per chamber.
                    'links': links}
        
        return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp], details

    return channels, [coords, s_srf.simp, merged_cavities, s_clr.simp]


def calcPoresFromChannels(channels, details, min_end_to_end=None, max_end_to_end=None,
    min_bottleneck=None, max_bottleneck=None, min_length=None, max_length=None,
    min_volume=None, max_volume=None, output_path=None, separate=False,
    output_format='pqr', atoms=None):
    """Construct potential pores from previously identified channels using 
    :func:`calcChannels`. This function performs a post-processing analysis of 
    channels and requires ``return_details`` set to ``True`` in :func:`calcChannels`.
    The `separate` parameter controls whether each pore is additionally saved to a 
    separate file.
    
    The pore-construction procedure consists of the following steps:

    1. Group channels according to their starting tetrahedron.
    2. Generate all unique pairs of channels within each group.
    3. Identify the common initial segment and the last tetrahedron shared by
       each pair of channel paths.
    4. Join the non-overlapping parts of the two channels at their branching
       tetrahedron to obtain a surface-to-surface path.
    5. Reject paths containing loops or discontinuities between neighboring
       tetrahedra.
    6. Remove identical paths and paths differing only in direction.
    7. Recalculate the centerline spline, radius profile, length, bottleneck,
       and volume of each resulting pore using approach implemented for channels
       identification and visualization.
    8. Pores are filtered based on the given criteria (``min_end_to_end``, 
       ``max_end_to_end``, ``min_bottleneck``, ``max_bottleneck``, ``min_length``, 
       ``max_length``, ``min_volume``, ``max_volume``). 
    
    :arg channels: A list of channel objects or a single channel object. Each
        channel should have a `getSplines()` method that returns two
        interpolators over one parameter domain: one for the centerline and one
        for the radii.
    :type channels: list or single channel object

    :arg details: Additional calculation data returned by
        :func:`calcChannels` with ``return_details=True``. The dictionary must
        contain ``calculator``, ``simplices``, ``neighbors``, ``vertices``,
        ``coords``, and ``vdw_radii``.
    :type details: dict

    :arg min_end_to_end: Minimum allowed distance between the two pore
        openings. Pores with a smaller end-to-end distance will be excluded.
        Default is None.
    :type min_end_to_end: int, float

    :arg max_end_to_end: Maximum allowed distance between the two pore
        openings. Pores with a larger end-to-end distance will be excluded.
        Default is None.
    :type max_end_to_end: int, float

    :arg min_bottleneck: Minimum allowed bottleneck radius of the pore.
        Pores with a smaller bottleneck will be excluded.
        Default is None.
    :type min_bottleneck: int, float

    :arg max_bottleneck: Maximum allowed bottleneck radius of the pore.
        Pores with a larger bottleneck will be excluded.
        Default is None.
    :type max_bottleneck: int, float

    :arg min_length: Minimum allowed length of the pore. Pores shorter than this 
        value will be excluded. Default is None.
    :type min_length: int, float

    :arg max_length: Maximum allowed length of the pore. Pores longer than this 
        value will be excluded. Default is None.
    :type max_length: int, float

    :arg min_volume: Minimum allowed volume of the pore. Pores with a smaller volume 
        will be excluded. The value is given in cubic Angstroms. Default is None.
    :type min_volume: int, float

    :arg max_volume: Maximum allowed volume of the pore. Pores with a larger volume 
        will be excluded. The value is given in cubic Angstroms. Default is None.
    :type max_volume: int, float

    :arg output_path: Optional path to save the resulting pores and
        associated data in PQR (or PDB) format. If None, results are not saved.
        Default is None.
    :type output_path: str or None

    :arg output_format: What ``output_path`` is written as; ``"pqr"`` (the
        default) or ``"mmcif"``, as in :func:`calcChannels`. A directory takes
        ``pores.cif``, named apart from ``channels.cif`` so that a run writing
        both into one folder does not have the second overwrite the first.
        ``separate`` is ignored for mmCIF.
    :type output_format: str

    :arg atoms: Structure to measure the pore lining against, used only for
        ``output_format="mmcif"``. This function is given channels and their
        tessellation rather than a structure, so without it the mmCIF holds the
        geometry alone - every category that describes what surrounds the route
        is left out rather than guessed at.
    :type atoms: :class:`.Atomic` or None

    :returns: Potential pores constructed from compatible channel pairs.
    :rtype: list of Channel
    
    Usage:
    channels, surface, details = calcChannels(protein, return_details=True)
    pores = calcPoresFromChannels(channels, details, output_path='pores', separate=True)   
    """

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path
    
    calculator = details['calculator']
    simplices = details['simplices']
    vertices = details['vertices']
    coords = details['coords']
    vdw_radii = details['vdw_radii']
    neighbors = details['neighbors']
    
    pores = []
    pore_paths = []
    seen_paths = set()
    channel_groups = {}
    
    # Group channels using their starting tetrahedron
    for channel_index, channel in enumerate(channels):
        path = np.asarray(channel.tetrahedra, dtype=np.intp)
        # If channel is smaller than two tetrahedra (probably very rare)
        # Those channels should be excluded because they can not be connected with others
        if len(path) < 2:
            continue

        start_tetrahedron = int(path[0])
        channel_groups.setdefault(start_tetrahedron, []).append((channel_index, channel, path))
        
    from itertools import combinations
    # Generate all channel pairs within each group
    for start_tetrahedron, group in channel_groups.items():
        if len(group) < 2:
            continue

        for (channel1_index, channel1, path1), (channel2_index, channel2, path2) in combinations(group, 2):
            common_length = 0

            for tetrahedron1, tetrahedron2 in zip(path1, path2):
                if tetrahedron1 != tetrahedron2:
                    break
                common_length += 1

            if common_length == 0:
                continue

            # If we have for example: path1: start → A → B → C → mouth 1 and path2: start → A → B → D → mouth 2
            # it will create mouth 1 → C → B → D → mouth 2
            branch_index = common_length - 1
            pore_path = np.concatenate((path1[branch_index:][::-1], path2[branch_index + 1:]))
            
            # Reject paths containing loops
            if len(np.unique(pore_path)) != len(pore_path):
                continue

            # Continulity check of the pores (neighbours)
            is_continuous = True
            for tetrahedron1, tetrahedron2 in zip(pore_path[:-1], pore_path[1:]):
                if tetrahedron2 not in neighbors[tetrahedron1]:
                    is_continuous = False
                    break
            if not is_continuous:
                continue
                
            # Remove identical paths            
            path_key = tuple(int(tetrahedron) for tetrahedron in pore_path)
            canonical_key = min(path_key, path_key[::-1])

            if canonical_key in seen_paths:
                continue

            seen_paths.add(canonical_key)
            pore_paths.append(pore_path)
    
    # Pores reconstruction
    for pore_path in pore_paths:
        # Filters - Distance between two ends of the pore
        end_to_end = np.linalg.norm(vertices[pore_path[0]] - vertices[pore_path[-1]])
        
        if min_end_to_end is not None and end_to_end < min_end_to_end:
            continue
        if max_end_to_end is not None and end_to_end > max_end_to_end:
            continue
    
        centerline_spline, radius_spline, length, bottleneck, volume, centers, radii, gates = calculator.processChannel(
                                                pore_path, vertices, coords, vdw_radii, simplices)
        
        # Filters - bottleneck, length, volume
        if min_bottleneck is not None and bottleneck < min_bottleneck:
            continue
        if max_bottleneck is not None and bottleneck > max_bottleneck:
            continue
            
        if min_length is not None and length < min_length:
            continue
        if max_length is not None and length > max_length:
            continue

        if min_volume is not None and volume < min_volume:
            continue
        if max_volume is not None and volume > max_volume:
            continue
        
        pore = Channel(pore_path, centerline_spline, 
                        radius_spline, length, bottleneck,
                         volume, centers, radii, gates, 0.0)
        pores.append(pore)
    
    if output_path and not pores:
        # As in calcChannels: no pores, no file - and a word about whatever an
        # earlier run left at the same place. Pores are named apart from
        # channels, so the tag globbed for and the mmCIF name are the pore ones.
        # Through both resolvers, in the order the writing path applies them:
        # _poreCifPath picks the name inside a directory, and writeChannelsCIF
        # then adds the suffix. Checking only the first would look for `out`
        # where `out.cif` was written.
        _warnStaleOutputs(output_path, output_format, separate, tags=('pore',),
                          cif_path=_cifOutputPath(_poreCifPath(output_path)))

    elif output_path and _isMmcifFormat(output_format, separate):
        # A directory takes pores.cif rather than channels.cif, for the same
        # reason the PQR path names them apart: a run writing both into one
        # folder would otherwise have the second overwrite the first.
        writeChannelsCIF(_poreCifPath(output_path), pores, atoms,
                         object_type='pore')

    elif output_path:
        output_path = Path(output_path)
        # As in calcChannels: a directory names no run, so its placeholder file
        # name is kept out of the per-pore ones.
        separate_stem = None
        if output_path.is_dir():
            output_path = output_path / "pores.pqr"
            separate_stem = ''
        elif output_path.suffix not in (".pdb", ".pqr"):
            output_path = output_path.with_suffix(".pqr")

        # Named as pores, not as channels. The writer's defaults were the only
        # thing naming them, so a pore came out as <stem>_chl0.pqr with a REMARK
        # reading "channel 0" - harmless while the stem told them apart, and a
        # collision once a run told only a directory drops the stem and writes
        # both into one folder.
        calculator.saveChannelsToPdb(pores, output_path, separate=separate,
                                     tag='pore', label='pore',
                                     separate_stem=separate_stem)

    return pores
            
def connectChannelsToSurfaceCavities(channels, channel_details, cavities,
    cavity_surface, tolerance=1.0, min_contact_points=2, cavity_margin=2.0,
    output_path=None, separate=False):
    """Connect independently calculated channels with surface cavities.

    Channels and surface cavities may be calculated using different parameters.
    For every connected channel-cavity pair, the overlapping surface-end part
    of the channel is removed and only the local connected region of the
    surface cavity associated with that channel is retained.

    :arg channels: Channels returned by :func:`calcChannels`.
    :type channels: list

    :arg channel_details: Additional data returned by
        ``calcChannels(..., return_details=True)``.
    :type channel_details: dict

    :arg cavities: Surface cavities returned by
        :func:`calcSurfaceCavities`.
    :type cavities: list

    :arg cavity_surface: Surface information returned by
        :func:`calcSurfaceCavities`.
    :type cavity_surface: list

    :arg tolerance: Maximum distance in Angstrom between channel and
        surface-cavity Voronoi vertices used to define direct contact.
        Default is 1.0.
    :type tolerance: float

    :arg min_contact_points: Minimum number of consecutive channel
        tetrahedra required to identify a channel-cavity contact.
        Default is 2.
    :type min_contact_points: int

    :arg cavity_margin: Additional distance in Angstrom beyond the local
        channel radius used when selecting the cavity region associated
        with the channel. Default is 2.0.
    :type cavity_margin: float

    :arg output_path: Optional path for the resulting PQR/PDB file.
    :type output_path: str or pathlib.Path or None

    :arg separate: If True, each connected cavity-channel pair is also
        saved to a separate file. Default is False.
    :type separate: bool

    :returns: Connected channel-cavity results.
    :rtype: list of dict 
    
    Example:
    protein = parsePDB('1tqn').select('protein')

    channels, channel_surface, channel_details = calcChannels(protein, inner_radius=0.8,
        min_depth=3, return_details=True, output_path='channels', separate=True)

    cavities, cavity_surface = calcSurfaceCavities(protein, surf_radius=3.8,
        inner_radius=1.1, min_depth=5, min_volume=500,
        output_path='surface_cavities', separate=True)

    connected = connectChannelsToSurfaceCavities(channels, channel_details,
        cavities, cavity_surface,
        tolerance=1.0, min_contact_points=3, cavity_margin=4.0,
        output_path='connected_cavities_channels.pqr', separate=True) """

    if tolerance <= 0:
        raise ValueError("tolerance must be greater than zero")

    if min_contact_points < 1:
        raise ValueError("min_contact_points must be at least 1")

    if cavity_margin is not None and cavity_margin < 0:
        raise ValueError("cavity_margin must be non-negative or None")

    if not isinstance(channel_details, dict):
        raise TypeError("channel_details must be returned by "
            "calcChannels(..., return_details=True)")

    required = ('calculator', 'simplices', 'vertices', 'coords', 'vdw_radii')
    missing = [key for key in required if key not in channel_details]

    if missing:
        raise ValueError("channel_details is missing: {0}".format(", ".join(missing)))

    if cavity_surface is None or len(cavity_surface) < 5:
        raise ValueError("cavity_surface must be returned by calcSurfaceCavities()")

    calculator = channel_details['calculator']
    simplices = channel_details['simplices']
    vertices = channel_details['vertices']
    coords = channel_details['coords']
    vdw_radii = channel_details['vdw_radii']

    cavity_vertices = cavity_surface[4]
    connected = []
    cavity_trees = []

    for cavity in cavities:
        tetrahedra = np.asarray(cavity.tetrahedra, dtype=np.intp)

        if len(tetrahedra) == 0:
            cavity_trees.append(None)
            continue

        cavity_trees.append(_kdTree(cavity_vertices[tetrahedra]))

    for channel_index, channel in enumerate(channels):
        path = np.asarray(channel.tetrahedra, dtype=np.intp)

        if len(path) < 2:
            continue

        for cavity_index, (cavity, cavity_tree) in enumerate(zip(cavities, cavity_trees)):

            if cavity_tree is None:
                continue

            channel_xyz = vertices[path]
            distances, _ = cavity_tree.query(channel_xyz)

            # Orient the path from protein interior toward the cavity.
            # A channel returned by calcChannels normally already follows
            # seed -> surface, but this keeps the post-processing independent
            # of path orientation.
            if distances[0] < distances[-1]:
                oriented_path = path[::-1].copy()
                oriented_distances = distances[::-1].copy()
            else:
                oriented_path = path.copy()
                oriented_distances = distances.copy()

            contact_mask = oriented_distances <= tolerance
            connection_index = _findContactRun(contact_mask, min_contact_points)

            if connection_index is None:
                continue

            # At least two vertices are required by CubicSpline.
            if connection_index < 1:
                continue

            # Retain the first contact vertex so the reconstructed channel
            # terminates directly at the cavity.
            trimmed_path = oriented_path[:connection_index + 1]

            if len(trimmed_path) < 2:
                continue

            centerline_spline, radius_spline, length, bottleneck, volume, centers, radii, gates = \
                calculator.processChannel(trimmed_path, vertices, coords, vdw_radii, simplices)

            trimmed_channel = Channel(trimmed_path, centerline_spline, radius_spline,
                                        length, bottleneck, volume,
                                         centers, radii, gates, cost=None)

            trimmed_channel.origin = getattr(channel, 'origin', None)
            trimmed_channel.destination = getattr(channel, 'destination', None)

            connection_point = vertices[trimmed_path[-1]].copy()

            local_cavity_tetrahedra = _selectLocalSurfaceCavity(cavity, cavity_surface,
                    trimmed_channel, tolerance=tolerance, cavity_margin=cavity_margin)

            if len(local_cavity_tetrahedra) == 0:
                continue

            connected.append({'cavity_index': cavity_index, 
                'channel_index': channel_index,
                'cavity': cavity, 'channel': channel, 
                'trimmed_channel': trimmed_channel,
                'cavity_tetrahedra': local_cavity_tetrahedra, 
                'connection_point': connection_point,
                'minimum_distance': float(np.min(distances))})

    LOGGER.info(
        "Detected {0} connected surface cavity-channel pair(s).".format(len(connected)))

    connected_cavities = set()
    connected_channels = set()

    if connected:
        LOGGER.info("Connected surface cavities and channels:")

        for result in connected:
            cavity_index = result['cavity_index']
            channel_index = result['channel_index']
            channel = result['channel']
            origin = getattr(channel, 'origin', None)

            connected_cavities.add(cavity_index)
            connected_channels.add(channel_index)

            if origin is None:
                channel_label = "channel {0}".format(channel_index)
            else:
                channel_label = "channel {0} (sp{1})".format(channel_index, origin)

            LOGGER.info("    cavity {0} <-> {1}, minimum distance {2:.2f} A, "
                    "local cavity {3}/{4} tetrahedra".format(
                    cavity_index, 
                    channel_label, 
                    result['minimum_distance'],
                    len(result['cavity_tetrahedra']), 
                    len(result['cavity'].tetrahedra)))

    unconnected_cavities = [i for i in range(len(cavities))
        if i not in connected_cavities]

    unconnected_channels = [i for i in range(len(channels))
        if i not in connected_channels]

    if unconnected_cavities:
        LOGGER.info("Surface cavities without connected channels: {0}.".format(
                ", ".join("cavity {0}".format(i)
                    for i in unconnected_cavities)))
    else:
        LOGGER.info("All surface cavities have at least one connected channel.")

    if unconnected_channels:
        channel_labels = []

        for i in unconnected_channels:
            origin = getattr(channels[i], 'origin', None)

            if origin is None:
                channel_labels.append("channel {0}".format(i))
            else:
                channel_labels.append("channel {0} (sp{1})".format(i, origin))

        LOGGER.info("Channels without connected surface cavities: {0}.".format(", ".join(channel_labels)))
    
    else:
        LOGGER.info("All channels have at least one connected surface cavity.")

    if output_path is not None:
        _saveConnectedCavityChannels(connected, cavity_surface, output_path, separate=separate)

    return connected

def calcChannelsMultipleAtomGroups(atomgroups, output_path=None, 
    separate=False, max_proc=2,start_point = None, mp_context=None,filenames = None,**kwargs):
    """Calculate channels for multiple atom groups."""
    import multiprocessing
    if PY3K:
        if not checkAndImport('pathlib'):
            errorMsg = 'To run showChannels, please install open3d.'
            raise ImportError(errorMsg)
                
        from pathlib import Path
    else:
        if not checkAndImport('pathlib2'):
            errorMsg = 'To run showChannels, please install pathlib2 for Python 2.7.'
            raise ImportError(errorMsg)
        
        from pathlib2 import Path

    return_details = kwargs.pop('return_details', False)

    for atoms in atomgroups:
        _requireCoords(atoms)

    channels_all = []
    surfaces_all = []
    details_all = []
    tasks = []

    from multiprocessing import Pool, cpu_count
    filenames = kwargs.pop('filenames',None)
    for i, ag in enumerate(atomgroups):
        if output_path:
            if filenames:
                frame_output_path = _frameOutputPath(output_path, filenames[i], "channels")
            else:
                frame_output_path = _frameOutputPath(output_path, i, "channels")
        else:
            frame_output_path = None
        tasks.append((i, ag, np.array(ag.getCoords(), copy=True),
                frame_output_path, separate, start_point, return_details, kwargs))

    if len(tasks) == 0:
        if return_details:
            return channels_all, surfaces_all, details_all
        return channels_all, surfaces_all
    
    max_proc = kwargs.pop('max_proc',None)
    
    if max_proc is None:
        max_proc = max(1, cpu_count()//2)
    chunksize = max(1, len(atomgroups)//(max_proc*4))
    if max_proc == 1:
        results = [_calcChannelsMultipleFramesWorker(task)
                    for task in tasks]
    else:
        if mp_context is None:
            ctx = multiprocessing.get_context()
        else:
            ctx = multiprocessing.get_context(mp_context)
        
        with ctx.Pool(processes=max_proc) as pool:
            results_unordered = pool.imap_unordered(_calcChannelsMultipleFramesWorker, tasks,chunksize=chunksize )
            results = [result for result in results_unordered]
            results.sort(key=lambda x: x[0])

    failed_idx = []
    for result in results:
        if return_details:
            idx, channels, surfaces, details = result
            details_all.append(details)
        else:
            idx, channels, surfaces = result
        # if channels is None and surfaces is None:
        #     LOGGER.warning(f"WARNING: {idx} failed or No Channels Detected")
        #     failed_idx.append(idx)
        channels_all.append(channels)
        surfaces_all.append(surfaces)

    for channels in channels_all:
        for channel in channels:
            channel._buildSplines()

    if return_details:
        return channels_all, surfaces_all, details_all

    return channels_all, surfaces_all

              
def calcChannelsMultipleFrames(atoms, trajectory=None, output_path=None, 
    separate=False, start_point=None, max_proc=2, mp_context=None, **kwargs):
    """Compute channels for each frame in a given trajectory or multi-model 
    PDB file.

    This function calculates the channels for each frame in a trajectory or for
     each model in a multi-model PDB file. The `kwargs` can include parameters 
     necessary for channel calculation. If the `separate` parameter is set to 
     True, each detected channel will be saved in a separate PDB file.

    :arg atoms: Atomic data or object containing atomic coordinates and methods 
        for accessing them.
    :type atoms: object

    :arg trajectory: Trajectory object containing multiple frames or a 
        multi-model PDB file.
    :type trajectory: Atomic or Ensemble object

    :arg output_path: Optional path to save the resulting channels and 
        associated data in PDB format. If a directory is specified, each 
        frame/model will have its results saved in separate files. If None, 
        results are not saved. Default is None.
    :type output_path: str or None

    :arg separate: If True, each detected channel is saved to a separate PDB 
        file for each frame/model.
        If False, all channels for each frame/model are saved in a single file. 
        Default is False.
    :type separate: bool

    :arg start_point: Optional starting point for channel search, applied to every
        frame. If provided, the search is restricted to the cavity holding the
        tetrahedron nearest the point and is seeded there, overriding the default
        automatic seed selection; see :func:`calcChannels` for how the seed is
        picked and for ``start_point_search``, which bounds how far from the point
        it may lie. Coordinates must be given in Å.
    :type start_point: list, tuple, or ndarray (length 3), or None

    :arg max_proc: Maximum number of parallel processes used for calculation. 
        If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg mp_context: Multiprocessing start method used for parallel pore
        calculations. If `None`, the default method for the operating system
        is used. Windows and macOS use the ``'spawn'`` method by default,
        whereas Linux typically uses ``'fork'``. Setting
        ``mp_context='spawn'`` can be potentially used on Linux, but might 
        be slower. Available values may include ``'spawn'``, ``'fork'``,
        and ``'forkserver'``, depending on the operating system. 
        Default is `None`.
    :type mp_context: str or None

    :arg kwargs: Additional parameters required for channel calculation. This can 
        include parameters such as radius values (surf_radius, inner_radius), minimum depth (min_depth), 
        bottleneck values, etc. 
        See the available parameters in calcChannels().
    :type kwargs: dict

    :returns: List of channels and surfaces computed for each frame or model. 
        Each entry in the list corresponds to a specific frame or model.
    :rtype: list of lists

    Example usage:
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, 
                                    output_path="channels.pdb", separate=False, surf_radius=15, 
                                    inner_radius=1.2, min_depth=5, bottleneck=1, sparsity=6)
                                  
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms, trajectory=traj, 
                                    output_path="channels.pdb", separate=False, 
                                    start_point=[-10.353, -0.133, 5.608]) """
    
    if PY3K:
        if not checkAndImport('pathlib'):
            errorMsg = 'To run showChannels, please install open3d.'
            raise ImportError(errorMsg)
                
        from pathlib import Path
    else:
        if not checkAndImport('pathlib2'):
            errorMsg = 'To run showChannels, please install pathlib2 for Python 2.7.'
            raise ImportError(errorMsg)
        
        from pathlib2 import Path
        
    _requireCoords(atoms)

    channels_all = []
    surfaces_all = []
    details_all = []
    tasks = []

    return_details = kwargs.pop('return_details', False)
    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    # Read rather than popped: the worker forwards the rest of kwargs to
    # calcChannels, which is where the format takes effect. It is wanted here only
    # to name the per-frame files. The schema describes one structure and has no
    # frame of its own, so a frame per file is what keeps each written file
    # something the schema can describe - the same arrangement the PQR path uses.
    frame_suffix = '.cif' if _isMmcifFormat(
        kwargs.get('output_format', 'pqr'), separate) else '.pqr'
    if output_path:
        output_path = Path(output_path)
        if output_path.suffix in ('.pqr', '.cif'):
            output_path = output_path.with_suffix('')

    if trajectory is not None:
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)

        # `_nfi` is a DCD read cursor, so only a file-backed trajectory has one;
        # an Ensemble holds its coordinates in memory and has nothing to rewind.
        # Reading it unguarded made every in-memory input raise, the conversion
        # just above included - the one branch written to accept an Atomic turned
        # it into the very type that could not survive the next line.
        nfi = getattr(trajectory, '_nfi', None)
        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        first, last = _frameBounds(None, start_frame, stop_frame)
        traj = trajectory[first:last]

        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=first):
            if output_path:
                frame_output_path = _frameOutputPath(output_path, j0, "channels",
                                                     frame_suffix)
            else:
                frame_output_path = None
            
            tasks.append((j0, atoms_copy, np.array(frame0.getCoords(), copy=True),
                            frame_output_path, separate, start_point, return_details, kwargs))
        if nfi is not None:
            trajectory._nfi = nfi

    else:
        if atoms.numCoordsets() > 1:
            coordsets = atoms.getCoordsets()
            first, last = _frameBounds(len(coordsets), start_frame, stop_frame)
            for model_nr in range(first, last):

                if output_path:
                    frame_output_path = _frameOutputPath(output_path, model_nr,
                                                         "channels",
                                                         frame_suffix)
                else:
                    frame_output_path = None
                
                tasks.append((model_nr, atoms, np.array(coordsets[model_nr], copy=True),
                                frame_output_path, separate, start_point, return_details, kwargs))
                
        else:
            LOGGER.info("Include trajectory or use multi-model PDB file.")

    import multiprocessing

    if len(tasks) == 0:
        if return_details:
            return channels_all, surfaces_all, details_all
        return channels_all, surfaces_all

    if max_proc is None:
        max_proc = multiprocessing.cpu_count()

    max_proc = max(1, min(int(max_proc), len(tasks)))
    chunksize = max(1, len(tasks)//(max_proc*4))
    if max_proc == 1:
        results = [_calcChannelsMultipleFramesWorker(task)
                    for task in tasks]
    else:
        if mp_context is None:
            ctx = multiprocessing.get_context()
        else:
            ctx = multiprocessing.get_context(mp_context)

        with ctx.Pool(processes=max_proc) as pool:
            results = pool.map(_calcChannelsMultipleFramesWorker, tasks)
            

    for result in results:
        if return_details:
            _, channels, surfaces, details = result
            details_all.append(details)
        else:
            _, channels, surfaces = result

        channels_all.append(channels)
        surfaces_all.append(surfaces)

    for channels in channels_all:
        for channel in channels:
            channel._buildSplines()

    output_format = kwargs.get('output_format', 'pqr')
    if output_path is not None and Path(output_path).is_dir() and (
            output_format == 'mmcif' or separate):

        channel_pattern = ('channels.cif' if output_format == 'mmcif' else 'channels*_sp*_chl*.pqr')

        _writeVisTrajScript(output_path,pattern=channel_pattern)
    elif output_path is not None and not separate and output_format != 'mmcif':
        LOGGER.info("Wrote combined channels.pqr per frame; re-run with "
                    "separate=True for a multi-state PyMOL viewer, which "
                    "needs one file per channel to group across frames.")

    if return_details:
        return channels_all, surfaces_all, details_all

    return channels_all, surfaces_all


def calcSurfaceCavitiesMultipleFrames(atoms, trajectory=None, output_path=None, 
    separate=False, max_proc=2, mp_context=None, **kwargs):
    """Compute surface cavities for each frame in a trajectory or multi-model PDB.

    This function calculates surface cavities for each frame of a trajectory or
    for each model of a multi-model PDB structure. For every frame/model, it
    calls :func:`calcSurfaceCavities` and stores the detected cavities together
    with the corresponding surface representation. The `kwargs` argument is
    passed directly to :func:`calcSurfaceCavities` and can include parameters
    controlling cavity detection, filtering, and output generation.

    :arg atoms: Atomic object containing the molecular structure. For trajectory 
        analysis, this object provides the reference topology and is updated with 
        coordinates from each frame. For multi-model PDB files, the individual 
        coordinate sets are analyzed one by one.
    :type atoms: :class:`.Atomic`

    :arg trajectory: Optional trajectory or ensemble object containing multiple
        coordinate frames. If provided, surface cavities are calculated for each
        selected trajectory frame. If not provided, the function attempts to use
        multiple coordinate sets stored in `atoms`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg output_path: Optional filename used to save detected surface cavities.
        If provided, one output file is generated for each frame/model by
        appending the frame/model index to the file name. If `None`, results are
        returned but not written in the folder. Default is `None`.
    :type output_path: str or None

    :arg separate: If `True`, each detected surface cavity is saved as a separate 
        PQR/PDB file for each frame/model. If `False`, all cavities detected 
        in a given frame/model are saved in a single file. Default is `False`.
    :type separate: bool

    :arg max_proc: Maximum number of parallel processes used for calculation. 
        If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg mp_context: Multiprocessing start method used for parallel pore
        calculations. If `None`, the default method for the operating system
        is used. Windows and macOS use the ``'spawn'`` method by default,
        whereas Linux typically uses ``'fork'``. Setting
        ``mp_context='spawn'`` can be potentially used on Linux, but might 
        be slower. Available values may include ``'spawn'``, ``'fork'``,
        and ``'forkserver'``, depending on the operating system. 
        Default is `None`.
    :type mp_context: str or None

    :arg kwargs: Additional parameters passed to :func:`calcSurfaceCavities`.
        These can include `surf_radius`, `inner_radius`, `min_depth`, `max_depth`,
        `min_tetrahedra`, `max_tetrahedra`, `min_volume`, `max_volume`,
        `start_frame`, and `stop_frame`.
    :type kwargs: dict

    :returns: Two lists:
        - `cavities_all`: a list containing detected surface cavities for each
          analyzed frame/model,
        - `surfaces_all`: a list containing the corresponding surface representations 
          for each analyzed frame/model.
    :rtype: tuple (list, list)

    Example usage:
    protein = parsePDB('1tqn').select('protein')
    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, 
                                trajectory=traj, output_path="surface_cavities",
                                surf_radius=4.5, inner_radius=2.0, min_depth=1.5, max_depth=2.5, min_volume=50)

    cavities_all, surfaces_all = calcSurfaceCavitiesMultipleFrames(protein, start_frame=0, 
                                stop_frame=10, surf_radius=4.5, inner_radius=2.0) """

    if PY3K:
        if not checkAndImport('pathlib'):
            raise ImportError('To run calcSurfaceCavitiesMultipleFrames, please install pathlib.')
        from pathlib import Path
    else:
        if not checkAndImport('pathlib2'):
            raise ImportError('To run calcSurfaceCavitiesMultipleFrames, please install pathlib2 for Python 2.7.')
        from pathlib2 import Path

    _requireCoords(atoms)

    cavities_all = []
    surfaces_all = []
    tasks = []

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)

    if output_path:
        output_path = Path(output_path)
        if output_path.suffix == ".pqr":
            output_path = output_path.with_suffix('')

    if trajectory is not None:
        if isinstance(trajectory, Atomic):
            trajectory = Ensemble(trajectory)

        # As in calcChannelsMultipleFrames: only a file-backed trajectory carries
        # a read cursor, and an in-memory one has nothing to save or rewind.
        nfi = getattr(trajectory, '_nfi', None)
        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        first, last = _frameBounds(None, start_frame, stop_frame)
        traj = trajectory[first:last]

        atoms_copy = atoms.copy()
        for j0, frame0 in enumerate(traj, start=first):
            if output_path:
                frame_output_path = _frameOutputPath(output_path, j0,
                                                     "cavities")
            else:
                frame_output_path = None

            tasks.append((j0, atoms_copy, np.array(frame0.getCoords(), copy=True),
                          frame_output_path, separate, kwargs))

        if nfi is not None:
            trajectory._nfi = nfi

    else:
        if atoms.numCoordsets() > 1:
            coordsets = atoms.getCoordsets()

            first, last = _frameBounds(len(coordsets), start_frame, stop_frame)

            for i in range(first, last):
                if output_path:
                    frame_output_path = _frameOutputPath(output_path, i,
                                                         "cavities")
                else:
                    frame_output_path = None

                tasks.append((i, atoms, np.array(coordsets[i], copy=True),
                              frame_output_path, separate, kwargs))

        else:
            LOGGER.info("Include trajectory or use multi-model PDB file.")
    
    import multiprocessing

    if len(tasks) == 0:
        _warn("No frames or models were found for surface cavity calculation.")
        return cavities_all, surfaces_all

    if max_proc is None:
        max_proc = max(1, multiprocessing.cpu_count() // 2)

    max_proc = max(1, min(int(max_proc), len(tasks)))

    if max_proc == 1:
        results = [_calcSurfaceCavitiesMultipleFramesWorker(task) for task in tasks]
    else:
        if mp_context is None:
            ctx = multiprocessing.get_context()
        else:
            ctx = multiprocessing.get_context(mp_context)

        with ctx.Pool(processes=max_proc) as pool:
            results = pool.map(_calcSurfaceCavitiesMultipleFramesWorker, tasks)

    for cavities, surface in results:
        cavities_all.append(cavities)
        surfaces_all.append(surface)

    return cavities_all, surfaces_all


def calcPoresFromChannelsMultipleFrames(channels_all, details_all, output_path=None, 
    separate=False, max_proc=2, mp_context=None, **kwargs):
    """Construct pores for multiple trajectory frames or multi-model PDBs from 
    channels previously calculated with :func:`calcChannelsMultipleFrames`.

    This function applies :func:`calcPoresFromChannels` independently to each
    frame or model. The channel list and calculation details at the same index
    must correspond to the same frame/model.
    
    When using parallel calculations on Windows or macOS, the call to this
    function should be placed inside an ``if __name__ == '__main__':`` block
    to prevent child processes from executing the main script again.
    
    :arg channels_all: Lists of channels returned by :func:`calcChannelsMultipleFrames`. 
        Each element contains channels calculated for one trajectory frame or model.
    :type channels_all: list of lists

    :arg details_all: Calculation details returned by :func:`calcChannelsMultipleFrames` 
        with ``return_details=True``. Each dictionary must contain ``calculator``, 
        ``simplices``, ``neighbors``, ``vertices``, ``coords``, and ``vdw_radii`` 
        for the corresponding frame/model.
    :type details_all: list of dict

    :arg max_proc: Maximum number of parallel processes used for calculation. 
        If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg mp_context: Multiprocessing start method used for parallel pore
        calculations. If `None`, the default method for the operating system
        is used. Windows and macOS use the ``'spawn'`` method by default,
        whereas Linux typically uses ``'fork'``. Setting
        ``mp_context='spawn'`` can be potentially used on Linux, but might 
        be slower. Available values may include ``'spawn'``, ``'fork'``,
        and ``'forkserver'``, depending on the operating system. 
        Default is `None`.
    :type mp_context: str or None
    
    :arg kwargs: Pore-filtering parameters passed to
        :func:`calcPoresFromChannels`, including ``min_end_to_end``,
        ``max_end_to_end``, ``min_bottleneck``, ``max_bottleneck``,
        ``min_length``, ``max_length``, ``min_volume``, and ``max_volume``.
    :type kwargs: dict

    :returns: A list containing pores constructed for each frame/model.
        Each element corresponds to the channels and calculation details at
        the same index in ``channels_all`` and ``details_all``.
    :rtype: list of lists

    Example usage:
    channels_all, surfaces_all, details_all = calcChannelsMultipleFrames(
        protein, trajectory=dcd, return_details=True)

    pores_all = calcPoresFromChannelsMultipleFrames(channels_all, details_all, 
        min_end_to_end=40, min_bottleneck=0.7, output_path='poresALL_', separate=True)"""

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path
    
    import os
    import multiprocessing

    if len(channels_all) != len(details_all):
        raise ValueError("channels_all and details_all must contain the same number of frames")

    # A directory takes the frames inside it, as everywhere else; anything else
    # is the name they are numbered from, as before.
    into_directory = output_path is not None and os.path.isdir(str(output_path))
    if output_path is not None and not into_directory:
        output_path = Path(output_path)
        if output_path.suffix not in ('.pqr', '.pdb'):
            output_path = output_path.with_suffix('.pqr')

    tasks = []
    for position, (channels, details) in enumerate(zip(channels_all, details_all)):
        # The frame the channels came from, as calcChannelsMultipleFrames records
        # it, so that a run started past frame 0 keeps its numbers; details put
        # together otherwise are numbered by position, as they always were.
        frame_nr = details.get('frame', position)
        if into_directory:
            frame_output_path = _frameOutputPath(output_path, frame_nr, "pores")
        elif output_path is not None:
            frame_output_path = output_path.with_name(
                    "{0}_frame{1}{2}".format(output_path.stem, frame_nr, output_path.suffix))
        else:
            frame_output_path = None
        
        tasks.append((frame_nr, channels, details, frame_output_path, separate, kwargs))

    if max_proc is None:
        max_proc = multiprocessing.cpu_count()

    max_proc = max(1, min(int(max_proc), len(tasks)))

    if max_proc == 1:
        pores_all = [_calcPoresFromChannelsWorker(task) for task in tasks]
    else:
        if mp_context is None:
            ctx = multiprocessing.get_context()
        else:
            ctx = multiprocessing.get_context(mp_context)
        
        with ctx.Pool(processes=max_proc) as pool:
            pores_all = pool.map(_calcPoresFromChannelsWorker, tasks)

    return pores_all


def mergeFramesPQR(inputs, filename, prefix=None, frames=None, links=False):
    """Merge the per-frame PQR files of a multi-frame run into one multi-model
    PQR, a ``MODEL`` block per frame.

    :func:`calcChannelsMultipleFrames` and the other multi-frame functions write
    a file per frame, the frame number ending its name: ``channels12.pqr`` in a
    directory, ``<output_path>12.pqr`` otherwise. This gathers them after the
    run into the one file the multi-model viewer scripts read.

    The ``MODEL`` number is the frame number rather than a count, so a frame
    keeps its number whatever else was merged, and the files of two frame
    ranges concatenate without renumbering.

    A frame that found nothing left no file, but still gets its ``MODEL``, so
    that every frame keeps its place in a viewer that steps through the models
    in order rather than by their numbers. The block holds a ``REMARK`` saying
    so and a single ``NIL`` atom of radius 0 at the previous frame's first
    sphere: an empty ``MODEL`` is dropped altogether by ChimeraX, and ``NIL``,
    unlike the ``FIL`` of every sphere, is never read as a channel. Such frames
    are found between the first and the last file; an empty first or last
    frame leaves no trace in the files, so *frames* names them.

    Each frame's records are copied as they stand, blank lines aside. Its serial
    numbers start from 1 as they did in its own file, which keeps them within
    the five columns a PQR serial has, however many frames there are.

    :arg inputs: the per-frame files: a directory, a glob pattern, or a list
        of paths.
    :type inputs: str or list

    :arg filename: the multi-model PQR to write.
    :type filename: str

    :arg prefix: the part of the names before the frame number, such as
        ``'channels'``. Only files named ``<prefix><frame>.pqr`` are taken. Needed
        only where the files given are of more than one kind, as in a directory
        written with ``separate=True``, which holds ``channels12_chl3.pqr`` beside
        ``channels12.pqr``. Default is the one prefix all the names share.
    :type prefix: str

    :arg frames: every frame of the run, such as ``range(len(channels_all))``
        for one over a whole trajectory. Those with no file are written with
        a placeholder, first and last included; a file of a frame not listed
        is still merged. Default is the frames between the first and the last
        file found.
    :type frames: list or range

    :arg links: merge the chamber links a run writes beside each frame's
        channels, ``<prefix><frame>_links.pqr``, rather than the channels
        themselves, into a multi-model file of their own. Default is **False**.
    :type links: bool

    :returns: the filename written
    :rtype: str

    Usage:
    channels_all, surfaces_all = calcChannelsMultipleFrames(atoms,
        trajectory=dcd, output_path='frames')
    mergeFramesPQR('frames', 'channels_frames.pqr',
                   frames=range(len(channels_all)))
    mergeFramesPQR('frames', 'links_frames.pqr',
                   frames=range(len(channels_all)), links=True)"""

    import glob
    import os
    import re

    if isinstance(inputs, (list, tuple)):
        paths = [str(path) for path in inputs]
    elif os.path.isdir(str(inputs)):
        paths = glob.glob(os.path.join(str(inputs), '*.pqr'))
    else:
        paths = glob.glob(str(inputs))

    # The frame is the number ending the name - or ending it before _links, for
    # the links written beside a frame's channels - and what comes before it
    # says what kind of file it is. A name ending otherwise is no frame's file
    # and is passed over, as the links are when the channels are merged; one of
    # another kind, a separate channel's <name>12_chl3.pqr, would be read as
    # frame 3 and is refused below instead.
    tail = '_links.pqr' if links else '.pqr'
    pattern = re.compile(r'^(.*?)(\d+)' + re.escape(tail) + '$')
    found = {}
    for path in paths:
        match = pattern.match(os.path.basename(path))
        if match is None or (prefix is not None and match.group(1) != prefix):
            continue
        found.setdefault(match.group(1), {}).setdefault(
            int(match.group(2)), []).append(path)

    if not found:
        raise ValueError('no per-frame PQR files found in {0!r} named {1}<frame>'
                         '{2}'.format(inputs, prefix or '<prefix>', tail))
    if len(found) > 1:
        # Shortest first, which puts the frames' own files ahead of the
        # per-channel ones, of which there are as many kinds as frames.
        named = sorted(found, key=lambda name: (len(name), name))
        raise ValueError('the files are of more than one kind, named {0}{1}; '
                         'choose the one to merge with prefix='.format(
                             ', '.join(name + '<frame>' + tail for name in named[:4]),
                             ', ...' if len(named) > 4 else ''))

    (stem, files), = found.items()
    for frame, given in sorted(files.items()):
        if len(given) > 1:
            raise ValueError('frame {0} is given more than once: {1}'.format(
                frame, ', '.join(given)))
    files = dict((frame, given[0]) for frame, given in files.items())

    if frames is None:
        numbers = list(range(min(files), max(files) + 1))
    else:
        numbers = sorted(set(int(frame) for frame in frames) | set(files))
    empty = [frame for frame in numbers if frame not in files]

    # Where a placeholder sits: the first sphere of the frame before it, or of
    # the first frame with a file for frames that come before any.
    with open(files[min(files)]) as handle:
        anchor = next((line[30:54] for line in handle
                       if line.startswith(('ATOM', 'HETATM'))), None)

    with open(str(filename), 'w') as out:
        out.write('REMARK   {0} frames of {1}<frame>{2}, one MODEL each\n'.format(
            len(numbers), stem, tail))
        out.write('REMARK   MODEL numbers are frame numbers\n')
        if empty:
            out.write('REMARK   {0} frame{1} found nothing: a NIL atom of radius 0 '
                      'holds each place\n'.format(
                          len(empty), '' if len(empty) == 1 else 's'))
        for frame in numbers:
            out.write('MODEL%9d\n' % frame)
            if frame in files:
                with open(files[frame]) as handle:
                    lines = [line for line in handle
                             if line.strip() and line.split()[0] != 'END']
                out.writelines(lines)
                anchor = next((line[30:54] for line in lines
                               if line.startswith(('ATOM', 'HETATM'))), anchor)
            else:
                out.write('REMARK   frame {0} found nothing\n'.format(frame))
                out.write('HETATM    1  H   NIL X   1    {0}  1.00  0.00\n'.format(
                    anchor or '%8.3f%8.3f%8.3f' % (0, 0, 0)))
            out.write('ENDMDL\n')
        out.write('END\n')

    LOGGER.info('{0} frames merged into {1}{2}.'.format(
        len(numbers), filename,
        '' if not empty else ', {0} of them empty'.format(len(empty))))
    return str(filename)


def parseParameters(channels, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each
    channel in a given list of channels.

    ``object_name`` names the objects in the rows and in the file the rows are
    written to, so that a list of chamber links is not filed as channels.
    Default ``'channel'``. """

    lengths = []
    bottlenecks = []
    volumes = []
    param_file_name = kwargs.pop('param_file_name', None)
    object_name = kwargs.pop('object_name', 'channel')

    for nr_ch, channel in enumerate(channels):
        lengths.append(channel.length)
        bottlenecks.append(channel.bottleneck)
        volumes.append(channel.volume)

        if param_file_name is not None:
            with open('{0}_Parameters_All_{1}s.txt'.format(param_file_name,
                                                           object_name), "a") as f_par:
                f_par.write(("{0}_{1}{2}: {3} {4} {5}\n".format(param_file_name, object_name, nr_ch, channel.length, channel.bottleneck, channel.volume)))

    return lengths, bottlenecks, volumes


def getChannelParameters(channels, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each 
    channel in a given list of channels.

    This function iterates through a list of channel objects, extracting the
    length, bottleneck, and volume of each channel. These values are collected
    into separate lists, which are returned as a tuple for further use.

    :arg channels: A list of channel objects, where each channel has attributes
      `length`, `bottleneck`,and `volume`. These attributes represent the 
      length of the channel, the minimum radius (bottleneck) along its path, 
      and the total volume of the channel, respectively.
    :type channels: list

    :arg param_file_name: The files with parameters will be saved in a text 
        file with the provided name. Use one word which will be added to
        '_Parameters_All_channels.txt' suffix. If further analysis will be
        performed with selectChannelBySelection() function, the preferable 
        param_file_name is PDB+chain for example: '1bbhA'.
    :type param_file_name: str 

    :returns: Three lists containing the lengths, bottlenecks, and volumes of 
        the channels.
    :rtype: tuple (list, list, list)

    Example usage:
    lengths, bottlenecks, volumes = getChannelParameters(channels) """
    
    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_L_B_V = parseParameters(channels, **kwargs)
        lengths, bottlenecks, volumes = results_L_B_V
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]',
                                                             'Length [Å]', 
                                                             'Bottleneck [Å]'))
        for i in range(len(lengths)):
            LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return results_L_B_V

    except:
        for nr_i,i in enumerate(channels):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseParameters(channels[nr_i], param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results) 
            
        LOGGER.info("Channel {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                            'Length [Å]', 
                                                            'Bottleneck [Å]'))
        for frame_nr, frame in enumerate(multi_model_param):
            lengths, bottlenecks, volumes = frame
            LOGGER.info("Frame {0}".format(frame_nr))
            for i in range(len(lengths)):
                LOGGER.info("channel {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return multi_model_param


def getPoreParameters(pores, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each 
    pore in a given list of pores identified using :func:`calcPoresFromChannels`.

    This function iterates through a list of pore objects, extracting the
    length, bottleneck, and volume of each pore. These values are collected
    into separate lists, which are returned as a tuple for further use.

    :arg pores: A list of pores objects, where each pore has attributes
      `length`, `bottleneck`,and `volume`. These attributes represent the 
      length of the pore, the minimum radius (bottleneck) along its path, 
      and the total volume of the pore, respectively.
    :type pores: list

    :arg param_file_name: The files with parameters will be saved in a text 
        file with the provided name. Use one word which will be added to
        '_Parameters_All_pores.txt' suffix.
    :type param_file_name: str 

    :returns: Three lists containing the lengths, bottlenecks, and volumes of 
        the pores.
    :rtype: tuple (list, list, list)

    Example usage:
    lengths, bottlenecks, volumes = getPoreParameters(pores) """
    
    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_L_B_V = parseParameters(pores, object_name='pore', **kwargs)
        lengths, bottlenecks, volumes = results_L_B_V
        LOGGER.info("Pore {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]',
                                                             'Length [Å]', 
                                                             'Bottleneck [Å]'))
        for i in range(len(lengths)):
            LOGGER.info("pore {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return results_L_B_V

    except:
        for nr_i,i in enumerate(pores):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseParameters(pores[nr_i], object_name='pore',
                                      param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results)
            
        LOGGER.info("Pore {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                            'Length [Å]', 
                                                            'Bottleneck [Å]'))
        for frame_nr, frame in enumerate(multi_model_param):
            lengths, bottlenecks, volumes = frame
            LOGGER.info("Frame {0}".format(frame_nr))
            for i in range(len(lengths)):
                LOGGER.info("pore {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i],2), np.round(lengths[i], 2), np.round(bottlenecks[i], 2)))
        return multi_model_param


def getLinkParameters(links, **kwargs):
    """Extracts and returns the lengths, bottlenecks, and volumes of each
    chamber link returned in ``details['links']`` by :func:`calcChannels`.

    A link is the route from one chamber of a cavity into a shallower one, cut
    where it joins that chamber, so its bottleneck is the neck between the two -
    the width that decides whether the deeper site can be reached at all. Its
    length and volume measure that neck only, not the way out to the solvent,
    which is a channel of the chamber it joins.

    :arg links: A list of link objects, each with `length`, `bottleneck` and
        `volume` attributes.
    :type links: list

    :arg param_file_name: The files with parameters will be saved in a text
        file with the provided name. Use one word which will be added to
        '_Parameters_All_links.txt' suffix.
    :type param_file_name: str

    :returns: Three lists containing the lengths, bottlenecks, and volumes of
        the links.
    :rtype: tuple (list, list, list)

    Example usage:
    channels, surface, details = calcChannels(atoms, return_details=True)
    lengths, bottlenecks, volumes = getLinkParameters(details['links']) """

    results = parseParameters(links, object_name='link', **kwargs)
    lengths, bottlenecks, volumes = results
    LOGGER.info("Link {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]',
                                                     'Length [Å]',
                                                     'Bottleneck [Å]'))
    for i, link in enumerate(links):
        # Both ends named, since a link is the one object with two of them, and
        # a list of lengths says nothing about what is connected to what.
        ends = ''
        if getattr(link, 'origin', None) is not None:
            ends = "\tsp{0}".format(link.origin)
            if getattr(link, 'destination', None) is not None:
                ends += " -> sp{0}".format(link.destination)
        LOGGER.info("link {0}: \t{1} \t\t{2} \t\t{3}{4}".format(
            i, np.round(volumes[i], 2), np.round(lengths[i], 2),
            np.round(bottlenecks[i], 2), ends))
    return results


def _frameParamFileName(param_file_name, frame_nr, trajectory):
    """``2kid/2kid`` becomes ``2kid/2kid_model3`` (or ``_frame3``).

    Every frame needs its own file: the rows are keyed by object number alone,
    so writing all the frames under one name either overwrites the file each
    time - the surface cavity writer opens it with 'w', and only the last frame
    survived - or piles them up under colliding keys, which the channel and pore
    writers did by opening it with 'a'. The word matches the one the residue
    wrappers use, so the parameters and the lining of one frame stay findable
    under the same prefix."""

    if param_file_name is None:
        return None
    return '{0}_{1}{2}'.format(param_file_name,
                               'frame' if trajectory is not None else 'model',
                               frame_nr)


def getChannelParametersMultipleFrames(channels_all, trajectory=None, **kwargs):
    """Extract channel parameters for multiple frames or models.

    This function is a multi-frame wrapper for :func:`getChannelParameters`.
    It extracts channel parameters for each model or trajectory frame separately.
    Each element of ``channels_all`` is treated as the list of channels calculated
    for one frame/model.

    This function should be used with channels returned by
    :func:`calcChannelsMultipleFrames`.

    :arg channels_all: list of channel lists returned by
        :func:`calcChannelsMultipleFrames`. Each element corresponds to one
        model or trajectory frame.
    :type channels_all: list

    :arg trajectory: The trajectory the frames came from, if any. Only its
        presence is used, to name the files ``_frame<i>`` rather than
         ``_model<i>``, matching :func:`getChannelResidueNamesMultipleFrames`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list of parameter tuples for each model/frame. Each tuple contains
        channel lengths, bottlenecks, and volumes.
    :rtype: list  """

    param_file_name = kwargs.pop('param_file_name', None)
    parameters_all = []
    for frame_nr, channels in enumerate(channels_all):
        LOGGER.info("Frame/model: {0}".format(frame_nr))
        frame_name = _frameParamFileName(param_file_name, frame_nr, trajectory)
        if frame_name is not None:
            kwargs['param_file_name'] = frame_name
        params = getChannelParameters(channels, **kwargs)
        parameters_all.append(params)

    return parameters_all


def getPoreParametersMultipleFrames(pores_all, trajectory=None, **kwargs):
    """Extract pore parameters for multiple frames or models.

    This function is a multi-frame wrapper for :func:`getPoreParameters`.
    It extracts pore parameters for each model or trajectory frame separately.
    Each element of ``pores_all`` is treated as the list of pores calculated
    for one frame/model.

    This function should be used with pores returned by
    :func:`calcPoresMultipleFrames`.

    :arg pores_all: list of pore lists returned by
        :func:`calcChannelsMultipleFrames`. Each element corresponds to one
        model or trajectory frame.
    :type pores_all: list

    :arg trajectory: The trajectory the frames came from, if any. Only its
        presence is used, to name the files ``_frame<i>`` rather than
        ``_model<i>``, matching :func:`getPoreResidueNamesMultipleFrames`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list of parameter tuples for each model/frame. Each tuple contains
        pore lengths, bottlenecks, and volumes.
    :rtype: list  """

    param_file_name = kwargs.pop('param_file_name', None)
    parameters_all = []
    for frame_nr, pores in enumerate(pores_all):
        LOGGER.info("Frame/model: {0}".format(frame_nr))
        frame_name = _frameParamFileName(param_file_name, frame_nr, trajectory)
        if frame_name is not None:
            kwargs['param_file_name'] = frame_name
        params = getPoreParameters(pores, **kwargs)
        parameters_all.append(params)

    return parameters_all


def getLinkParametersMultipleFrames(links_all, trajectory=None, **kwargs):
    """Extract chamber link parameters for multiple frames or models.

    This function is a multi-frame wrapper for :func:`getLinkParameters`. Each
    element of ``links_all`` is the ``details['links']`` of one frame, which
    :func:`calcChannelsMultipleFrames` returns in its third value when called
    with ``return_details=True``; the links of a frame are not part of its
    channel list.

    :arg links_all: list of link lists, one per model or trajectory frame.
    :type links_all: list

    :arg trajectory: The trajectory the frames came from, if any. Only its
        presence is used, to name the files ``_frame<i>`` rather than
        ``_model<i>``, matching :func:`getLinkResidueNamesMultipleFrames`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list of parameter tuples for each model/frame. Each tuple contains
        link lengths, bottlenecks, and volumes.
    :rtype: list

    Example usage:
    channels_all, surfaces_all, details_all = calcChannelsMultipleFrames(
        atoms, trajectory=dcd, return_details=True)
    params = getLinkParametersMultipleFrames([d['links'] for d in details_all],
                                             trajectory=dcd)  """

    param_file_name = kwargs.pop('param_file_name', None)
    parameters_all = []
    for frame_nr, links in enumerate(links_all):
        LOGGER.info("Frame/model: {0}".format(frame_nr))
        frame_name = _frameParamFileName(param_file_name, frame_nr, trajectory)
        if frame_name is not None:
            kwargs['param_file_name'] = frame_name
        parameters_all.append(getLinkParameters(links, **kwargs))

    return parameters_all


def parseSurfaceCavityParameters(cavities, **kwargs):
    """Extract depths, volumes, and tetrahedra counts for surface cavities."""

    depths = []
    volumes = []
    tetrahedra_counts = []
    param_file_name = kwargs.pop('param_file_name', None)
    lines = []
    
    if param_file_name is not None:
        lines.append("# Cavity_id Volume [Å³] Depth [Å] Tetrahedra_count\n")
    
    for nr_cav, cavity in enumerate(cavities):
        depth = cavity.depth
        volume = cavity.volume
        tetrahedra_count = len(cavity.tetrahedra)
        depths.append(depth)
        volumes.append(volume)
        tetrahedra_counts.append(tetrahedra_count)

        if param_file_name is not None:
            lines.append("{0}_cavity{1}: {2:.3f} {3:.2f} {4}\n".format(
                param_file_name, nr_cav, volume, depth, tetrahedra_count))

    if param_file_name is not None:
        with open(param_file_name + '_Parameters_All_surface_cavities.txt', "w") as f_par:
            f_par.writelines(lines)
    
    return volumes, depths, tetrahedra_counts


def getSurfaceCavityParameters(cavities, **kwargs):
    """Extract volumes, depths, and tetrahedra counts of surface cavities.

    This function iterates through a list of surface cavity objects and extracts
    the volume, depth, and number of tetrahedra assigned to each cavity. These
    values are returned as separate lists and can optionally be saved to a text file.

    :arg cavities: A list of surface cavity objects returned by :func:`calcSurfaceCavities`.
    :type cavities: list

    :arg param_file_name: Optional name used to save cavity parameters to a text file. 
        The suffix '_Parameters_All_surface_cavities.txt' will be added.
    :type param_file_name: str

    :returns: Three lists containing volumes, depths, and tetrahedra counts.
    :rtype: tuple (list, list, list)

    Example usage:
    volumes, depths, tetrahedra_counts = getSurfaceCavityParameters(cavities)
    """

    multi_model_param = []
    param_file_name = kwargs.get('param_file_name', None)

    try:
        results_V_D_T = parseSurfaceCavityParameters(cavities, **kwargs)
        volumes, depths, tetrahedra_counts = results_V_D_T

        LOGGER.info("Cavity {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                           'Depth [Å]', 
                                                           'Tetrahedra count'))

        for i in range(len(volumes)):
            LOGGER.info("cavity {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i], 2), np.round(depths[i], 2),
                tetrahedra_counts[i]))

        return results_V_D_T

    except:
        for nr_i, i in enumerate(cavities):
            safe_param_file_name = param_file_name if param_file_name is not None else ""
            results = parseSurfaceCavityParameters(cavities[nr_i],
                param_file_name=safe_param_file_name + str(nr_i))
            multi_model_param.append(results)

        LOGGER.info("Cavity {0}: \t{1} \t{2} \t{3}".format('ID', 'Volume [Å³]', 
                                                           'Depth [Å]', 
                                                           'Tetrahedra count'))

        for frame_nr, frame in enumerate(multi_model_param):
            volumes, depths, tetrahedra_counts = frame
            LOGGER.info("Frame {0}".format(frame_nr))

            for i in range(len(volumes)):
                LOGGER.info("cavity {0}: \t{1} \t\t{2} \t\t{3}".format(i, np.round(volumes[i], 2), np.round(depths[i], 2),
                    tetrahedra_counts[i]))

        return multi_model_param


def getSurfaceCavityParametersMultipleFrames(cavities_all, trajectory=None, **kwargs):
    """Provides surface cavity parameters for multiple frames or models.

    It analyzes surface cavities calculated for multi-model PDB files or
    trajectories and returns cavity parameters for each model/frame.

    :arg cavities_all: list of surface cavity lists returned by
        :func:`calcSurfaceCavitiesMultipleFrames`.
    :type cavities_all: list
    
    :arg trajectory: The trajectory the frames came from, if any. Only its
        presence is used, to name the files ``_frame<i>`` rather than
        ``_model<i>``, matching
        :func:`getSurfaceCavityResidueNamesMultipleFrames`.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg param_file_name: base name for the output parameter files. If provided,
        one file will be written for each model/frame with the frame/model index
        added to the file name.
    :type param_file_name: str

    :returns: A list with surface cavity parameters for each frame/model.
    :rtype: list """

    param_file_name = kwargs.pop('param_file_name', None)
    parameters_all = []

    for i, cavities in enumerate(cavities_all):
        LOGGER.info("Model/frame: {0}".format(i))
        frame_name = _frameParamFileName(param_file_name, i, trajectory)
        if frame_name is not None:
            kwargs['param_file_name'] = frame_name
        params = getSurfaceCavityParameters(cavities, **kwargs)
        parameters_all.append(params)

    return parameters_all


def _sampleObjectSpheres(object, num_samples=5):
    """``(centres, radii)`` of the probe spheres along a channel or a pore.

    ``num_samples`` points per tetrahedron of the route, evenly spaced in the
    spline parameter - which, the centerline being parameterized by the square
    root of the step between circumcenters, puts them roughly evenly along the
    route rather than crowding wherever circumcenters happen to cluster.
    This is the one definition of "the spheres of an object":
    :func:`getChannelAtoms` and
    :meth:`~ChannelCalculator.saveChannelsToPdb` write these very spheres out as
    FIL atoms, and the lining queries take them straight from here, since the
    radius they need survives no better in a PDB radius column than in memory."""

    centerline_spline, radius_spline = object.getSplines()
    t = _objectSampleParameters(object, num_samples)
    if len(t) == 0:
        return np.empty((0, 3)), np.empty(0)

    return centerline_spline(t), radius_spline(t)


def _objectSampleParameters(object, num_samples=5):
    """The spline parameters :func:`_sampleObjectSpheres` samples at.

    Split out because the mmCIF export needs the parameters themselves and not only
    the spheres read off them: arc length along the route has to be evaluated at the
    same points the spheres came from, or the profile's ``distance`` column would
    describe a slightly different sampling than its coordinates do. Recomputing the
    ``linspace`` at the call site would work until one of the two definitions was
    edited and the other was not."""

    samples = len(object.tetrahedra) * num_samples
    if samples < 1:
        return np.empty(0)

    spline = object.getSplines()[0]
    return np.linspace(spline.x[0], spline.x[-1], samples)

def _sampleSplineByArcLength(object, n_samples=1000):
    """
    Sample a Channel spline approximately uniformly in physical arc length.

    Returns
    -------
    t : (n_samples,) ndarray
        Spline parameters.
    points : (n_samples, 3) ndarray
        Centerline coordinates.
    radii : (n_samples,) ndarray
        Channel radii at those coordinates.
    s : (n_samples,) ndarray
        Cumulative arc length corresponding to the samples.
    """
    centerline_spline, radius_spline = object.getSplines()

    if centerline_spline is None or radius_spline is None:
        object._buildSplines()
        centerline_spline, radius_spline = object.getSplines()

    t0 = centerline_spline.x[0]
    t1 = centerline_spline.x[-1]

    if n_samples < 2:
        n_samples = 2
    n_dense = max(10 * n_samples, 1000)
    td = np.linspace(t0, t1, n_dense)
    xd = np.asarray(centerline_spline(td))
    dx = np.diff(xd, axis=0)
    ds = np.linalg.norm(dx, axis=1)
    sd = np.empty(n_dense, dtype=float)
    sd[0] = 0.0
    sd[1:] = np.cumsum(ds)
    total_length = sd[-1]

    if total_length <= 1e-12:
        t = np.full(n_samples, t0)
        points = np.repeat(xd[:1], n_samples, axis=0)
        radii = np.full(n_samples, float(radius_spline(t0)))
        s = np.zeros(n_samples)
        return t, points, radii, s

    s = np.linspace(0.0, total_length, n_samples)
    t = np.interp(s, sd, td)

    points = np.asarray(centerline_spline(t))
    radii = np.asarray(radius_spline(t), dtype=float)

    return t, points, radii, s


def getChannelAtoms(channels, protein=None, num_samples=5):
    """Generates an AtomGroup object representing the atoms along the paths of 
    the given channels and optionally combines them with an existing protein 
    structure.

    This function takes a list of channel objects and generates atomic 
    representations of the channels based on their centerline splines and 
    radius splines. The function samples points along each channel's centerline
     and assigns atom positions at these points with corresponding radii, 
    creating a list of PDB-formatted lines. These lines are then converted 
    into an AtomGroup object using the ProDy library. If a protein structure is
     provided, it is combined with the generated channel atoms by merging their
     respective PDB streams.

    :arg channels: A list of channel objects. Each channel has a method 
        `getSplines()` that
        returns the centerline spline and radius spline of the channel.
    :type channels: list

    :arg protein: An optional AtomGroup object representing a protein structure.
        If provided, it will be combined with the generated channel atoms.
    :type protein: prody.atomic.AtomGroup or None

    :arg num_samples: The number of atom samples to generate along each segment
         of the channel. More samples result in a finer representation of the
         channel. Default is 5.
    :type num_samples: int

    :returns: An AtomGroup object representing the combined atoms of the 
        channels and the protein, if a protein is provided.
    :rtype: prody.atomic.AtomGroup

    Example usage:
    atomic_structure = getChannelAtoms(channels, protein) """
    
    if PY3K:
        import io
    else:
        import StringIO as io
    
    from prody import parsePDBStream, writePDBStream

    def convert_lines_to_atomic(atom_lines):
        pdb_text = "\n".join(atom_lines)
        pdb_stream = io.StringIO(pdb_text)
        structure = parsePDBStream(pdb_stream)
        return structure

    atom_index = 1
    pdb_lines = []

    if not isinstance(channels, list):
        channels = [channels]
    
    for channel in channels:
        centers, radii = _sampleObjectSpheres(channel, num_samples)

        for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1], 
                                                  centers[:, 2], radii), 
                                                  start=atom_index):
            pdb_lines.append("ATOM  %5d  H   FIL T   1    %8.3f%8.3f%8.3f%6.2f%6.2f\n" % (i, x, y, z, 1.00, radius))
    
    if protein is not None:
        protein_stream = io.StringIO()
        writePDBStream(protein_stream, protein, csets=[protein.getACSIndex()])
        
        protein_stream.seek(0)

        protein_lines = protein_stream.readlines()
        if protein_lines[-1].strip() == 'END':
            protein_lines = protein_lines[:-1]
        
        combined_pdb_text = "".join(protein_lines) + "\n".join(pdb_lines) + "\nEND\n"
        combined_stream = io.StringIO(combined_pdb_text)
        combined_structure = parsePDBStream(combined_stream)

        return combined_structure

    channels_atomic = convert_lines_to_atomic(pdb_lines)
    return channels_atomic

def _atomRadii(atoms, warn=True):
    """Van der Waals radius of every atom of *atoms*, from :data:`VDW_RADII`.

    :meth:`ChannelCalculator.getVdwRadii` reading an :class:`.Atomic` rather than
    a list of symbols, and covering the one case it cannot: a file that gave no
    element column at all is read as unknown throughout instead of as the string
    ``None``. Both reporting paths measure through here, so a structure that one
    of them can describe is never one the other dies on.

    *warn* is passed down. A cavity report reads the radii twice, over all the
    atoms and again over the dry ones, and the substitution is worth saying once
    rather than twice for what is nearly the same set."""

    elements = atoms.getElements()
    if elements is None:
        elements = np.zeros(atoms.numAtoms(), dtype='<U2')
    return ChannelCalculator.getVdwRadii(elements, warn=warn)

def _vertexRadiiSource(atoms):
    """``(tree, vdw_radii)`` over the atoms a tessellation would have used.

    Water is excluded, because :func:`calcChannels` drops it before tessellating:
    a cavity holding a water is still a cavity, and measuring the vertex against
    that water would report no room where the diagram found plenty. On 3A2M this is
    the difference between a vertex of radius 2.08 and one of radius -1.27, sitting
    0.25 A from the oxygen of HOH 610."""

    dry = atoms.select('not water')
    if dry is None:
        dry = atoms

    # Silent: getSurfaceCavityResidueNames has already measured the same atoms
    # through _liningSource, which said whatever there was to say about them.
    return _kdTree(dry), _atomRadii(dry, warn=False)


def _vertexRadii(points, source, k=24):
    """Inscribed-sphere radius at each Voronoi vertex, ``min(|v - x| - vdw)``.

    Surface cavities keep no radius of their own. The radius column of a cavity file
    is a marker size (see :meth:`ChannelCalculator.saveCavitiesToPdb`), not the room
    at the vertex, so the radius the diagram implies is recomputed here. Only the
    nearest *k* atoms are consulted, which is ample: van der Waals radii span barely
    half an Angstrom, so no atom outside the nearest handful can hold the minimum."""

    tree, vdw = source
    k = min(k, len(vdw))
    distances, indices = tree.query(np.asarray(points, dtype=float), k=k)
    if k == 1:
        distances, indices = distances[:, None], indices[:, None]
    return (distances - vdw[indices]).min(axis=1)

def _liningSource(atoms):
    """``(tree, coords, vdw_radii)`` for the lining queries over *atoms*.

    The radii come from :data:`VDW_RADII`, the table :func:`calcChannels`
    tessellates with, so the clearance the lining measures is a clearance from
    the surface the channel was actually carved against, and not from some other
    scale's idea of the same atom. Several such scales exist and they disagree by
    more than the distance being measured, so which one is used matters less than
    that it is the same one throughout.

    An element the table does not cover falls back on its ``UNKNOWN`` entry rather
    than raising, and says so. The lining is reported against whatever structure
    the caller passes, which is routinely wider than the selection that was traced
    - the ions, cofactors and metals the tessellation never saw - so it is the
    report rather than the trace that meets an uncovered element first."""

    return _kdTree(atoms), atoms.getCoords(), _atomRadii(atoms)


def _liningResidues(atoms, source, points, radii, distA, deep=None):
    """Whole residues of *atoms* whose surface comes within *distA* of the probe's.

    The probe spheres are given by *points* and *radii*, and a residue lines the
    object when one of its atoms satisfies ``|x_p - x_a| - r_p - r_a <= distA`` -
    the gap between the two surfaces. Both radii belong in it.

    The probe radius, because the probe is a sphere and not a point: measured from
    its centre a fixed *distA* reaches only ``distA - r_p`` past the surface, so it
    gathers a second shell where the object is narrow and misses the wall where it
    is wide.

    The atom radius, because without it the reach past an atom's *own* surface is
    ``distA - r_a``, which varies with the element. It is widest at a hydrogen,
    narrower at every heavy atom, and negative for the alkali and alkaline-earth
    ions, which then have to overlap the probe before being reported at all.
    Whole-residue completion hides that wherever a residue is polyatomic, one
    qualifying atom carrying the rest, and leaves it exposed for exactly the
    monatomic species a lining report exists to name. It also made the criterion
    depend on whether the depositor modelled hydrogens, a hydrogen being the
    element the missing term favours most.

    *source* is a :func:`_liningSource` bundle and not something derived from
    *atoms* here, because the reports call this once per object against one
    structure: on a large structure, building the tree and the radii again per
    object costs more than the query it serves.

    *deep* is an optional dict, updated in place with the atoms lying inside the
    route rather than beside it; :func:`_reportDeepLiningAtoms` reads it.

    Residues are completed but never widened past what the caller supplied, which
    is what the ``same residue as`` selection this replaces also did."""

    # The probe each atom answered for is what separates the two functions; this
    # one wants only the set of atoms, so the pairing is dropped here.
    _, candidates, _ = _liningContacts(source, points, radii, distA, deep)
    if len(candidates) == 0:
        return None

    resindices = atoms.getResindices()
    lining = np.unique(resindices[np.unique(candidates)])
    selected = np.flatnonzero(np.isin(resindices, lining))

    if hasattr(atoms, 'getAtomGroup'):   # a selection: index back into its group
        return atoms.getAtomGroup()[atoms.getIndices()[selected].tolist()]
    return atoms[selected.tolist()]

def _reportDeepLiningAtoms(atoms, deep):
    """Name the residues lying inside the route rather than beside it.

    *deep* is what :func:`_liningResidues` collected. These are atoms the
    tessellation never saw: an ion or a cofactor sitting in the route, which is a
    result worth having and is what measuring to the atom centre was most likely
    to miss, or a structure that is not the one the channels were traced on, which
    is not. Geometry cannot tell the two apart - a permeant ion and a mismatched
    structure overlap the route alike - so this states what was found and leaves
    the reading to the caller, rather than warning about a mismatch that is
    usually a ligand."""

    if not deep:
        return

    indices = np.fromiter(deep, dtype=int, count=len(deep))
    resnames, resnums = atoms.getResnames(), atoms.getResnums()
    _, first = np.unique(atoms.getResindices()[indices], return_index=True)
    labels = ['{0}{1}'.format(resnames[i], resnums[i])
              for i in indices[first].tolist()]

    LOGGER.info("{0} residue(s) lie inside the route rather than beside it "
                "({1}), reaching {2:.2f} A past its surface. These atoms were "
                "not part of the tessellation the channels came from.".format(
                    len(labels), ', '.join(labels[:8]) +
                    (', ...' if len(labels) > 8 else ''), -min(deep.values())))

def _oneLetterResname(residue):
    """One-letter code for an amino acid or a nucleotide, residue name for the rest.

    Only residue names that ProDy knows as amino acids or nucleic acids are
    translated. AAMAP holds both directions -- LYS to K, and also K to LYS -- so
    looking a ligand up in it renames it into an amino acid: the ion K becomes LYS
    and F becomes PHE, while the cofactors SAM and SAH become M and C, each
    indistinguishable in the report from the residue it now names. Ligands, ions,
    cofactors and sugars therefore keep their PDB chemical component ID, which
    identifies them and cannot be misread as a residue.

    The test is the name alone, against ``flags.AMINOACIDS`` and the nucleic
    definition, not the ``protein`` and ``nucleic`` flags: ProDy plants the
    ``protein`` flag through ``calpha``, so a histidine deposited as HETATM, or any
    residue whose CA is missing, is not flagged protein and would lose its one-letter
    code for want of an atom that has nothing to do with what the residue is."""

    from prody.atomic import flags
    from prody.atomic.atomic import AAMAP

    resname = residue.getResname()

    if resname in flags.DEFINITIONS['nucleic']:
        # NAMAP is ProDy's own nucleotide map and covers the modified bases. AAMAP
        # is no use here: it reads GUA as glutamate and CYT as tyrosine.
        return flags.NAMAP.get(resname, resname)

    if resname in flags.AMINOACIDS:
        # The CHARMM/AMBER histidine names are either absent from AAMAP or, for HSE,
        # mapped to serine, so they are resolved before the lookup.
        if resname in ('HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP'):
            return AAMAP['HIS']
        return AAMAP.get(resname, resname)

    return resname


_LiningOptions = namedtuple('_LiningOptions',
                            'distA residues_file_name one_letter_aa '
                            'include_water include_chain')


def _popLiningOptions(kwargs):
    """The reporting options shared by every lining function, popped from *kwargs*.

    :func:`getObjectResidueNames` and :func:`getSurfaceCavityResidueNames` both
    receive them as keywords forwarded by their wrappers, so a default is written
    once here instead of in two lists that have to be kept in step."""

    return _LiningOptions(
        distA=kwargs.pop('distA', 1.5),
        residues_file_name=kwargs.pop('residues_file_name', None),
        one_letter_aa=kwargs.pop('one_letter_aa', False),
        include_water=kwargs.pop('include_water', False),
        include_chain=kwargs.pop('include_chain', True))

def _infoAsListOfStrings(info,key='Residues',object_type = 'channel', sep=", "):
    """Format channel metadata as a list of human-readable strings.

    Parameters
    ----------
    info : dict or list
        Channel information to format. This can be a mapping from channel names to
        per-channel dictionaries, or a list of such mappings for multiple
        conformations. Each dictionary is expected to contain a key named by
        ``key`` whose value is a sequence of residue or atom labels.
    key : str, optional
        The dictionary key whose associated values should be converted to a string
        representation for each channel.
    sep : str, optional
        Separator used to join multiple entries within a single channel.

    Returns
    -------
    list
        A list of formatted strings. For a single structure, each entry is
        ``"{channel}: {values}"``. For multiple conformations, each element is a
        list of such strings, one per channel.

    Notes
    -----
    This helper is used to turn nested channel dictionaries into compact output
    suitable for logging or reporting while preserving the per-channel
    organization of the source data.
    """

    if not info:
        return []
      
    def _joinList(lst,sep=sep ):
        if not lst:
            return ""
        if len(lst) == 1:
            return str(lst[0])
        return sep.join(map(str, lst))

    result = []
    if isinstance(info[0],dict):
        for i in range(len(info)):
            lst = _joinList(info[i][key])
            result.append(f"{object_type}_{i}: {lst}")
        return result
    
    for j in range(len(info)):
        tmp = []
        for i in range(len(conf)):
            lst = _joinList(conf[i][key])
            tmp.append(f"Conformation_{j}_{object_type}_{i}: {lst}")
        result.append(tmp)
    return result

def getLiningResindexMask(resindices, tree, points, radii, buffer=2.5):
    """Return a boolean mask selecting residues lining a set of channel points.

    Parameters
    ----------
    resindices : array-like
        Residue indices for each coordinate in the structure. The values are used
        to map queried point hits back to the residues they belong to.
    tree : scipy.spatial.cKDTree or similar
        Spatial index queried with the channel point coordinates and radii.
    points : array-like
        Coordinates of points sampled along the channel centerline.
    radii : array-like
        Radius associated with each channel point. A point is considered a hit if
        a structure coordinate lies within the corresponding distance.
    buffer : float, optional
        Additional padding added to each radius before querying the spatial tree.
        Defaults to 2.5.

    Returns
    -------
    numpy.ndarray or None
        A boolean mask over the input residue indices, with ``True`` values for
        residues that intersect any of the point/radius neighborhoods. Returns
        ``None`` when no points are provided or when no hits are found.
    """
    
    from itertools import chain
    if len(points) == 0:
        return None
    
    r = radii + buffer
    
    hits = tree.query_ball_point(np.asarray(points, dtype=float),
                                 np.asarray(r,dtype=float))

    flat_idx = np.fromiter(chain.from_iterable(hits), dtype=np.int64)
    if flat_idx.size==0:
        return None

    selected = np.unique(flat_idx)
    selected_resindices = resindices[selected]
    mask = np.isin(resindices, selected_resindices)

    return mask


def getResidueOrientation(context):
    """Return the signed orientation of each lining residue relative to a channel center.

    This callback is intended for ``getChannelLiningInfo`` and evaluates whether a
    residue points toward or away from the channel interior. For each residue, the
    function combines the direction from the residue C-alpha to the side-chain
    center-of-mass with the direction from the residue C-alpha to the nearest
    channel center point, then returns the dot product of these normalized vectors.

    Parameters
    ----------
    context : dict
        A dictionary containing the filtered atom data for the channel context,
        including the atom mask, coordinates, residue indices, and channel center
        points. Expected keys are ``mask``, ``coords``, ``resindices``,
        ``all_ca``, ``all_hydrogen``, ``all_backbone``, and ``points``.

    Returns
    -------
    numpy.ndarray
        One value per residue with a positive sign indicating the residue points
        toward the channel center and a negative sign indicating it points away
        from the channel.
    """

    mask = context['mask']
    coords_sub =context['coords'][mask]
    resindices_sub = context['resindices'][mask]
    ca_mask = context['all_ca'][mask]
    hydrogen = context['all_hydrogen'][mask]
    backbone = context['all_backbone'][mask]
    centers = context['points']
    unique_res, inverse = np.unique(resindices_sub, return_inverse=True)
    n_res = len(unique_res)
    
    ca_coords = np.full((n_res, 3), np.nan)
    ca_coords[inverse[ca_mask]] = coords_sub[ca_mask]
    
    sc_mask = ~backbone & ~hydrogen
    sc_sum = np.zeros((n_res, 3))
    sc_count = np.zeros(n_res)
    np.add.at(sc_sum, inverse[sc_mask], coords_sub[sc_mask])
    np.add.at(sc_count, inverse[sc_mask], 1)

    sc_com = np.full((n_res, 3), np.nan)
    has_sc = sc_count > 0
    sc_com[has_sc] = sc_sum[has_sc] / sc_count[has_sc, None]

    valid = ~np.isnan(sc_com[:, 0]) & ~np.isnan(ca_coords[:, 0])
    ca_coords = ca_coords[valid]
    sc_com = sc_com[valid]

    center_tree = _kdTree(centers) 
    _, idx = center_tree.query(ca_coords)

    v = sc_com - ca_coords
    v /= np.linalg.norm(v, axis=1, keepdims=True)

    n = centers[idx] - ca_coords
    n /= np.linalg.norm(n, axis=1, keepdims=True)

    return np.einsum("ij,ij->i", v, n)


_LINING_SUFFIXES = ('.pkl', '.txt', '.csv', '.tsv')


def _jsonDefault(obj):
    """json.dumps fallback so callback results such as ndarrays and numpy
    scalars (e.g. getObjectResponseProfile) can be written to .txt."""
    if hasattr(obj, 'tolist'):
        return obj.tolist()
    return str(obj)


def saveLiningInfo(info, filepath, **kwargs):
    """
    Save a dict of dicts to disk in the requested format.

    Parameters
    ----------
    info : dict
        Dict of dicts, e.g. {record_id: {field: value, ...}, ...}
    filepath : str
        Output path. .pkl extension is added if no extension is passed.
    kwargs :
        Extra options:
        - for txt: none currently
        - for csv/tsv: fieldnames (list) to control/limit columns and order
    """
    import pickle
    import csv
    import json
    from pathlib import Path

    filepath = Path(filepath)
    if filepath.suffix == "":
        filepath = filepath.with_suffix(".pkl")
    fmt = filepath.suffix

    if fmt == ".pkl":
        with open(filepath, "wb") as f:
            pickle.dump(info, f)

    elif fmt == ".txt":
        with open(filepath, "w") as f:
            for key, subdict in info.items():
                line = json.dumps({key: subdict}, default=_jsonDefault)
                f.write(line + "\n")

    elif fmt in (".csv", ".tsv"):
        delimiter = "," if fmt == ".csv" else "\t"

        fieldnames = kwargs.get("fieldnames")
        if fieldnames is None:
            fieldnames = []
            for subdict in info.values():
                for k in subdict.keys():
                    if k not in fieldnames:
                        fieldnames.append(k)

        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(
                f, fieldnames=["id"] + fieldnames, delimiter=delimiter
            )
            writer.writeheader()
            for key, subdict in info.items():
                row = {"id": key}
                row.update(subdict)
                writer.writerow(row)

    else:
        raise ValueError("Unknown format '{0}'. Use .pkl, .txt, .csv or .tsv.".format(fmt))

    return filepath

def _liningContacts(source, points, radii, distA, deep=None):
    """
    Find atom/probe contacts using a surface-to-surface distance criterion.
    For a probe sphere p and atom a:
        gap = |x_p - x_a| - r_p - r_a
    gap <= distA  -> atom lines the object
    gap < -0.5    -> atom substantially overlaps the sampled object surface
    Returns
    -------
    probes : ndarray
        Probe indices for qualifying contacts.
    candidates : ndarray
        Atom indices for qualifying contacts.
    gaps : ndarray
        Surface-to-surface gaps for qualifying contacts.
    """
   
    tree, coords, atom_radii = source
    points = np.asarray(points, dtype=float)
    radii = np.asarray(radii, dtype=float)
    empty = (
        np.empty(0, dtype=int),
        np.empty(0, dtype=int),
        np.empty(0, dtype=float),
    )

    if len(points) == 0 or len(coords) == 0:
        return empty

    distA = float(distA)
    query_radii = radii + distA + float(atom_radii.max())
    hits = tree.query_ball_point(points, query_radii)
    counts = np.fromiter(
        (len(hit) for hit in hits),
        dtype=np.int64,
        count=len(hits))
    if not counts.any():
        return empty
    candidates = np.fromiter((atom for hit in hits for atom in hit),dtype=np.int64,count=int(counts.sum()))
    probes = np.repeat(np.arange(len(points), dtype=np.int64),counts)

    delta = points[probes] - coords[candidates]
    center_dist = np.sqrt(np.einsum('ij,ij->i', delta, delta))
    gaps = (center_dist- radii[probes]- atom_radii[candidates])

    if deep is not None:
        inside = gaps < -0.5

        if inside.any():
            for atom, gap in zip(
                candidates[inside].tolist(),
                gaps[inside].tolist(),
            ):
                if gap < deep.get(atom, 0.0):
                    deep[atom] = gap
    within = gaps <= distA

    if not within.any():
        return empty

    return (probes[within],candidates[within],gaps[within])


def _liningResindexMaskFromContacts(resindices,source,points,radii,distA,deep=None):
    """
    Return a boolean atom mask for residues lining an object.

    A residue lines the object when at least one of its atoms satisfies

        |xp - xa| - rp - ra <= distA

    where rp is the sampled object-sphere radius and ra is the atom's
    van der Waals radius.

    Also returns the minimum qualifying surface gap for every atom.
    Atoms with no qualifying contact have np.inf.

    Returns
    -------
    mask : ndarray or None
        Boolean mask over atoms.
    atom_gaps : ndarray or None
        Minimum qualifying gap for each atom.
    """
    _, candidates, gaps = _liningContacts(source,points,radii,distA,deep,)

    if len(candidates) == 0:
        return None, None

    atom_gaps = np.full(len(resindices),np.inf,dtype=float)

    np.minimum.at(atom_gaps,candidates,gaps)
    
    lining_resindices = np.unique(resindices[np.unique(candidates)])
    mask = np.isin(resindices,lining_resindices)
    return mask, atom_gaps


def getLiningResidues(context):
    """Return unique residue labels for residues lining the object."""
    from prody.atomic.atomic import AAMAP

    mask = context.get('mask')
    if mask is None or not np.any(mask):
        return []

    resindices = context['resindices'][mask]
    resnames = context['all_resnames'][mask]
    resnums = context['all_resnums'][mask]
    chains = context['all_chains'][mask]

    # One representative atom per residue.
    _, first_idx = np.unique(resindices,return_index=True,)
    resnames = resnames[first_idx]
    resnums = resnums[first_idx]
    chains = chains[first_idx]

    if context['options'].one_letter_aa:
        resnames = np.array([AAMAP['HIS'] if aa in ('HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP')
            else AAMAP.get(aa, aa) for aa in resnames])

    labels = np.char.add(resnames.astype(str),resnums.astype(str),)

    if context['options'].include_chain:
        labels = np.char.add(labels,np.char.add(':', chains.astype(str)),)
    return np.unique(labels).tolist()


def getSiteScores(context):
    kwargs = context['kwargs']
    mask = context['mask']
    scores = kwargs.get('scores')
    if isinstance(scores,dict):
        subscores = {}
        for key,values in scores.items():
            subscores[key] = values[mask]
        return subscores
    subscores = scores[mask]
    return subscores

def getSiteScoreStats(context):
    subscores = getSiteScores(context)
    if subscores is None:
        return {}
    if isinstance(subscores, dict):
        return {key:{'max': float(np.max(values)),
                    'median': float(np.median(values))}
                    for key, values in subscores.items()}

    return {'max': float(np.max(subscores)),
            'median': float(np.median(subscores))}

def getLiningResIndices(context):
    """Return the unique residue indices for the currently selected lining atoms.

    Parameters
    ----------
    context : dict
        Calculation context containing the atom mask and residue index array.

    Returns
    -------
    numpy.ndarray
        Unique residue numbers for the atoms that line the object.
    """
    mask = context['mask']
    resindices = context['resindices'][mask]
    unique_resindices = np.unique(resindices)
    return unique_resindices

def getObjectLiningInfo(atoms, objects, object_type='channel',
                        callbacks=None, **kwargs):
    """
    Provides residue information for atoms lining channels, pores,
    or links.

    Lining is defined by the surface-to-surface distance between a
    sampled object sphere and an atom's van der Waals sphere:

        gap = |xp - xa| - rp - ra

    A residue lines the object when at least one atom satisfies:

        gap <= distA

    ``distA`` therefore represents an actual surface-to-surface
    clearance, rather than a centre-to-centre buffer.

    Deep atoms are detected during the same contact calculation:

        gap < -0.5

    Deep-atom findings accumulate across all objects and are reported
    once at the end of the function.

    Parameters
    ----------
    atoms : AtomGroup
        Structure used for residue/atom identification and VDW radii.
    objects : object or list
        Channel/pore/link object or objects.
    object_type : {'channel', 'pore', 'link'}
        Type of object.
    callbacks : dict, optional
        Callback functions receiving the calculation context.

    Returns
    -------
    info : dict
        Lining information for each object.
    """

    aliases = {
        'channel': 'channels',
        'channels': 'channels',
        'pore': 'pores',
        'pores': 'pores',
        'link': 'links',
        'links': 'links'}

    info = []
    context = {}
    deep = {}

    options = _popLiningOptions(kwargs)
    context['options'] = options

    if options.include_water:
        atoms = atoms.select('not resname FIL')
    else:
        atoms = atoms.select('not water and not resname FIL')

    if object_type not in aliases:
        raise ValueError(f"object_type must be  in {list(aliases.keys())}")

    object_type = aliases[object_type]

    _requireCoords(atoms)

    source = _liningSource(atoms)
    tree = source[0]
    context['kwargs'] = kwargs
    context['source'] = source
    context['coords'] = source[1]
    context['vdw_radii'] = source[2]
    context['tree'] = tree

    context['resindices'] = atoms.getResindices()
    context['atom_names'] = atoms.getNames()
    context['elements'] = atoms.getElements()

    context['all_resnames'] = atoms.getResnames()
    context['all_resnums'] = atoms.getResnums()
    context['all_chains'] = atoms.getChids()

    context['all_ca'] = (context['atom_names'] == 'CA')

    context['all_hydrogen'] = ( context['elements'] == 'H')

    context['all_backbone'] = np.isin(context['atom_names'],['N', 'CA', 'C', 'O'],)

    object_list = (objects if isinstance(objects, list)else [objects])

    for i, obj in enumerate(object_list):
        info.append({})
        points, radii = _sampleObjectSpheres(obj)

        mask, atom_gaps = _liningResindexMaskFromContacts(context['resindices'],
                                                            source, points, radii,
                                                            options.distA, deep)
        context['points'] = points
        context['radii'] = radii
        context['mask'] = mask
        context['atom_gaps'] = atom_gaps

        for name, func in callbacks.items():
            info[i][name] = (func(context) if mask is not None else [])

    _reportDeepLiningAtoms(atoms,deep)

    if options.residues_file_name is not None:
        output_file = saveLiningInfo(info,options.residues_file_name,**kwargs)
        LOGGER.info("Lining Info was saved to: {0}".format(output_file))
        
    return info


def getObjectResidueNames(atoms, channels, **kwargs):
    '''Provides the resnames and resid of residues that are forming the object(s). 
    Residues are extracted based on distA, the clearance between the surface of
    the FIL atoms (object atoms) and the van der Waals surface of the residue's
    own atoms.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg objects: A list of objects. Each object has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the object.
    :type objects: list
    
    :arg object_type: Type of the object; "channel", "pore" or "link".
        Default is "channel".
    :type object_type: str

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. It is a clearance between surfaces: the local probe radius
        and the atom's van der Waals radius are both taken off, so a touching
        atom sits at 0 and the reach is the same in a wide part of the object as
        in a narrow one, and the same at a hydrogen as at a potassium ion.
        Default is 1.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_{object_type}.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''


    info = getObjectLiningInfo(atoms,
                                channels,
                                object_type='channel',
                                callbacks = {'Residues': getLiningResidues},
                                **kwargs )
    
    resnames = _infoAsListOfStrings(info, object_type = object_type, **kwargs)
    return resnames

def getObjectLiningInfoMultipleFrames(atoms, objects_all, trajectory=None, object_type='channel',
                                                callbacks={'Residues':getLiningResidues}, **kwargs):
    '''Provides the resnames and resid of residues that are forming the object(s) in
    multiple frames/models. 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (object atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg objects_all: A list of objects. Each object has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the object.
    :type objects_all: list
    
    :arg trajectory: optional trajectory object. If provided, coordinates are
        taken from trajectory frames. If None, a multi-model PDB is assumed and
        models are selected using ``setACSIndex``.
    :type trajectory: :class:`.Trajectory` or None

    :arg object_type: Type of the object; "channel", "pore" or "link".
        Default is "channel".
    :type object_type: str

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_{object_type}.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    residues_file_name = kwargs.pop('residues_file_name', None)
    info_all = []

    if object_type not in ('channel', 'pore', 'link'):
        raise ValueError("object_type must be 'channel', 'pore' or 'link'")

    if trajectory is None:
        # multi-model PDB
        for frame_pos, objects in enumerate(objects_all):
            model_index = start_frame + frame_pos

            if stop_frame != -1 and model_index > stop_frame:
                break

            LOGGER.info("Model: {0}".format(model_index))
            atoms.setACSIndex(model_index)

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_model{}".format(model_index)
            else:
                frame_residues_file_name = None
            
            # Straight to the general form: the per-type wrappers do nothing but
            # pass object_type on, and branching over them here meant every new
            # type had to be added in three places.
            lining_info = getObjectLiningInfo(atoms, objects,
                                    object_type=object_type,
                                    residues_file_name=frame_residues_file_name,callbacks=callbacks, **kwargs)

            info_all.append(lining_info)

    else:
        # trajectory / DCD
        nfi = getattr(trajectory, '_nfi', None)

        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        if stop_frame == -1:
            traj = trajectory[start_frame:]
        else:
            traj = trajectory[start_frame:stop_frame + 1]

        atoms_copy = atoms.copy()
        for frame_pos, frame in enumerate(traj):
            frame_index = start_frame + frame_pos

            if frame_pos >= len(objects_all):
                break

            LOGGER.info("Frame: {0}".format(frame_index))
            atoms_copy.setCoords(frame.getCoords())

            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_frame{}".format(frame_index)
            else:
                frame_residues_file_name = None

            lining_info = getObjectLiningInfo(atoms_copy, objects_all[frame_pos],
                                    object_type=object_type,
                                    residues_file_name=frame_residues_file_name,callbacks=callbacks, **kwargs)

            info_all.append(lining_info)

        if nfi is not None:
            trajectory._nfi = nfi

    return info_all

def getChannelResidueNames(atoms, channels, **kwargs):
    '''Provides the resnames and resid of residues that are forming the pore(s). 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (pore atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg pores: A list of pore objects. Each pore has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the pore.
    :type pores: list

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_pores.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''

    info = getObjectLiningInfo(atoms, 
                                channels, 
                                object_type='channel',
                                callbacks={'Residues':getLiningResidues},
                                **kwargs)
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames

def getPoreResidueNames(atoms, pores, **kwargs):
    '''Provides the resnames and resid of residues that are forming the pore(s). 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (pore atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg pores: A list of pore objects. Each pore has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the pore.
    :type pores: list

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_pores.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''

    info = getObjectLiningInfo(atoms, 
                                pores, 
                                object_type='pore',
                                callbacks={'Residues':getLiningResidues},
                                **kwargs)
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames


def getLinkResidueNames(atoms, links, **kwargs):
    '''Provides the resnames and resid of residues that are forming the chamber
    link(s) returned in ``details['links']`` by :func:`calcChannels`.
    Residues are extracted based on distA which is the distance between FIL atoms
    (link atoms) and protein residues. A link runs from one chamber into a
    shallower one and is cut where it joins it, so the residues reported are
    those lining the neck between the two sites, not a whole route to the solvent.
    Results could be save as txt file by providing the `residues_file_name` parameter.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg links: A list of link objects. Each link has a method `getSplines()`
        that returns the centerline spline and radius spline of the link.
    :type links: list

    :arg distA: Residues reaching within this distance of the cavity's surface
        are reported. It is a clearance between surfaces: the inscribed radius at
        each cavity vertex and the atom's van der Waals radius are both taken
        off, so a touching atom sits at 0 and the reach is the same in a wide
        cavity as in a narrow one, and the same at a hydrogen as at a potassium
        ion. Default is 1.5 Å
    :type distA: int, float

    :arg residues_file_name: The file with residues will be saved in a text
        file with the provided name. Use one word which will be added to
        '_Residues_All_links.txt' sufix.
    :type residues_file_name: str

    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''
    info = getObjectLiningInfo(atoms, 
                                links, 
                                object_type='link',
                                callbacks={'Residues':getLiningResidues},
                                **kwargs)
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames

def getChannelResidueNamesMultipleFrames(atoms, channels, trajectory=None, **kwargs):
    '''Provides the resnames and resid of residues that are forming the channel(s). 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (channel atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg channels: A list of channel objects. Each channel has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the channel.
    :type channels: list

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_channels.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''
    info = getObjectLiningInfoMultipleFrames(atoms,
                                     channels, 
                                     trajectory=trajectory, 
                                    object_type='channel',
                                      callbacks={'Residues':getLiningResidues},
                                        **kwargs)
    
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames

def getPoreResidueNamesMultipleFrames(atoms, pores, trajectory=None, **kwargs):
    '''Provides the resnames and resid of residues that are forming the pore(s). 
    Residues are extracted based on distA which is the distance between FIL atoms 
    (pore atoms) and protein residues.
    Results could be save as txt file by providing the `residues_file_name` parameter.
    
    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`

    :arg pores: A list of pore objects. Each pore has a method 
        `getSplines()` that returns the centerline spline and radius spline of 
        the pore.
    :type pores: list

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float
    
    :arg residues_file_name: The file with residues will be saved in a text 
        file with the provided name. Use one word which will be added to 
        '_Residues_All_pores.txt' sufix. If further analysis will be 
        performed with selectChannelBySelection() function, the preferable 
        residues_file_name is PDB+chain for example: '1bbhA'.
    :type residues_file_name: str  
    
    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''
    info = getObjectLiningInfoMultipleFrames(atoms,
                                     pores, 
                                     trajectory=trajectory, 
                                    object_type='pore',
                                      callbacks={'Residues':getLiningResidues},
                                        **kwargs)
    
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames


def getLinkResidueNamesMultipleFrames(atoms, links, trajectory=None, callbacks={'Residues':getLiningResidues}, **kwargs):
    '''Provides the resnames and resid of residues that are forming the chamber
    link(s) in multiple frames/models.

    Each element of ``links`` is the ``details['links']`` of one frame, which
    :func:`calcChannelsMultipleFrames` returns in its third value when called
    with ``return_details=True``; the links of a frame are not part of its
    channel list.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`

    :arg links: list of link lists, one per model or trajectory frame.
    :type links: list

    :arg trajectory: Trajectory object containing multiple frames. If not given,
        the coordinate sets of *atoms* are used and the files are named
        ``_model<i>`` rather than ``_frame<i>``.
    :type trajectory: :class:`.Atomic`, :class:`.Ensemble`, or trajectory-like object

    :arg distA: Residues reaching within this distance of the object's surface
        are reported. The local probe radius is added to it, so the reach past
        the surface is the same in a wide part of the object as in a narrow one.
        The distance runs to the atom centre, so a touching atom sits at about
        one van der Waals radius. Default is 2.5 [Ang]
    :type distA: int, float

    :arg residues_file_name: The file with residues will be saved in a text
        file with the provided name, one per frame/model, with
        '_Residues_All_links.txt' added.
    :type residues_file_name: str

    :arg one_letter_aa: Whether to apply 1-latter code to residue name
        by defult is False. Only amino acids and nucleotides are translated;
        ligands, cofactors and ions keep their residue name, so that the ion K
        stays K rather than being read as a lysine.
    :type one_letter_aa: bool

    :arg include_water: Whether to list water molecules among the lining
        residues. They are reported one entry per molecule, so a solvated
        structure gives hundreds of them. Default is False.
    :type include_water: bool

    :arg include_chain: Whether to append the chain identifier to each residue,
        as ``ASP108:A``. Default is True: without it a residue number does not
        identify a residue, since an oligomer lines a channel with residues of
        the same number from several chains. Pass False for the plain
        ``ASP108`` labels written by earlier versions.
    :type include_chain: bool  '''
    info = getObjectLiningInfoMultipleFrames(atoms,
                                     links, 
                                     trajectory=trajectory, 
                                    object_type='link',
                                      callbacks={'Residues':getLiningResidues},
                                        **kwargs)
    
    resnames = _infoAsListOfStrings(info,object_type=object_type,**kwargs)
    return resnames


def getSurfaceCavityLiningInfo(atoms, cavities, surface, callbacks=None,
                               **kwargs):
    '''Provides lining information for surface cavities, in the same form as
    :func:`getObjectLiningInfo`.

    Cavity points are the Voronoi vertices assigned to each cavity
    (``surface[4][cavity.tetrahedra]``), each with the probe radius that
    :func:`_vertexRadii` fits to it. Lining is then the same surface-to-surface
    criterion used for channels, pores and links::

        gap = |xp - xa| - rp - ra  <=  distA

    A residue lines the cavity when at least one of its atoms satisfies it.

    :arg atoms: an Atomic object from which residues are selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg cavities: surface cavity objects returned by :func:`calcSurfaceCavities`.
    :type cavities: list

    :arg surface: surface data returned by :func:`calcSurfaceCavities`;
        ``surface[4]`` must hold the Voronoi vertices.
    :type surface: list

    :arg callbacks: dict ``{name: func}``. Each ``func`` receives the
        calculation context and its return value is stored under ``name``.
        Default is ``{'Residues': getLiningResidues}``. Besides the keys that
        :func:`getObjectLiningInfo` provides (``mask``, ``points``, ``radii``,
        ``atom_gaps``, ``resindices`` ...) the context carries ``cavity`` and
        ``object_index``. Cavities without vertices give ``[]`` for every
        callback.
    :type callbacks: dict

    The remaining keywords (``distA``, ``residues_file_name``,
    ``one_letter_aa``, ``include_water``, ``include_chain``) behave as in
    :func:`getObjectResidueNames`. ``residues_file_name`` is passed to
    :func:`saveLiningInfo`, so its extension picks the format (.pkl default,
    .txt, .csv, .tsv). Any other keywords reach the callbacks through
    ``context['kwargs']``.

    :returns: ``{'cavity0': {callback_name: result, ...}, ...}``
    :rtype: dict
    '''

    if callbacks is None:
        callbacks = {'Residues': getLiningResidues}

    _requireCoords(atoms)
    options = _popLiningOptions(kwargs)

    if surface is None or len(surface) < 5:
        raise ValueError('surface must contain Voronoi vertices in surface[4]')

    vertices = surface[4]
    if not isinstance(cavities, list):
        cavities = [cavities]

    # Probe radii come from the full structure the diagram was built on, not
    # from the water-filtered one used to decide what lines the cavity.
    radii_source = _vertexRadiiSource(atoms)
    lining_atoms = None
    if options.include_water:
        lining_atoms = atoms.select('not resname FIL')
    else:
        lining_atoms = atoms.select('not water and not resname FIL')

    source = _liningSource(lining_atoms)
    context = {'options': options,
                'kwargs': kwargs,
                'source': source,
                'tree': source[0],
                'coords': source[1],
                'vdw_radii': source[2],
                'resindices': lining_atoms.getResindices(),
                'atom_names': lining_atoms.getNames(),
                'elements': lining_atoms.getElements(),
                'all_resnames': lining_atoms.getResnames(),
                'all_resnums': lining_atoms.getResnums(),
                'all_chains': lining_atoms.getChids()}

    context['all_ca'] = (context['atom_names'] == 'CA')
    context['all_hydrogen'] = (context['elements'] == 'H')
    context['all_backbone'] = np.isin(context['atom_names'], ['N', 'CA', 'C', 'O'])
    source = context['source']

    info = []
    intruding = 0

    for i, cavity in enumerate(cavities):
        info.append({})

        if cavity.tetrahedra is None or len(cavity.tetrahedra) == 0:
            for name in callbacks:
                info[key][name] = []
            continue

        points = vertices[cavity.tetrahedra]
        # No deep-atom accounting: a cavity keeps no radius of its own, so the
        # probe is the min over the very atoms being reported and the gap to
        # the atom holding that minimum is identically zero. The intruding
        # count below is the signal that survives.
        radii = _vertexRadii(points, radii_source)
        intruding += int((radii < 0).sum())

        mask, atom_gaps = _liningResindexMaskFromContacts(context['resindices'], source, points, radii, options.distA)

        context['cavity'] = cavity
        context['object_index'] = i
        context['points'] = points
        context['radii'] = radii
        context['mask'] = mask
        context['atom_gaps'] = atom_gaps

        for name, func in callbacks.items():
            info[i][name] = func(context) if mask is not None else []

    # A vertex of the diagram cannot lie inside an atom the diagram was built
    # from, so a negative radius means these atoms are not the ones the
    # cavities came from.
    if intruding:
        _warn("{0} cavity vertices lie inside an atom of the supplied structure. "
              "The cavities were calculated on a different set of atoms than the "
              "one given here, so their lining is reported against the wrong "
              "structure.".format(intruding))

    if options.residues_file_name is not None:
        output_file = saveLiningInfo(info, options.residues_file_name, **kwargs)
        LOGGER.info("Lining Info was saved to: {0}".format(output_file))

    return info


def getSurfaceCavityResidueNames(atoms, cavities, surface, **kwargs):
    '''Provides the resnames and resid of residues that form surface cavities.

    Thin wrapper over :func:`getSurfaceCavityLiningInfo` with the default
    ``Residues`` callback, returning one string per cavity as before. See that
    function for the arguments (``distA``, ``residues_file_name``,
    ``one_letter_aa``, ``include_water``, ``include_chain``).

    :returns: A list of residue names and residue numbers for each surface cavity.
    :rtype: list
    '''

    info = getSurfaceCavityLiningInfo(atoms, cavities, surface,
                                      callbacks={'Residues': getLiningResidues},
                                      **kwargs)
    return _infoAsListOfStrings(info,object_type='SurfaceCavity',**kwargs)



def getSurfaceCavityLiningInfoMultipleFrames(atoms, cavities_all, surfaces_all,
                                             trajectory=None, callbacks=None,
                                             **kwargs):
    """Multi-frame version of :func:`getSurfaceCavityLiningInfo`, for results
    of :func:`calcSurfaceCavitiesMultipleFrames`.

    Cavities calculated for frame/model ``i`` are analysed against the
    coordinates of frame/model ``i``.

    :arg trajectory: optional trajectory object. If None, a multi-model PDB
        is assumed and models are selected with ``setACSIndex``.
    :arg start_frame: first frame/model index to analyze. Default is 0.
    :arg stop_frame: last frame/model index to analyze. Default is -1 (all
        available).
    :arg callbacks: as in :func:`getSurfaceCavityLiningInfo`.
    :arg residues_file_name: base name; ``_modelX`` / ``_frameX`` is inserted
        before the extension, giving one file per frame.

    :returns: list of ``getSurfaceCavityLiningInfo`` dicts, one per frame/model.
    :rtype: list
    """

    if callbacks is None:
        callbacks = {'Residues': getLiningResidues}

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    residues_file_name = kwargs.pop('residues_file_name', None)

    info_all = []

    if trajectory is None:
        # multi-model PDB: the lists hold one entry per analysed frame,
        # numbered from start_frame.
        available = min(len(cavities_all), len(surfaces_all))
        first, last = _frameBounds(start_frame + available, start_frame,
                                   stop_frame)
        for model_index in range(first, last):
            cavities = cavities_all[model_index - first]
            surface = surfaces_all[model_index - first]
            atoms.setACSIndex(model_index)
            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_frame{}".format(model_index)
            else:
                frame_residues_file_name = None
            info = getSurfaceCavityLiningInfo(
                atoms, cavities, surface, callbacks=callbacks,
                residues_file_name=frame_residues_file_name,
                **kwargs)
            info_all.append(info)

    else:
        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        first, last = _frameBounds(None, start_frame, stop_frame)
        atoms_copy = atoms.copy()
        for frame_pos, frame in enumerate(trajectory[first:last]):
            # frames outlast the results when the trajectory is longer than
            # the run that produced them
            if frame_pos >= min(len(cavities_all), len(surfaces_all)):
                break

            frame_index = first + frame_pos
            atoms_copy.setCoords(frame.getCoords())
            if residues_file_name is not None:
                frame_residues_file_name = residues_file_name + "_frame{}".format(frame_index)
            else:
                frame_residues_file_name = None

            info = getSurfaceCavityLiningInfo(
                atoms_copy, cavities_all[frame_pos], surfaces_all[frame_pos],
                callbacks=callbacks,
                residues_file_name=frame_residues_file_name,
                **kwargs)
            info_all.append(info)

    return info_all


def _getContext(context):
    """Return a copy of the current lining context."""
    return context.copy()


def getConnectedCavityChannelLiningInfo(atoms, connected, surface,
                                        callbacks=None, merge=True, **kwargs):
    """
    Calculate lining information for connected surface-cavity/channel pairs.

    Each unique cavity and channel is processed only once.

    If merge=False, cavity and channel lining information are reported
    separately as ``cavity_info`` and ``channel_info``.

    If merge=True, the cavity and channel contexts for each connection are
    merged and callbacks are evaluated on the combined lining region.
    """
    if surface is None:
        raise ValueError("surface is required")

    if callbacks is None:
        callbacks = {'Residues': getLiningResidues}

    # Collect unique objects by their original indices.
    cavities = {}
    channels = {}

    for connection in connected:
        cavity_index = connection['cavity_index']
        channel_index = connection['channel_index']

        cavities.setdefault(cavity_index, connection['cavity'])
        channels.setdefault(channel_index, connection['trimmed_channel'])

    cavity_indices = list(cavities)
    channel_indices = list(channels)

    if merge:
        # Get the calculation context for every unique object.
        context_callback = {'context': _getContext}

        cavity_results = getSurfaceCavityLiningInfo(atoms,
                                                    list(cavities.values()),
                                                    surface, index=False,
                                                    callbacks=context_callback,
                                                    **kwargs)

        channel_results = getObjectLiningInfo(atoms,
                                            list(channels.values()),
                                            object_type='channel', index=False,
                                            callbacks=context_callback,
                                            **kwargs)

        # Map original indices -> contexts.
        cavity_contexts = {index: result['context'] for index, result in zip(cavity_indices, cavity_results)}

        channel_contexts = {index: result['context'] for index, result in zip(channel_indices, channel_results)}

        # Merge the appropriate cavity/channel context for each connection.
        for connection in connected:
            cavity_context = cavity_contexts[connection['cavity_index']]
            channel_context = channel_contexts[connection['channel_index']]

            combined_context = cavity_context.copy()
            combined_context['points'] = np.concatenate((cavity_context['points'], channel_context['points']))
            combined_context['radii'] = np.concatenate(( cavity_context['radii'], channel_context['radii']))
            combined_context['atom_gaps'] = np.minimum( cavity_context['atom_gaps'], channel_context['atom_gaps'] )
            combined_context['mask'] = (cavity_context['mask'] | channel_context['mask'])
            connection['merged_info'] = {key: callback(combined_context) for key, callback in callbacks.items()}

        return connected

    # Normal separate cavity/channel reporting.
    cavity_info = getSurfaceCavityLiningInfo(atoms,
                                    list(cavities.values()),
                                    surface, index=False,
                                    callbacks=callbacks,
                                    **kwargs)

    channel_info = getObjectLiningInfo(atoms,
                                list(channels.values()),
                                object_type='channel',
                                callbacks=callbacks, index=False,
                                **kwargs)

    cavity_info = {index: result for index, result in zip(cavity_indices, cavity_info)}
    channel_info = {index: result for index, result in zip(channel_indices, channel_info)}
    for connection in connected:
        connection['cavity_info'] = (cavity_info[connection['cavity_index']])
        connection['channel_info'] = (channel_info[connection['channel_index']])

    return connected

def getConnectedCavityChannelLiningInfoMultipleFrames(atoms, connected, surface, trajectory=None,
                                                      callbacks=None, merge=True, **kwargs):
    """
    Calculate connected cavity/channel lining information over multiple frames.

    Parameters
    ----------
    atoms : AtomGroup
        Structure used for residue/atom identification and VDW radii.
    connected : list
        Per-frame lists of connected cavity/channel dictionaries.
    surface : list
        Surface object for each analysed frame.
    trajectory : Trajectory, optional
        Trajectory supplying coordinates. If None, coordinates are taken
        from the coordinate sets in ``atoms``.
    callbacks : dict, optional
        Callback functions receiving the lining calculation context.
    merge : bool, default=True
        If True, report merged cavity/channel lining information.

    Returns
    -------
    info_all : list
        Lining information for each analysed frame.
    """
    if callbacks is None:
        callbacks = {'Residues': getLiningResidues}

    start_frame = kwargs.pop('start_frame', 0)
    stop_frame = kwargs.pop('stop_frame', -1)
    residues_file_name = kwargs.pop('residues_file_name', None)

    connected_all = connected
    surfaces_all = surface

    if len(connected_all) != len(surfaces_all):
        LOGGER.warning("The number of connected-cavity/channel results ({0}) does not "
            "match the number of surfaces ({1}). Only frames for which both "
            "are available will be analysed.".format(len(connected_all), len(surfaces_all)))

    available = min(len(connected_all), len(surfaces_all))

    if trajectory is None:
        # Multi-model PDB: results are assumed to correspond to coordinate
        # sets beginning at start_frame.
        first, last = _frameBounds(start_frame + available, start_frame, stop_frame)

        for model_index in range(first, last):
            result_index = model_index - first

            if result_index >= available:
                break

            atoms.setACSIndex(model_index)

            if residues_file_name is not None:
                frame_residues_file_name = (residues_file_name + "_frame{}".format(model_index))
            else:
                frame_residues_file_name = None

            connected_all[result_index] = getConnectedCavityChannelLiningInfo(atoms,
                                    connected_all[result_index],
                                    surfaces_all[result_index],
                                    callbacks=callbacks,
                                    residues_file_name=frame_residues_file_name,
                                    merge=merge,
                                    **kwargs)

    else:
        if hasattr(trajectory, 'reset'):
            trajectory.reset()

        first, last = _frameBounds(None, start_frame, stop_frame)

        atoms_copy = atoms.copy()

        for frame_pos, frame in enumerate(trajectory[first:last]):
            # The trajectory may contain more frames than the channel/cavity
            # calculation that produced connected_all and surfaces_all.
            if frame_pos >= available:
                break

            frame_index = first + frame_pos
            atoms_copy.setCoords(frame.getCoords())

            if residues_file_name is not None:
                frame_residues_file_name = (residues_file_name + "_frame{}".format(frame_index))
            else:
                frame_residues_file_name = None

            connected_all[frame_pos]= getConnectedCavityChannelLiningInfo(atoms_copy,
                                                connected_all[frame_pos],
                                                surfaces_all[frame_pos],
                                                callbacks=callbacks,
                                                residues_file_name=frame_residues_file_name,
                                                merge=merge,
                                                **kwargs)


    return connected_all

def getSurfaceCavityResidueNamesMultipleFrames(atoms, cavities_all,
                                               surfaces_all, trajectory=None,
                                               **kwargs):
    """Per-frame lists of cavity residue strings; wrapper over
    :func:`getSurfaceCavityLiningInfoMultipleFrames` with the ``Residues``
    callback. Arguments as in that function."""

    info_all = getSurfaceCavityLiningInfoMultipleFrames(
        atoms, cavities_all, surfaces_all, trajectory=trajectory,
        callbacks={'Residues': getLiningResidues}, **kwargs)
    return [_infoAsListOfStrings(info,object_type='SurfaceCavity',**kwargs) for info in info_all]


def selectChannelBySelection(atoms, residue_sele, **kwargs):
    """Select PQR files with channels that are having FIL residues within 
    certain distance (distA) from selected residue (temporarily one residue).
    If not all files should be included use pqr_files to provide the new list. 
    For example:
    pqr_files = [file for file in os.listdir('.') if file.startswith('7lafA_') and file.endswith('.pqr')]
    pqr_files = [file for file in os.listdir('.') if '5kbd' in file and file.endswith('.pqr')]

    :arg atoms: an Atomic object from which residues are selected 
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg residue_sele: selection string
                        for example: 'resid 377 and chain A', 'resid 10 to 20'
    :type residue_sele: str
   
    :arg pqr_files: list of PQR files to analyze
                    default is False, which means that all .pqr files from the 
                    current directory will be analyzed.
    :type pqr_files: bool or list
    
    :arg folder_name: The name of the folder to which PDBs will be extracted
    :type folder_name: str

    :arg distA: non-zero value, maximal distance from selected region to 
        channel (FIL atoms). Default is 5.
    :type distA: int, float 
        
    :arg residues_file: File with residues forming the channel created by 
        getChannelResidues(), default is False 
    :type residues_file: bool

    :arg param_file: File with residues forming the channel created by
        getChannelParameters(). Default is False.
    :type param_file: bool

    :arg file_prefix: The name the residue and parameter files were written
        under - what was passed as ``residues_file_name`` / ``param_file_name``
        to :func:`getChannelResidueNames` and :func:`getChannelParameters`. By
        default it is read from the PQR file names, which works when they share
        that name (``1bbhA_sp0_chl3.pqr`` beside ``1bbhA_Residues_All_channels.txt``).
        Give it where the two differ, and where the PQR files carry no name at
        all because the run was told only a directory (``sp0_chl3.pqr``).
    :type file_prefix: str

    :arg object_name: The kind of object to extract rows for: ``'channel'``,
        ``'link'``, ``'pore'`` or ``'cavity'``. Files of the other kinds are
        copied as before but contribute no rows, so a directory holding channels
        and links can be searched for either. Default is ``'channel'``; the
        matching ``residues_suffix``, ``parameters_suffix``,
        ``selected_residues_output`` and ``selected_parameters_output`` are
        accepted as keywords beside it.
    :type object_name: str  """

    _requireCoords(atoms)
    
    import os, shutil
    
    pqr_files = kwargs.pop('pqr_files', False)
    distA = kwargs.pop('distA', 5)
    folder_name = kwargs.pop('folder_name', 'selected_files')
    residues_file = kwargs.pop('residues_file', False)
    param_file = kwargs.pop('param_file', False)

    object_name = kwargs.pop('object_name', 'channel')
    file_prefix = kwargs.pop('file_prefix', None)
    residues_suffix = kwargs.pop('residues_suffix', '_Residues_All_channels.txt')
    parameters_suffix = kwargs.pop('parameters_suffix', '_Parameters_All_channels.txt')
    selected_residues_output = kwargs.pop('selected_residues_output', 
                                          'Selected_channel_residues.txt')
    selected_parameters_output = kwargs.pop('selected_parameters_output', 
                                            'Selected_channel_parameters.txt')

    copied_files_list = []
    
    if pqr_files == False:
        # take all PDBs from the current dir
        pqr_files = [file for file in os.listdir('.') if file.endswith('.pqr')]

    residue_sele = atoms.select(residue_sele)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    for i in pqr_files:
        channel = parsePQR(i)
        if 'FIL' in np.unique(channel.getResnames()):
            sele_FIL = channel.select('same residue as exwithin '+str(distA)+' of center', center=residue_sele.getCoords())

            if sele_FIL is not None:
                shutil.copy(i, folder_name)
                copied_files_list.append(i)
            else:
                pass

    if copied_files_list:
        LOGGER.info("Filtered files are now in: {0}".format(folder_name))

    # Extract parameters and/or residues for the selected objects. The rows are
    # found by reading each file name back into (prefix, object, number) rather
    # than by splitting on the object name, which never matched: the name in the
    # file is the short tag - chl, lnk - while the row inside is written under
    # the full word.
    if residues_file == True:
        selected_residues = _selectObjectRows(copied_files_list, residues_suffix,
                                              object_name, file_prefix)
        with open(selected_residues_output, 'w') as f_out:
            f_out.writelines(selected_residues)
        LOGGER.info("{0} residue row(s) saved to: {1}".format(
            len(selected_residues), selected_residues_output))

    if param_file == True:
        selected_param = _selectObjectRows(copied_files_list, parameters_suffix,
                                           object_name, file_prefix)
        with open(selected_parameters_output, 'w') as f_out:
            f_out.writelines(selected_param)
        LOGGER.info("{0} parameter row(s) saved to: {1}".format(
            len(selected_param), selected_parameters_output))

    LOGGER.info("Selected files: ")
    LOGGER.info(' '.join(copied_files_list))


def selectSurfaceCavityBySelection(atoms, residue_sele, **kwargs):
    """Select PQR files with surface cavities located close to a selected region.

    This function is a surface-cavity wrapper for :func:`selectChannelBySelection`.
    It selects PQR files containing surface cavities represented by FIL pseudoatoms
    that are located within a user-defined distance from a selected protein region.

    Surface cavity files should be generated by :func:`calcSurfaceCavities` or
    :func:`calcSurfaceCavitiesMultipleFrames`, preferably with ``separate=True``,
    so that each cavity is saved in an individual PQR file.

    :arg atoms: an Atomic object from which the reference region is selected
    :type atoms: :class:`.Atomic`, :class:`.LigandInteractionsTrajectory`

    :arg residue_sele: selection string defining the reference region.
        For example: ``'resid 173'``, ``'resid 173 and chain A'``, ``'resid 170 to 180'``.
    :type residue_sele: str

    :arg pqr_files: list of PQR files to analyze. If not provided, all PQR files
        from the current directory will be analyzed.
    :type pqr_files: bool or list

    :arg folder_name: name of the folder to which selected PQR files will be copied.
        Default is ``'selected_surface_cavities'``.
    :type folder_name: str

    :arg distA: maximal distance between the selected region and surface cavity
        FIL atoms. Default is 5 Å.
    :type distA: int, float

    :arg residues_file: if True, residue information for selected surface cavities
        will be extracted from ``*_Residues_All_surface_cavities.txt`` and saved
        to ``Selected_surface_cavity_residues.txt``.
        Default is False.
    :type residues_file: bool

    :arg param_file: if True, parameter information for selected surface cavities
        will be extracted from ``*_Parameters_All_surface_cavities.txt`` and saved
        to ``Selected_surface_cavity_parameters.txt``.
        Default is False.
    :type param_file: bool """

    kwargs.setdefault('object_name', 'cavity')
    kwargs.setdefault('residues_suffix', '_Residues_All_surface_cavities.txt')
    kwargs.setdefault('parameters_suffix', '_Parameters_All_surface_cavities.txt')
    kwargs.setdefault('selected_residues_output', 
                      'Selected_surface_cavity_residues.txt')
    kwargs.setdefault('selected_parameters_output', 
                      'Selected_surface_cavity_parameters.txt')

    return selectChannelBySelection(atoms, residue_sele, **kwargs)


def calcChannelSurfaceOverlaps(**kwargs):
    """Calculate overlapping parts of the predicted channels, tunnels, and 
    pores denote as 'FIL' atoms. Results are normalized within [0,1].

    :arg resolution: Surface sampling resolution.
        default is 0.5
    :type resolution: float
    
    :arg max_proc: Maximum number of parallel processes used to voxelize individual
        PQR files. If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None
    
    :arg mp_context: Multiprocessing start method used for parallel pore
        calculations. If `None`, the default method for the operating system
        is used. Windows and macOS use the ``'spawn'`` method by default,
        whereas Linux typically uses ``'fork'``. Setting
        ``mp_context='spawn'`` can be potentially used on Linux, but might 
        be slower. Available values may include ``'spawn'``, ``'fork'``,
        and ``'forkserver'``, depending on the operating system. 
        Default is `None`.
    :type mp_context: str or None

    :arg output_file_name: The name of the PDB file with overlapping surfaces.
    :type output_file_name: str

    :arg pqr_files: File with residues forming the channel created by
        getChannelResidues(). Default is False, then all the files from the
        current directory will be analyzed. When providing a list, only the
        PDBs from the list will be analyzed. When providing str, it will be
        treated as a folder path.
    :type pqr_files: bool, list or str
    
    Example usage:
    calcChannelSurfaceOverlaps() - all the files in the current directory will 
    be analyzed
    
    from pathlib import Path
    pqr_files = [str(f) for f in Path(".").glob("channels_*.pqr")]
    calcChannelSurfaceOverlaps(pqr_files=pqr_files, 
                output_file_name='results.pdb', max_proc=4)
    - files with the "channels_" prefix will be selected from the current folder 
    and analyzed using four parallel processes.
    
    calcChannelSurfaceOverlaps(pqr_files='./DATA', output_file_name='results.pdb') 
    - only files from the DATA folder will be analyzed and results will be saved 
    as results.pdb
    
    list_of_files = ['file1.pqr', 'file2.pqr', 'file3.pqr', ..]
    calcChannelSurfaceOverlaps(pqr_files=list_of_files, output_file_name='results.pdb') 
    - files from the list will be analyzed and results will be saved as results.pdb
    """
    
    import os
    import multiprocessing
    from collections import Counter

    resolution = kwargs.pop('resolution', 0.5)
    max_proc = kwargs.pop('max_proc', 2)
    mp_context = kwargs.pop('mp_context', None)
     
    pqr_files = kwargs.pop('pqr_files', False)
    if pqr_files == False or pqr_files is None:
        # take all PQRs from the current dir
        pqr_files = [file for file in os.listdir('.') if file.endswith('.pqr')]
    elif isinstance(pqr_files, str):
        # folder path
        if os.path.isdir(pqr_files):
            pqr_files = [os.path.join(pqr_files, file) for file in os.listdir(pqr_files) if file.endswith('.pqr')]
    elif isinstance(pqr_files, list):
        # list of PQRs
        pqr_files = [file for file in pqr_files if file.endswith('.pqr')]
    else:
        raise ValueError('Please provide list with PQR files, folder path, or nothing to analyze PQRs in the current folder')

    output_file_name = kwargs.pop('output_file_name','overlap_regions.pdb')

    # PRQ files might be empty    
    valid_pqr_files = []
    for pqr_file in pqr_files:
        if not os.path.isfile(pqr_file) or os.path.getsize(pqr_file) == 0:
            LOGGER.warn("Skipping empty PQR file: {0}".format(pqr_file))
            continue

        valid_pqr_files.append(pqr_file)

    pqr_files = valid_pqr_files

    if len(pqr_files) == 0:
        LOGGER.info("No PQR files found.")
        return None

    if os.path.exists(output_file_name):
        os.rename(output_file_name, output_file_name + '-old')

    if max_proc is None:
        max_proc = multiprocessing.cpu_count()

    max_proc = int(max_proc)
    if max_proc < 1:
        max_proc = 1

    max_proc = min(max_proc, len(pqr_files))

    LOGGER.info("Number of PQR files: {0}".format(len(pqr_files)))
    LOGGER.info("Resolution: {0}".format(resolution))
    LOGGER.info("max_proc: {0}".format(max_proc))

    merged_surface = Counter()
    tasks = [(pqr_file, resolution) for pqr_file in pqr_files]

    if max_proc > 1:
        LOGGER.info("Calculating overlaps using {0} processes.".format(max_proc))
        chunksize = max(1, len(tasks) // (max_proc * 4))
        if mp_context is None:
            ctx = multiprocessing.get_context()
        else:
            ctx = multiprocessing.get_context(mp_context)
        
        with ctx.Pool(processes=max_proc) as pool:
            for surface in pool.imap_unordered(_surfaceFromPqrWorker, tasks,
                                               chunksize=chunksize):
                merged_surface.update(surface)

    else:
        for pqr_file in pqr_files:
            LOGGER.info("Processing file: {0}".format(pqr_file))
            surface = _surfaceFromPqrWorker((pqr_file, resolution))
            merged_surface.update(surface)

    with open(output_file_name, 'w') as out:
        atom_id = 1

        for (ix, iy, iz), count in merged_surface.items():
            x = ix * resolution
            y = iy * resolution
            z = iz * resolution

            norm_count = float(count) / float(len(pqr_files))

            out.write("ATOM  {:5d}  H   FIL T   1    {:8.3f}{:8.3f}{:8.3f}{:6.2f}  1.00\n"
                .format(atom_id, x, y, z, norm_count))

            atom_id += 1

    LOGGER.info("Overlap written to: {0}".format(output_file_name))
    LOGGER.info("Number of occupied overlap voxels: {0}".format(len(merged_surface)))

    return output_file_name
    

def calcSurfaceCavityOverlaps(**kwargs):
    """Calculate overlapping regions of surface cavities represented as FIL atoms.

    It calculates spatial overlap between surface cavities saved as PQR files with
    FIL pseudoatoms, as generated by :func:`calcSurfaceCavities` or
    :func:`calcSurfaceCavitiesMultipleFrames`.

    Results are normalized within [0, 1], where the value corresponds to the
    fraction of analyzed PQR files contributing to a given spatial region.

    :arg resolution: surface sampling resolution. Default is 0.5.
    :type resolution: float
    
    :arg max_proc: Maximum number of parallel processes used to voxelize individual
        PQR files. If 1, files are processed serially. If None, all available CPU
        cores are used. Default is 2.
    :type max_proc: int or None

    :arg output_file_name: name of the output PDB file with overlapping cavity
        regions. Default is ``'surface_cavity_overlap_regions.pdb'``.
    :type output_file_name: str

    :arg pqr_files: PQR files with surface cavities represented as FIL atoms.
        If not provided, all PQR files from the current directory will be analyzed.
        A list of PQR files can also be provided.
    :type pqr_files: bool or list """

    kwargs.setdefault('output_file_name', 'surface_cavity_overlap_regions.pdb')

    return calcChannelSurfaceOverlaps(**kwargs)


def calcSurfaceCavities(atoms, output_path=None, surf_radius=4.5, inner_radius=2.0, min_depth=1.5,
                        max_depth=2.5, min_tetrahedra=None, max_tetrahedra=None,
                        min_volume=50, max_volume=None, sparsity=None,
                        separate=False):
    """Calculate surface cavities (pockets) on protein surface using CaviTracer 
    approach.

    :arg atoms: An object representing the molecular structure, typically 
        containing atomic coordinates and element types.
    :type atoms: `Atoms` object

    :arg output_path: Optional path to save the resulting cavities and
        associated data in PQR (or PDB) format. If None, results are not saved.
         Default is None.

        A cavity is written as one FIL pseudoatom per Voronoi vertex, preceded by
        a REMARK reporting its volume, depth and tetrahedron count. The cloud
        marks where the cavity is; every atom carries the same marker radius, so
        the extent it draws is not the extent of the cavity and the radius column
        is not a measurement. The REMARK is where the size is.
    :type output_path: str or None

    :arg separate: If True, each detected cavity is saved to a separate PQR
        file. If False, all cavities are saved in a single PQR file. Default is
        False.
    :type separate: bool

    :arg surf_radius: The first radius threshold used during the deletion of simplices, 
        which is used to define the outer surface of the cavities. Default is 4.5.
    :type surf_radius: float

    :arg inner_radius: The second radius threshold used to define the inner surface of 
        the cavities. Default is 2.
    :type inner_radius: float

    :arg min_depth: The minimum depth, in Angstrom, a cavity must reach to be
        considered. Depth is the geodesic distance from the surface opening along 
        the Voronoi network, a physical length independent of tessellation density.
        Default is 1.5.
    :type min_depth: float

    :arg max_depth: Maximum cavity depth, in Angstrom. Portions of a cavity deeper
        than this value are trimmed away, keeping the shallow surface shell that
        defines a pocket. Default is 2.5.
    :type max_depth: float

    :arg sparsity: Deprecated and ignored; accepted only so that existing calls
        keep working. It never affected surface cavities. In :func:`calcChannels`
        it is the smallest separation at which two channel openings still count
        as separate, applied when the finished channels are deduplicated, and no
        cavity property reads it. Cavity extent, depth, volume and filtering are
        all derived from the exit tetrahedra, so passing 1 or 15 returns the same
        cavities.
    :type sparsity: int

    :arg min_tetrahedra: Minimum number of tetrahedra required for a cavity to
        be retained. Smaller cavities are discarded. Default is None.
    :type min_tetrahedra: int

    :arg max_tetrahedra: Maximum number of tetrahedra allowed for a cavity to 
        be retained. Larger cavities are discarded. Default is None.
    :type max_tetrahedra: int

    :arg min_volume: Minimum volume, in cubic Angstrom, required for a cavity to
        be retained. Default is 50.

        The volume is the cavity's Delaunay tetrahedra summed. Their corners are
        atom centres, so they lie against the wall of the pocket rather than
        filling it, and the number is not the room inside. Compare a threshold
        against other cavities, then, not against a volume measured some other
        way - and note that a channel volume is on another scale entirely, the
        probe swept along the centerline.
    :type min_volume: float

    :arg max_volume: Maximum volume allowed for a cavity to be retained, in cubic
        Angstrom and on the same scale as ``min_volume``. Default is None.
    :type max_volume: float

    :returns: A tuple containing two elements:
        - `cavities`: A list of detected cavities, where each channel is an 
            object containing information about its path and geometry.
        - `surface`: A list containing additional information for further 
            visualization, including the atomic coordinates, simplices defining
            the surface, and merged cavities.
    :rtype: tuple (list, list)

    This function performs the following steps:
    1. **Selection and Filtering:** Selects non-hetero atoms from the protein, 
        calculates van der Waals radii, and performs 3D Delaunay triangulation 
        and Voronoi tessellation on the coordinates.
    2. **Surface and Interior Filtering:** Iteratively removes simplices based 
        on the user-defined radii (`surf_radius` and `inner_radius`) to distinguish the molecular 
        surface from the internal void space.
    3. **Surface Cavity Identification:** Detects connected void regions and 
        identifies those that remain connected to the protein surface, 
        corresponding to surface-accessible cavities and pockets.
    4. **Depth Calculation and Filtering:** Estimates cavity depth using a 
        graph-based traversal from the cavity openings, identifies the deepest
         tetrahedra, and filters cavities according to the specified depth criteria.
    5. **Output Generation:** Optionally trims cavities exceeding the specified
    	maximum depth, saves detected cavities to PDB/PQR files, and returns 
        cavity objects together with the surface representation for further 
        analysis and visualization.
       
    Example usage:
    p = parsePDB('1tqn')
    protein = p.select('protein')
    cavities, surface = calcSurfaceCavities(protein, output_path='test_surf_cav.pqr')   """

    if sparsity is not None:
        _warn("sparsity is deprecated in calcSurfaceCavities and is "
              "ignored. It separates the openings of finished channels in "
              "calcChannels; cavities are built from the exit tetrahedra "
              "before that, so it never changed them.")

    # No peel (min_enclosure=0). The enclosure peel strips the shell of true
    # exterior that a large surf_radius probe bridges over instead of entering, because it
    # offers a channel wide, low-cost routes along the outside of the protein. A
    # surface cavity *is* that shell: a pocket is shallow and open by definition,
    # so the peel deletes these cavities 
    cavities, surface = calcChannels(
            atoms,
            output_path=output_path,
            separate=separate,
            surf_radius=surf_radius, inner_radius=inner_radius,
            min_depth=min_depth, max_depth=max_depth,
            min_volume=min_volume, max_volume=max_volume,
            min_tetrahedra=min_tetrahedra, max_tetrahedra=max_tetrahedra,
            min_enclosure=0.0, cavities_only=True)
    
    return cavities, surface

def scanChannelParameters(atoms, inner_radius_values=(1.2, 1.4, 1.6),
    sparsity_values=(2.0, 6.0, 10.0), min_depth_values=(3.0, 5.0, 10.0),
    output_path='channel_parameter_grid', resolution=0.5, max_proc=2,
    start_point=None, **kwargs):
    """Calculate channels over a combination grid of parameters.

    This function evaluates every combination of ``inner_radius``, ``sparsity``, and
    ``min_depth`` for one molecular structure. All channels obtained for one
    parameter combination are saved together in one PQR file, so every grid
    point contributes one equally weighted result to the final spatial
    occupancy map. The occupancy value written by
    :func:`calcChannelSurfaceOverlaps` is the fraction of grid combinations 
    in which a voxel belongs to at least one predicted channel.

    :arg atoms: Atomic structure analyzed with :func:`calcChannels`.
    :type atoms: :class:`.Atomic`

    :arg inner_radius_values: Probe radii defining the internal void space. Default is
        ``(1.2, 1.4, 1.6)``.
    :type inner_radius_values: float or sequence of float

    :arg sparsity_values: Mouth-separation values, in Angstrom, tested in the
        grid. Default is ``(2.0, 6.0, 10.0)``.
    :type sparsity_values: float or sequence of float

    :arg min_depth_values: Minimum cavity depths tested in the grid. Default is
        ``(3.0, 5.0, 10.0)``.
    :type min_depth_values: float or sequence of float

    :arg output_path: Directory in which individual PQR files, summaries, and
        the final occupancy map are saved. Default is
        ``'channel_parameter_grid'``.
    :type output_path: str or pathlib.Path

    :arg resolution: Voxel resolution used by
        :func:`calcChannelSurfaceOverlaps`. Default is 0.5 Angstrom.
    :type resolution: float

    :arg max_proc: Maximum number of processes used for overlap calculation.
        Default is 2; ``None`` uses all available CPU cores.
    :type max_proc: int or None

    :arg start_point: Optional starting point or atomic selection passed to
        :func:`calcChannels` for every grid combination.
    :type start_point: array-like, :class:`.Atomic`, or None

    :arg kwargs: Additional parameters passed unchanged to :func:`calcChannels`.
        Grid-controlled parameters ``inner_radius``, ``sparsity``, and ``min_depth`` must
        not be supplied here.
    :type kwargs: dict

    :returns: Lists of channels and parameter dictionaries in matching order,
        followed by the path to the occupancy PDB file.
    :rtype: tuple

    Example usage:
    channels_all, parameter_sets, occupancy_file = scanChannelParameters(
        protein, inner_radius_values=[1.2, 1.4, 1.6], sparsity_values=[2, 6, 10],
        min_depth_values=[3, 5, 10], output_path='channel_parameter_grid') """

    from itertools import product
    
    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    if not isinstance(atoms, Atomic):
        raise TypeError("atoms must be a ProDy Atomic object")

    def prepareValues(values, name, allow_zero=False):
        # Checking the input data format
        if np.isscalar(values):
            values = [values]
        try:
            values = [float(value) for value in values]
        except (TypeError, ValueError):
            raise TypeError("{0} must be a number or a sequence of numbers".format(name))
        
        if len(values) == 0:
            raise ValueError("{0} must contain at least one value".format(name))
        if not np.all(np.isfinite(values)):
            raise ValueError("{0} must contain finite values".format(name))
        if any(value < 0 if allow_zero else value <= 0 for value in values):
            relation = "non-negative" if allow_zero else "greater than zero"
            raise ValueError("{0} values must be {1}".format(name, relation))
        
        return list(dict.fromkeys(values))

    forbidden_params = sorted(set(kwargs).intersection(
        {'inner_radius', 'sparsity', 'min_depth', 'output_path', 
         'separate', 'cavities_only', 'return_details'}))
    if forbidden_params:
        raise ValueError("Grid-controlled arguments must not be passed in kwargs: {0}".format(
            ', '.join(forbidden_params)))

    # Parameters for checkup
    inner_radius_values = prepareValues(inner_radius_values, 'inner_radius_values')
    sparsity_values = prepareValues(sparsity_values, 'sparsity_values', allow_zero=True)
    min_depth_values = prepareValues(min_depth_values, 'min_depth_values', allow_zero=True)
    
    resolution = float(resolution)
    if resolution <= 0:
        raise ValueError("resolution must be greater than zero")

    output_dir = Path(output_path)
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError("output_path must be a directory")
    output_dir.mkdir(parents=True, exist_ok=True)

    parameter_grid = list(product(inner_radius_values, sparsity_values, min_depth_values))
    channels_all = []
    parameter_sets = []
    pqr_files = []
    summary_file = output_dir / 'channel_parameter_grid_summary.txt'
    details_file = output_dir / 'channel_parameter_grid_channels.txt'
    occupancy_file = output_dir / 'channel_parameter_occupancy.pdb'

    LOGGER.timeit('_prody_scanChannelParameters')
    LOGGER.info("Calculating channels for {0} parameter combinations.".format(len(parameter_grid)))

    with open(str(summary_file), 'w') as summary, open(str(details_file), 'w') as details:
        summary.write("# Run inner_radius [Å] sparsity [Å] min_depth [Å] Number_of_channels PQR_file\n")
        details.write("# Run Channel_id inner_radius [Å] sparsity [Å] min_depth [Å] Length [Å] Bottleneck [Å] Volume [Å^3] Curvature Cost\n")

        for run_index, (inner_radius, sparsity, min_depth) in enumerate(parameter_grid):
            tag = "run{0:03d}_inner_radius_{1}_sparsity_{2}_depth_{3}".format(
                run_index, *("{0:g}".format(value).replace('-', 'm').replace('.', 'p')
                             for value in (inner_radius, sparsity, min_depth)))
            pqr_file = output_dir / ('channels_' + tag + '.pqr')

            LOGGER.info("Grid run {0}/{1}: inner_radius={2:g}, sparsity={3:g}, min_depth={4:g}".format(
                run_index + 1, len(parameter_grid), inner_radius, sparsity, min_depth))

            channels, _ = calcChannels(atoms, output_path=str(pqr_file), separate=False,
                start_point=start_point, inner_radius=inner_radius, sparsity=sparsity,
                min_depth=min_depth, **kwargs)

            params = {'run': run_index, 'inner_radius': inner_radius, 'sparsity': sparsity,
                      'min_depth': min_depth, 'pqr_file': str(pqr_file)}

            channels_all.append(channels)
            parameter_sets.append(params)
            pqr_files.append(str(pqr_file))

            summary.write("{0} {1:.3f} {2:.3f} {3:.3f} {4} {5}\n".format(
                run_index, inner_radius, sparsity, min_depth, len(channels), pqr_file.name))

            for channel_index, channel in enumerate(channels):
                curvature = channel.curvature if np.isfinite(channel.curvature) else float('nan')
                cost = channel.cost if channel.cost is not None else float('nan')
                details.write("{0} {1} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f} {8:.3f} {9:.6g}\n".format(
                    run_index, channel_index, inner_radius, sparsity, min_depth,
                    channel.length, channel.bottleneck, channel.volume,
                    curvature, cost))

    calcChannelSurfaceOverlaps(pqr_files=pqr_files,
                               output_file_name=str(occupancy_file),
                               resolution=resolution,
                               max_proc=max_proc)
    LOGGER.report('Channel parameters scan completed in %.2fs.', '_prody_scanChannelParameters')

    return channels_all, parameter_sets, str(occupancy_file)


def scanSurfaceCavityParameters(atoms, surf_radius_values=(4.0, 4.5, 5.0),
    inner_radius_values=(1.5, 2.0), min_depth_values=(1.5, 2.0),
    max_depth_values=(2.5, 3.0), min_volume_values=(None, 50),
    output_path='surface_cavity_parameter_grid', resolution=0.5, max_proc=2, **kwargs):
    """Calculate surface cavities over a combination grid of parameters.

    This function evaluates every combination of ``surf_radius``,
    ``inner_radius``, ``min_depth``, ``max_depth``, and ``min_volume`` for one
    molecular structure. All surface cavities obtained for one parameter
    combination are saved together in one PQR file, so every grid point
    contributes one equally weighted result to the final spatial occupancy map 
    (obtained using :func:`calcSurfaceCavityOverlaps`).
    The occupancy value written by :func:`calcSurfaceCavityOverlaps` is the
    fraction of grid combinations in which a voxel belongs to at least one
    predicted surface cavity.

    :arg atoms: Atomic structure analyzed with :func:`calcSurfaceCavities`.
    :type atoms: :class:`.Atomic`

    :arg surf_radius_values: Probe radii used to define the outer molecular
        surface. Default is ``(4.0, 4.5, 5.0)``.
    :type surf_radius_values: float or sequence of float

    :arg inner_radius_values: Probe radii used to define the inner accessible
        cavity space. Default is ``(1.5, 2.0)``.
    :type inner_radius_values: float or sequence of float

    :arg min_depth_values: Minimum cavity depths tested in the grid. Default is
        ``(1.5, 2.0)``.
    :type min_depth_values: float or sequence of float

    :arg max_depth_values: Maximum cavity depths tested in the grid. Portions
        deeper than this value are trimmed. Default is ``(2.5, 3.0)``.
    :type max_depth_values: float or sequence of float

    :arg min_volume_values: Minimum cavity volumes tested in the grid. ``None``
        disables volume filtering for a given run. Default is ``(None, 50)``.
    :type min_volume_values: float, None, or sequence

    :arg output_path: Directory in which individual PQR files, summaries, and
        the final occupancy map are saved. Default is
        ``'surface_cavity_parameter_grid'``.
    :type output_path: str or pathlib.Path

    :arg resolution: Voxel resolution used by
        :func:`calcSurfaceCavityOverlaps`. Default is 0.5 Angstrom.
    :type resolution: float

    :arg max_proc: Maximum number of processes used for overlap calculation.
        Default is 2; ``None`` uses all available CPU cores.
    :type max_proc: int or None

    :arg kwargs: Additional parameters passed unchanged to
        :func:`calcSurfaceCavities`. Grid-controlled parameters
        ``surf_radius``, ``inner_radius``, ``min_depth``, ``max_depth``,
        ``min_volume``, ``output_path`` and ``separate`` must not be supplied
        here.
    :type kwargs: dict

    :returns: Lists of surface cavities and parameter dictionaries in matching
        order, followed by the path to the occupancy PDB file.
    :rtype: tuple

    Example usage:
    cavities_all, parameter_sets, occupancy_file = scanSurfaceCavityParameters(
        protein, surf_radius_values=[4.0, 4.5, 5.0],
        inner_radius_values=[1.5, 2.0], min_depth_values=[1.5, 2.0],
        max_depth_values=[2.5, 3.0], min_volume_values=[None, 50],
        output_path='surface_cavity_parameter_grid')  """

    from itertools import product

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    if not isinstance(atoms, Atomic):
        raise TypeError("atoms must be a ProDy Atomic object")

    def prepareValues(values, name, allow_zero=False, allow_none=False):
        if values is None:
            values = [None]
        elif np.isscalar(values):
            values = [values]

        checked = []
        for value in values:
            if value is None:
                if allow_none:
                    checked.append(None)
                    continue
                else:
                    raise ValueError("{0} values must not contain None".format(name))

            try:
                value = float(value)
            except (TypeError, ValueError):
                raise TypeError("{0} must be a number or a sequence of numbers".format(name))

            if not np.isfinite(value):
                raise ValueError("{0} must contain finite values".format(name))

            if value < 0 if allow_zero else value <= 0:
                relation = "non-negative" if allow_zero else "greater than zero"
                raise ValueError("{0} values must be {1}".format(name, relation))

            checked.append(value)

        if len(checked) == 0:
            raise ValueError("{0} must contain at least one value".format(name))

        unique_values = []
        for value in checked:
            if value not in unique_values:
                unique_values.append(value)

        return unique_values

    forbidden_params = sorted(set(kwargs).intersection(
        {'surf_radius', 'inner_radius', 'min_depth', 'max_depth',
         'min_volume', 'output_path', 'separate', 'cavities_only',
         'sparsity'}))

    if forbidden_params:
        raise ValueError(
        "Grid-controlled arguments must not be passed in kwargs: {0}".format(
        ', '.join(forbidden_params)))

    surf_radius_values = prepareValues(surf_radius_values, 'surf_radius_values')
    inner_radius_values = prepareValues(inner_radius_values, 'inner_radius_values')
    min_depth_values = prepareValues(min_depth_values, 'min_depth_values', allow_zero=True)
    max_depth_values = prepareValues(max_depth_values, 'max_depth_values', allow_zero=True)
    min_volume_values = prepareValues(min_volume_values, 'min_volume_values', allow_zero=True, allow_none=True)

    resolution = float(resolution)
    if resolution <= 0:
        raise ValueError("resolution must be greater than zero")

    output_dir = Path(output_path)
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError("output_path must be a directory")
    output_dir.mkdir(parents=True, exist_ok=True)

    parameter_grid = list(product(surf_radius_values, inner_radius_values,
                                  min_depth_values, max_depth_values,
                                  min_volume_values))

    cavities_all = []
    parameter_sets = []
    pqr_files = []
    summary_file = output_dir / 'surface_cavity_param_grid_summary.txt'
    details_file = output_dir / 'surface_cavity_param_grid_cavities.txt'
    occupancy_file = output_dir / 'surface_cavity_param_occupancy.pdb'

    LOGGER.timeit('_prody_scanSurfaceCavityParameters')
    LOGGER.info("Calculating surface cavities for {0} parameter combinations.".format(
        len(parameter_grid)))

    with open(str(summary_file), 'w') as summary, open(str(details_file), 'w') as details:
        summary.write("# Run surf_radius [A] inner_radius [A] min_depth [A] max_depth [A] min_volume [A^3] Number_of_cavities PQR_file\n")
        details.write("# Run Cavity_id surf_radius [A] inner_radius [A] min_depth [A] max_depth [A] min_volume [A^3] Volume [A^3] Depth [A] Tetrahedra_count\n")

        for run_index, (surf_radius, inner_radius, min_depth, max_depth,
                        min_volume) in enumerate(parameter_grid):

            values_for_tag = (surf_radius, inner_radius, min_depth, max_depth,
                              'none' if min_volume is None else min_volume)

            tag = "run{0:03d}_surf_radius_{1}_inner_radius_{2}_mindepth_{3}_maxdepth_{4}_minvol_{5}".format(
                run_index, *("{0:g}".format(value).replace('-', 'm').replace('.', 'p')
                             if isinstance(value, float) else str(value)
                             for value in values_for_tag))

            pqr_file = output_dir / ('surface_cavities_' + tag + '.pqr')

            LOGGER.info(
                "Grid run {0}/{1}: surf_radius={2:g}, inner_radius={3:g}, "
                "min_depth={4:g}, max_depth={5:g}, min_volume={6}".format(
                    run_index + 1, len(parameter_grid), surf_radius,
                    inner_radius, min_depth, max_depth, min_volume))

            cavities, _ = calcSurfaceCavities(atoms, output_path=str(pqr_file), separate=False,
                surf_radius=surf_radius, inner_radius=inner_radius, min_depth=min_depth,
                max_depth=max_depth, min_volume=min_volume, **kwargs)

            params = {'run': run_index,
                      'surf_radius': surf_radius,
                      'inner_radius': inner_radius,
                      'min_depth': min_depth,
                      'max_depth': max_depth,
                      'min_volume': min_volume,
                      'pqr_file': str(pqr_file)}

            cavities_all.append(cavities)
            parameter_sets.append(params)
            pqr_files.append(str(pqr_file))

            min_volume_text = 'None' if min_volume is None else "{0:.3f}".format(min_volume)

            summary.write("{0} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5} {6} {7}\n".format(
                run_index, surf_radius, inner_radius, min_depth, max_depth,
                min_volume_text, len(cavities), pqr_file.name))

            for cavity_index, cavity in enumerate(cavities):
                volume = cavity.volume
                depth = cavity.depth
                tetrahedra_count = len(cavity.tetrahedra)

                details.write("{0} {1} {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6} {7:.3f} {8:.3f} {9}\n".format(
                    run_index, cavity_index, surf_radius, inner_radius,
                    min_depth, max_depth, min_volume_text, volume, depth,
                    tetrahedra_count))

    calcSurfaceCavityOverlaps(pqr_files=pqr_files,
                              output_file_name=str(occupancy_file),
                              resolution=resolution,
                              max_proc=max_proc)

    LOGGER.report('Surface cavity parameters scan completed in %.2fs.',
                  '_prody_scanSurfaceCavityParameters')

    return cavities_all, parameter_sets, str(occupancy_file)


def calcFrequentObjectResidues(residues_all, count_residue_names=False, 
                        count_once_per_frame=True, output_file_name=None):
    """Count residues lining channels, pores, or surface cavities by chain.

    This function analyzes the output returned by:
    - getChannelResidueNamesMultipleFrames()
    - getPoreResidueNamesMultipleFrames()
    - getSurfaceCavityResidueNamesMultipleFrames()

    :arg residues_all: Residue lists returned by one of the multiple-frame
        residue-reporting functions. The expected input is a list of frame/model
        entries, where each entry contains strings such as
        ``'channel0: ASP108:A, LYS245:A'`` or ``'cavity1: GLY20:B, SER55:B'``.
    :type residues_all: list

    :arg count_residue_names: If **False**, individual residues are counted,
        for example ``ASP108`` or ``LYS245``. If **True**, residue names are
        counted instead, for example ``ASP`` or ``LYS``. Default is **False**.
    :type count_residue_names: bool

    :arg count_once_per_frame: If **True**, the same residue is counted only
        once per frame/model, even if it appears in multiple objects in that
        frame. If **False**, every occurrence is counted. Default is **True**.
    :type count_once_per_frame: bool

    :arg output_file_name: Optional base name for saving residue counts to a text
        file. The suffix ``'_Residue_counts.txt'`` will be added. If
        **None**, no file is written. Default is **None**.
    :type output_file_name: str or None

    Examples:
    residues_all = getChannelResidueNamesMultipleFrames(protein, channels_all, trajectory=dcd)

    counts = countObjectResiduesByChain(residues_all)

    counts = countObjectResiduesByChain(residues_all, count_residue_names=True,
        output_file_name='channel_res_counts') """

    from collections import Counter, defaultdict

    try:
        string_types = (basestring,)
    except NameError:
        string_types = (str,)

    def asFrameList(residues):
        if isinstance(residues, string_types):
            return [[residues]]
        if len(residues) == 0:
            return []
        if all(isinstance(item, string_types) for item in residues):
            return [residues]
        return residues

    def residueNameFromLabel(label):
        name = []
        for char in label:
            if char.isdigit() or char in ('-', '+'):
                break
            name.append(char)
        return ''.join(name) if name else label

    counts_by_chain = defaultdict(Counter)
    frames = asFrameList(residues_all)

    for frame in frames:
        frame_seen = set()

        for object_line in frame:
            if object_line is None:
                continue

            if ': ' in object_line:
                residues_part = object_line.split(': ', 1)[1]
            else:
                residues_part = object_line

            if residues_part == 'None':
                continue

            residues = [res.strip() for res in residues_part.split(',')]

            for residue in residues:
                if residue == '' or residue == 'None':
                    continue

                if ':' in residue:
                    residue_id, chain = residue.rsplit(':', 1)
                else:
                    residue_id = residue
                    chain = ''

                if count_residue_names:
                    key = residueNameFromLabel(residue_id)
                else:
                    key = residue_id

                if count_once_per_frame:
                    frame_seen.add((chain, key))
                else:
                    counts_by_chain[chain][key] += 1

        if count_once_per_frame:
            for chain, key in frame_seen:
                counts_by_chain[chain][key] += 1

    counts_by_chain = dict(counts_by_chain)

    if output_file_name is not None:
        output_file = output_file_name + '_ResCounts.txt'

        with open(output_file, 'w') as out:
            out.write('# Chain Residue Count\n')

            for chain in sorted(counts_by_chain):
                chain_label = chain if chain else 'no_chain'

                for residue, count in counts_by_chain[chain].most_common():
                    out.write('{0} {1} {2} \n'.format(chain_label, residue, count))

        LOGGER.info("Residue counts by chain were saved to: {0}".format(output_file))

    return counts_by_chain


def showFrequentObjectResidues(counts_by_chain, top=50):
    """Show residue-count bar plots.

    This function visualizes residue counts returned by
    :func:`calcFrequentObjectResidues`. Counts can be shown for one selected
    chain or, if ``chain`` is **None**, separately for all chains. Residues are
    displayed on the x-axis and their counts on the y-axis.

    :arg counts_by_chain: Dictionary returned by
        :func:`countObjectResiduesByChain`, with chain identifiers as keys and
        :class:`collections.Counter` objects as values.
    :type counts_by_chain: dict

    :arg chain: Chain identifier to plot. If **None**, residue counts are plotted
        separately for all chains. Default is **None**.
    :type chain: str or None

    :arg top: Maximum number of most frequent residues to display. If **None**,
        all residues are shown. Default is 50.
    :type top: int or None

    Example usage:
    import matplotlib.pyplot as plt
    counts = countObjectResiduesByChain(residues_all)
    showObjectResidueCountsByChain(counts, chain='A', top=None)
    plt.show()
    
    or instead:
    showObjectResidueCountsByChain(counts, top=20) """

    if not isinstance(counts_by_chain, dict):
        raise TypeError("counts_by_chain must be a dictionary returned by "
                    "calcFrequentObjectResidues().")
    
    import matplotlib.pyplot as plt
    
    axes = []
    for chain in sorted(counts_by_chain):
        items = counts_by_chain[chain].most_common()

        if top is not None:
            items = items[:top]

        labels = [item[0] for item in items]
        values = [item[1] for item in items]

        fig, ax = plt.subplots(figsize=(14, 4))
        ax.bar(labels, values)
        ax.set_xlabel("Residue")
        ax.set_ylabel("Count")
        ax.set_title("Residue counts for chain {0}".format(chain if chain else 'no_chain'))
        ax.tick_params(axis='x', rotation=90)
        fig.tight_layout()
        axes.append(ax)
        
    return axes[0] if len(axes) == 1 else axes

# The dictionary the written file declares conformance to, quoted from the schema's
# own README so that a reader can check the file against the same document we wrote
# it from. https://github.com/sb-ncbr/tunnels-schema
_CIF_DICT_NAME = 'mmcif_tunnels.dic'
_CIF_DICT_VERSION = '1.0'
_CIF_DICT_LOCATION = \
    'https://sb-ncbr.github.io/tunnels-schema/schemas/mmcif_tunnels_v10.dic'

# The schema's vocabulary for _sb_ncbr_channel.type against this module's. The
# schema names no enumeration - the item is free text and its description says
# "Pore, Path, etc." - so these are MOLE's words for the same three things, chosen
# because a file nobody else's reader recognises is not worth the format.
_CIF_OBJECT_TYPES = {'channel': 'Tunnel', 'pore': 'Pore', 'link': 'Path'}

def _cifValue(value, precision=3):
    """One data item, formatted for a CIF loop.

    ``None`` and NaN become ``?``, the CIF's own "value not given", which is what
    every unmeasured item in this export writes; the schema marks most items
    mandatory, and a mandatory item is satisfied by an explicit ``?`` but not by a
    guess.

    Text is quoted only where it has to be - whitespace, or a leading character that
    would otherwise start a comment, a data block or a quoted string. Nothing this
    module writes needs quoting today, since residue and chain identifiers are bare
    codes, but a structure with a chain named ``'`` should produce an unreadable file
    rather than a silently wrong one."""

    if value is None:
        return '?'
    if isinstance(value, bool):
        # Ahead of the numeric branch: bool is an int in Python, and the schema's
        # own examples for its boolean items read True/False.
        return 'True' if value else 'False'
    if isinstance(value, float):
        if not np.isfinite(value):
            return '?'
        return '{0:.{1}f}'.format(value, precision)
    if isinstance(value, (int, np.integer)):
        return str(int(value))

    text = str(value)
    if text == '':
        return "''"
    if any(c.isspace() for c in text) or text[0] in '_#$[];\'"':
        return "'{0}'".format(text.replace("'", "\\'"))
    return text


def _cifLoop(out, category, columns, rows, precision=None):
    """Write one ``loop_`` of *rows* under *category*, or nothing when it is empty.

    *columns* are the item names without their category prefix. *precision* is an
    optional per-column dict of float precisions; anything unnamed uses the default.

    An empty loop is skipped rather than written headerless, because a ``loop_`` with
    no rows is a syntax error in CIF and a category with nothing to say is better
    absent - the schema marks every category ``mandatory_code no`` precisely so that
    a producer can leave out what it does not compute."""

    if not rows:
        return

    precision = precision or {}

    out.write('#\n')
    out.write('loop_\n')
    for column in columns:
        out.write('_{0}.{1}\n'.format(category, column))

    for row in rows:
        out.write(' '.join(
            _cifValue(row[column], precision.get(column, 3))
            for column in columns) + '\n')


def _pqrOutputPaths(output_path):
    """``(channels_path, links_path, into_directory, separate_stem)`` for a PQR run.

    A directory names no run, so nothing here is named after one: the files are
    named after what they hold, ``channels.pqr`` beside ``links.pqr``, and the
    per-object files carry no stem either - ``sp0_chl3.pqr`` rather than
    ``output_sp0_chl3.pqr``, a stem that would be the same word for every run and
    so tell the reader nothing.

    Links go beside the channels rather than among them: a link does not reach the
    surface, so a viewer loading the channel file should not find one in it.

    Split out so that the run which writes these files and the check for stale ones
    left by an earlier run resolve them the same way. Two copies of this would
    disagree eventually, and the disagreement would be a warning about the wrong
    file, or no warning at all."""

    from pathlib import Path

    output_path = Path(output_path)

    into_directory = output_path.is_dir()
    separate_stem = None
    if into_directory:
        output_path = output_path / "channels.pqr"
        separate_stem = ''

    elif not (output_path.suffix == ".pdb" or output_path.suffix == ".pqr"):
        output_path = output_path.with_suffix(".pqr")

    links_path = output_path.with_name(
        ('links' if into_directory else output_path.stem + '_links')
        + output_path.suffix)

    return output_path, links_path, into_directory, separate_stem


def _warnStaleOutputs(output_path, output_format='pqr', separate=False,
                      tags=('chl', 'lnk'), cif_path=None):
    """Name the files an earlier run left where this one has just written nothing.

    A run that finds nothing now leaves no file at all, which is right - an empty
    one cannot be told from a write that failed. It does mean an older file at the
    same place survives and goes on looking like this run's output, and that is
    worse than an empty file, because it is wrong rather than merely uninformative.
    So the paths are resolved exactly as a writing run would resolve them, and
    whatever is already there is named.

    The per-object files cannot be listed from a run that produced none, so they are
    globbed by the tag they carry in their names instead."""

    from pathlib import Path

    if not output_path:
        return

    stale = []
    if _isMmcifFormat(output_format, separate=False):
        path = _cifOutputPath(output_path) if cif_path is None else Path(cif_path)
        if path.exists():
            stale.append(path)
    else:
        channels_path, links_path, into_directory, _ = \
            _pqrOutputPaths(output_path)
        stale.extend(p for p in (channels_path, links_path) if p.exists())

        if separate:
            prefix = '' if into_directory else channels_path.stem + '_'
            for tag in tags:
                stale.extend(sorted(channels_path.parent.glob(
                    '{0}*{1}*{2}'.format(prefix, tag, channels_path.suffix))))

    if not stale:
        return

    stale = sorted(set(str(p) for p in stale))
    _warn("Nothing was found, so nothing was written, but {0} file(s) from an "
          "earlier run are still there and now describe no run at all: {1}{2}. "
          "Delete them, or write this run somewhere else.".format(
              len(stale), ', '.join(stale[:6]),
              ', ...' if len(stale) > 6 else ''))


def _cifOutputPath(path):
    """Resolve *path* to the file the mmCIF is written to.

    An existing directory takes ``channels.cif`` inside it, which is what
    :func:`calcChannels` already does with a directory given as ``output_path``.
    There is no ``links.cif`` beside it: links share the file and are told apart by
    ``_sb_ncbr_channel.type``, the column the schema provides for exactly that.

    One trap this exists to keep in one place: a path that does not exist yet is not
    a directory as far as :meth:`~pathlib.Path.is_dir` is concerned, so a suffixless
    ``results`` becomes ``results.cif`` beside the working directory rather than a
    file within it - surprising, but it is what ``output_path`` has always done, and
    two rules would be worse than one."""

    from pathlib import Path

    path = Path(path)
    if path.is_dir():
        return path / 'channels.cif'

    if path.suffix.lower() == '.cif':
        return path
    return path.with_name(path.name + '.cif')


def _isMmcifFormat(output_format, separate=False):
    """Whether *output_format* asks for mmCIF, rejecting anything else outright.

    One flag rather than a second output path, so that a run has one place its
    results go and the format is a property of that place. An unknown value raises
    rather than falling back on PQR: a misspelled format that quietly wrote the old
    one would be found only by opening the file.

    *separate* is ignored for mmCIF, since the value of the format is the joins
    between its categories and a file per channel would cut exactly those. Said once
    here, so that a run asking for both does not leave the reader wondering which
    won."""

    if output_format is None:
        return False

    fmt = str(output_format).lower().lstrip('.')
    if fmt == 'pqr':
        return False
    if fmt not in ('mmcif', 'cif'):
        # 'pdb' is not among them on purpose. The flag picks the format family
        # and the path picks the extension, which is how PDB has always been
        # chosen here - output_path ending in .pdb. A 'pdb' accepted alongside
        # would name an extension it could not promise, since the path would
        # still decide it. It is the one wrong value worth answering rather
        # than only rejecting, being the one a reader would reasonably try.
        raise ValueError(
            "output_format must be 'pqr' or 'mmcif', not {0!r}.{1}".format(
                output_format,
                " For PDB rather than PQR, keep output_format='pqr' and give "
                "output_path a .pdb suffix." if fmt == 'pdb' else ''))

    if separate:
        LOGGER.info('separate is ignored for mmCIF output: every channel goes '
                    'into one file, which is what lets its categories refer to '
                    'one another.')

    return True


def _poreCifPath(output_path):
    """``pores.cif`` where *output_path* is a directory, otherwise unchanged.

    :func:`_cifOutputPath` would call it ``channels.cif``, which is right for a
    channel run and wrong for a pore one told the same folder - the second would
    overwrite the first. The PQR path names them apart for the same reason."""

    from pathlib import Path

    path = Path(output_path)
    return path / 'pores.cif' if path.is_dir() else path


def _cumulativeArcLength(centerline_spline, t):
    """Arc length from the start of the route to each parameter in *t*.

    Measured on the same polyline
    :meth:`~ChannelCalculator.calculateChannelLength` uses - the same ten points per
    knot - so that the last value of the profile's ``distance`` column equals the
    channel's reported length by construction rather than by coincidence. Evaluating
    the integral more accurately here would be worse, not better: it would put a
    length in the file that disagrees with the one the same run reports everywhere
    else."""

    dense = np.linspace(centerline_spline.x[0], centerline_spline.x[-1],
                        max(2, len(centerline_spline.x) * 10))
    points = centerline_spline(dense)
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(steps)))

    return np.interp(t, dense, cumulative)


def _backboneRadiiSource(atoms):
    """``(tree, vdw_radii)`` over backbone atoms, or ``None`` where there are none.

    The schema's ``free_radius`` is the radius a probe would have if only main-chain
    atoms confined it, which is MOLE's way of saying how much room a side chain could
    make by moving. Only the source differs from an ordinary clearance query, so
    :func:`_vertexRadii` does the measuring for both.

    Main-chain atoms being a subset of all of them, the free radius reads as though
    it must be at least the radius, and mostly it is. It is not an invariant, though,
    and the exception is intrinsic rather than a sign of trouble: the spheres are read
    off the splines, not off the Voronoi vertices they were fitted to, so a sample
    overshoots its inscribed sphere by a little and clips the wall. The shortfall is a
    few hundredths of an Angstrom - the same overshoot the ``-0.5`` floor in
    :func:`_liningContacts` is set well above. Nothing is clamped, because a clamp
    would hide the one case that does mean something, a report measured against a
    structure other than the traced one."""

    backbone = atoms.select('backbone')
    if backbone is None or len(backbone) == 0:
        return None

    return _kdTree(backbone), _atomRadii(backbone, warn=False)


def _residueProperties(resnames):
    """The schema's physicochemical items for one set of lining residues.

    *resnames* is one entry per residue, repeated names included where a lining holds
    several residues of a kind. Returns a dict keyed by the schema's own item names.

    The averages are taken over the residues the scales cover, and a lining that
    covers none of them - a channel through a nucleic acid, or between ligands - gets
    ``None`` rather than a mean over nothing. The counts and the charge are taken over
    every residue, since a residue absent from a scale is genuinely uncharged rather
    than unmeasured.

    The averages are taken over lining residues, which is how MOLE defines every one
    of these, and over side-chain values, since that is what the scales measure. MOLE
    additionally gives a residue touching the route only through its backbone a
    per-scale mainchain value rather than its own; that is not reproduced here, being
    undocumented for four of the seven scales and inferable only by fitting. The
    numbers are therefore MOLE's scales computed over this module's lining, and agree
    with MOLEonline's to the extent the two linings agree - not bit for bit."""

    def mean(scale):
        values = [scale[name] for name in resnames if name in scale]
        return float(np.mean(values)) if values else None

    # Everything is computed over the standard residues alone - the ones the scales
    # and the schema's own formulas are defined for - and never over a default
    # substituted for the rest. Calling an unscored residue neutral would put a
    # figure in the file that is untrue rather than merely partial: a route walled
    # by nucleotides is not uncharged, whatever its charged amino acids sum to.
    #
    # What makes the numbers honest is therefore not their arithmetic but the
    # fraction they cover, which the caller records and reports; see
    # _reportUncoveredLining.
    scored = [name for name in resnames if name in _SCALE_RESIDUES]

    charges = [_CHARGED_RESIDUES.get(name, 0) for name in scored]

    # A lining with nothing scorable in it yields no count either, not a zero.
    # Summing an empty set gives 0 and averaging one is undefined, so left alone
    # the two halves of this dict would describe the same empty set differently -
    # `?` for the averages and a confident 0 for the counts. Zero charged amino
    # acids is also the truthful reading of a wall made entirely of nucleotides,
    # and the least useful one.
    return {
        'charge': int(sum(charges)) if scored else None,
        'numPositives': (int(sum(1 for c in charges if c > 0))
                         if scored else None),
        'numNegatives': (int(sum(1 for c in charges if c < 0))
                         if scored else None),
        'hydropathy': mean(_KYTE_DOOLITTLE),
        'hydrophobicity': mean(_CID_HYDROPHOBICITY),
        'polarity': mean(_ZIMMERMAN_POLARITY),
        'ionizable': (int(sum(1 for name in scored
                              if name in _IONIZABLE_RESIDUES))
                      if scored else None),
        'mutability': mean(_JTT_MUTABILITY),
        'logD': mean(_MOLE_LOGD),
        'logP': mean(_MOLE_LOGP),
        'logS': mean(_MOLE_LOGS),
    }


def _sampleLining(atoms, source, centers, radii, distA, deep=None):
    """Which residues line each probe sphere, and which of their atoms did it.

    Returns ``(residues_per_sample, backbone_residues)``: a list with one set of
    residue indices per sphere, and the set of residue indices that touched the route
    through a main-chain atom anywhere along it.

    One query for the whole route rather than one per sphere - :func:`_liningContacts`
    is already batched, and the per-sphere split is a regrouping of its output."""

    probes, candidates, _ = _liningContacts(source, centers, radii, distA, deep)

    per_sample = [set() for _ in range(len(centers))]
    backbone = set()
    if len(candidates) == 0:
        return per_sample, backbone

    resindices = atoms.getResindices()
    names = atoms.getNames()
    bb_names = ('N', 'CA', 'C', 'O')

    for probe, atom in zip(probes.tolist(), candidates.tolist()):
        residue = int(resindices[atom])
        per_sample[probe].add(residue)
        if names[atom] in bb_names:
            backbone.add(residue)

    return per_sample, backbone


def _channelLayers(per_sample):
    """Group consecutive spheres sharing one lining into layers.

    Returns a list of ``(first, last, residues)``, both bounds inclusive.

    This is MOLE's definition of a layer: the route is cut wherever the set of
    residues around it changes, so a layer is a stretch of channel with one wall.
    MOLE additionally merges layers below a length it does not publish; that is not
    reproduced here, because a merge rule guessed at would produce a different
    partition wearing the same name. The consequence is that our layers are shorter
    and more numerous than MOLEonline's for the same structure."""

    layers = []
    if not per_sample:
        return layers

    start = 0
    for i in range(1, len(per_sample) + 1):
        if i == len(per_sample) or per_sample[i] != per_sample[start]:
            layers.append((start, i - 1, per_sample[start]))
            start = i

    return layers


def _objectCifRows(object_id, obj, type_name, context, num_samples):
    """Every schema row one channel, pore or link contributes.

    *context* bundles what is shared across the objects of one file: the structure
    the lining is measured against, its query sources, the clearance, and the
    residue lookups. Returns a dict keyed by the schema's category names.

    Without a structure only the geometry is produced - the lining categories and
    everything derived from them are simply absent, which the schema permits and a
    guess would not."""

    from collections import Counter

    atoms = context['atoms']
    rows = {'channel': [], 'profile': [], 'props': [], 'layer': [],
            'residue': [], 'layer_residue': [], 'weighted': []}

    centers, radii = _sampleObjectSpheres(obj, num_samples)
    if len(centers) == 0:
        return rows

    t = _objectSampleParameters(obj, num_samples)
    distance = _cumulativeArcLength(obj.centerline_spline, t)
    total = float(distance[-1])
    # A route of no length has no normalized position along it; 0 everywhere keeps
    # the column present and the key unique, there being only one sphere to key.
    fraction = distance / total if total > 0 else np.zeros_like(distance)

    free = (_vertexRadii(centers, context['backbone'])
            if context['backbone'] is not None else None)

    rows['channel'].append({
        'id': object_id,
        'type': type_name,
        'method': 'CaviTracer',
        'software': context['software'],
        'auto': context['auto'],
        # The item is typed int and described as a count of the structure's
        # cavities, yet it sits on a per-channel row, where only an index makes
        # sense. Read as an index, and filled with the search site (sp<n>) rather
        # than with the index of the Delaunay cavity the site lies in.
        #
        # Deliberately, because the site is this module's operative void: chambers
        # are carved at a probe radius of their own, and each seeded chamber is
        # searched independently, so two chambers of one cavity are two separate
        # sub-cavities here and their channels have no more in common than any
        # other pair. The cavity index would group them; that grouping is the one
        # the trace does not use.
        #
        # The consequence to know is that channels sharing a Delaunay cavity can
        # carry different values, two chambers of one cavity being two sites. The
        # log's site table gives the cavity for every site, and the
        # schema has nowhere to record the chamber, its `cavity` being an int.
        'cavity': 0 if obj.origin is None else int(obj.origin),
    })

    if atoms is None:
        for i in range(len(centers)):
            rows['profile'].append({
                'channel_id': object_id, 'radius': float(radii[i]),
                'free_radius': None, 'distance': float(distance[i]),
                'T': float(fraction[i]), 'x': float(centers[i, 0]),
                'y': float(centers[i, 1]), 'z': float(centers[i, 2]),
                'charge': None})
        return rows

    per_sample, backbone_residues = _sampleLining(
        atoms, context['lining'], centers, radii, context['distA'],
        context['deep'])

    allowed = context['allowed']
    per_sample = [residues & allowed for residues in per_sample]

    resname_of = context['resname_of']

    for i in range(len(centers)):
        names = [resname_of[r] for r in per_sample[i]]
        rows['profile'].append({
            'channel_id': object_id,
            'radius': float(radii[i]),
            'free_radius': None if free is None else float(free[i]),
            'distance': float(distance[i]),
            'T': float(fraction[i]),
            'x': float(centers[i, 0]),
            'y': float(centers[i, 1]),
            'z': float(centers[i, 2]),
            'charge': int(sum(_CHARGED_RESIDUES.get(n, 0) for n in names)),
        })

    # Residues in the order they are first met along the route, which is what the
    # schema's `order` means and what `layer_residue.residue_id` points back to.
    first_seen = {}
    for i, residues in enumerate(per_sample):
        for residue in residues:
            first_seen.setdefault(residue, i)

    ordered = sorted(first_seen, key=lambda r: (first_seen[r], r))
    residue_order = {}
    for position, residue in enumerate(ordered, start=1):
        residue_order[residue] = position
        rows['residue'].append({
            'channel_id': object_id,
            'order': position,
            'sequence_number': int(context['resnum_of'][residue]),
            'chain_id': context['chid_of'][residue] or '.',
            'backbone': residue in backbone_residues,
        })

    layers = _channelLayers(per_sample)
    narrowest = int(np.argmin(radii))
    layer_minima = [float(radii[first:last + 1].min())
                    for first, last, _ in layers]

    # Kept beside the rows rather than read back out of them: the weighted props
    # average logD, logP and logS over the layers, and sb_ncbr_channel_layer has no
    # column for those three, so taking the weighting input from the written rows
    # would silently weight nothing and write `?`.
    layer_props = []

    for index, (first, last, residues) in enumerate(layers):
        props = _residueProperties([resname_of[r] for r in residues])
        layer_props.append(props)
        order = index + 1

        local_minimum = (
            (index == 0 or layer_minima[index - 1] > layer_minima[index]) and
            (index == len(layers) - 1 or
             layer_minima[index + 1] > layer_minima[index]))

        row = {
            'channel_id': object_id,
            'order': order,
            'min_radius': layer_minima[index],
            'min_free_radius': (None if free is None
                                else float(free[first:last + 1].min())),
            'start_distance': float(distance[first]),
            'end_distance': float(distance[last]),
            'local_minimum': bool(local_minimum),
            'bottleneck': first <= narrowest <= last,
        }
        for item in ('charge', 'numPositives', 'numNegatives',
                     'hydrophobicity', 'hydropathy', 'polarity', 'mutability'):
            row[item] = props[item]
        rows['layer'].append(row)

        for position, residue in enumerate(
                sorted(residues, key=lambda r: residue_order[r]), start=1):
            rows['layer_residue'].append({
                'channel_id': object_id,
                'layer_id': order,
                'order': position,
                'residue_id': residue_order[residue],
            })

    lining = [resname_of[r] for r in ordered]

    # Counted once per object, on the whole lining, rather than per layer, where
    # the same residue would be tallied for every layer it touches.
    #
    # The schema defines every one of these items over amino acids - "a sum of
    # charged amino acid residues", "amino acid polarities" - so a nucleotide, a
    # cofactor or a modified residue has no value under it, and MOLE excludes them
    # too. Nothing here departs from that. What is worth saying out loud is how much
    # of the wall it leaves out, because the residues left out are rarely a random
    # sample of it: a D-peptide loses its D-residues, a ribosomal route its
    # nucleotides, and either can be half of a lining.
    missing = [name for name in lining if name not in _SCALE_RESIDUES]
    if missing and lining:
        context['uncovered'].update(missing)
        context['coverage'][object_id] = (len(lining) - len(missing), len(lining),
                                          Counter(missing))

    props = _residueProperties(lining)
    props['channel_id'] = object_id
    props['bRadius'] = _bRadius(atoms, context, layers, layer_minima, narrowest)
    rows['props'].append(props)

    weights = [row['end_distance'] - row['start_distance']
               for row in rows['layer']]
    weighted = {'channel_id': object_id}
    for item in ('hydropathy', 'hydrophobicity', 'mutability', 'polarity',
                 'logD', 'logP', 'logS'):
        weighted[item] = _weightedMean(
            [props[item] for props in layer_props], weights)
    rows['weighted'].append(weighted)

    return rows


def _weightedMean(values, weights):
    """Mean of *values* weighted by *weights*, ignoring the entries that are ``None``.

    A layer whose lining no scale covers contributes nothing rather than a zero, and
    an item no layer could measure stays ``None`` all the way up."""

    pairs = [(v, w) for v, w in zip(values, weights)
             if v is not None and w > 0]
    if not pairs:
        return None

    total = sum(w for _, w in pairs)
    if total <= 0:
        return None

    return float(sum(v * w for v, w in pairs) / total)


def _bRadius(atoms, context, layers, layer_minima, narrowest):
    """``min_radius + RMSF`` at the bottleneck, or ``None`` without B-factors.

    The schema describes the item as a radius widened by the mobility of the residues
    around it, with the mobility read off the B-factors as
    ``RMSF = sqrt(3B / 8 pi^2)``. A structure that carries no B-factors - an MD frame,
    or a model - gets ``None``; the schema's own example value for this item is
    ``null``, so an absent value is expected here rather than exceptional."""

    betas = context['betas']
    if betas is None or not layers:
        return None

    for index, (first, last, residues) in enumerate(layers):
        if not (first <= narrowest <= last):
            continue
        if not residues:
            return None
        mask = np.isin(atoms.getResindices(), list(residues))
        values = betas[mask]
        if len(values) == 0:
            return None
        rmsf = np.sqrt(3.0 * np.maximum(values, 0.0) / (8.0 * np.pi ** 2))
        return float(layer_minima[index] + rmsf.mean())

    return None


def writeChannelsCIF(filename, channels, atoms=None, structure=None, links=None,
                     object_type='channel', num_samples=5, auto=True,
                     autoext=True, **kwargs):
    """Write *channels* as an mmCIF conforming to the sb-ncbr tunnels schema.

    The schema is the one MOLEonline and ChannelsDB publish their tunnels under
    (https://github.com/sb-ncbr/tunnels-schema, v1.0). Writing it puts the profile,
    the layers and the lining of every channel in one file with keys joining them,
    where this module otherwise spreads them over a PQR and two text tables joined
    only by filename.

    :arg filename: file to write. ``.cif`` is appended when *autoext* and the name
        does not carry it; an existing directory takes ``channels.cif`` within it.
    :type filename: str

    :arg channels: the objects to write, from :func:`calcChannels` or
        :func:`calcPoresFromChannels`.
    :type channels: list

    :arg atoms: structure the lining is measured against. Without it only
        ``sb_ncbr_channel`` and the geometry of ``sb_ncbr_channel_profile`` are
        written, every other category depending on knowing what surrounds the route.
    :type atoms: :class:`.Atomic`

    :arg structure: structure to write as ``_atom_site`` ahead of the channel
        categories, so that the file opens on its own in a viewer. Needs Biopython,
        as :func:`.writeMMCIF` does the writing; without it the file holds channels
        alone and a reader supplies the structure.
    :type structure: :class:`.Atomic`

    :arg links: chamber links, written into the same file as
        ``_sb_ncbr_channel.type`` ``Path``. Unlike the PQR output, which keeps links
        in a file of their own so that a viewer loading channels does not find one
        among them, the schema has a type column and one file is what a consumer of
        it wants.
    :type links: list

    :arg object_type: what *channels* holds; ``"channel"`` or ``"pore"``.
    :type object_type: str

    :arg num_samples: samples per tetrahedron of the route, as elsewhere in this
        module. The profile is the very spheres the PQR writes, not a second
        sampling of the same splines.
    :type num_samples: int

    :arg auto: whether the start point was chosen automatically. A channel does not
        record this, so :func:`calcChannels` passes ``start_point is None``.
    :type auto: bool

    Accepts ``distA`` (default 1.5) and ``include_water`` as
    :func:`getChannelResidueNames` does.

    Three mappings are this module's reading of the schema rather than the schema's
    own words, since it fixes no vocabulary for them. ``type`` follows MOLE, a
    channel being a ``Tunnel``, a pore a ``Pore`` and a link a ``Path``. ``cavity``
    is typed as an int and described as a count, but sits on a per-channel row where
    only an index makes sense, so it carries the search site (``sp<n>``) the object
    was traced from - the void this module actually works in, a chamber being carved
    at its own probe radius and searched on its own. Channels in one Delaunay cavity
    but different chambers therefore differ here; the log's site table maps every
    site to its cavity. ``method`` is ``CaviTracer`` and ``software`` is ProDy and
    its version.

    Every physicochemical item is written, on the scales MOLE's method page names:
    hydropathy from Kyte and Doolittle, hydrophobicity from Cid *et al.*, polarity
    from Zimmerman *et al.*, mutability from Jones, Taylor and Thornton, and the
    ``logP``, ``logD`` and ``logS`` fragment values MOLE publishes as prose. Only
    ``bRadius`` can be absent, where the structure carries no B-factors.

    They are those scales averaged over *this* module's lining, not a reproduction
    of MOLEonline's numbers: the two disagree wherever the two linings do, and MOLE
    also applies an undocumented per-scale mainchain value to backbone-only contacts
    that is not imitated here.

    Unlike the PQR writers in this module, and like :func:`.writePDB` and
    :func:`.writePQR`, this takes the filename first, goes through
    :func:`.openFile` - so the backup setting is honoured - and returns the name it
    wrote. ``separate`` has no counterpart here: the value of the format is the joins
    between its categories, and a file per channel would cut exactly those.

    :returns: the filename written
    :rtype: str"""

    from pathlib import Path
    from prody.utilities import openFile

    if object_type not in ('channel', 'pore'):
        raise ValueError("object_type must be 'channel' or 'pore'")

    options = _popLiningOptions(kwargs)
    if kwargs:
        # As calcChannels does: a misspelled option that was quietly dropped would
        # be found only by noticing the file was written with the default.
        raise TypeError('writeChannelsCIF() got an unexpected keyword argument '
                        '{0}.'.format(', '.join(repr(k) for k in sorted(kwargs))))

    if not isinstance(channels, list):
        channels = [channels]
    links = list(links) if links else []

    path = _cifOutputPath(filename) if autoext else Path(filename)

    context = _cifContext(atoms, options, auto)

    rows = {'channel': [], 'profile': [], 'props': [], 'layer': [],
            'residue': [], 'layer_residue': [], 'weighted': []}

    for objects, kind in ((channels, object_type), (links, 'link')):
        for index, obj in enumerate(objects):
            produced = _objectCifRows(
                '{0}{1}'.format(kind, index), obj, _CIF_OBJECT_TYPES[kind],
                context, num_samples)
            for category, values in produced.items():
                rows[category].extend(values)

    if not rows['channel']:
        # Nothing traced, so no file: a block holding only its own audit record
        # says a run happened and nothing else, and a reader cannot tell it from
        # a file whose channels failed to write. The count is already reported by
        # the caller. Returns None rather than a name nothing is under.
        return None

    notes = []
    if atoms is not None:
        # Once for the file, not once per object: a cofactor sitting in several
        # channels of one protein is one finding.
        _reportDeepLiningAtoms(atoms, context['deep'])
        notes = _reportUncoveredLining(context)

    if structure is not None:
        from prody.proteins.ciffile import writeMMCIF
        writeMMCIF(str(path), structure, autoext=False)
        out = openFile(str(path), 'a')
    else:
        out = openFile(str(path), 'w')
        out.write('data_{0}\n'.format(context['block']))

    try:
        _writeCifCategories(out, rows, notes)
    finally:
        out.close()

    LOGGER.info("{0} channel(s) written to {1}{2}.".format(
        len(rows['channel']), path,
        '' if atoms is not None else
        ' (geometry only - no structure given to measure the lining against)'))

    return str(path)


def _cifContext(atoms, options, auto):
    """What every object of one file shares: the structure and its lookups.

    Built once because the queries are per structure and not per channel - the tree
    and the radii cost more to rebuild than the query they serve, which is the same
    reason :func:`_liningSource` is passed in rather than derived per object."""

    from prody import __version__

    from collections import Counter

    context = {'atoms': atoms, 'auto': bool(auto), 'deep': {},
               'uncovered': Counter(), 'coverage': {},
               'distA': options.distA, 'block': 'channels',
               'software': 'ProDy {0}'.format(__version__),
               'lining': None, 'backbone': None, 'betas': None,
               'allowed': set(), 'resname_of': {}, 'resnum_of': {},
               'chid_of': {}}

    if atoms is None:
        return context

    title = (atoms.getTitle() or '').strip()
    if title:
        context['block'] = '_'.join(title.split())

    context['lining'] = _liningSource(atoms)
    context['backbone'] = _backboneRadiiSource(atoms)

    betas = atoms.getBetas()
    context['betas'] = betas if betas is not None and np.any(betas) else None

    # The residues a lining report may name, as _formatLiningResidues decides it:
    # waters shape no channel, calcChannels having dropped them before tessellating,
    # and FIL pseudoatoms are this module's own output.
    reportable = atoms.select('not resname FIL' if options.include_water
                              else 'not water and not resname FIL')
    if reportable is not None:
        context['allowed'] = set(reportable.getResindices().tolist())

    resindices = atoms.getResindices()
    unique, first = np.unique(resindices, return_index=True)
    resnames, resnums = atoms.getResnames(), atoms.getResnums()
    chids, icodes = atoms.getChids(), atoms.getIcodes()

    for residue, index in zip(unique.tolist(), first.tolist()):
        context['resname_of'][residue] = resnames[index]
        context['resnum_of'][residue] = resnums[index]
        context['chid_of'][residue] = chids[index].strip()
        if icodes[index].strip():
            # The schema has no home for an insertion code: sequence_number is an
            # int linked to _atom_site.auth_seq_id. Two residues differing only by
            # icode therefore write the same number, which is worth saying once.
            context.setdefault('icodes', []).append(residue)

    if context.get('icodes'):
        _warn("{0} residue(s) carry insertion codes, which the tunnels schema "
              "cannot record - _sb_ncbr_channel_residue.sequence_number is an "
              "integer. Those residues are written under their number "
              "alone.".format(len(context['icodes'])))

    return context


def _objectSortKey(object_id):
    """``('channel', 9)`` for ``channel9``, so ids sort numerically not lexically."""

    import re

    match = re.match(r'^(\D*)(\d+)$', object_id)
    return (match.group(1), int(match.group(2))) if match else (object_id, -1)


def _reportUncoveredLining(context):
    """Say how much of each lining the physicochemical properties actually cover.

    Returns the lines to write into the file as comments, having already logged the
    same thing.

    The properties are computed over standard residues alone, which is how the
    schema defines them and what MOLE does. That is a true statement about those
    residues and a misleading one about the channel, wherever the rest of the wall
    is a nucleotide, a cofactor or a modified residue - and the rest is rarely a
    random sample, a D-peptide losing its D-residues and a ribosomal route its
    nucleotides. So the fraction covered travels with the numbers.

    It goes in as a CIF comment because the schema has nowhere to put it: no
    category carries a coverage or a het-residue field, and inventing an item would
    make the file non-conformant. A comment is ignored by every parser and read by
    every person, which is the right audience for a caveat. The residues themselves
    are named in ``sb_ncbr_channel_residue`` regardless, so a reader can always
    recompute this."""

    coverage = context.get('coverage') or {}
    if not coverage:
        return []

    def named(counts):
        """``HEM`` for one, ``DLE(4), DVA(2)`` where a name recurs."""
        return ', '.join(name if count == 1 else '{0}({1})'.format(name, count)
                         for name, count in counts.most_common())

    worst = min(scored / float(total) for scored, total, _ in coverage.values())

    _warn("{0} channel(s) are lined in part by residues no scale covers, so their "
          "properties describe only the standard residues - as little as {1:.0%} of "
          "the lining. Not covered: {2}.".format(
              len(coverage), worst, named(context['uncovered'])))

    lines = ['# Physicochemical properties are computed over standard residues only,',
             '# as the schema defines them. Where a lining holds anything else, the',
             '# properties describe the standard part of it alone:']
    # Numerically, so channel9 precedes channel11 as it does everywhere else.
    for object_id in sorted(coverage, key=_objectSortKey):
        scored, total, counts = coverage[object_id]
        lines.append('#   {0}: {1} of {2} residues ({3:.0%}); not covered: {4}'.format(
            object_id, scored, total, scored / float(total), named(counts)))
    return lines


def _writeCifCategories(out, rows, notes=()):
    """Write the audit block and every non-empty category, in schema order."""

    out.write('#\n')
    out.write('loop_\n')
    out.write('_audit_conform.dict_name\n')
    out.write('_audit_conform.dict_version\n')
    out.write('_audit_conform.dict_location\n')
    out.write('{0} {1} {2}\n'.format(_CIF_DICT_NAME, _CIF_DICT_VERSION,
                                     _CIF_DICT_LOCATION))

    _cifLoop(out, 'sb_ncbr_channel', ['id', 'type', 'method', 'software',
                                      'auto', 'cavity'], rows['channel'])

    _cifLoop(out, 'sb_ncbr_channel_profile',
             ['channel_id', 'radius', 'free_radius', 'distance', 'T',
              'x', 'y', 'z', 'charge'], rows['profile'],
             # T keys the category together with channel_id, and three decimals
             # collide on a route sampled at more than a thousand points.
             precision={'T': 5})

    if notes:
        out.write('#\n')
        for line in notes:
            out.write(line + '\n')

    _cifLoop(out, 'sb_ncbr_channel_props',
             ['channel_id', 'charge', 'hydropathy', 'hydrophobicity',
              'mutability', 'numNegatives', 'numPositives', 'polarity',
              'logD', 'logP', 'logS', 'ionizable', 'bRadius'], rows['props'])

    _cifLoop(out, 'sb_ncbr_channel_layer',
             ['channel_id', 'order', 'min_radius', 'min_free_radius',
              'start_distance', 'end_distance', 'local_minimum', 'bottleneck',
              'charge', 'numPositives', 'numNegatives', 'hydrophobicity',
              'hydropathy', 'polarity', 'mutability'], rows['layer'])

    _cifLoop(out, 'sb_ncbr_channel_layer_weighted_props',
             ['channel_id', 'hydropathy', 'hydrophobicity', 'mutability',
              'polarity', 'logD', 'logP', 'logS'], rows['weighted'])

    _cifLoop(out, 'sb_ncbr_channel_residue',
             ['channel_id', 'order', 'sequence_number', 'chain_id',
              'backbone'], rows['residue'])

    _cifLoop(out, 'sb_ncbr_channel_layer_residue',
             ['channel_id', 'layer_id', 'order', 'residue_id'],
             rows['layer_residue'])

    out.write('#\n')


class Channel:
    def __init__(self, tetrahedra, centerline_spline, radius_spline, length, 
                 bottleneck, volume, centers, radii, gates, cost=None):
        self.tetrahedra = tetrahedra
        self.centerline_spline = centerline_spline
        self.radius_spline = radius_spline
        self.length = length
        self.bottleneck = bottleneck
        self.volume = volume
        self.centers = centers
        self.radii = radii
        self.gates = gates
        # cost: accumulated Dijkstra path weight (sum of l / (d**2 + b) edge
        # costs) from the seed to the exit. Lower is better - a short path
        # through wide tetrahedra. Set by dijkstra(); None when not computed.
        self.cost = cost
        # curvature: path length / straight-line end-to-end distance
        # (dimensionless, >= 1; 1.0 == perfectly straight).
        self.curvature = self._computeCurvature()
        # origin: index of the start point this object was traced from, so that
        # the objects of one void can be analysed together. One label is enough
        # because the voids nest strictly - a chamber never spans two cavities,
        # since the chamber carve refines the cavity components - so the start
        # point names a chamber where chambers are defined and the whole cavity
        # where none qualified. Set by calcChannels; None when it did not run.
        self.origin = None
        # destination: for a chamber link, the start point of the chamber it
        # joins - the far end of the route, which a channel does not have,
        # since a channel ends at the solvent. None on channels and pores.
        self.destination = None
        # joined_chamber: the chamber label behind `destination`, set where the
        # link is built and kept only until calcChannels can translate it into
        # a start point, the labels themselves being internal to one run.
        self.joined_chamber = None

    def __getstate__(self):
        # Return a copy of the object state, but without the splines
        # which are not pickleable
        state = self.__dict__.copy()
        state["centerline_spline"] = None
        state["radius_spline"] = None
        return state
    
    def _computeCurvature(self):
        """Path length divided by straight-line end-to-end distance."""
        x = self.centerline_spline.x
        start = np.asarray(self.centerline_spline(x[0]))
        end = np.asarray(self.centerline_spline(x[-1]))
        straight = float(np.linalg.norm(end - start))
        if straight <= 1e-9:
            return float('nan')
        return float(self.length / straight)
    
    def _buildSplines(self):
        # computes splines for channels
        from scipy.interpolate import CubicSpline
        centers = self.centers
        radii = self.radii
        t = np.arange(len(centers))
        self.centerline_spline = CubicSpline(t, centers, bc_type='natural')
        gates = self.gates
        if len(gates):
            knot_t = np.empty(2 * len(centers) - 1)
            knot_t[0::2] = t
            knot_t[1::2] = t[:-1] + 0.5
            knot_r = np.empty_like(knot_t)
            knot_r[0::2] = radii
            knot_r[1::2] = gates
            self.radius_spline = CubicSpline(knot_t, knot_r, bc_type='natural')
        else:
            self.radius_spline = CubicSpline(t, radii, bc_type='natural')

    def getSplines(self):
        return self.centerline_spline, self.radius_spline
    
    
class State:
    def __init__(self, simplices, neighbors, vertices):
        self.simp = simplices
        self.neigh = neighbors
        self.verti = vertices
        
    def __eq__(self, other):
        if not isinstance(other, State):
            return False
        return (np.array_equal(self.simp, other.simp) and
                np.array_equal(self.neigh, other.neigh) and
                np.array_equal(self.verti, other.verti))
        
    def setState(self, simplices, neighbors, vertices):
        self.simp = simplices
        self.neigh = neighbors
        self.verti = vertices
        
    def getState(self):
        return self.simp, self.neigh, self.verti

class Cavity:
    def __init__(self, tetrahedra, is_connected_to_surface):
        self.tetrahedra = tetrahedra
        self.is_connected_to_surface = is_connected_to_surface
        self.starting_tetrahedron = None
        self.channels = []
        # Routes from one chamber of this cavity to a shallower one; see
        # ChannelCalculator.dijkstra. Empty unless the cavity was seeded per
        # chamber and chamber_links is on.
        self.links = []
        # Chamber of each seed, {starting tetrahedron: chamber label}, and its
        # volume, both filled by setStartingTetrahedraFromChambers.
        self.seed_chambers = {}
        self.seed_volumes = {}
        # How close each seeded chamber comes to the surface, {chamber label:
        # smallest depth over its tetrahedra}. This is what orders a pair of
        # chambers for the links, not the depth of either seed; see dijkstra.
        self.chamber_depths = {}
        # How many chambers lie in this cavity at all, seeded or not.
        self.chambers_found = 0
        self.depth = 0
        self.tetrahedra_depths = {}
        self.volume = 0.0

    def makeSurface(self):
        self.is_connected_to_surface = True
        
    def setExitTetrahedra(self, exit_tetrahedra):
        self.exit_tetrahedra = exit_tetrahedra
        
        
    def setStartingTetrahedron(self, tetrahedron):
        self.starting_tetrahedron = tetrahedron
        
    def setDepth(self, depth):
        self.depth = depth
        
    def addChannel(self, channel):
        self.channels.append(channel)

    def addLink(self, link):
        self.links.append(link)
    
        
def _rowsIsin(a, b):
    """Boolean mask marking which rows of 2D integer array ``a`` occur as a row
    in 2D array ``b`` (exact, order-sensitive match).

    Uses a void-dtype view so each row is treated as a single hashable scalar,
    turning an O(len(a) x len(b)) row-by-row scan into an O(len(a) + len(b))
    hashed membership test.
    """
    a = np.ascontiguousarray(a)
    b = np.ascontiguousarray(b)
    if a.size == 0 or b.size == 0:
        return np.zeros(a.shape[0], dtype=bool)
    if a.dtype != b.dtype:
        b = b.astype(a.dtype)
    va = a.view(np.dtype((np.void, a.dtype.itemsize * a.shape[1]))).ravel()
    vb = b.view(np.dtype((np.void, b.dtype.itemsize * b.shape[1]))).ravel()
    return np.isin(va, vb)


class ChannelCalculator:

    # Every cavity atom is written with this radius. It is a marker size, not a
    # measurement - see saveCavitiesToPdb for why a cavity cannot be drawn at its
    # true width.
    CAVITY_MARKER_RADIUS = 1.00
    def __init__(self, atoms, inner_radius=1.2, sparsity=6,
                 edge_cost='integral'):
        # Only the parameters the class actually consults are held here. surf_radius,
        # min_depth and bottleneck are stages of the pipeline, applied to the
        # tessellation and to the finished channels by calcChannels; keeping copies
        # of them on the calculator suggested it filtered by them, which it does
        # not.
        self.atoms = atoms
        self.inner_radius = inner_radius
        self.sparsity = sparsity
        # 'integral' (clearance-profile integral) or 'bottleneck' (l/(d^2+b));
        # the Dijkstra edge weight in buildSparseGraph. Resolved per diagram by
        # calcChannels (weighted defaults to 'bottleneck').
        self.edge_cost = edge_cost
        # Filled once by buildSparseGraph and read by the channel geometry:
        # the per-simplex Voronoi-vertex clearance (the spline knots) and the
        # per-edge gate clearance on each shared Delaunay face (the reported
        # bottleneck and, later, the Dijkstra cost). Cached so the two consumers
        # share one definition of width instead of recomputing it apart.
        self._vertex_clearance = None
        self._edge_bottleneck = None
        # The real atoms and a tree over them, built by _enclosureAtoms only
        # when an exchange of exits has to be judged.
        self._enclosure_atoms = None

    def sphereFit(self, points, simplices, vertices, vdw_radii, r, rows=None):
        """Sum-based clearance test: for each tetrahedron, decide whether a probe
        of radius ``r`` fits at its Voronoi vertex.

        Compares the sum of the distances from the Voronoi vertex to the four
        atom centres against the sum of ``r + vdw_radius`` over the same four
        atoms. Summing over the four atoms instead of testing each one is the
        tangent-sphere test exactly when the vertex is equidistant from all four,
        and a mild relaxation of it otherwise; the erosion defaults are
        calibrated against that behaviour.

        :arg points: coordinates of all atoms, shape ``(n_atoms, 3)``.
        :type points: :class:`~numpy.ndarray`

        :arg simplices: atom indices of each tetrahedron, shape ``(n, 4)``.
        :type simplices: :class:`~numpy.ndarray`

        :arg vertices: Voronoi vertex of each tetrahedron, shape ``(n, 3)``.
        :type vertices: :class:`~numpy.ndarray`

        :arg vdw_radii: van der Waals radius of each atom.
        :type vdw_radii: :class:`~numpy.ndarray`

        :arg r: probe radius in Angstrom.
        :type r: float

        :arg rows: optional boolean mask over the tetrahedra. Rows outside it are
            reported as ``False`` without paying for the distance computation.
            :meth:`deleteSimplices3d` uses it to restrict the surface pass to the
            boundary shell.
        :type rows: :class:`~numpy.ndarray`, optional

        :returns: boolean array of length ``len(simplices)``, ``True`` where the
            probe fits.
        :rtype: :class:`~numpy.ndarray`
        """
        fits = np.zeros(len(simplices), dtype=bool)
        if rows is None:
            rows = slice(None)
        elif not rows.any():
            return fits

        atom_coords = points[simplices[rows]]                          # (m, 4, 3)
        d_sum = np.linalg.norm(
            atom_coords - vertices[rows][:, None, :], axis=2).sum(axis=1)
        r_sum = (r + vdw_radii[simplices[rows]]).sum(axis=1)
        fits[rows] = d_sum >= r_sum

        return fits

    def deleteSimplices3d(self, points, simplices, neighbors, vertices,
                          vdw_radii, r, surface):
        """Delete the tetrahedra that fail the :meth:`sphereFit` probe test and
        return the compacted tessellation.

        Which side of the test is deleted depends on the pass. The ``surface``
        pass erodes from the outside in: a boundary tetrahedron the probe fits
        into is open to the solvent, hence exterior, and goes. Only tetrahedra on
        the boundary (those with a ``-1`` neighbour) are reachable from outside,
        so the test is restricted to that shell (~n^(2/3) rows) instead of being
        evaluated over every tetrahedron on each erosion iteration, and the
        caller iterates to convergence. The inner pass instead drops every
        tetrahedron too tight for the probe, leaving the space it can actually
        occupy.

        :arg points: coordinates of all atoms, shape ``(n_atoms, 3)``.
        :type points: :class:`~numpy.ndarray`

        :arg simplices: atom indices of each tetrahedron, shape ``(n, 4)``.
        :type simplices: :class:`~numpy.ndarray`

        :arg neighbors: index of the tetrahedron opposite each vertex, ``-1`` on
            the boundary, shape ``(n, 4)``.
        :type neighbors: :class:`~numpy.ndarray`

        :arg vertices: Voronoi vertex of each tetrahedron, shape ``(n, 3)``.
        :type vertices: :class:`~numpy.ndarray`

        :arg vdw_radii: van der Waals radius of each atom.
        :type vdw_radii: :class:`~numpy.ndarray`

        :arg r: probe radius in Angstrom.
        :type r: float

        :arg surface: ``True`` for one erosion step of the surface pass, ``False``
            for the inner pass.
        :type surface: bool

        :returns: the surviving ``(simplices, neighbors, vertices)``, with
            neighbour indices remapped to the new numbering and deleted
            neighbours set to ``-1``.
        :rtype: tuple
        """
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        n = len(simplices)
        if n == 0:
            return simplices, neighbors, vertices

        if surface:
            # Erode: a boundary tetrahedron the probe fits into is open to the
            # solvent, hence exterior. Only the boundary shell is reachable from
            # outside, so only it is tested.
            boundary = (neighbors == -1).any(axis=1)
            should_delete = self.sphereFit(points, simplices, vertices,
                                           vdw_radii, r, rows=boundary)
        else:
            # Carve: drop every tetrahedron too tight for the probe, anywhere,
            # leaving the space it can actually occupy.
            should_delete = ~self.sphereFit(points, simplices, vertices,
                                            vdw_radii, r)

        keep = ~should_delete
        simp = simplices[keep]
        neigh = neighbors[keep].copy()
        verti = vertices[keep]

        # Remap neighbour indices from the old numbering to the compacted one in a
        # single pass (deleted neighbours -> -1), replacing the previous
        # O(len(deleted) x len(neigh)) decrement loop.
        new_index = np.full(n, -1, dtype=neigh.dtype)
        new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
        neigh = np.where(neigh == -1, -1, new_index[neigh])

        return simp, neigh, verti

    def calcEnclosure(self, query, centers, tree=None, reach=25.0, rays=32,
                      step=0.75, radius=1.7):
        """Fraction of the directions seen from each point of ``query`` that are
        blocked by an atom within ``reach`` Angstrom.

        A local, probe-independent measure of burial: a point in the open solvent
        sees sky in most directions, a point inside a channel is surrounded
        whatever the channel's width.

        Rays are marched outwards and a ray is dropped as soon as it is blocked,
        which is what keeps this affordable: in a buried region most directions
        hit protein within the first few Angstrom, and only the few that escape
        are followed the whole way out. Marching costs ``rays x steps`` tree
        queries per point and so is all but insensitive to how many atoms there
        are, whereas testing every atom in range against every ray costs a
        multiple of the atom count, and with rays this sparse nearly all of that
        work is wasted on atoms that lie near no ray at all.

        Pass the real atoms here, not the balls of a homogenized diagram: burial
        is a property of the protein, not of the tessellation. Atoms are all given
        the same ``radius``, so one tree and one plain nearest-neighbour test
        suffice. This is a burial heuristic and not a surface calculation, and the
        alternative -- a per-atom radius, which no nearest-neighbour query can
        express -- buys nothing that ``min_enclosure`` cannot absorb.

        The sampling is part of the definition of the quantity and not an accuracy
        knob, which is why every caller that compares enclosures against a
        threshold leaves it alone. Adding rays is not a free refinement: more
        directions discover more of the thin escape routes out of a channel, so
        the enclosure of every point drifts downwards and a threshold calibrated
        at one ray count does not carry over to another.

        :arg query: points to evaluate, ``(n, 3)``.
        :arg centers: atom centres.
        :arg tree: optional prebuilt :class:`~scipy.spatial.cKDTree` over
            ``centers``, to avoid rebuilding it on every call.
        :arg reach: how far a ray is followed before the direction counts as
            open, in Angstrom. The default asks whether a point is inside the
            protein at all; a few Angstrom asks the much more local question of
            whether it sits in a niche.
        :arg rays: how many directions are sampled.
        :arg step: spacing of the samples along a ray, in Angstrom.
        :arg radius: one radius for every atom, on the scale of a heavy-atom vdW
            radius.
        :returns: ``n`` fractions in ``[0, 1]``."""
        query = np.asarray(query, dtype=float)
        if len(query) == 0:
            return np.empty(0)
        if tree is None:
            tree = _kdTree(centers)

        i = np.arange(rays) + 0.5
        phi = np.arccos(1 - 2 * i / rays)
        theta = np.pi * (1 + 5 ** 0.5) * i              # Fibonacci sphere
        directions = np.stack([np.sin(phi) * np.cos(theta),
                               np.sin(phi) * np.sin(theta),
                               np.cos(phi)], axis=1)

        blocked = np.zeros((len(query), rays), dtype=bool)
        live = np.ones_like(blocked)
        for distance in np.arange(step, reach + step, step):
            point, ray = np.nonzero(live)
            if not len(point):
                break
            samples = query[point] + directions[ray] * distance
            hit = tree.query(samples, distance_upper_bound=radius,
                             workers=-1)[0] <= radius
            blocked[point[hit], ray[hit]] = True
            live[point[hit], ray[hit]] = False

        return blocked.mean(axis=1)

    def peelSurfaceByEnclosure(self, points, simplices, neighbors, vertices,
                               vdw_radii, r, atom_coords, min_enclosure,
                               max_depth=None):
        """Erode the surface inward with a probe of radius ``r``, stopping where
        the tetrahedra stop being open to the solvent.

        This removes the "moat": the shell of true exterior that lies inside the
        ``surf_radius`` surface, because an ``surf_radius`` probe cannot enter the concavities it
        bridges over. Left in place the moat joins the cavity and offers wide,
        cheap routes along the outside of the protein.

        Neither obvious way of bounding the erosion works. A count of tetrahedron
        layers is not mesh-invariant, since a layer is one tetrahedron thick and
        tetrahedra shrink as the tessellation is refined. A depth in Angstrom is
        not ``surf_radius``-invariant, since the moat has no constant thickness: it is as
        deep as ``surf_radius - r`` inside a concavity and vanishes on a flat face, so a
        depth large enough to clear it where it is thick also marches down the
        channel mouths and erodes the channels themselves. At ``surf_radius = 10``, a
        reasonable setting for a porin or a ribosome, that leaves almost nothing.

        The rule used here is local instead. A boundary tetrahedron is stripped
        only while it is *open*, that is while its enclosure is below
        ``min_enclosure`` (see :meth:`calcEnclosure`). The moat is open by
        construction and goes; erosion halts by itself at the first buried layer.
        ``surf_radius`` then decides only where the erosion starts, not where it stops, so
        the result no longer depends on it, and ``surf_radius`` is left doing the one job
        it should: capping the mouths.

        Since enclosure is a static field, the peel is really "delete the
        outside-connected component of ``{enclosure < min_enclosure}``". A
        threshold above the enclosure of a channel interior (empirically about
        0.93, as a channel is itself an escape direction) therefore percolates
        along the channels and erodes the cavity away entirely. ``max_depth`` is
        available as a hard backstop, but the default threshold leaves a wide
        margin and the failure mode is loud - no channels at all - rather than a
        plausible-looking result with the real channels missing.

        Note that the probe test and the enclosure test deliberately run against
        different spheres. ``points`` and ``vdw_radii`` are the balls the diagram
        is built on, and the probe test has to use them or its geometry stops
        agreeing with :meth:`deleteSimplices3d`. ``atom_coords`` are the real
        atoms, and the enclosure test has to use those, or burial would depend on
        the tessellation. Under ``diagram="homogenized"`` the two are not the same
        set, as some 4700 atoms become some 33000 equal balls, which would make
        the enclosure test both slower and a function of ``max_deviation``. Under
        ``"simple"`` and ``"weighted"`` they coincide.

        :arg atom_coords: the real atoms, for the enclosure test.
        :arg min_enclosure: fraction of directions that must be blocked for a
            tetrahedron to count as interior and stop the erosion. ``<= 0``
            returns the state unchanged.
        :arg max_depth: optional cap, in Angstrom, on how far the front may
            advance from the initial surface. ``None`` (default) is uncapped.
        :returns: ``(simplices, neighbors, vertices)``, compacted."""
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        if min_enclosure <= 0 or len(simplices) == 0:
            return simplices, neighbors, vertices

        boundary = (neighbors == -1).any(axis=1)
        if not boundary.any():
            return simplices, neighbors, vertices

        # Fixed for the whole peel, so the cap bounds the total advance of the
        # front rather than its advance per pass.
        surface = _kdTree(vertices[boundary]) if max_depth is not None else None
        atoms = _kdTree(atom_coords)
        # Enclosure is a property of a point, not of the shrinking mesh, so a
        # tetrahedron re-examined on a later pass is never re-traced. Tetrahedra
        # are renumbered by the compaction below, but the four balls they are
        # built on are not, so those index the cache.
        traced = {}

        while True:
            n = len(simplices)
            if n == 0:
                break
            boundary = np.nonzero((neighbors == -1).any(axis=1))[0]
            if not len(boundary):
                break

            # The cheap tests first: does the probe fit (the same sum-based test
            # as deleteSimplices3d), and are we still inside the optional cap?
            ball_coords = points[simplices[boundary]]
            d_sum = np.linalg.norm(
                ball_coords - vertices[boundary][:, None, :], axis=2).sum(axis=1)
            r_sum = (r + vdw_radii[simplices[boundary]]).sum(axis=1)
            candidate = d_sum >= r_sum
            if max_depth is not None:
                candidate &= surface.query(vertices[boundary])[0] <= max_depth
            candidates = boundary[candidate]
            if not len(candidates):
                break

            # Ray tracing runs only on what survived those, and only once each.
            keys = [tuple(key) for key in simplices[candidates]]
            fresh = [i for i, key in enumerate(keys) if key not in traced]
            if fresh:
                values = self.calcEnclosure(vertices[candidates[fresh]],
                                            atom_coords, tree=atoms)
                for i, value in zip(fresh, values):
                    traced[keys[i]] = value
            enclosure = np.array([traced[key] for key in keys])

            should_delete = np.zeros(n, dtype=bool)
            should_delete[candidates[enclosure < min_enclosure]] = True
            if not should_delete.any():
                break

            keep = ~should_delete
            simplices = simplices[keep]
            neigh = neighbors[keep].copy()
            vertices = vertices[keep]

            new_index = np.full(n, -1, dtype=neigh.dtype)
            new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
            neighbors = np.where(neigh == -1, -1, new_index[neigh])

        return simplices, neighbors, vertices

    def deleteSection(self, simplices_subset, simplices, neighbors, vertices,
                      reverse=False):
        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        vertices = np.asarray(vertices)

        n = len(simplices)
        if n == 0:
            return simplices, neighbors, vertices

        # Which rows of `simplices` also appear in `simplices_subset` (exact,
        # order-sensitive row match via hashed membership instead of the former
        # O(n x len(subset)) per-row scan).
        matches = _rowsIsin(simplices, np.asarray(simplices_subset))
        keep = matches if reverse else ~matches

        simp = simplices[keep]
        neigh = neighbors[keep].copy()
        verti = vertices[keep]

        new_index = np.full(n, -1, dtype=neigh.dtype)
        new_index[keep] = np.arange(keep.sum(), dtype=neigh.dtype)
        neigh = np.where(neigh == -1, -1, new_index[neigh])

        return simp, neigh, verti

    @staticmethod
    def getVdwRadii(atoms, warn=True):
        """Van der Waals radius of each element in *atoms*, from :data:`VDW_RADII`.

        *atoms* is a sequence of element symbols, matched case-insensitively. An
        element the table does not cover takes its ``UNKNOWN`` entry: a structure
        is tessellated as the caller supplied it, and only water is dropped
        first, so any ion or cofactor kept in the selection reaches here - a
        metal outside the table would otherwise stop the trace rather than the
        report, which is the wrong end to fail at. The substitution is announced
        unless *warn* is false, for callers that have already said it themselves.

        A static method, so the radii a diagram would have been built on can be
        had without a calculator to build one (see :func:`_atomRadii`)."""

        elements = np.char.upper(np.asarray(atoms, dtype=str))
        if warn:
            unknown = sorted(set(elements.tolist()) - set(VDW_RADII))
            if unknown:
                _warn("No van der Waals radius for {0}; {1} atom(s) are "
                      "measured with the UNKNOWN radius of {2} A instead."
                      .format(', '.join(repr(element) for element in unknown),
                              int(np.isin(elements, unknown).sum()),
                              VDW_RADII['UNKNOWN']))

        return np.array([VDW_RADII.get(element, VDW_RADII['UNKNOWN'])
                         for element in elements])

    def _fibonacciSphere(self, n):
        """Return ``n`` roughly evenly distributed unit vectors on a sphere using
        the Fibonacci (golden spiral) lattice."""
        n = int(np.maximum(1, n))
        indices = np.arange(n) + 0.5
        phi = np.arccos(1.0 - 2.0 * indices / n)
        theta = np.pi * (1.0 + 5.0 ** 0.5) * indices
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        return np.stack([x, y, z], axis=1)

    def _shellPointCount(self, rad, rho, max_deviation):
        """Number of equal balls of radius ``rho`` to place on a shell of radius
        ``rad`` so that the outer envelope of their union stays within
        ``max_deviation`` of ``rad + rho``.

        A ball centered at radius ``rad`` only touches the target sphere of radius
        ``rad + rho`` at a single point, so the shell must be sampled densely
        enough that the "valleys" between neighbouring balls do not dip more than
        ``max_deviation``. Each ball covers a spherical cap of half-angle
        ``alpha`` on the ``rad + rho - max_deviation`` sphere; the count is the
        number of such caps needed to tile the sphere (with an overlap factor).
        """
        r = rad + rho - max_deviation
        cos_alpha = (rad * rad + r * r - rho * rho) / (2.0 * rad * r)
        cos_alpha = float(np.clip(cos_alpha, -1.0, 1.0))
        if cos_alpha >= 1.0:
            return 1
        # The exact number of caps to tile the sphere is 2 / (1 - cos_alpha); we use
        # a factor of 4 (a ~2x overlap margin). This is NOT slack to be trimmed: the
        # Fibonacci lattice is not an optimal packing and its coverage efficiency
        # degrades as the shell (and count) grows, so the factor needed to actually
        # hold the max_deviation bound increases with atom size. A single constant
        # must therefore be sized for the largest atoms (e.g. metals); 4 keeps the
        # measured dip below max_deviation across the whole range, whereas 3 already
        # fails for anything larger than the thinnest shell. Lowering it silently
        # breaks large-atom accuracy - tune max_deviation instead to change cost.
        return int(np.ceil(4.0 / (1.0 - cos_alpha)))

    def homogenizeAtoms(self, coords, vdw_radii, max_deviation=0.2):
        """Substitute every atom by a set of homogeneous balls whose common radius
        equals the smallest van der Waals radius present in the structure.

        Each atom of radius ``R`` is replaced by a collection of overlapping balls
        of radius ``rho = min(vdw_radii)`` arranged on concentric shells (plus a
        central ball) so that their union approximates the original atomic sphere
        to within ``max_deviation``. Because all resulting balls share the same
        radius, an ordinary Voronoi / Delaunay tessellation of their centers yields
        an accurate estimate of the additively weighted (power) Voronoi diagram of
        the original atoms. This is the approach used by MolAxis and CAVER 3 and
        avoids simply discarding the smaller (e.g. hydrogen) atoms.

        Atoms whose radius is within ``max_deviation`` of ``rho`` are kept as a
        single ball, so a structure of similarly sized atoms is left essentially
        unchanged while a structure containing hydrogens (small ``rho``) fills its
        larger atoms with several balls.

        :arg coords: atomic coordinates, shape ``(N, 3)``
        :arg vdw_radii: per-atom van der Waals radii, shape ``(N,)``
        :arg max_deviation: maximum tolerated deviation (in Angstrom) between the
            union surface of the substitute balls and the original atomic surface.
            Smaller values are more accurate but generate more balls. Default 0.2.
        :returns: a tuple ``(new_coords, new_radii)`` where every entry of
            ``new_radii`` equals ``rho``.
        """
        coords = np.asarray(coords, dtype=float)
        vdw_radii = np.asarray(vdw_radii, dtype=float)

        rho = float(np.min(vdw_radii))
        tol = 1e-6
        new_points = []

        for center, R in zip(coords, vdw_radii):
            # Atoms within max_deviation of the smallest radius stay a single ball.
            if R - rho <= max_deviation + tol:
                new_points.append(center)
                continue

            # Central ball plus concentric shells stepped by rho, with the
            # outermost shell at (R - rho) so the union surface reaches R.
            new_points.append(center)
            shell_radii = list(np.arange(rho, R - rho, rho))
            if not shell_radii or (R - rho) - shell_radii[-1] > tol:
                shell_radii.append(R - rho)

            for rad in shell_radii:
                if rad <= tol:
                    continue
                n = self._shellPointCount(rad, rho, max_deviation)
                new_points.extend(center + rad * self._fibonacciSphere(n))

        new_points = np.array(new_points)
        new_radii = np.full(len(new_points), rho)

        return new_points, new_radii

    def buildSurfaceDepthOracle(self, coords, vdw_radii, surf_radius, max_deviation, max_depth):
        """Interior/exterior depth oracle for relabeling the additively-weighted
        diagram's surface mouths (``diagram="weighted"``).

        The AW tessellation is not a clean simplicial complex: many interior 3-ball
        faces are left unpaired and masquerade as surface boundaries, so channels
        truncate to stubs. This builds a *homogenized* Voronoi diagram of the same
        atoms (a clean simplicial complex), erodes it with an ``surf_radius`` probe to separate
        solvent (exterior) from protein (interior), and labels every tetrahedron by its
        geodesic distance (A) below the molecular surface (0 = exterior/solvent).
        :meth:`getSurfaceCavities` then keeps only AW exit tetrahedra whose Voronoi
        vertex maps (via ``find_simplex``) to depth ``<= max_depth`` Angstrom.

        :returns: ``(delaunay, depth, max_depth)`` -- the homogenized
            :class:`~scipy.spatial.Delaunay`, its per-tetrahedron geodesic depth in
            Angstrom (points outside the hull are treated as depth 0), and the
            passed-through threshold.
        """
        from scipy.spatial import Delaunay

        hp, hrho = self.homogenizeAtoms(coords, vdw_radii, max_deviation)
        delaunay = Delaunay(hp)
        centers = self.calcCircumcenters(delaunay)
        clearance = (np.linalg.norm(hp[delaunay.simplices] - centers[:, None, :], axis=2)
                     - hrho[delaunay.simplices]).min(axis=1)
        neighbors = delaunay.neighbors
        n = len(delaunay.simplices)

        # surf_radius surface erosion: peel boundary tetrahedra wide enough for the probe, from
        # the hull inward, until nothing more can be removed. Survivors == interior.
        alive = np.ones(n, dtype=bool)
        while True:
            dead = np.zeros(n, dtype=bool)
            for k in range(4):
                col = neighbors[:, k]
                dead |= (col == -1) | ((col >= 0) & ~alive[col])
            peel = alive & dead & (clearance >= surf_radius)
            if not peel.any():
                break
            alive[peel] = False

        # Geodesic depth (A) below the surface: shortest path from the exterior
        # (solvent) tetrahedra inward along Voronoi edges. A physical distance, not a
        # tetrahedron-layer count, so weighted_mouth_depth is a real Angstrom threshold
        # that does not drift with the homogenization density.
        degenerate = self._degenerateTetrahedra(delaunay.simplices, centers, hp)
        scratch = np.full(n, -1, dtype=np.intp)
        depth = self._geodesicDepth(np.arange(n), np.nonzero(~alive)[0], neighbors,
                                    centers, degenerate, scratch)
        # Enclosed pockets never reached from the exterior are deep (never a mouth).
        reached = depth[np.isfinite(depth)]
        depth[~np.isfinite(depth)] = (float(reached.max()) + 5.0) if reached.size \
            else float(max_depth + 1)

        return delaunay, depth, max_depth
    
    def surfaceLayer(self, shape_simplices, filtered_simplices, shape_neighbors):
        shape_simplices = np.asarray(shape_simplices)
        shape_neighbors = np.asarray(shape_neighbors)
        filtered_simplices = np.asarray(filtered_simplices)

        # Split simplices into those touching the boundary (a -1 neighbour) and
        # the interior ones, preserving order.
        boundary = (shape_neighbors == -1).any(axis=1)
        surface_simplices = shape_simplices[boundary]
        surface_neighbors = shape_neighbors[boundary]
        interior_simplices = shape_simplices[~boundary]

        # Row-membership tests replace the former (N, M, 4) broadcast temporaries.
        surf_keep = _rowsIsin(surface_simplices, filtered_simplices)
        filtered_surface_simplices = surface_simplices[surf_keep]
        filtered_surface_neighbors = surface_neighbors[surf_keep]

        filtered_surface_neighbors = np.unique(filtered_surface_neighbors)
        filtered_surface_neighbors = filtered_surface_neighbors[filtered_surface_neighbors != 0]

        filtered_interior_simplices = interior_simplices[
            _rowsIsin(interior_simplices, filtered_simplices)]

        surface_layer_neighbor_simplices = shape_simplices[filtered_surface_neighbors]

        second_layer = filtered_interior_simplices[
            _rowsIsin(filtered_interior_simplices, surface_layer_neighbor_simplices)]

        return filtered_surface_simplices, second_layer

            
    def findGroups(self, neigh, is_cavity=True):
        x = neigh.shape[0]
        visited = np.zeros(x, dtype=bool)
        groups = []

        def dfs(tetra_index):
            stack = [tetra_index]
            current_group = []
            while stack:
                index = stack.pop()
                if not visited[index]:
                    visited[index] = True
                    current_group.append(index)
                    stack.extend(neighbor for neighbor in neigh[index] if neighbor != -1 and not visited[neighbor])
            return np.array(current_group)

        for i in range(x):
            if not visited[i]:
                current_group = dfs(i)
                if is_cavity:
                    groups.append(Cavity(current_group, False))
                else:
                    groups.append(current_group)

        return groups

    def findChambers(self, simplices, neighbors, vertices, points, vdw_radii,
                     seed_radius, depths=None, min_depth=0.0):
        """Label the tetrahedra of the cleared state by *chamber*: the connected
        components of the sub-network a probe of ``seed_radius`` fits into.

        A chamber is a place a channel should start from; ``inner_radius`` says
        what it may then squeeze through. One probe doing both jobs is what makes
        a single seed arbitrary: the interstices that merely connect the real
        pockets are themselves cavity under a small probe, so the deepest point
        of a cavity can be an accident of that network rather than a site.
        Carving again with a larger probe prunes those necks and leaves the
        pockets standing, which is where the seeds belong.

        The carve runs on the already-cleared state, so the labels index it
        directly - unlike :meth:`deleteSimplices3d`, which compacts the arrays
        and renumbers. No second tessellation is involved: one :meth:`sphereFit`
        pass and one connected-components call.

        :arg seed_radius: chamber-defining probe radius in Angstrom. It has to
            exceed ``inner_radius`` to prune anything: every tetrahedron of the
            cleared state passes ``sphereFit(inner_radius)`` by construction, so
            an equal radius keeps each cavity whole as a single chamber.
        :type seed_radius: float

        :arg depths: geodesic depth of every tetrahedron, in Angstrom, as
            :meth:`findDeepestTetrahedra` measured it, and ``-inf`` outside the
            cavities. Given, the carve keeps only what lies at least
            ``min_depth`` below the surface, so that a chamber is a site rather
            than a lobe holding its own way out. ``None`` carves on width alone.
        :type depths: :class:`~numpy.ndarray` or None

        :arg min_depth: depth floor in Angstrom, applied only when ``depths`` is
            given.
        :type min_depth: float

        :returns: ``(labels, volumes)`` - the chamber of each tetrahedron, ``-1``
            where the probe does not fit, and the approximate volume of each
            chamber indexed by that label, on the same scale
            :meth:`calculate_cavity_volumes` measures a cavity on.
        :rtype: tuple"""

        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components

        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        n = len(simplices)
        labels = np.full(n, -1, dtype=np.intp)

        fits = self.sphereFit(points, simplices, vertices, vdw_radii,
                              seed_radius)
        if depths is not None:
            # Stop the carve where everything else stops. A cavity ends at its
            # mouths and a channel is absorbed at one, but the carve knows only
            # about width, so a chamber runs right up to the mouths and takes
            # them in.The floor is min_depth, the same "buried enough to be a site" the
            # cavities are filtered on and the seeds are picked under, so the
            # boundary is one notion throughout. Nothing is lost from the
            # seeding: a seed had to clear min_depth already, so the fringe this
            # drops could never have supplied one.
            fits = fits & (np.asarray(depths) >= min_depth)
        if not fits.any():
            return labels, np.zeros(0, dtype=np.intp)

        # An edge is kept only when the probe fits at both of its ends, so every
        # component is either wholly inside the carve or a single tetrahedron
        # outside it. That is what lets the sizes below be counted over the
        # fitting tetrahedra alone: a label appearing among them can never be
        # shared with a tetrahedron the probe does not fit into.
        rows = np.repeat(np.arange(n), neighbors.shape[1])
        cols = neighbors.ravel()
        keep = (cols != -1) & fits[rows] & fits[np.where(cols < 0, 0, cols)]
        graph = csr_matrix((np.ones(int(keep.sum()), dtype=bool),
                            (rows[keep], cols[keep])), shape=(n, n))
        _, components = connected_components(graph, directed=False)

        labels[fits] = components[fits]

        return labels, self.calculateChamberVolumes(labels, simplices, points)

    def getSurfaceCavities(self, cavities, interior_simplices, second_layer,
                           state, mouth_oracle=None):
        surface_cavities = []
        
        for cavity in cavities:
            tetrahedra = cavity.tetrahedra
            second_layer_mask = np.isin(interior_simplices[tetrahedra], second_layer).all(axis=1)
            
            if np.any(second_layer_mask):
                exit_tetrahedra = tetrahedra[second_layer_mask]
                if mouth_oracle is not None:
                    # diagram="weighted": drop the false (buried) mouths the leaky AW
                    # diagram produces. Keep an exit tetrahedron only if its Voronoi
                    # vertex lies within max_depth Angstrom (geodesic) of the true
                    # molecular surface, per a homogenized interior/exterior oracle.
                    delaunay, depth, max_depth = mouth_oracle
                    located = delaunay.find_simplex(state.verti[exit_tetrahedra])
                    surface_depth = np.where(located >= 0, depth[located.clip(0)], 0.0)
                    exit_tetrahedra = exit_tetrahedra[surface_depth <= max_depth]
                    if len(exit_tetrahedra) == 0:
                        continue
                cavity.makeSurface()
                cavity.setExitTetrahedra(exit_tetrahedra)
                surface_cavities.append(cavity)
                
        return surface_cavities


    def mergeCavities(self, cavities, simplices):
        if not cavities:
            # No cavities survived filtering (e.g. a start_point selected a single
            # cavity shallower than min_depth, or none within start_point_search).
            # Return an empty (0, 4) slice so the pipeline yields zero channels
            # instead of crashing in np.concatenate on an empty list.
            return simplices[np.empty(0, dtype=np.intp)]
        merged_tetrahedra = np.concatenate([cavity.tetrahedra for cavity in cavities])
        return simplices[merged_tetrahedra]

    def _degenerateTetrahedra(self, simplices, vertices, points):
        # Scale-invariant flatness flag for the geodesic depth graph. A near-flat
        # tetrahedron has a runaway circumcenter, so its incident Voronoi edges can be
        # astronomically long (measured up to ~1e14 A) and would corrupt a shortest-path
        # depth. Flatness is the ratio of circumradius R to the tetrahedron's own atom
        # span L: well-shaped cells sit at R/L < ~4 whatever their absolute size - the
        # large fat cells that span a wide pore under a big probe included - while flat
        # slivers diverge to R/L -> infinity. Using the dimensionless R/L, never an
        # absolute length, keeps this correct for porins, ribosome tunnels and
        # large-probe (e.g. surf_radius=20) runs alike. Edges touching a flagged tetrahedron are
        # dropped from the graph in _geodesicDepth.
        apex = points[simplices]                                  # (n, 4, 3)
        R = np.linalg.norm(apex[:, 0] - vertices, axis=1)         # circumradius
        L = np.zeros(len(simplices))
        for a in range(4):
            for b in range(a + 1, 4):
                L = np.maximum(L, np.linalg.norm(apex[:, a] - apex[:, b], axis=1))
        # Well-shaped cells measure R/L up to ~3.7 (p99.9); flat slivers reach ~1e14.
        # Flag above 5 - the ~14-order gap makes the exact cut irrelevant across [5, 100].
        return R / np.maximum(L, 1e-9) > 5.0

    def _geodesicDepth(self, tetrahedra, sources, neighbors, vertices, degenerate, scratch):
        # Shortest-path distance (A) from the source set to every tetrahedron in
        # `tetrahedra`, along Voronoi edges weighted by circumcenter-to-circumcenter
        # distance (a multi-source Dijkstra over the induced subgraph). Edges touching a
        # degenerate (runaway-circumcenter) tetrahedron are dropped. `scratch` is a
        # reusable global->local index buffer (-1 outside `tetrahedra`), reset before
        # return so the caller can pass it again. Returns the distance aligned to
        # `tetrahedra`, np.inf where a tetrahedron is unreachable from every source.
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import dijkstra

        T = np.asarray(tetrahedra, dtype=np.intp)
        m = len(T)
        if m == 0:
            return np.empty(0)
        scratch[T] = np.arange(m)
        nb = neighbors[T]                                         # (m, deg) global
        safe_nb = np.where(nb < 0, 0, nb)
        local = np.where(nb < 0, -1, scratch[safe_nb])
        keep = local >= 0
        keep &= ~degenerate[T][:, None]                          # drop from a flat tetra
        keep &= ~np.where(nb < 0, True, degenerate[safe_nb])     # drop into a flat tetra
        keep = keep.ravel()
        row = np.repeat(np.arange(m), nb.shape[1])[keep]
        col = local.ravel()[keep]
        w = np.linalg.norm(vertices[T[row]] - vertices[T[col]], axis=1)
        graph = csr_matrix((w, (row, col)), shape=(m, m))
        src = scratch[np.asarray(sources, dtype=np.intp)]
        src = src[src >= 0]
        scratch[T] = -1                                          # reset for reuse
        if len(src) == 0:
            return np.full(m, np.inf)
        return dijkstra(graph, indices=src, min_only=True)

    def findDeepestTetrahedra(self, cavities, neighbors, vertices, points, simplices):
        # Cavity depth = the geodesic distance (A) from the cavity's openings (its exit
        # tetrahedra) to its farthest point, as a shortest path along Voronoi edges. This
        # is a physical length independent of tetrahedron size, so min_depth is mesh
        # invariant; the old +1-per-tetrahedron layer count grew as the mesh refined.
        degenerate = self._degenerateTetrahedra(simplices, vertices, points)
        scratch = np.full(neighbors.shape[0], -1, dtype=np.intp)
        for cavity in cavities:
            tetra = np.asarray(cavity.tetrahedra, dtype=np.intp)
            dist = self._geodesicDepth(tetra, cavity.exit_tetrahedra, neighbors,
                                       vertices, degenerate, scratch)
            finite = np.isfinite(dist)
            if not finite.any():
                # No exit reached any tetrahedron (degenerate cavity); keep it minimal.
                cavity.setStartingTetrahedron(np.array([int(tetra[0])]))
                cavity.setDepth(0.0)
                cavity.tetrahedra_depths = {int(tetra[0]): 0.0}
                continue
            deepest = int(np.argmax(np.where(finite, dist, -np.inf)))
            cavity.setStartingTetrahedron(np.array([int(tetra[deepest])]))
            cavity.setDepth(float(dist[deepest]))
            cavity.tetrahedra_depths = {int(tetra[k]): float(dist[k])
                                        for k in np.nonzero(finite)[0]}
            
    def calcCircumcenters(self, dela):
        # per-simplex circumcenters recovered analytically from the Delaunay 
        # paraboloid lifting, avoiding a second Qhull pass. Identical to scipy 
        # Voronoi vertices in general position.
        eq = dela.equations
        scale = dela.paraboloid_scale
        centers = -eq[:, :-2] / (2 * scale * eq[:, -2][:, None])
        return centers

    def _edgeBottleneck(self, ci, cj, shared_atoms, points, vdw_radii):
        # Minimum clearance along the Voronoi edge - the segment between the two
        # circumcenters ci, cj, dual to the Delaunay face the two tetrahedra
        # share - measured against that face's atoms. This is the edge bottleneck
        # radius: the circumcenters are local clearance maxima, so the tightest
        # point of the segment (the gate) generally lies between them and is
        # narrower than either endpoint.
        #
        # min over t in [0, 1] and over the shared atoms of |p(t) - a| - vdw(a),
        # with p(t) = ci + t (cj - ci). Because the min over the (t, atom)
        # product equals the min of the per-atom minima, each atom reduces to an
        # independent clamped point-to-segment distance - a closed form, no
        # sampling. The foot clamps to an endpoint when it falls outside the
        # segment, so the gate value is always <= both vertex clearances. Exact
        # for the straight edges of the homogenized/simple diagrams; for the
        # weighted (Apollonius) diagram the true edge is a slight arc and the
        # chord is a local approximation.
        a = points[shared_atoms]
        r = vdw_radii[shared_atoms]
        u = cj - ci
        uu = float(u @ u)
        if uu <= 1e-12:
            # Twin tetrahedra: the circumcenters coincide, the edge is a point,
            # so return the shared vertex clearance directly. This guard is
            # load-bearing, not defensive boilerplate - do not remove it. Two
            # near-cospherical tetrahedra (the "T5" degeneracy CAVER's paper
            # jitters against) produce coincident circumcenters; measured on real
            # structures these are absent at the default max_deviation but do
            # appear - a handful per structure - as the diagram is refined
            # (max_deviation -> 0.02). scipy/Qhull's cospherical facet merging
            # keeps the geometry finite (no NaN/inf circumcenters) so we do not
            # need CAVER's jitter/random-rotation precautions, but the coincident
            # circumcenters would still divide by ~0 here. Just above the
            # threshold the clip below handles it: for a tiny but non-zero edge
            # the foot clamps to an endpoint and the gate degrades continuously
            # to the endpoint clearance.
            return float(np.min(np.linalg.norm(a - ci, axis=1) - r))
        t = np.clip((a - ci) @ u / uu, 0.0, 1.0)
        p = ci + t[:, None] * u
        return float(np.min(np.linalg.norm(p - a, axis=1) - r))

    def _edgeBottleneckBatch(self, ci, cj, shared_atoms, points, vdw_radii):
        # Vectorized _edgeBottleneck over F edges at once. ``ci``, ``cj`` are
        # (F, 3) circumcenters (``ci`` the lower-index endpoint, so the reverse
        # edge yields a bitwise-identical gate) and ``shared_atoms`` is (F, 3)
        # atom indices of the shared Delaunay face. Same closed-form clamped
        # point-to-segment as the scalar version, with the twin-tetrahedron
        # (coincident circumcenter) guard applied row-wise. See _edgeBottleneck.
        # Returns ``(gate, tstar)``: the min clearance and the parameter t in
        # [0, 1] where it is attained (the binding atom's clamped foot), the
        # latter used by _edgeCostIntegralBatch to force a quadrature node on
        # the pinch.
        a = points[shared_atoms]                          # (F, 3, 3)
        r = vdw_radii[shared_atoms]                       # (F, 3)
        u = cj - ci                                       # (F, 3)
        uu = np.einsum('fj,fj->f', u, u)                  # (F,)
        twin = uu <= 1e-12
        uu_safe = np.where(twin, 1.0, uu)
        diff = a - ci[:, None, :]                         # (F, 3, 3)
        t = np.clip(np.einsum('faj,fj->fa', diff, u) / uu_safe[:, None],
                    0.0, 1.0)                             # (F, 3)
        p = ci[:, None, :] + t[:, :, None] * u[:, None, :]
        clr = np.linalg.norm(p - a, axis=2) - r          # (F, 3)
        amin = clr.argmin(axis=1)
        rows = np.arange(len(u))
        gate = clr[rows, amin]
        tstar = t[rows, amin]
        if twin.any():
            # coincident circumcenters: the edge is a point, so the gate is the
            # shared vertex clearance measured at ci (== cj) and tstar is moot.
            twin_gate = (np.linalg.norm(diff, axis=2) - r).min(axis=1)
            gate = np.where(twin, twin_gate, gate)
            tstar = np.where(twin, 0.0, tstar)
        return gate, tstar

    def _edgeCostIntegralBatch(self, ci, cj, tstar, shared_atoms, points,
                               vdw_radii, z=2.0, delta=0.3, r_floor=1e-2):
        # Price the edge by the integral of its clearance profile r(t)^-z, not by
        # the whole length charged at its single narrowest point (the l/gate^2 MOLE
        # cost). This profile-integral idea follows CAVER (TCBB'15, Eq. 1), but this
        # is NOT a CAVER reimplementation - the quadrature deliberately differs (see
        # below). r(t) = min over the 3 shared-face balls of |ci + t (cj-ci) - a| -
        # vdw is the exact clearance profile of the straight (homogenized/simple)
        # edge. The integral is additive under subdivision, so unlike l/gate^2 the
        # cost is mesh-invariant; the vertex-only formula's error grows with edge
        # length (measured ~19% -> ~1500% p95 across length bins) and drifts the
        # routing as max_deviation coarsens.
        #
        # The quadrature differs from CAVER's: where CAVER samples a plain uniform
        # grid (and takes its grid-minimum as the bottleneck), we take the uniform
        # grid linspace(0, 1, K) AND force a node at the exact analytic gate t*, so
        # the narrowest point is never missed - and the reported edge bottleneck is
        # that exact gate (see _edgeBottleneckBatch), not a grid sample. K =
        # ceil(L/delta) is a fixed ARCLENGTH step in Angstrom, so the sample count
        # scales with physical length (what makes it mesh-invariant). A short edge
        # (L <= delta) collapses to {0, t*, 1}, the three clearances already in
        # hand. r(t) is floored at r_floor so the integrand cannot diverge on a
        # sub-inner_radius edge that dips through an atom (the integral's analog of the
        # l/(d^2+b) regularizer); traversable edges have r(t) >= inner_radius >> r_floor and
        # are untouched. Exact for straight edges; a chord approximation for the
        # weighted (Apollonius) diagram, whose edges are arcs.
        chunk = 20000
        u = cj - ci
        L = np.linalg.norm(u, axis=1)
        cost = np.zeros(len(ci))
        alive = np.nonzero(L > 1e-6)[0]                   # twin/zero-length -> 0
        K = np.maximum(2, np.ceil(L / delta).astype(int) + 1)
        for k in np.unique(K[alive]):
            bucket = alive[K[alive] == k]
            grid = np.linspace(0.0, 1.0, k)
            for s in range(0, len(bucket), chunk):
                ii = bucket[s:s + chunk]
                nodes = np.sort(np.concatenate(
                    [np.broadcast_to(grid, (len(ii), k)), tstar[ii, None]],
                    axis=1), axis=1)                      # (n, k+1)
                p = ci[ii][:, None, :] + nodes[:, :, None] * u[ii][:, None, :]
                a = points[shared_atoms[ii]]              # (n, 3, 3)
                rr = np.linalg.norm(p[:, :, None, :] - a[:, None, :, :], axis=3) \
                    - vdw_radii[shared_atoms[ii]][:, None, :]
                rr = np.maximum(rr.min(axis=2), r_floor)  # (n, k+1) profile
                cost[ii] = np.trapz(rr ** (-z), nodes, axis=1) * L[ii]
        return cost

    def buildSparseGraph(self, simplices, neighbors, vertices, points, vdw_radii,
                         gate_floor=0.0):
        # One weighted CSR adjacency matrix for the whole cleared state, built
        # from array ops over the (N, deg) neighbour table - no Python loop.
        # Edges narrower than ``gate_floor`` are left out of it; see the pruning
        # at the end.
        # Edge (tetra -> neigh) weight is l / (d**2 + b), where l is the
        # vertex-to-vertex distance and d is the gate clearance on the shared
        # Delaunay face (min clearance along the connecting Voronoi edge). The
        # gate is the width the cost should see: the clearance between the two
        # circumcenters, not the entered node's own vertex clearance, which is a
        # local maximum and lets the search prefer a route that is actually
        # narrower at a face it never measures. Because the gate is symmetric,
        # the cost no longer depends on traversal direction the way the
        # entered-node vertex clearance did.
        #
        # Two edge-cost modes (self.edge_cost): 'bottleneck' is the l/(d**2 + b)
        # above (d = gate); 'integral' replaces it on face edges with a clearance-
        # profile integral (_edgeCostIntegralBatch), which is mesh-invariant.
        # Either way the gate is still cached as the reported edge bottleneck.
        from scipy.sparse import csr_matrix

        simplices = np.asarray(simplices)
        neighbors = np.asarray(neighbors)
        N, deg = neighbors.shape

        tetra_points = points[simplices]
        distances = np.linalg.norm(tetra_points - vertices[:, None, :], axis=2)
        bottleneck = np.min(distances - vdw_radii[simplices], axis=1)
        # Cache the per-simplex vertex clearance: the channel geometry otherwise
        # recomputes this identical min over each path's tetrahedra.
        self._vertex_clearance = bottleneck

        # Directed edge list straight off the neighbour table.
        rows = np.repeat(np.arange(N), deg)
        cols = neighbors.ravel()
        keep = cols != -1
        rows = rows[keep]
        cols = cols[keep]

        # Drive the gate geometry from the lower-index endpoint so an edge's two
        # directed copies see identical inputs and get a bitwise identical,
        # direction-symmetric gate - i.e. compute each undirected edge once.
        lo = np.minimum(rows, cols)
        hi = np.maximum(rows, cols)

        # Shared 3 atoms of each Delaunay face, convention-agnostic (set
        # intersection, not scipy's opposite-vertex rule, so it holds for the
        # weighted diagram too): ``present`` marks which of the lower simplex's
        # four atoms also appear in the higher one; a face-adjacent edge has 3.
        slo = simplices[lo]
        shi = simplices[hi]
        present = (slo[:, :, None] == shi[:, None, :]).any(axis=2)   # (M, 4)
        face = present.sum(axis=1) == 3

        # Rare non-face links fall back to the entered node's vertex clearance;
        # face links (essentially all of them) overwrite it with the gate.
        d = bottleneck[cols].astype(float, copy=True)
        fi = tstar = shared = None
        if face.any():
            fi = np.nonzero(face)[0]
            shared = slo[fi][present[fi]].reshape(-1, 3)
            d[fi], tstar = self._edgeBottleneckBatch(vertices[lo[fi]],
                                                     vertices[hi[fi]], shared,
                                                     points, vdw_radii)

        l = np.linalg.norm(vertices[rows] - vertices[cols], axis=1)
        b = 1e-3
        weight = l / (d * d + b)          # 'bottleneck' cost; non-face fallback
        if self.edge_cost == 'integral' and fi is not None:
            # 'integral' cost: replace the l/(d^2+b) of every face edge (the
            # bottleneck-only MOLE cost) with the clearance-profile integral. The
            # l/(d^2+b) fallback stays on the rare non-face links.
            #
            # No R/L flatness guard is needed here (unlike buildSurfaceDepthOracle,
            # which runs on the full diagram): this graph is the *cleared interior*
            # state, and a runaway circumcenter means a huge clearance, so the surf_radius
            # erosion has already stripped every flat boundary tetra - measured 0
            # degenerate tetra and max edge ~4 A in the cleared graph. That matters
            # because l/(d^2+b) *over*-prices a runaway edge (huge l) so the search
            # avoids it for free, whereas the integral would *under*-price it (r(t)
            # is huge over most of a runaway edge, so r(t)^-z ~ 0 there).
            weight[fi] = self._edgeCostIntegralBatch(
                vertices[lo[fi]], vertices[hi[fi]], tstar, shared, points, vdw_radii)

        # Per-edge gate cache (unordered key) read by _pathGates: each face edge
        # stored once, keyed (lo, hi). The reported bottleneck reads the same map.
        # Filled before the pruning below, so a gate is available whether or not
        # its edge is left in the graph.
        undirected = face & (rows < cols)
        self._edge_bottleneck = {
            (int(i), int(j)): float(v)
            for i, j, v in zip(rows[undirected], cols[undirected],
                               d[undirected])
        }

        # Edges the floor's probe does not fit through are dropped, so the search
        # routes around them instead of crossing one and having the channel
        # discarded afterwards by that same floor. calcChannels passes
        # ``bottleneck`` here under prune_narrow_edges and 0 otherwise. The floor
        # is that and not inner_radius: set below inner_radius it asks for
        # channels narrower than the traversal probe, and pruning at the probe
        # would delete exactly those.
        #
        # No channel can be lost this way. A path over a dropped edge pinches at
        # or below that edge's gate, so the floor filters it either way, while
        # every path that survives the filter uses only kept edges at unchanged
        # cost and so remains the cheapest route to its exit. What can appear is
        # a wider route to an exit whose cheapest one was narrow.
        #
        # ``d`` is the gate on face edges and the entered node's clearance on the
        # rare non-face ones, where _pathGates reports the tighter of the two
        # endpoints instead; d is then the larger of them, so an edge is dropped
        # only when the gate it would report is below the floor as well.
        if gate_floor > 0:
            open_enough = d >= gate_floor
            rows = rows[open_enough]
            cols = cols[open_enough]
            weight = weight[open_enough]

        return csr_matrix((weight, (rows, cols)), shape=(N, N))

    def dijkstra(self, cavity, graph, simplices, neighbors, vertices, points,
                 vdw_radii, divergence=0.2, chamber_labels=None):
        # a single multi-target Dijkstra from the seed over the cavity subgraph,
        # then every exit path reconstructed from the predecessor tree - 
        # instead of one heap search per (seed, exit) pair.
        # Channel geometry still goes through the current 
        # process_channel/Channel (Simpson-based volume).
        from scipy.sparse.csgraph import dijkstra
        from collections import defaultdict

        cavity_tetra = np.asarray(cavity.tetrahedra)
        if len(cavity_tetra) == 0:
            return
        global_to_local = {tetra: i for i, tetra in enumerate(cavity_tetra)}
        cavity_graph = graph[np.ix_(cavity_tetra, cavity_tetra)]
        # The same subgraph before the mouths stop conducting. The search out to
        # the surface needs them absorbing, but the short hop from an opening to
        # its own mouths does not: there the constraint only forces the path to
        # detour around every neighbouring mouth. Used by _addOpeningChannels.
        transit_graph = cavity_graph

        # A tunnel ends at the surface, but the Dijkstra cost has no such term:
        # it rewards width, and the widest places are the surface grooves. Left
        # free, the cheapest path to a far exit leaves the pocket at one mouth,
        # runs along the outside and re-enters at another - which is not a
        # tunnel. So a mouth must not *conduct*, only *absorb*: we drop the
        # outgoing edges of every mouth before the search, making "a channel
        # ends at its first surface contact" a hard constraint of the search
        # rather than a cut applied afterwards to the winning path.
        # The ordering is the whole point. Truncating after selection cuts a
        # path that was itself chosen *because* it ran along the surface, while
        # a genuine narrow interior corridor to the same mouth loses the argmin
        # to that groove and is never enumerated at all - it vanishes from the
        # output even though it is open. Mesh refinement makes the groove
        # cheaper, so interior tunnels drop out one by one as max_deviation
        # shrinks; forbidding transit removes that dependence entirely.
        # A mouth is a surface (exit) tetrahedron a probe of the traversal
        # radius inner_radius can leave through. Note the gate is inner_radius, not bottleneck:
        # bottleneck is a reporting filter, and letting it decide which mouths
        # absorb would let it silently re-open narrow mouths as transit nodes,
        # i.e. change the routes rather than filter them. For the homogenized
        # and weighted diagrams every surviving tetrahedron already clears inner_radius by
        # construction (equal radii + equidistant circumcenter collapse the
        # sum-based test in deleteSimplices3d to the per-atom clearance), so the
        # test is a no-op there; it earns its keep for diagram="simple", where
        # unequal radii break that identity.
        # Local indices of the tetrahedra a channel is allowed to end at; filled
        # from the mouths below. A surface cavity always has exit tetrahedra (it
        # is classified as one by having them), so this stays empty only for a
        # cavity that should produce no channels at all.
        terminals_local = []
        exit_tetra = np.asarray(getattr(cavity, 'exit_tetrahedra',
                                        np.empty(0, dtype=np.intp)))
        if len(exit_tetra):
            verts = vertices[exit_tetra]
            atom_pos = points[simplices[exit_tetra]]
            atom_rad = vdw_radii[simplices[exit_tetra]]
            clearance = (np.linalg.norm(atom_pos - verts[:, None, :],
                                        axis=2) - atom_rad).min(axis=1)
            seeds = set(int(s) for s in cavity.starting_tetrahedron)

            # Only the mouths themselves absorb. A tetrahedron that merely lies
            # inside a mouth's inscribed ball must NOT be absorbed: the ball's
            # radius is the clearance (up to ~2 A) and it reaches inward as
            # well as outward, so absorbing on it eats the corridors that
            # approach the surface and truncates real tunnels before they
            # arrive - measured to delete both known side tunnels at
            # max_deviation=0.1 while keeping them at 0.02, i.e. exactly the
            # silent, mesh-dependent tunnel loss this design exists to prevent.
            # A path can consequently still slip *past* a mouth through a twin
            # tetrahedron - a neighbour sharing almost the same circumcenter,
            # not itself in the second layer and so still conducting - and
            # surface again somewhere else. That leak is real but narrow (the
            # twins sit 0.1-0.7 A from a mouth, in the surface shell at depth
            # 1-4). It is closed downstream rather than by walling the graph off
            # against every mouth: a route is judged to have reached the surface
            # when it enters a mouth's inscribed *ball*, so slipping past the
            # tetrahedron does not slip past the opening. See the arrival test
            # below.
            absorbing = [global_to_local[int(t)]
                         for t, c in zip(exit_tetra, clearance)
                         if c >= self.inner_radius and int(t) in global_to_local
                         and int(t) not in seeds]
            absorbing_radius = np.array(
                [c for t, c in zip(exit_tetra, clearance)
                 if c >= self.inner_radius and int(t) in global_to_local
                 and int(t) not in seeds])
            if absorbing:
                # Zero the mouths' rows: edges *into* a mouth survive (a
                # channel may end there), edges *out of* it are gone.
                cavity_graph = cavity_graph.tolil()
                for i in absorbing:
                    cavity_graph.rows[i] = []
                    cavity_graph.data[i] = []
                cavity_graph = cavity_graph.tocsr()
            # Every mouth is a terminus; the dedup decides which of them are
            # one opening. See the comment at the target loop below.
            terminals_local = absorbing

        # The mouth layer is one tetrahedron thick, so zeroing its outgoing edges
        # does not actually stop a route from walking *around* a mouth through
        # the shell behind it and surfacing elsewhere - measured on 16 of 139
        # channels over a set of eight structures. An opening is therefore taken
        # to be a mouth's inscribed ball rather than the single tetrahedron at
        # its centre: a route that enters one has reached the surface there,
        # whatever tetrahedron it happens to stand in.
        opening_tree = None
        opening_count = None
        # A cavity with no mouth at all reports no channel, and the machinery
        # below has nothing to build its trees from, so everything it guards is
        # skipped rather than run on empty arrays.
        has_mouths = bool(terminals_local)
        if has_mouths:
            mouth_local = np.asarray(terminals_local, dtype=np.intp)
            mouth_xyz = vertices[cavity_tetra[mouth_local]]
            mouth_radius = np.asarray(absorbing_radius, dtype=float)
            # How many openings cover each tetrahedron. Queried once per mouth
            # over a tree of the cavity's vertices, not once per tetrahedron:
            # mouths are an order of magnitude fewer than tetrahedra, so this
            # costs O(total marked) instead of a sweep over the whole cavity.
            # A count rather than a flag because a seed can lie inside an opening
            # itself, and each seed then takes its own openings back out of the
            # mask - see the arrival test in the search below.
            node_tree = _kdTree(vertices[cavity_tetra])
            opening_count = np.zeros(len(cavity_tetra), dtype=np.int32)
            for centre, radius in zip(mouth_xyz, mouth_radius):
                hit = node_tree.query_ball_point(centre, radius)
                if hit:
                    opening_count[np.asarray(hit, dtype=np.intp)] += 1
            # For the few arrival nodes, which mouths' balls cover them.
            opening_tree = (_kdTree(mouth_xyz), mouth_xyz, mouth_radius,
                            float(mouth_radius.max()))
            # NOTE: the openings are deliberately *not* made absorbing in the
            # search itself. Doing so was measured and is worse: Dijkstra then
            # reroutes around every opening, discovers alternative interior
            # corridors to regions it previously reached through one, and each of
            # those becomes a new arrival - channel counts rise sharply. Leaving
            # the graph alone and stopping only the *candidate* at the first
            # arrival keeps the routes the search would have taken anyway.

        def mouthsAt(node_local):
            """Indices into ``mouth_local`` of the openings covering this node."""
            tree, centres, radii, reach = opening_tree
            here = vertices[cavity_tetra[node_local]]
            near = tree.query_ball_point(here, reach)
            if not near:
                return frozenset()
            near = np.asarray(near, dtype=np.intp)
            covers = np.linalg.norm(centres[near] - here, axis=1) < radii[near]
            return frozenset(int(k) for k in near[covers])

        # (path_local, cost, arrival_local, opening, seed_index) per arrival.
        opening_candidates = []
        # Per seed, (start_local, distances, predecessors) of its search, so that
        # _addOpeningChannels can read off the seed's own path to whichever mouth
        # phase two settles on. See the note there.
        seed_trees = []
        # Chamber links: the deep chamber of a cavity often has no way out of its
        # own, and reaches the surface only by joining a shallower chamber and
        # using that one's channels. Once every chamber is seeded, the dedup
        # rightly drops such a route as a duplicate of the shallower chamber's
        # shorter one, so the deep chamber's own access would go unreported. The
        # link is that access: the route from the deep seed, cut where it first
        # joins a shallower seeded chamber, so the full way out reads as
        # link(deep -> shallow) + channel(shallow -> surface) and the link's
        # bottleneck is the neck between them, which is the number that governs.
        link_candidates = []
        chamber_of_seed = cavity.seed_chambers if chamber_labels is not None else {}
        seed_of_chamber = {label: seed for seed, label in chamber_of_seed.items()}
        chamber_depth = cavity.chamber_depths

        for start_global in cavity.starting_tetrahedron:
            if start_global not in global_to_local:
                continue
            start_local = global_to_local[start_global]

            # A seed can lie inside an opening itself: a mouth's ball reaches
            # inward as well as outward, and on a wide mouth that is far enough
            # to swallow the widest buried tetrahedron of a shallow pocket. Two
            # things follow, and both are needed - with neither, every branch
            # arrives at once at the root and the cavity reports nothing at all.
            # The openings the seed already sits in are not this seed's arrival,
            # so they come out of the mask; the route has to get somewhere else
            # before it counts as having reached the surface. But the way out
            # through them is still a channel - the shortest one there is - so
            # the seed is emitted as an arrival in its own right below and phase
            # two walks it to the cheapest of those mouths.
            seed_openings = mouthsAt(start_local) if has_mouths else frozenset()
            inside_opening = None
            if has_mouths:
                if seed_openings:
                    own = np.zeros_like(opening_count)
                    for k in seed_openings:
                        hit = node_tree.query_ball_point(mouth_xyz[k],
                                                         mouth_radius[k])
                        if hit:
                            own[np.asarray(hit, dtype=np.intp)] += 1
                    inside_opening = (opening_count - own) > 0
                else:
                    inside_opening = opening_count > 0

            # directed=True: edge (u -> v) keeps weight l / (d_v**2 + b), i.e.
            # clearance of the node being *entered* - exactly the current heap
            # Dijkstra's cost model. (directed=False would symmetrize each edge 
            # to l / (max(d_u, d_v)**2 + b) and pick slightly different paths)
            distances, predecessors = dijkstra(
                cavity_graph, directed=True, indices=start_local,
                return_predecessors=True)
            seed_index = len(seed_trees)
            seed_trees.append((start_local, distances, predecessors,
                               seed_openings))
            parent_to_children = defaultdict(list)

            for node, parent in enumerate(predecessors):
                if parent >= 0:
                    parent_to_children[parent].append(node)

            own_chamber = chamber_of_seed.get(int(start_global))

            def firstForeignChamber(path_local):
                """Where along `path_local` the route first joins a chamber other
                than the one it started in, as ``(index, chamber)``.

                This is the point a route stops being this chamber's business.
                Both the channels and the links are cut here: past it the route is
                the joined chamber's own way out, which that chamber reports from
                its own seed, so continuing would report one corridor twice - once
                as the shallow chamber's channel and once as a longer, narrower
                copy owned by every chamber behind it."""
                for i, node in enumerate(path_local):
                    chamber = int(chamber_labels[cavity_tetra[node]])
                    if (chamber >= 0 and chamber != own_chamber
                            and chamber in seed_of_chamber):
                        return i, chamber
                return None, None

            paths = {}
            arrivals = []
            # `arrived` rides down the tree so that a route arrives once: the
            # node where it *first* enters an opening. A local test - inside here
            # and not inside at the parent - would fire again at every re-entry,
            # and a long route crossing in and out of the surface region would
            # spawn a candidate at each crossing rather than stopping at the
            # first, which is the whole point.
            stack = [(start_local, [start_local], False)]
            while stack:
                node, path, arrived = stack.pop()
                paths[node] = path
                if has_mouths and not arrived and inside_opening[node]:
                    arrivals.append(node)
                    arrived = True
                for child in parent_to_children.get(node, []):
                    stack.append((child, path + [child], arrived))

            # One candidate per *arrival*: the node where a route first enters an
            # opening, collected during the walk above so that each branch
            # contributes exactly one. The seed leads the list when it sits in an
            # opening of its own, standing in for the way straight out through it;
            # its route is the single node, and phase two supplies the whole of
            # the path.
            #
            # Note what is NOT done here: thinning the mouths by `sparsity` before
            # the search, as an earlier version did. That can pick a target
            # sitting behind another mouth - the path is absorbed at the nearer
            # one and never reaches the sampled target, so the tunnel is reported
            # nowhere at all. Which mouth shadows which is a tessellation
            # accident, so real tunnels vanished at some meshes and not others.
            # Every mouth is a legitimate terminus; exit identity is left to the
            # dedup, where `sparsity` merges the mouths that share one opening.
            for arrival in ([start_local] if seed_openings else []) + arrivals:
                arrival = int(arrival)
                path_local = paths.get(arrival)
                if path_local is None:
                    continue
                if chamber_labels is not None:
                    joined_at, _ = firstForeignChamber(path_local)
                    if joined_at is not None:
                        continue
                # The seed's own openings are the seed's business alone: for any
                # other arrival they are behind it, and letting one serve as the
                # exit would send phase two back the way it came.
                opening = (seed_openings if arrival == start_local
                           else mouthsAt(arrival) - seed_openings)
                if not opening:
                    continue
                opening_candidates.append((path_local,
                                           float(distances[arrival]),
                                           arrival, opening, seed_index))

            # One link per other seed reachable from this one, cut at the first
            # chamber the route joins. The cut chamber need not be the chamber of
            # the seed aimed at: a route that crosses chamber X on its way to seed
            # T is really S -> X, and X -> T is a link of its own, so cutting at
            # the first foreign chamber decomposes the chain instead of reporting
            # the composite. Emitted only in the deep -> shallow direction, so a
            # pair yields one link and it runs the way the molecule does: the
            # shallow chamber talks to the surface, the deep one talks to it.
            for target_global in chamber_of_seed:
                target_local = global_to_local.get(int(target_global))
                if target_local is None or target_local == start_local:
                    continue
                if np.isinf(distances[target_local]):
                    continue
                path_local = paths.get(target_local)
                if path_local is None:
                    continue

                cut, joined = firstForeignChamber(path_local)
                if not cut:            # never joined one, or started inside it
                    continue

                # Deeper source only, where "deeper" is how close each chamber
                # comes to the surface anywhere along itself, not how deep its
                # seed sits. The seed is the widest tetrahedron at least
                # min_depth down, so it is buried in every chamber and says
                # nothing about whether the chamber has its own way out. An exact
                # tie (two chambers both reaching the surface) is broken by seed
                # index so that the pair still yields exactly one link.
                here = chamber_depth.get(own_chamber, 0.0)
                there = chamber_depth.get(joined, 0.0)
                if here < there or (here == there and
                                    int(start_global) > int(seed_of_chamber[joined])):
                    continue

                path_global = cavity_tetra[path_local[:cut + 1]]
                link = Channel(path_global, *self.processChannel(
                    path_global, vertices, points, vdw_radii, simplices),
                    cost=float(distances[path_local[cut]]))
                link.joined_chamber = joined
                link_candidates.append((link, joined))

        # Channels first: a link is judged against the openings they report, so
        # they have to exist by the time the links are deduped.
        openings = np.empty((0, 3)), np.empty(0)
        if has_mouths:
            openings = self._addOpeningChannels(
                cavity, opening_candidates, divergence, transit_graph,
                cavity_tetra, vertices, points, vdw_radii, simplices,
                mouth_local, mouth_xyz, mouth_radius, seed_trees)
        self._addDedupedLinks(cavity, link_candidates, divergence, openings)

    def _exitImproves(self, home, candidate, seed_trees):
        """Is ``candidate`` a better way out than the ``home`` already reported?

        Both are routes the corridor test has already called one channel, so this
        is not asking which corridor to keep but which of two endings to show for
        it. A bare "is the mouth wider" rule is not enough on its own; the tests
        below are what separates an exchange worth making from one that trades a
        longer route for a number.

        An exchange has to be an improvement in one of two ways, and may not be a
        regression in either.

        *Enabler: a wider mouth.* By at least ``min_gain``, because clearances a
        few hundredths apart are one hole sampled by two tetrahedra, and buying
        that costs a longer route for nothing.

        *Enabler: a wider bottleneck.* The channel's own headline number, which a
        mouth width can move without touching. It matters on exactly the channels
        whose bottleneck *is* their exit - there, widening the mouth widens the
        channel, and holding those to the mouth threshold would refuse a real
        improvement over a few hundredths of an Angstrom.

        *Guard: the bottleneck may never go backwards.* The two routes are
        identical up to their fork, so a candidate whose bottleneck is lower must
        be pinching within the stretch it adds - it is reaching a better exit by
        threading a tighter gap than anything on the route it replaces.

        *Guard: nor may the mouth, while the bottleneck is unchanged.* Only
        while: the bottleneck is the tightest gate on the whole route, the
        terminal one into the mouth included, so a narrower exit cannot be hiding
        a tighter approach behind it - a neck just inside a wide mouth is a gate
        like any other and is already in the number. Where the candidate is
        genuinely wider at its tightest point, that is the passage improving, and
        a wider mouth on the route being replaced buys nothing when its own pinch
        lies upstream of it: the wide hole cannot be reached through the narrow
        one. Where the bottleneck does not improve, the mouth is the only thing
        left to compare and this guard is the whole of it. What the exchange may
        still not do is bury the exit, which clearance cannot see and the
        enclosure test below can.

        *Not much dearer, measured on the tails.* Only what follows the fork is
        being chosen; the shared head is common to both by construction. Charging
        the difference against the whole cost would make the verdict depend on how
        much identical prefix happens to precede the fork, so that one and the
        same decision reads as a small increase on a long channel and a large one
        on a short channel, and is refused only on the short one. The head is
        discounted only when both routes come from one seed's tree and so really
        do share it; otherwise the whole cost is compared, as before.

        *Not the same ending continued.* Both routes are paths in one Dijkstra
        tree, so they share a prefix and part at a fork. When the fork is at the
        last node or two the candidate does not go anywhere else, it carries on
        past the reported mouth and surfaces a little further out - which is the
        thing ending a channel at its first opening exists to prevent, and is
        refused on that ground rather than on any measurement.

        *No more buried than what it replaces.* This is the test that catches what
        clearance cannot. A wider mouth can sit deeper in a groove than the
        narrower one beside it, and then the channel gains a number and loses its
        ending. Burial is read at two scales because neither works alone - inside
        a wide lumen the closer scale saturates near zero and stops
        discriminating, while the further one misses a mouth that is pinched only
        locally. A tolerance of one ray keeps the sampling's own quantum from
        reading as a change."""
        # Least widening of the mouth, in Angstrom, that enables an exchange on
        # its own: above the width two tetrahedra sampling one hole differ by,
        # below the width at which an exit is really a different one.
        min_gain = 0.15
        # How much dearer the candidate's tail may be than the tail it replaces.
        max_tail_cost = 0.30
        scales = (6.0, 8.0)          # Angstrom; see the paragraph on burial above
        # How much of its own path a candidate must give up to be going somewhere
        # else rather than continuing past the reported mouth.
        min_fork = 1.0
        rays = 32
        eps = 1e-9

        if candidate['bottleneck'] < home['bottleneck'] - eps:
            return False
        if (candidate['bottleneck'] <= home['bottleneck'] + eps
                and candidate['clear'] < home['clear'] - eps):
            return False
        if not (candidate['clear'] >= home['clear'] + min_gain
                or candidate['bottleneck'] > home['bottleneck'] + eps):
            return False

        kept, offered = home['path'], candidate['path']
        shared = 0
        while (shared < min(len(kept), len(offered))
               and kept[shared] == offered[shared]):
            shared += 1
        fork = max(shared - 1, 0)

        forsaken = home['route'][fork:]
        if len(forsaken) < 2 or np.linalg.norm(
                np.diff(forsaken, axis=0), axis=1).sum() < min_fork:
            return False

        head = 0.0
        if (shared and home['seed'] == candidate['seed']
                and home['from_tree'] and candidate['from_tree']):
            head = float(seed_trees[home['seed']][1][kept[fork]])
        if not (0.0 <= head < home['cost']):
            head = 0.0               # a spliced route, or a degenerate prefix
        if (candidate['cost'] - head > (home['cost'] - head)
                * (1.0 + max_tail_cost)):
            return False

        coords, tree = self._enclosureAtoms()
        pair = np.vstack((home['xyz'], candidate['xyz']))
        for reach in scales:
            here, there = self.calcEnclosure(pair, coords, tree=tree,
                                             reach=reach, rays=rays)
            if there > here + 1.0 / rays:
                return False
        return True

    def _enclosureAtoms(self):
        """The real atoms and a tree over them, built once and only if asked for.

        Burial is a property of the protein, so this is the structure itself and
        not the balls the tessellation was built on, which may be homogenized."""
        if self._enclosure_atoms is None:
            coords = getCoords(self.atoms)
            self._enclosure_atoms = (coords, _kdTree(coords))
        return self._enclosure_atoms

    def _addOpeningChannels(self, cavity, candidates, divergence, cavity_graph,
                            cavity_tetra, vertices, points, vdw_radii, simplices,
                            mouth_local, mouth_xyz, mouth_radius, seed_trees):
        """Deduplicate routes at the opening they arrive at, then carry the
        survivors out through a mouth of that same opening.

        Three questions, answered separately because they are separate.

        *Which* channels exist is decided on the interior corridors, before the
        routes fan out across the mouth - that fan is the splay
        :meth:`_routeCoverage` has to discount when the comparison is made at the
        exits instead, and comparing before it removes the need. Only routes whose
        arrivals share a covering mouth are compared here; the rest meet later, at
        their exits, where the ordinary opening-and-corridor test still runs.

        *Which mouth* of that opening a survivor leaves by is the cheapest one
        from the seed, among those its arrival can reach.

        *What the channel looks like* is then the seed's own path to that mouth,
        read off the search that is already done - not seed->arrival spliced onto
        a fresh arrival->mouth search, whose two halves are each cheapest for
        their own endpoints and meet at an angle wherever the mouth sits on
        another branch of the tree. The splice remains as the fallback for the
        two cases the tree cannot answer: a mouth the seed's search never reached,
        being absorbed behind another, and a path to it that surfaces at some
        different opening on the way.

        Deduplicating before all of this is what keeps the arrival search cheap:
        one per reported channel rather than one per arrival, which on a large
        structure is the difference between tens and hundreds of searches."""

        from scipy.sparse.csgraph import dijkstra as sparse_dijkstra

        mouth_tree = _kdTree(mouth_xyz)
        mouth_reach = float(mouth_radius.max()) if len(mouth_radius) else 0.0

        def coveringMouths(node_local):
            """Indices into ``mouth_local`` of the openings covering this node."""
            here = vertices[cavity_tetra[node_local]]
            near = mouth_tree.query_ball_point(here, mouth_reach)
            if not near:
                return frozenset()
            near = np.asarray(near, dtype=np.intp)
            covers = np.linalg.norm(mouth_xyz[near] - here,
                                    axis=1) < mouth_radius[near]
            return frozenset(int(k) for k in near[covers])

        def surfacesElsewhere(route_local, seed_openings):
            """Does this route leave by an opening other than the one it first
            reached?

            The seed's path to a mouth is the cheapest one, but nothing stops it
            crossing a different opening on the way - which is the very leak the
            arrival test exists to close, surfacing at one opening and leaving by
            another. The test is the one the whole design is judged by: the
            opening a route first enters must be the opening it ends in.
            Neighbouring openings overlap, so this asks whether the two sets of
            mouths meet, not whether they are the same. Openings the seed already
            sits in do not count - the route starts inside those, and has not
            gone anywhere by being there."""
            exits_in = coveringMouths(int(route_local[-1]))
            for node in route_local[:-1]:
                here = coveringMouths(int(node)) - seed_openings
                if here:
                    return not (here & exits_in)
            return False

        kept = []    # (route, clearance, opening, path_local, arrival, cost,
                     #  seed_index)
        for path_local, cost, arrival, opening, seed in sorted(
                candidates, key=lambda c: c[1]):
            route = vertices[cavity_tetra[path_local]]
            route_clear = self._vertex_clearance[cavity_tetra[path_local]]
            duplicate = False
            for kept_route, kept_clear, kept_opening, _p, _a, _c, _s in kept:
                # Two routes share an opening when the balls covering their
                # arrivals overlap at all. Mouth circumcenters sit a fraction of
                # an Angstrom apart, so one opening is sampled by many nearly
                # coincident tetrahedra, and asking whether two routes reached the
                # *same* one is a question the tessellation answers arbitrarily;
                # asking whether their covering sets meet is stable.
                if not (opening & kept_opening):
                    continue
                # No opening discount here: these routes stop at the opening, so
                # the splay the discount exists for has not happened yet.
                if self._routeDivergence(route, kept_route, route_clear,
                                         kept_clear) <= divergence:
                    duplicate = True
                    break
            if not duplicate:
                kept.append((route, route_clear, opening, path_local, arrival,
                             cost, seed))

        centres = np.empty((0, 3))
        radii = np.empty(0)
        extended = []               # (cost, channel, route, xyz, radius, clear)
        reported = []               # one dict per channel kept, see below
        if not kept:
            return centres, radii

        # Phase two, batched over the survivors: the cheapest way out of each
        # one's opening, costed from its arrival rather than from the seed. The
        # distance from the seed would pick whichever mouth is cheapest to reach
        # overall, by a path that need not pass through this arrival at all.
        sources = [arrival for _r, _rc, _o, _p, arrival, _c, _s in kept]
        distances, predecessors = sparse_dijkstra(
            cavity_graph, directed=True, indices=sources,
            return_predecessors=True)
        distances = np.atleast_2d(distances)
        predecessors = np.atleast_2d(predecessors)

        for i, (_route, _clear, opening, path_local, arrival, cost, seed) in \
                enumerate(kept):
            reachable = [k for k in sorted(opening)
                         if np.isfinite(distances[i][mouth_local[k]])]
            path_full = np.asarray(path_local, dtype=np.intp)
            exit_k = None
            # Whether the reported route is the seed's own tree path. Only then
            # do two routes of one seed share a prefix whose cost is the same for
            # both, which is what _exitImproves discounts before comparing them.
            from_tree = False
            if reachable:
                # The route reported is the SEED's own path to the mouth, not
                # seed->arrival spliced onto a fresh arrival->mouth search. Each
                # half of such a splice is cheapest for its own endpoints, and
                # where the chosen mouth sits on a different branch of the seed's
                # tree the two meet at an angle: measured a 92 degree turn that
                # added 17% to a channel's length and 0.20 to its curvature,
                # against a straight path to the same exit the tree already held.
                # The arrival settles *which* opening a channel may end at, and
                # is where routes are deduplicated; it is not a waypoint the
                # geometry has to pass through.
                #
                # Which mouth of that opening is therefore chosen on the seed's
                # cost too. Choosing on the arrival's cost picks a mouth that is
                # cheap to reach from the arrival and then reports a route that
                # was never optimised for it, so the channel wanders where a
                # straighter one to a neighbouring mouth of the same opening
                # existed. The candidates stay gated on being reachable from the
                # arrival: that is what makes them mouths this route can leave
                # by at all.
                (start_local, seed_distances, seed_predecessors,
                 seed_openings) = seed_trees[seed]
                from_seed = [k for k in reachable
                             if np.isfinite(seed_distances[mouth_local[k]])]
                route_local = []
                if from_seed:
                    exit_k = min(from_seed,
                                 key=lambda k: seed_distances[mouth_local[k]])
                    walk = self._tracePath(seed_predecessors, start_local,
                                           int(mouth_local[exit_k]))
                    # Only if it does not surface somewhere else on the way: the
                    # cheapest path to a mouth is free to cross a neighbouring
                    # opening, and taking it then would reinstate exactly the
                    # leak this whole search exists to close.
                    if len(walk) > 1 and not surfacesElsewhere(
                            walk, seed_openings):
                        route_local = walk
                if route_local:
                    path_full = np.asarray(route_local, dtype=np.intp)
                    cost = float(seed_distances[int(mouth_local[exit_k])])
                    from_tree = True
                else:
                    # Either the seed's search never reached a mouth of this
                    # opening (all of them absorbed behind another), or its path
                    # to the cheapest one surfaces elsewhere first. Fall back on
                    # the splice, whose second half starts at the arrival and so
                    # cannot surface before it.
                    exit_k = min(reachable,
                                 key=lambda k: distances[i][mouth_local[k]])
                    tail = self._tracePath(predecessors[i], arrival,
                                           int(mouth_local[exit_k]))
                    if len(tail) > 1:
                        path_full = np.concatenate(
                            [path_full, np.asarray(tail[1:], dtype=np.intp)])
                        cost = cost + float(distances[i][mouth_local[exit_k]])
                    else:
                        exit_k = None
            # exit_k stays None when the opening covers this node but none of its
            # mouths can be walked to, every route there being absorbed first. The
            # channel is then reported as it stands, at the arrival - unless that
            # leaves nothing to report, which happens for a seed that sits in an
            # opening whose mouths are all unreachable: its route is the seed
            # alone, and a one-tetrahedron channel has no centerline.
            if len(path_full) < 2:
                continue

            path_global = cavity_tetra[path_full]
            channel = Channel(path_global, *self.processChannel(
                path_global, vertices, points, vdw_radii, simplices),
                cost=float(cost))
            exit_xyz = vertices[path_global[-1]]
            exit_radius = max(float(mouth_radius[exit_k]), self.sparsity / 2.0) \
                if exit_k is not None else self.sparsity / 2.0
            # The mouth's own clearance as well as the sphere it registers: the
            # sphere carries the sparsity floor, so it cannot say which of two
            # exits is the wider hole.
            exit_clear = float(mouth_radius[exit_k]) if exit_k is not None \
                else 0.0
            extended.append(dict(
                cost=float(cost), channel=channel, route=vertices[path_global],
                route_clear=self._vertex_clearance[path_global],
                xyz=exit_xyz, radius=exit_radius, clear=exit_clear,
                bottleneck=float(channel.bottleneck), path=path_full,
                seed=seed, from_tree=from_tree))

        # Second pass, on the finished channels: the arrival test above settles
        # only whether two routes came out of the *same* opening, which is a much
        # narrower question than `sparsity` asks. Two arrivals a few Angstrom
        # apart have disjoint mouth sets and are never compared there, yet
        # `sparsity` may well call their openings one. So the ordinary
        # opening-and-corridor identity still runs, now on exits that are
        # guaranteed to lie in the opening the route actually reached.
        # A duplicate is not simply dropped: if _exitImproves finds it a better
        # way out, it becomes the one reported. Two routes that the corridor test
        # calls the same channel can still end at very different mouths -
        # measured 2.30 A against 4.16 A on a 0.9% cost difference - and cost,
        # which is dominated by the interior, is close to blind to the difference.
        for offer in sorted(extended, key=lambda c: c['cost']):
            route, exit_xyz = offer['route'], offer['xyz']
            home = None
            for entry in reported:
                kept_xyz, kept_radius = entry['xyz'], entry['radius']
                if np.linalg.norm(exit_xyz - kept_xyz) >= (kept_radius
                                                           + offer['radius']):
                    continue
                # No opening discount, unlike the one-pass dedup. That discount
                # exists to stop the fan two routes make as they splay across a
                # shared mouth from reading as divergence, and that fan has
                # already been settled at the arrival, where the routes were
                # compared before they reached a mouth at all. Applying it again
                # here would discount the very stretch that tells apart two
                # corridors arriving at neighbouring openings.
                if self._routeDivergence(route, entry['route'],
                                         offer['route_clear'],
                                         entry['route_clear']) <= divergence:
                    home = entry
                    break
            if home is None:
                reported.append(offer)
            elif self._exitImproves(home, offer, seed_trees):
                founding_cost = home['cost']
                home.update(offer)
                # The founder's cost stays, so every later candidate is judged
                # against the cheapest route into this opening and the tolerance
                # cannot ratchet upwards through a chain of exchanges.
                home['cost'] = founding_cost

        # Added only now, so that an exchange replaces a representative rather
        # than leaving the superseded channel in the cavity. The list is still in
        # cost order, which is the order the channels are numbered in.
        for entry in reported:
            cavity.addChannel(entry['channel'])
            centres = np.vstack((centres, entry['xyz']))
            radii = np.append(radii, entry['radius'])
        return centres, radii

    def _tracePath(self, predecessors, source, target):
        """Local indices along ``source`` -> ``target`` in a predecessor array."""
        path = [int(target)]
        node = int(target)
        while node != int(source):
            node = int(predecessors[node])
            if node < 0:
                return []
            path.append(node)
        path.reverse()
        return path
      
    def _addDedupedLinks(self, cavity, candidates, divergence, openings=None):
        """Keep one link per (chamber joined, corridor taken), cheapest first.

        The same two-part identity the channels use, with the chamber standing in
        for the opening: two routes into one chamber along one corridor are one
        link, while two corridors into the same chamber are two, as are two
        chambers reached along what starts out as the same corridor. Cheapest
        first, and each candidate judged only against what is already kept, so
        the result does not depend on the order the seeds were searched in.

        A link that runs through a reported opening is dropped first, on the same
        ground the channels themselves are held to - a route passing through the
        exit sphere of a reported channel has left the protein there - and links
        were the one object never held to it, the
        chamber cut having been assumed to stop them soon enough. It does not: a
        route aimed at another seed has no reason to enter a mouth, so absorption
        never fires, and it steps around one through a neighbouring tetrahedron
        and travels the outside instead. 

        Where a channel is *cut* at the opening and keeps what came before, a
        link is dropped outright: what came before is the seed's route to the
        surface, which is a channel, and is already reported as one."""

        kept = []                     # (points, clearance, chamber joined)
        centres, radii = openings if openings is not None else (None, None)
        for link, joined in sorted(candidates, key=lambda c: c[0].cost):
            points_on_route = np.asarray(link.centerline_spline(
                link.centerline_spline.x))
            # The radius profile shares the centerline's parameter, so reading it
            # at the same knots gives the clearance at each of those points.
            clear_on_route = np.asarray(link.radius_spline(
                link.centerline_spline.x))
            if centres is not None and len(centres):
                if (np.linalg.norm(points_on_route[:, None, :] - centres,
                                   axis=2) < radii).any():
                    continue
            duplicate = False
            for kept_points, kept_clear, kept_chamber in kept:
                if kept_chamber != joined:
                    continue
                if self._routeDivergence(points_on_route, kept_points,
                                         clear_on_route,
                                         kept_clear) <= divergence:
                    duplicate = True
                    break
            if not duplicate:
                kept.append((points_on_route, clear_on_route, joined))
                cavity.addLink(link)


    def _routeDivergence(self, a, b, ra, rb):
        """How far two routes part company, per Angstrom of corridor they share.

        ``ra`` and ``rb`` are the clearances along ``a`` and ``b`` - the radius
        of a ball centred on the centerline that holds no atom - so each route
        is a tube, and the question asked at every vertex is how far the two
        tube *surfaces* are apart::

            g = |x - nearest vertex of the other route| - (r_here + r_there)

        Negative means the tubes overlap there: the two balls share a point,
        both are free of atoms, and a path from one centerline to the other runs
        through free space, so no wall separates them. Positive is the width of
        the gap between them. Both routes are measured, each against the other,
        and every vertex is weighted by the arc length it stands for. The result
        is the divergence integral per Angstrom of channel::

            I = sum of max(0, g) * weight              how far, times how long
            L = (length of a + length of b) / 2        the average channel
            return I / L                               Angstrom

        Two things have to enter the answer and neither is enough alone. *How
        far* apart they get says nothing about whether it is a brief excursion
        or a parting of the ways; *how much* of the route diverges says nothing
        about whether it strays by half an Angstrom or by eight. The integral
        holds both, and dividing by the length is what makes a two Angstrom arm
        off a ninety Angstrom trunk read differently from the same arm off a ten
        Angstrom one.

        Scaling by the clearance rather than by a fixed distance is what lets one
        number serve everywhere: the question has no absolute scale, since two
        paths a couple of Angstrom apart in a wide chamber have nothing between
        them while the same distance in a narrow throat spans a wall.

        Both routes are measured, rather than only the shorter one against the
        longer, because a route that carries on past where the other ended is
        the whole difference between them in the contained case - and the far
        stretch then registers exactly as much as it deserves: a continuation
        that stays inside the tube it came from contributes nothing, while one
        that leaves contributes its full length.

        Three limits worth knowing, all of them in the pairing rather than the
        formula. The nearest vertex of the other route is not a correspondence:
        where two routes fork, the fork stays the nearest point on the other
        one, so a diverging tail is measured against the fork rather than
        against anything it runs beside. The clearance grows towards a mouth, so
        the reach grows exactly where sibling routes fan out, and a divergence
        there can go unseen from one side - measured 7.2 Angstrom of separation
        still reading as touching tubes, caught only because the other route's
        arm registered from its own side. And a pair whose routes both fork into
        wide mouths could therefore be under-reported from both sides at once."""

        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        ra = np.asarray(ra, dtype=float)
        rb = np.asarray(rb, dtype=float)

        def gaps(here, r_here, there, r_there):
            """Surface gap at each vertex of ``here``, and its arc weight."""
            distance, nearest = _kdTree(there).query(here)
            gap = distance - (r_here + r_there[nearest])
            steps = np.linalg.norm(np.diff(here, axis=0), axis=1)
            weight = np.zeros(len(here))
            weight[:-1] += steps / 2.0
            weight[1:] += steps / 2.0
            return gap, weight

        if len(a) < 2 or len(b) < 2:
            # A route of one tetrahedron has no arc length to weigh, so the
            # integral is undefined - but the question still has an answer: is
            # that point inside the other's tube or outside it?
            point, r_point, other, r_other = (a, ra, b, rb) if len(a) < 2 \
                else (b, rb, a, ra)
            if not len(point) or len(other) < 1:
                return float('inf')
            distance, nearest = _kdTree(other).query(point)
            inside = (distance - (r_point + r_other[nearest])) <= 0
            return 0.0 if bool(np.all(inside)) else float('inf')

        gap_a, weight_a = gaps(a, ra, b, rb)
        gap_b, weight_b = gaps(b, rb, a, ra)
        gap = np.concatenate((gap_a, gap_b))
        weight = np.concatenate((weight_a, weight_b))

        # Each route's weights sum to its own arc length, so the two together
        # are the total and half of that is the average channel.
        span = float(weight.sum()) / 2.0
        if span <= 0:
            return float('inf')
        return float((np.maximum(gap, 0.0) * weight).sum() / span)
    
    def calculateMaxRadius(self, vertice, points, vdw_radii, simp):
        atom_positions = points[simp]
        radii = vdw_radii[simp]
        distances = np.linalg.norm(atom_positions - vertice, axis=1) - radii
        return np.min(distances)

    def _pathGates(self, tetrahedra, voronoi_vertices, points, vdw_radii,
                   simp, vertex_radii):
        # Per-edge gate clearance along the path: the minimum clearance on each
        # shared Delaunay face between consecutive circumcenters (the edge
        # bottleneck radius), read from the cache and recomputed for any edge the
        # map lacks. Each gate is <= the clearance at both its endpoints, so the
        # path minimum is the reported bottleneck and the gates are where the
        # tube pinches for the volume. Length len(tetrahedra) - 1; empty for a
        # single-tetrahedron path.
        n = len(tetrahedra)
        if n < 2:
            return np.empty(0)
        eb = self._edge_bottleneck
        gates = np.empty(n - 1)
        for k in range(n - 1):
            i, j = int(tetrahedra[k]), int(tetrahedra[k + 1])
            key = (i, j) if i < j else (j, i)
            g = eb.get(key) if eb is not None else None
            if g is None:
                shared = np.intersect1d(simp[i], simp[j], assume_unique=True)
                if len(shared) != 3:
                    # not face-adjacent (should not happen on a graph path);
                    # fall back to the tighter of the two endpoints
                    g = float(min(vertex_radii[k], vertex_radii[k + 1]))
                else:
                    g = self._edgeBottleneck(voronoi_vertices[i],
                                             voronoi_vertices[j], shared,
                                             points, vdw_radii)
            gates[k] = g
        return gates

    def calculateRadiusSpline(self, tetrahedra, voronoi_vertices, points,
                              vdw_radii, simp):
        tetrahedra = np.asarray(tetrahedra)
        # The per-vertex clearance is the same min buildSparseGraph already took
        # over every simplex; read it back instead of recomputing it per path.
        if self._vertex_clearance is not None:
            radii = self._vertex_clearance[tetrahedra]
        else:
            vertices = voronoi_vertices[tetrahedra]
            radii = np.array([self.calculateMaxRadius(v, points, vdw_radii, s)
                              for v, s in zip(vertices, simp[tetrahedra])])
        gates = self._pathGates(tetrahedra, voronoi_vertices, points,
                                vdw_radii, simp, radii)
        return radii, gates

    def processChannel(self, tetrahedra, voronoi_vertices, points, vdw_radii,
                       simp):
        """The geometry of one route: ``(centerline_spline, radius_spline,
        length, bottleneck, volume)``.

        The centerline runs through the circumcenters of ``tetrahedra`` and the
        radius profile through their clearances and the gates between them, over
        one shared parameter domain - so a value of the parameter names the same
        place on both, and the endpoints (hence the cap radii) are the route's
        own ends. The parameter measures distance, not tetrahedra; see the
        comments below for why, and what it costs to get that wrong.

        ``length`` and ``volume`` are the centerline's arc length and the volume
        of the tube it sweeps; ``bottleneck`` is the tightest gate, taken from
        the measurements rather than from the interpolated profile."""

        from scipy.interpolate import CubicSpline, PchipInterpolator

        centers = voronoi_vertices[tetrahedra]
        radii, gates = self.calculateRadiusSpline(tetrahedra,
                                                  voronoi_vertices,
                                                  points, vdw_radii, simp)
        bottleneck = float(np.min(gates)) if len(gates) else float(np.min(radii))

        # Coincident circumcenters - the twin tetrahedra _edgeBottleneck guards
        # against - are a repeated knot to a spline, and a parameter that
        # advances by distance would stand still at one. Collapse each run of
        # them onto its first point, at the same 1e-6 A separation that test
        # uses. Nothing is lost: the gate of a collapsed span is taken as the
        # tightest of the edges it covers, and the bottleneck above is over
        # every edge in any case.
        keep = [0]
        for i in range(1, len(centers)):
            if np.linalg.norm(centers[i] - centers[keep[-1]]) > 1e-6:
                keep.append(i)
        # The spline parameter advances by sqrt(step) - the centripetal
        # parameterization - rather than by one per tetrahedron. Circumcenters
        # are spaced anything but evenly along a route: neighbouring steps of
        # 0.1 and 3.7 A occur, and giving each of them one unit of parameter
        # makes the cubic overshoot the long edge and swing back, a bend the
        # route does not have. Measured over a set of proteins, the index
        # parameterization runs 6-8% (up to 22%) longer than the polyline
        # through the same circumcenters and wanders up to 0.6 A off it, so it
        # inflates every reported length, volume and curvature; centripetal
        # stays within 1% and roughly halves the wander. Plain chord length
        # fixes the length but wanders further still at abrupt turns, which is
        # the familiar Catmull-Rom result and holds here as well.
        if len(keep) < 2:
            # every circumcenter of the route sits in one place, so there is no
            # distance to parameterize by; fall back on the index.
            keep = np.arange(len(centers), dtype=np.intp)
            t = keep.astype(float)
        else:
            keep = np.asarray(keep, dtype=np.intp)
            step = np.linalg.norm(np.diff(centers[keep], axis=0), axis=1)
            t = np.concatenate([[0.0], np.cumsum(np.sqrt(step))])
        centerline_spline = CubicSpline(t, centers[keep], bc_type='natural')
        # The tube pinches at the gates, not at the wide circumcenters, so give
        # the radius profile a knot at each gate (at the parameter midpoint
        # between its two vertices) carrying the gate clearance. The centerline
        # keeps only the vertex knots; both splines share the same t domain, so
        # the volume integral samples them consistently and the endpoints (hence
        # the cap radii) are unchanged.
        # Shape-preserving (PCHIP) rather than a cubic spline: the profile is a
        # sequence of measured clearances, and an interpolant that overshoots
        # them writes spheres narrower than the bottleneck it reports - measured
        # at up to 0.14 A below, and 0.18 A above the widest gate. PCHIP is
        # monotone between knots, so the sampled tube is bounded by the numbers
        # the knots carry.
        if len(gates):
            knot_t = np.empty(2 * len(keep) - 1)
            knot_t[0::2] = t
            knot_t[1::2] = 0.5 * (t[:-1] + t[1:])
            knot_r = np.empty_like(knot_t)
            knot_r[0::2] = radii[keep]
            knot_r[1::2] = [float(gates[keep[k]:keep[k + 1]].min())
                            for k in range(len(keep) - 1)]
            radius_spline = PchipInterpolator(knot_t, knot_r)
        else:
            radius_spline = PchipInterpolator(t, radii[keep])

        length = self.calculateChannelLength(centerline_spline)
        volume = self.calculateChannelVolume(centerline_spline, radius_spline)
        
        return centerline_spline, radius_spline, length, bottleneck, volume, centers, radii, gates

    def filterCavities(self, cavities, min_depth):
        return [cavity for cavity in cavities if cavity.depth >= min_depth]

    def filterChannelsByBottleneck(self, cavities, bottleneck):
        for cavity in cavities:
            cavity.channels = [channel for channel in cavity.channels if channel.bottleneck >= bottleneck]
    
    def filterChannelsByVolume(self, cavities, min_volume=None, max_volume=None):
        """Filter channels by volume."""
    
        for cavity in cavities:
            filtered_channels = []
            for channel in cavity.channels:
                if min_volume is not None and channel.volume < min_volume:
                    continue
                if max_volume is not None and channel.volume > max_volume:
                    continue
                filtered_channels.append(channel)
            cavity.channels = filtered_channels

    def filterCavitiesByTetrahedra(self, cavities, min_tetrahedra=None, 
                                   max_tetrahedra=None):
        """Filter cavities by cavity volume."""
    
        filtered = []
        for cavity in cavities:
            n = len(cavity.tetrahedra)
            if min_tetrahedra is not None and n < min_tetrahedra:
                continue
            if max_tetrahedra is not None and n > max_tetrahedra:
                continue
            filtered.append(cavity)
        return filtered

    def calculateTetrahedronVolume(self, a, b, c, d):
        return abs(np.dot(a - d, np.cross(b - d, c - d))) / 6.0

    def calculate_cavity_volumes(self, cavities, simplices, coords):
        """Calculate approximate cavity volumes from Delaunay tetrahedra."""

        for cavity in cavities:
            volume = 0.0
            for tetra in cavity.tetrahedra:
                atom_ids = simplices[tetra]
                a, b, c, d = coords[atom_ids]
                volume += self.calculateTetrahedronVolume(a, b, c, d)
            cavity.volume = volume

    @staticmethod
    def orderCavitiesByVolume(cavities):
        """The cavities largest first, so that cavity 0 is the biggest void.

        Ties are broken by the lowest tetrahedron index a cavity holds, which is
        the order they arrived in, so the numbering is reproducible. Requires
        :meth:`calculate_cavity_volumes` to have run."""

        return sorted(cavities, key=lambda cavity: (
            -cavity.volume, int(np.min(cavity.tetrahedra))))

    def calculateChamberVolumes(self, labels, simplices, coords):
        """Approximate volume of every chamber, indexed by chamber label.

        The same measure :meth:`calculate_cavity_volumes` gives a cavity - its
        Delaunay tetrahedra summed - so that a chamber and a cavity are sized on
        one scale and can share a threshold."""

        if labels.size == 0 or labels.max() < 0:
            return np.zeros(0)

        volumes = np.zeros(int(labels.max()) + 1)
        for tetra, label in enumerate(labels):
            if label >= 0:
                a, b, c, d = coords[simplices[tetra]]
                volumes[label] += self.calculateTetrahedronVolume(a, b, c, d)
        return volumes

    def filterCavitiesByVolume(self, cavities, min_volume=None, max_volume=None):
        """Filter cavities by approximate volume."""

        filtered_cavities = []
        for cavity in cavities:
            if min_volume is not None and cavity.volume < min_volume:
                continue
            if max_volume is not None and cavity.volume > max_volume:
                continue
            filtered_cavities.append(cavity)
        return filtered_cavities

    @staticmethod
    def _channelRecords(channel_index, channel, atom_index, num_samples,
                        label='channel', name_sites=True):
        """The FIL records of one channel: its REMARK, one ATOM per sampled
        sphere and the CONECT bonds between them.

        Written the same way into the combined file and into the per-channel
        one, so the two differ only in which channels they hold and where the
        serial numbers start."""
        centers, radii = _sampleObjectSpheres(channel, num_samples)

        # Each channel gets its own residue number so the channels stay
        # separable at the record level, matching saveCavitiesToPdb.
        lines = [ChannelCalculator._channelRemark(channel_index, channel, label,
                                                  name_sites)]
        for i, (x, y, z, radius) in enumerate(zip(centers[:, 0], centers[:, 1],
                                                  centers[:, 2], radii),
                                              start=atom_index):
            lines.append("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                         % (i, channel_index + 1, x, y, z, 1.00, radius))

        # Bond consecutive samples of THIS channel only, using the global
        # atom serial numbers (start=atom_index). No CONECT spans two
        # channels, so each one is a separate strand in the viewer.
        for i in range(atom_index, atom_index + len(centers) - 1):
            lines.append("CONECT%5d%5d\n" % (i, i + 1))

        return lines, len(centers)

    def saveChannelsToPdb(self, channels, filename, separate=False, num_samples=5,
                          tag='chl', label='channel', separate_path=None,
                          separate_stem=None, name_sites=True):
        # ``channels`` is a flat list, already ordered by cost - that order is
        # the order they are written and numbered here. Each channel is preceded
        # by a REMARK reporting its length, bottleneck radius, curvature and cost.
        # ``tag``/``label`` name the objects in the per-object file names and in
        # that REMARK, so chamber links go through this same writer rather than a
        # copy of it. ``separate_path`` is the name the per-object files are
        # numbered from when it differs from the combined one - the links are
        # collected in ``<stem>_links.pqr`` but numbered from ``<stem>.pqr``, so
        # they come out as ``<stem>_lnk0.pqr`` and not ``<stem>_links_lnk0.pqr``.
        # ``separate_stem=''`` drops the stem from those names altogether, which
        # is what a run told only a directory does.
        #
        # ``name_sites=False`` drops the ``sp<n>`` tag from the names and the
        # REMARKs. It tells the objects of one search site apart from another's,
        # so it says nothing at all when the run had a single site: there every
        # channel would carry the same sp0, which is clutter and not a label.
        # Nothing found, nothing written - not even the empty file that opening
        # for writing would leave, which says only that something ran and cannot
        # be told apart from a write that failed halfway. The count the run found
        # is already reported by the caller.
        if not channels:
            return

        filename = str(filename)
        separate_path = str(separate_path) if separate_path else filename

        # All channels in a single file, one after another in list (cost) order.
        with open(filename, 'w') as pqr_file:
            atom_index = 1
            for channel_index, channel in enumerate(channels):
                lines, samples = self._channelRecords(channel_index, channel,
                                                      atom_index, num_samples,
                                                      label, name_sites)
                pqr_file.writelines(lines)
                pqr_file.write("\n")
                atom_index += samples

        # When separate is set to True also separate PDB/PQR files will be
        # created, one per channel, numbered by the same cost order.
        if separate:
            for channel_index, channel in enumerate(channels):
                # The start point goes in the name, so that the files belonging
                # to one void sort and glob together: <stem>_sp0_chl3.pqr. A
                # link carries the void it arrives at as well, after the number:
                # <stem>_sp13_lnk0_sp9.pqr. Both ends are then greppable - the
                # prefix still gathers everything traced out of sp13, and sp9
                # now matches everything touching sp9 from either side.
                origin = getattr(channel, 'origin', None) if name_sites else None
                destination = (getattr(channel, 'destination', None)
                               if name_sites else None)
                object_tag = (tag if origin is None
                              else 'sp{0}_{1}'.format(origin, tag))
                object_suffix = ('' if origin is None or destination is None
                                 else '_sp{0}'.format(destination))
                channel_filename = _numberedPath(separate_path, object_tag,
                                                 channel_index, object_suffix,
                                                 separate_stem)

                with open(str(channel_filename), 'w') as pqr_file:
                    lines, _ = self._channelRecords(channel_index, channel,
                                                    1, num_samples, label,
                                                    name_sites)
                    pqr_file.writelines(lines)

    @staticmethod
    def _channelRemark(channel_index, channel, label='channel', name_sites=True):
        """One-line PQR/PDB REMARK with a channel's basic geometry:
        length, bottleneck radius, curvature and Dijkstra cost."""
        curv = 'n/a' if np.isnan(channel.curvature) else "%.3f" % channel.curvature
        cost = 'n/a' if channel.cost is None else "%.4g" % channel.cost
        # "from sp13 to sp9" rather than a pair of assignments: both ends of a
        # link are start points, and reading them as a direction is the point.
        # A channel has only the near end - the far end is the solvent. Left off
        # entirely when the run had one site to name (see saveChannelsToPdb).
        origin = getattr(channel, 'origin', None) if name_sites else None
        destination = getattr(channel, 'destination', None) if name_sites else None
        where = '' if origin is None else "  from sp%d" % origin
        if where and destination is not None:
            where += " to sp%d" % destination
        return ("REMARK   %s %d  length=%.3f A  bottleneck=%.3f A  "
                "curvature=%s  cost=%s%s\n" % (
                    label, channel_index, channel.length, channel.bottleneck,
                    curv, cost, where))


    @staticmethod
    def _cavityRemark(cavity_index, cavity):
        """One-line PQR/PDB REMARK with a cavity's basic geometry: volume, depth
        and the number of tetrahedra it holds.

        The counterpart of :meth:`_channelRemark`, and the only place a cavity
        file states its size: the records below carry markers of a fixed radius,
        not a measurement (see :meth:`saveCavitiesToPdb`). The volume is the one
        :meth:`calculate_cavity_volumes` gives, the cavity's Delaunay tetrahedra
        summed."""

        return ("REMARK   cavity %d  volume=%.1f A^3  depth=%.1f A  "
                "tetrahedra=%d\n" % (cavity_index, cavity.volume, cavity.depth,
                                     len(cavity.tetrahedra)))

    @staticmethod
    def _cavityRecords(cavity_index, cavity, centers, atom_index):
        """The FIL records of one cavity: its REMARK and one ATOM per Voronoi
        vertex.

        Written the same way into the combined file and into the per-cavity one,
        so the two differ only in which cavities they hold and where the serial
        numbers start - the arrangement :meth:`_channelRecords` already uses."""

        radius = ChannelCalculator.CAVITY_MARKER_RADIUS
        # Each cavity gets its own residue number so the cavities stay separable
        # at the record level, matching saveChannelsToPdb.
        lines = [ChannelCalculator._cavityRemark(cavity_index, cavity)]
        for i, (x, y, z) in enumerate(centers, start=atom_index):
            lines.append("ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                         % (i, cavity_index + 1, x, y, z, 1.00, radius))

        return lines, len(centers)

    def saveCavitiesToPdb(self, cavities, vertices, filename, separate=False,
                          separate_stem=None):
        """Save surface cavities to a PDB/PQR file as dummy atoms.

        One atom per Voronoi vertex of the cavity, all of them of
        ``CAVITY_MARKER_RADIUS``. The cloud marks where the cavity is; it is not
        drawn at its true width, and the radius column must not be read as one.
        The size of a cavity is in the REMARK.

        That the radius is a constant looks like an oversight, so here is why it
        stays one. The obvious fix is the one the channels use - write the probe
        that fits at each vertex, as :meth:`_channelRecords` does - and it fails
        here. ``max_depth`` keeps only
        the shell of a pocket, so on protonated 3A2M the vertices sit a median 0.2 to
        0.5 A from the exterior the ``surf_radius`` probe carved away, while the
        probe that fits at them is 2.7 A. Those spheres put two thirds to nine
        tenths of their volume outside that surface, half of it outside the convex
        hull, reaching 7.7 A clear of every atom: a pocket drawn as balloons
        hanging off the protein. Clipping them to the surface is no better, since
        the mouth of an open pocket is continuous with the outside - it leaves
        radii of 0.2 to 0.5 A. Sizing each atom by its own tetrahedron instead
        lands within a few percent of the marker on every measure. What would work
        is filling the pocket - the probe spheres intersected with the inside of
        the ``surf_radius`` surface, voxelized - and that is a different file, of
        some 17000 atoms rather than 869.

        ``separate_stem=''`` names the per-cavity files ``cavity0.pqr`` rather
        than ``<stem>_cavity0.pqr``, which is what a run told only a directory
        does."""

        filename = str(filename)

        # (cavity index, cavity, vertices) of every cavity that has tetrahedra, so
        # that the combined file and the per-cavity ones number them alike.
        drawn = []
        for cavity in cavities:
            tetrahedra = cavity.tetrahedra
            if tetrahedra is None or len(tetrahedra) == 0:
                continue
            drawn.append((len(drawn), cavity, vertices[tetrahedra]))

        # As in saveChannelsToPdb: nothing to draw, nothing written. Tested on
        # what is actually drawn rather than on the cavities given, since a
        # cavity with no tetrahedra contributes no atoms either.
        if not drawn:
            return

        with open(filename, 'w') as pqr_file:
            atom_index = 1
            for cavity_index, cavity, centers in drawn:
                lines, written = self._cavityRecords(cavity_index, cavity,
                                                     centers, atom_index)
                pqr_file.writelines(lines)
                atom_index += written

        if separate:
            for cavity_index, cavity, centers in drawn:
                cavity_filename = _numberedPath(filename, 'cavity', cavity_index,
                                                stem=separate_stem)

                with open(str(cavity_filename), 'w') as pqr_file:
                    lines, _ = self._cavityRecords(cavity_index, cavity,
                                                   centers, 1)
                    pqr_file.writelines(lines)
            
    def calculateChannelLength(self, centerline_spline):
        t_values = np.linspace(centerline_spline.x[0], centerline_spline.x[-1], len(centerline_spline.x) * 10)
        points = centerline_spline(t_values)
        diffs = np.diff(points, axis=0)
        lengths = np.linalg.norm(diffs, axis=1)
        return np.sum(lengths)
    
    def calculateChannelVolume(self, centerline_spline, radius_spline):
        # Tube volume  V = \int pi r(t)^2 |x'(t)| dt  evaluated by vectorized
        # composite Simpson instead of an adaptive scipy.quad that called the
        # integrand thousands of times per channel. The centerline/radius are
        # piecewise cubic, so a uniform grid scaled to the number of spline
        # segments (~128 points/segment) reaches ~1e-8 relative accuracy - better
        # than quad's default tolerance - at a fraction of the cost.
        from scipy.integrate import simpson

        t_min = centerline_spline.x[0]
        t_max = centerline_spline.x[-1]
        n_segments = max(1, len(centerline_spline.x) - 1)
        n = n_segments * 128 + 1

        t = np.linspace(t_min, t_max, n)
        r = radius_spline(t)
        speed = np.linalg.norm(centerline_spline(t, 1), axis=1)
        volume = simpson(np.pi * r ** 2 * speed, x=t)

        r_start = radius_spline(t_min)
        r_end = radius_spline(t_max)

        hemisphere_volume_start = (2/3) * np.pi * r_start**3
        hemisphere_volume_end = (2/3) * np.pi * r_end**3

        total_volume = volume + hemisphere_volume_start + hemisphere_volume_end

        return total_volume
            
    def selectSeedTetrahedron(self, cavity, vertices, points, vdw_radii, simp,
                              neighbors, sp, search_radius, min_depth):
        '''Map `sp` to the seed tetrahedron of one cavity.

        The tetrahedron nearest `sp` (the anchor) is frequently a tight one, and every
        channel of the cavity leaves through it: its inscribed radius then caps all of
        their bottlenecks, and the shared first links show up as one common bottleneck
        at the joint beginning of the bundle. So the anchor only says where to look.
        The seed is the widest (largest inscribed radius) tetrahedron of this cavity
        that lies within `search_radius` of `sp`, is buried deeply enough to still be
        inside the site, and is reachable from the anchor through the tetrahedra within
        `search_radius`. Reachability is over the adjacency of the cleared tetrahedra,
        which is free space, so the seed can only move through the void the start point
        sits in and never hops across a wall into a lobe that merely passes nearby; the
        depth floor keeps it from sliding outward towards the mouth, where tetrahedra
        are wide but no longer inside the site. Note that the floor filters the seed,
        not the walk: a marginally shallower cell in between must not wall off the
        wider region behind it.

        The floor is `min_depth`, raised to the widest mouth of the cavity wherever an
        opening is wider than that -- the same rule, and for the same reason, that
        :meth:`setStartingTetrahedraFromChambers` applies to an automatically placed
        seed. It is deliberately not the anchor's own depth. A start point is placed by
        hand, usually on a bound ligand or a catalytic residue, and such a point
        routinely sits deeper than the widest part of the pocket around it; held to the
        anchor's depth the search then finds nothing eligible and leaves the seed on
        the narrow cell nearest the point, capping every channel of the site at that
        radius, which is the outcome the search exists to prevent -- and it fails
        silently, since a pocket whose every route reports the same bottleneck looks
        like a narrow pocket. Anything above the floor is inside the site by the
        measure the module uses everywhere else, so a seed somewhat shallower than the
        point is a seed in the same pocket. How much of the width the search can
        actually reach from there is then a question of `search_radius` alone.

        Where nothing in reach clears the raised floor the plain `min_depth` is tried,
        for a site lying wholly under a wide opening, and failing that the anchor's own
        depth, which always leaves at least the anchor itself.

        The anchor itself competes whatever its depth: the floor bounds where the seed
        may move to, not whether the tetrahedron at the point counts, so the seed is
        only ever moved to a tetrahedron wider than the anchor. Held to the floor, an
        anchor lying just under it would drop out and be traded for the widest
        tetrahedron clearing it, however much narrower that one is. And the anchor's
        depth moves with the probe: a smaller `inner_radius` opens mouths nearer the
        point, so the anchor crosses the floor as the probe shrinks, the seed jumps
        with it, and a site could read as sealed at one `inner_radius` and open at
        the next.

        `search_radius` <= 0 restores the plain nearest-vertex seed.

        :returns: dict of the seed and anchor properties (`seed`, `anchor`, and their
            `_vertex`, `_distance` from `sp`, inscribed `_radius` and `_depth`), plus
            the number of tetrahedra `searched`, how many of them were `eligible`, and
            the depth `floor` they had to clear'''

        from collections import deque

        depths = cavity.tetrahedra_depths
        tet = np.asarray(cavity.tetrahedra)
        d2 = np.sum((vertices[tet] - sp) ** 2, axis=1)
        anchor = int(tet[int(np.argmin(d2))])

        def properties(tetra):
            return dict(
                vertex=vertices[tetra],
                distance=float(np.linalg.norm(vertices[tetra] - sp)),
                radius=float(self.calculateMaxRadius(
                    vertices[tetra], points, vdw_radii, simp[tetra])),
                depth=float(depths.get(tetra, 0.0)))

        def report(seed, searched, eligible, floor):
            info = {'seed': seed, 'anchor': anchor, 'searched': searched,
                    'eligible': eligible, 'floor': float(floor)}
            for name, tetra in (('seed', seed), ('anchor', anchor)):
                for key, value in properties(tetra).items():
                    info['{0}_{1}'.format(name, key)] = value
            return info

        anchor_depth = float(depths.get(anchor, 0.0))

        if not search_radius or search_radius <= 0:
            return report(anchor, 1, 1, anchor_depth)

        def clearances(tetrahedra):
            atoms = simp[tetrahedra]
            return np.min(np.linalg.norm(points[atoms] - vertices[tetrahedra][:, None, :],
                                         axis=2) - vdw_radii[atoms], axis=1)

        near = set(int(t) for t, close in zip(tet, d2 <= search_radius ** 2) if close)

        # BFS from the anchor, staying inside the sphere and inside the cavity.
        reachable = [anchor]
        seen = {anchor}
        queue = deque([anchor])
        while queue:
            current = queue.popleft()
            for neighbor in neighbors[current]:
                neighbor = int(neighbor)
                if neighbor in near and neighbor not in seen:
                    seen.add(neighbor)
                    reachable.append(neighbor)
                    queue.append(neighbor)

        # How far down a seed has to be to be past the openings. The widest mouth
        # of the cavity stands for all of them, since a route out of the seed is
        # cut against whichever opening it reaches. _vertex_clearance is not built
        # yet at seeding time, so the radii are measured here.
        mouth = 0.0
        exits = getattr(cavity, 'exit_tetrahedra', None)
        if exits is not None and len(exits):
            mouth = float(np.max(clearances(np.asarray(exits, dtype=np.intp))))

        reach = np.array(reachable, dtype=np.intp)
        radii = clearances(reach)
        reach_depths = np.array([depths.get(int(t), 0.0) for t in reach])

        # The anchor's depth closes the list because it always keeps the anchor,
        # so some candidate always survives.
        for floor in (max(float(min_depth), mouth), float(min_depth), anchor_depth):
            eligible = np.flatnonzero(reach_depths >= floor)
            if len(eligible):
                break

        # The anchor competes whatever its depth (see above), so the seed only ever
        # moves to a wider tetrahedron; `eligible` still counts only what cleared
        # the floor, which is what the report states. The anchor comes first (BFS
        # order), so a tie goes to it.
        candidates = eligible if eligible[0] == 0 else np.concatenate(([0], eligible))
        best = int(reach[candidates[int(np.argmax(radii[candidates]))]])

        return report(best, len(reachable), len(eligible), floor)


    def setStartingTetrahedraFromPoint(self, cavities, vertices, start_point,
                                       points, vdw_radii, simp, neighbors,
                                       search_radius=5.0, min_depth=5.0):
        '''Restrict the search to the cavity `start_point` names, and seed it there.

        The starting tetrahedron is the widest one `selectSeedTetrahedron` finds in
        the neighbourhood of `start_point`; with ``search_radius=0`` it is simply the
        one whose Voronoi vertex is closest to `start_point`. The cavity holding it is
        the only one returned: a start point says which site to search, and channels
        are computed for that site instead of one bundle per cavity in the structure.

        `search_radius` is a requirement and not a hint. A cavity is a candidate only
        if it holds a tetrahedron within that distance of `start_point`, and where no
        cavity does, none is seeded and no channels are computed. The alternative is
        to seed the nearest tetrahedron however far away it is, which is precisely the
        seed the caller did not ask for: a start point that misses the void -- the
        centroid of a residue selection routinely lands inside an atom -- is then
        answered with channels through whatever cavity happens to lie nearest, and
        nothing in the result says so. Refusing instead, and reporting where the
        nearest tetrahedron actually is, lets the point be corrected. ``search_radius``
        of 0 asks for no neighbourhood at all and so imposes no such distance either.

        :arg cavities: list of cavity objects
        :arg vertices: Voronoi vertices (array of shape (n, 3))
        :arg start_point: point [x, y, z] in Å (list/tuple/ndarray of length 3)
        :arg points: atom coordinates (array of shape (n_atoms, 3)), used to compute
            the inscribed radius of candidate tetrahedra
        :arg vdw_radii: per-atom van der Waals radii (array of shape (n_atoms,))
        :arg simp: simplices (tetrahedron -> its 4 atom indices)
        :arg neighbors: tetrahedron adjacency (tetrahedron -> its 4 neighbours, -1 none)
        :arg search_radius: radius, in Angstrom, of the neighbourhood of `start_point`
            a seed must lie in. 0 disables both the search for a wider seed and the
            distance requirement.
        :type search_radius: float
        :arg min_depth: depth floor, in Angstrom, a seed must clear, so that the search
            for a wider one cannot drift out of the site; raised per cavity to its
            widest mouth. The same value that filters the cavities themselves.
        :type min_depth: float
        :returns: a single-element list holding the selected cavity, or an empty list
            if no cavity has a tetrahedron within `search_radius` of `start_point`'''

        sp = np.asarray(start_point, dtype=float).reshape(3,)

        # Anchor every cavity by distance alone. The cavity the point names is the
        # one holding the nearest tetrahedron, and that is settled here rather than
        # by the seed search: widening moves the seed inside a cavity, it must never
        # decide between cavities. `min` returns the first of any tie, so cavities
        # equidistant from the point resolve in their original order.
        anchored = []
        for cavity in cavities:
            tet = cavity.tetrahedra
            if tet is None or len(tet) == 0:
                continue
            tet = np.asarray(tet)
            distances = np.linalg.norm(vertices[tet] - sp, axis=1)
            nearest = int(np.argmin(distances))
            anchored.append((float(distances[nearest]), cavity, int(tet[nearest])))

        if not anchored:
            _warn("start_point was provided but no cavity contains any "
                "tetrahedron; no channels will be computed.")
            return []

        distance, cavity, tetra = min(anchored, key=lambda entry: entry[0])

        if search_radius and search_radius > 0 and distance > search_radius:
            vertex = vertices[tetra]
            radius = float(self.calculateMaxRadius(vertex, points, vdw_radii,
                                                   simp[tetra]))
            _warn("start_point [{0:.3f}, {1:.3f}, {2:.3f}] has no cavity tetrahedron "
                "within start_point_search={3:.1f} Å; the nearest one is {4:.3f} Å away "
                "at [{5:.3f}, {6:.3f}, {7:.3f}] (inscribed radius {8:.3f} Å). No "
                "channels will be computed. Move start_point towards that vertex, or "
                "raise start_point_search above {4:.3f} Å to accept it."
                .format(sp[0], sp[1], sp[2], float(search_radius), distance,
                        vertex[0], vertex[1], vertex[2], radius))
            return []

        info = self.selectSeedTetrahedron(
            cavity, vertices, points, vdw_radii, simp, neighbors, sp,
            search_radius, min_depth)

        cavity.setStartingTetrahedron(np.array([info['seed']]))
        self.reportSeedTetrahedron(info, search_radius)
        LOGGER.info("    restricting the channel search to the cavity that contains it "
            "({0} tetrahedra, depth {1:.1f} Å).".format(len(cavity.tetrahedra),
                                                        float(cavity.depth)))

        return [cavity]

    def setStartingTetrahedraFromChambers(self, cavities, labels, volumes,
                                          min_depth, seed_volume=50.0,
                                          max_seeds=None):
        '''Seed every cavity at each of its chambers instead of at its single
        deepest tetrahedron.

        :func:`findDeepestTetrahedra` gives a cavity one seed, the far end of its
        geodesic depth field. For a compact pocket that is the site. For a large
        or branched cavity it is arbitrary, and two things follow: every channel
        is forced to be at least as long as the cavity is deep, and the mouths of
        the other lobes go unreported, because a route from the far seed to them
        is absorbed at some nearer mouth on the way and never arrives. Seeding
        each chamber instead lets every lobe reach its own openings, and the
        openings that were already reported get traced from the chamber that
        actually feeds them - on 1tqn at ``inner_radius=1.4`` the four cavity-0
        channels kept their exits but went from 46-58 A to 16-25 A, and their
        bottlenecks widened, because the route no longer has to cross a neck.

        The seed of a chamber is its widest tetrahedron *among those no shallower
        than* `min_depth`. The order matters: a chamber can be deep and still
        have its widest tetrahedron sitting in the mouth, and seeding there
        produces a channel a fraction of an Angstrom long instead of the tunnel
        that chamber really has. Filtering by depth first and taking the widest
        of what remains keeps the seed inside the site.

        How deep is deep enough is not an absolute length, though, because how
        far below a mouth a point has to be to be past it is set by how wide that
        mouth is. `min_depth` is measured in Angstrom while "still in the
        opening" is measured by the opening's own radius, and under a wide one
        the two part company: in a porin at ``surf_radius=10`` the widest mouth
        has a 7.5 A radius, so a tetrahedron 5.3 A down clears `min_depth=5` and
        is still inside the mouth. Seeding there costs the whole site - every
        route out of it leaves by that one opening, the dedup cuts them all
        against it, and the chamber reports a single stub instead of the pore.
        So the depth floor is raised to the widest opening the cavity has, and
        the seed has to clear that. Ordinary mouths measure 1-3 A and sit below
        `min_depth`, which leaves the floor where it was: this only acts once an
        opening is wider than the depth demanded of a seed, which is where the
        absolute floor stops meaning what it says.

        Cavities are left with the seed `findDeepestTetrahedra` gave them when no
        chamber of theirs qualifies, so a pocket too narrow to hold the chamber
        probe still reports its channels.

        :arg labels: chamber of each tetrahedron, from :meth:`findChambers`.
        :type labels: :class:`~numpy.ndarray`

        :arg volumes: approximate volume of each chamber, indexed by label.
        :type volumes: :class:`~numpy.ndarray`

        :arg min_depth: depth floor, in Angstrom, a seed must clear. The same
            value that filters cavities, so "buried enough to be a site" means
            one thing throughout.
        :type min_depth: float

        :arg seed_volume: chambers smaller than this, in cubic Angstrom, are
            ignored as tessellation debris rather than sites, measured on the
            cavity scale (Delaunay tetrahedra summed). It is a floor on where a
            search may start, independent of the ``min_volume`` the caller
            applies to the objects it reports, which is measured on the
            swept-sphere scale instead. ``None`` applies no floor. Default 50.

            The caller applies the same floor to the cavities before handing them
            over, so every cavity here has already cleared it: the whole-cavity
            fallback below decides only *where* in a cavity to start, never
            whether it is a site at all.
        :type seed_volume: float or None

        :arg max_seeds: cap on the seeds of one cavity, widest chamber first.
            ``None`` for no cap. A percolated interior can otherwise offer
            dozens of chambers, and every seed multiplies the candidate routes
            the dedup has to compare.
        :type max_seeds: int or None

        :returns: the number of cavities that were re-seeded, and one report line
            per cavity that has chambers at all. The lines are returned rather
            than logged here so that the caller can print them under the summary
            they belong to; logged from inside this loop they would arrive before
            the summary, the detail of a heading not yet written.
        :rtype: (int, list)'''

        clearance = self._vertex_clearance
        floor = 0.0 if seed_volume is None else float(seed_volume)
        reseeded = 0
        notes = []

        for index, cavity in enumerate(cavities):
            tetrahedra = np.asarray(cavity.tetrahedra, dtype=np.intp)
            depths = cavity.tetrahedra_depths

            # Every chamber lying in this cavity, before either floor is applied,
            # so that the report can say how many of them a seed stands for.
            cavity.chambers_found = len(set(int(label)
                                            for label in labels[tetrahedra]
                                            if label >= 0))

            chambers = {}
            for tetra, label in zip(tetrahedra, labels[tetrahedra]):
                if label >= 0 and volumes[label] >= floor:
                    chambers.setdefault(int(label), []).append(int(tetra))

            # The depth a seed has to clear to be past the mouths, which is
            # min_depth until an opening is wider than that (see the note above).
            # The widest opening of the cavity stands for all of them: a seed is
            # cut against whichever opening its routes reach, not only the
            # nearest, so clearing the widest is what puts it past every one.
            # This is the same radius the dedup measures an opening by, so "below
            # the mouth" means one thing on both sides.
            mouth = 0.0
            exit_tetrahedra = getattr(cavity, 'exit_tetrahedra', None)
            if exit_tetrahedra is not None and len(exit_tetrahedra):
                mouth = float(np.max(
                    clearance[np.asarray(exit_tetrahedra, dtype=np.intp)]))
            seed_depth = max(float(min_depth), mouth)

            def widestDeepest(members):
                """The widest of `members` that lies below the mouths: no
                shallower than `seed_depth`, falling back to `min_depth` where a
                chamber holds nothing that deep, so that a site whose whole
                extent lies under a wide opening still reports its channels."""
                members = np.asarray(members, dtype=np.intp)
                member_depths = np.array([depths.get(int(t), 0.0)
                                          for t in members])
                deep = members[member_depths >= seed_depth]
                if len(deep) == 0:
                    deep = members[member_depths >= min_depth]
                if len(deep) == 0:
                    return None
                return int(deep[np.argmax(clearance[deep])])

            seeds = {}
            seed_volumes = {}
            for label, members in chambers.items():
                seed = widestDeepest(members)
                if seed is None:
                    continue
                seeds[seed] = label
                seed_volumes[seed] = float(volumes[label])

            found = cavity.chambers_found
            plural = 'chamber' if found == 1 else 'chambers'

            if not seeds:
                if found:
                    notes.append("cavity {0}: {1} {2}, none of them deep and "
                                 "large enough to seed; searched whole."
                                 .format(index, found, plural))
                # No chamber qualified, so the cavity stands as its own void and
                # is seeded by the same rule, over the whole of it. Left alone it
                # would keep the deepest tetrahedron findDeepestTetrahedra picked,
                # which is chosen without reference to width - and since every
                # channel of the cavity leaves through its seed, a narrow one caps
                # all of their bottlenecks. The cavity has already cleared
                # seed_volume in the caller, so this is a choice of seed and not a
                # way back in for a void the floor excluded.
                seed = widestDeepest(tetrahedra)
                if seed is not None:
                    cavity.setStartingTetrahedron(np.array([seed]))
                continue

            if max_seeds is not None and len(seeds) > max_seeds:
                # Largest first, on the volume the chambers were filtered by
                largest = sorted(seeds, key=lambda t: -seed_volumes[t])[:max_seeds]
                notes.append("cavity {0}: {1} {2}, {3} of them qualify as sites, "
                             "the {4} largest seeded (max_seeds={4})."
                             .format(index, found, plural, len(seeds), max_seeds))
                seeds = {t: seeds[t] for t in largest}
            elif len(seeds) == found:
                notes.append("cavity {0}: {1} {2}, {3}seeded."
                             .format(index, found, plural,
                                     '' if found == 1 else 'all '))
            else:
                notes.append("cavity {0}: {1} {2}, {3} of them seeded."
                             .format(index, found, plural, len(seeds)))

            # Sorted so that the seed order, and with it the order candidates
            # reach the dedup, does not depend on the dict.
            cavity.setStartingTetrahedron(np.array(sorted(seeds)))
            # Which chamber each seed speaks for; the chamber links read it.
            cavity.seed_chambers = seeds
            cavity.seed_volumes = {seed: seed_volumes[seed] for seed in seeds}
            # How close each of them comes to the surface anywhere along itself,
            # which is what says whether a chamber has its own way out. Its seed
            # cannot say: the seed is the widest tetrahedron at least min_depth
            # deep, so it is buried by construction, and a large chamber that
            # opens straight to the solvent can still hold the deepest seed of
            # the cavity.
            cavity.chamber_depths = {
                label: min(depths.get(int(t), 0.0) for t in chambers[label])
                for label in seeds.values()}
            reseeded += 1

        return reseeded, notes

    def reportSeedTetrahedron(self, info, search_radius):
        '''Log the seed tetrahedron `selectSeedTetrahedron` picked, and, when it is not
        the one nearest the start point, the anchor it replaced -- the two radii are what
        tell the user whether the seed was capping the bottlenecks of the cavity.'''

        LOGGER.info("start_point seeded at tetrahedron {0} (Voronoi vertex at "
            "[{1:.3f}, {2:.3f}, {3:.3f}], {4:.3f} Å from start_point, inscribed radius "
            "{5:.3f} Å, depth {6:.1f} Å)."
            .format(info['seed'], info['seed_vertex'][0], info['seed_vertex'][1],
                    info['seed_vertex'][2], info['seed_distance'], info['seed_radius'],
                    info['seed_depth']))

        if info['seed'] != info['anchor']:
            LOGGER.info("    widened from the nearest tetrahedron {0} ({1:.3f} Å away, "
                "inscribed radius {2:.3f} Å, depth {3:.1f} Å), the widest of the {4} tetrahedra "
                "at least {5:.1f} Å deep among the {6} reachable within {7:.1f} Å; seeding "
                "the narrow one would have capped every channel here at its radius."
                .format(info['anchor'], info['anchor_distance'], info['anchor_radius'],
                        info['anchor_depth'], info['eligible'], info['floor'],
                        info['searched'], float(search_radius)))
        elif search_radius and search_radius > 0:
            if info['anchor_depth'] < info['floor']:
                LOGGER.info("    kept although shallower than {0:.1f} Å: already as wide "
                    "as any of the {1} tetrahedra at least that deep among the {2} "
                    "reachable within {3:.1f} Å."
                    .format(info['floor'], info['eligible'], info['searched'],
                            float(search_radius)))
            else:
                LOGGER.info("    already the widest of the {0} tetrahedra at least {1:.1f} Å "
                    "deep among the {2} reachable within {3:.1f} Å."
                    .format(info['eligible'], info['floor'], info['searched'],
                            float(search_radius)))


    def trimCavitiesByDepth(self, cavities, max_depth):
        """Filtering cavities by max_depth."""
    
        for cavity in cavities:
            cavity.tetrahedra = np.array([
                tetra for tetra in cavity.tetrahedra
                if cavity.tetrahedra_depths.get(tetra, np.inf) <= max_depth])




def prepareChannelRecords(channels, atoms=None, trajectory=None, metric=None, voxelize=False,
                          max_proc=None, mp_context=None, **kwargs):
    """
    Build a flat list of ChannelRecords (one per channel, frame-major order).

    Parameters
    ----------
    channels : list[list[Channel]]
        One list of channels per frame.
    atoms, trajectory :
        If both are given, lining residues are computed and stored on each
        record (record.lining_residues).

    metric : str, callable, or iterable of str/callable, optional
        Metric or metrics the records should be prepared for.

        Recognized preparation requirements:

            chamfer
            frechet
            projection
                Basic ChannelRecord geometry only.
            centerline
                Additionally prepare representative points along the centerline.
            voxel
            voxel_overlap
                Additionally voxelize all records onto a common lattice.

            lining
            residue
            residue_lining
                Additionally calculate lining residues.

        Callables are accepted but do not imply any additional preparation,
        since their geometry requirements cannot be inferred. For example,
        if a callable requires voxels, pass:

            metric=['voxel', my_metric]
    voxelize : bool
        If True, every record is voxelized onto a shared lattice
        (common global_origin + voxel_size).
    max_proc : int or None
        Worker processes for record building. None -> cpu_count // 2.
    mp_context : str or None
        Multiprocessing context, e.g. 'fork' or 'spawn'.
    **kwargs :
        ClusterCalculator options (voxel_size, pad_factor, h, angle_deg, q,
        n_samples, alpha, beta) are consumed here; everything else is
        forwarded to getObjectLiningInfoMultipleFrames.

    Returns
    -------
    records : list[ClusterCalculator.ChannelRecord]
    """
    import multiprocessing

    if channels is None or len(channels) == 0:
        return []
    
    if metric is None:
        metrics = []

    elif isinstance(metric, str) or callable(metric):
        metrics = [metric]

    else:
        try:
            metrics = list(metric)
        except TypeError:
            raise TypeError("metric must be a string, callable, iterable of metrics, or None")

    sparse_centerline_metrics = {'centerline'} 
    dense_centerline_metrics = {'chamfer', 'frechet','projection'}
    voxel_metrics = {'voxel','voxel_overlap',}
    lining_metrics = { 'lining','residue','residue_lining', }
    known_metrics = ( sparse_centerline_metrics | dense_centerline_metrics | voxel_metrics | lining_metrics)
    string_metrics = set()

    for m in metrics:
        if callable(m):
            continue
        if not isinstance(m, str):
            raise TypeError( "each metric must be a string or callable")
        m = m.lower()
        if m not in known_metrics:
            raise ValueError(f"unknown metric {m!r}; expected one of {sorted(known_metrics)} or a callable")
        string_metrics.add(m)

    prepare_rep_points = bool(string_metrics & sparse_centerline_metrics)
    prepare_voxels = (voxelize or bool(string_metrics & voxel_metrics))
    prepare_lining = bool(string_metrics & lining_metrics)
    if metric is None and atoms is not None and trajectory is not None:
        prepare_lining = True
    if prepare_lining:
        if atoms is None or trajectory is None:
            raise ValueError("atoms and trajectory are required when preparing records for a residue-lining metric")

    calc_keys = ('voxel_size', 'pad_factor', 'h', 'angle_deg', 'q',
                 'n_samples', 'alpha', 'beta')
    calc_kwargs = {k: kwargs.pop(k) for k in calc_keys if k in kwargs}
    clustercalc = ClusterCalculator(atoms=atoms,**calc_kwargs)

    # Flatten channels, tracking which frame each belongs to
    flat_channels = []
    frame_ids = []
    for frame_id, channel_group in enumerate(channels):
        flat_channels.extend(channel_group)
        frame_ids.extend([frame_id] * len(channel_group))
    n = len(flat_channels)
    if n == 0:
        return []

    # Multiprocessing setup
    if max_proc is None:
        max_proc = max(1, multiprocessing.cpu_count() // 2)
    max_proc = max(1, int(max_proc))
    ctx = multiprocessing.get_context(mp_context) if mp_context else multiprocessing.get_context()

    # A shared lattice origin is required before any voxelization
    if prepare_voxels:
        raw_bboxes = [clustercalc.channelBoundBox(c) for c in flat_channels]
        clustercalc.global_origin = np.min([bb[0] for bb in raw_bboxes], axis=0)
        LOGGER.info("All Bounding Boxes Built")

    # Build records (optionally voxelized)
    if max_proc == 1 or n == 1:
        records = [clustercalc._buildChannelRecord(ch, fid, voxelize=prepare_voxels)
                   for fid, ch in zip(frame_ids, flat_channels)]
    else:
        tasks = [(ch, fid, prepare_voxels) for fid, ch in zip(frame_ids, flat_channels)]
        chunksize = max(1, n // (max_proc * 4))
        with ctx.Pool(processes=max_proc) as pool:
            records = pool.map(clustercalc._buildChannelRecordWorker, tasks,
                               chunksize=chunksize)
        # Splines may not survive pickling; rebuild if missing
        for record in records:
            if not record.channel.getSplines()[0]:
                record.channel._buildSplines()
    LOGGER.info("All Channels Voxelized" if prepare_voxels else "All Channel Records Built")

    # Optional lining residues

    if prepare_lining:
        resinds = getObjectLiningInfoMultipleFrames(
            atoms, channels, trajectory,
            callbacks={'ResIndices': getLiningResIndices},
            **kwargs)

        # Each frame's result is a dict keyed by whatever ID the routine uses
        # (not the Channel objects), so iterate the dict's own keys,
        # exactly as the original code did.
        residue_sets = [frame_info[key]['ResIndices']
                        for frame_info in resinds
                        for key in frame_info]

        if len(residue_sets) != len(records):
            raise ValueError(
                f"Lining info returned {len(residue_sets)} channels but "
                f"{len(records)} records were built; cannot match them up.")

        for residue_set, record in zip(residue_sets, records):
            record.lining_residues = residue_set
            record.packed_lining = clustercalc._pack(residue_set)
        LOGGER.info("Lining residues assigned")
    
    if prepare_rep_points:
        records = clustercalc.prepareCenterlineRecords(records)

    return records

def _scanGrid(scan, allowed, defaults=None):
    """Value tuples in `allowed` order. Missing keys take `defaults`; unknown keys raise."""
    from itertools import product
    unknown = set(scan) - set(allowed)
    if unknown:
        raise ValueError(f"unknown scan keys {sorted(unknown)}; allowed: {list(allowed)}")
    defaults = defaults or {}
    axes = []
    for k in allowed:
        if k in scan:
            vals = scan[k]
        elif k in defaults:
            vals = [defaults[k]]
        else:
            raise ValueError(f"scan is missing {k!r}")
        axes.append(np.atleast_1d(vals).astype(float).tolist())
    return list(product(*axes))

def buildCenterlineDistanceMatrix(channels=None, records=None, h=10, n_samples=5, angle_deg=5.0,
                                    q=1.0,dtype=np.float32,block_size=None,max_proc=None,
                                    normalize_radii=True,clip_distance=True,**kwargs ):
    """Build a pairwise distance matrix between channel centerlines.

    This routine accepts either a nested collection of Channel objects or a
    precomputed sequence of ChannelRecord instances, converts them into
    representative centerline records, and computes the geometric distances
    between all channels. The result is an NxN NumPy array where each entry is a
    distance between the corresponding pair of channel centerlines.

    Parameters
    ----------
    channels : sequence of sequence of Channel, optional
        Nested channel groups, typically one list per conformation. Each channel
        must have a valid cost and centerline spline before clustering.
    records : sequence of ChannelRecord, optional
        Prebuilt channel records. Pass either this or ``channels``, but not both.
    h : int, optional
        Number of representative points used along the centerline for distance
        evaluation.
    n_samples : int, optional
        Number of samples used when constructing the spline-based geometric
        representation.
    angle_deg : float, optional
        Angular tolerance used during centerline distance calculations.
    q : float, optional
        Distance weighting exponent.
    dtype : numpy.dtype, optional
        Data type of the output distance matrix.
    block_size : int, optional
        Number of channel pairs to process per block. If None, all channels are
        processed in one block.
    max_proc : int, optional
        Maximum number of worker processes used when evaluating the matrix.
    normalize_radii : bool, optional
        Whether to normalize radial contributions in the distance metric.
    clip_distance : bool, optional
        Whether to clip computed distances to a maximum value.
    **kwargs
        Additional keyword arguments passed through to lower-level routines.

    Returns
    -------
    numpy.ndarray or None
        An NxN symmetric distance matrix for the input channels, or None when no
        channels or records are provided.
    """

    import multiprocessing
    from multiprocessing import Pool
    scan = kwargs.pop('scan', None)
    clustercalc = ClusterCalculator(h=h, angle_deg=angle_deg, q=q, n_samples=n_samples)
    # Accept either nested channels or a flat list of ChannelRecords.
    if records is not None and channels is not None:
        raise ValueError("Pass either channels or records, not both.")

    if records is None:
        if channels is None or len(channels) == 0:
            return None

        # Flatten the list of lists of channels into records.
        records = []

        for conf_id, channel_group in enumerate(channels):
            for channel in channel_group:
                if channel.cost is None:
                    raise ValueError("Channel.cost is None -- required to rank channels.")

                if channel.centerline_spline is None:
                    raise ValueError("Channel.centerline_spline is None -- call " 
                                        "channel._buildSplines() before clustering.")

                records.append(clustercalc._buildChannelRecord(channel=channel, conformation_id=conf_id))

    else:
        # Materialize the input in case records is an iterator.
        records = list(records)
        if not records:
            return None
        if not np.all(isinstance(record, clustercalc.ChannelRecord) for record in records):
            raise TypeError("records must contain only ChannelRecord instances.")

    if not records:
        return None
    S = clustercalc.computeAverageStart(records)
    clustercalc.computeRepresentativePoints(records, S, h=h)
    n = len(records)
    if n == 0:
        return np.empty((0, 0), dtype=dtype)

    h = min(len(r.rep_points) for r in records)
    shell_data = clustercalc._prepareAllShells(records, h)   # independent of q
    centroids, radii_means = clustercalc._prepareChannelCentroids(records)

    if block_size is None:
        block_size = n
    block_size = max(1, int(block_size))
    blocks = [(start, min(start + block_size, n)) for start in range(0, n, block_size)]

    if max_proc is None:
        max_proc = max(1, multiprocessing.cpu_count() // 2)
    mp_context = kwargs.get('mp_context')
    ctx = multiprocessing.get_context(mp_context) if mp_context else multiprocessing.get_context()

    def _buildD(q_value):
        weights = clustercalc.computeDistanceWeights(h, q_value).astype(np.float32)
        initargs = (shell_data, weights, n, centroids, radii_means,
                    normalize_radii, clip_distance)
        D = np.empty((n, n), dtype=dtype)
        if max_proc == 1:
            clustercalc._initDistanceWorker(*initargs)
            for block in blocks:
                start, end, Db = clustercalc._distanceBlock(block)
                D[start:end] = Db
        else:
            with ctx.Pool(processes=max_proc, initializer=clustercalc._initDistanceWorker,
                          initargs=initargs) as pool:
                for start, end, Db in pool.imap(clustercalc._distanceBlock, blocks):
                    D[start:end] = Db
        D = np.minimum(D, D.T)          # enforce exact numerical symmetry
        np.fill_diagonal(D, 0.0)
        return D

    if scan:
        results = {(qv,): _buildD(qv) for (qv,) in _scanGrid(scan, ['q'])}
        return results, records
    
    return _buildD(q), records


def buildVoxelOverlapDistanceMatrix(channels=None, records=None, voxel_size=0.25,
                                    max_proc=1, mp_context=None, pad_factor=1.02,
                                    interframe_only=False, **kwargs):
    """
    Calculate pairwise spatial distance (1 - weighted voxel overlap score)
    between channels.

    Every channel is voxelized exactly once onto a shared lattice (common
    origin + voxel pitch). Pairwise overlap is then integer-indexed slicing
    and a bitwise AND between cached masks.

    Parameters
    ----------
    channels : list[list[Channel]] or None
        One list of channels per frame. Mutually exclusive with `records`.
    records : list[ChannelRecord] or None
        Pre-built records (e.g. from prepareChannelRecords). Frame membership
        is taken from each record's `conformation_id`. Unvoxelized records
        are voxelized here onto a shared lattice.
    voxel_size : float
        Voxel pitch. Ignored when records are already voxelized (their
        lattice is used instead).
    max_proc : int or None
        Worker processes for the per-channel voxelization step.
        None -> cpu_count // 2.
    mp_context : str or None
        Multiprocessing context, e.g. 'fork'.
    pad_factor : float
        Passed to ClusterCalculator.
    interframe_only : bool
        If True, pairs within the same frame are never compared and their
        entries in D are NaN. Cross-frame pairs that do not overlap are
        1.0 (max distance), never NaN.

    Returns
    -------
    D : ndarray (n, n)
        Symmetric distance matrix. Diagonal is 0.
    records : list[ChannelRecord]
        Voxelized records, in the same order as the rows/columns of D.
    """
    import multiprocessing

    if records is not None and channels is not None:
        raise ValueError("Pass either channels or records, not both.")

    alpha = kwargs.pop('alpha', 1.0)
    beta = kwargs.pop('beta', 2.0)
    scan = kwargs.pop('scan', None)
    clustercalc = ClusterCalculator(voxel_size=voxel_size, pad_factor=pad_factor,
                                    alpha=alpha, beta=beta)

    if max_proc is None:
        max_proc = max(1, multiprocessing.cpu_count() // 2)
    max_proc = max(1, int(max_proc))
    ctx = (multiprocessing.get_context(mp_context) if mp_context
           else multiprocessing.get_context())

    if records is None:
        if channels is None or len(channels) == 0:
            return None
        records = prepareChannelRecords(channels, voxelize=True,
                                        max_proc=max_proc, mp_context=mp_context,
                                        voxel_size=voxel_size, pad_factor=pad_factor,
                                        alpha=alpha, beta=beta)
    else:
        records = list(records)
        if not records:
            return None
        if not all(isinstance(rec, ClusterCalculator.ChannelRecord) for rec in records):
            raise TypeError("records must contain only ChannelRecord instances.")

        if all(rec.voxelized for rec in records):
            # Adopt the lattice the records were voxelized on, and make sure
            # they all share it - masks are only comparable on one lattice.
            ref = records[0].calculator
            for rec in records:
                cal = rec.calculator
                if (not np.isclose(cal.voxel_size, ref.voxel_size) or
                        not np.allclose(cal.global_origin, ref.global_origin)):
                    raise ValueError("Records were voxelized on different lattices "
                                     "(voxel_size / global_origin differ).")
            clustercalc.voxel_size = ref.voxel_size
            clustercalc.global_origin = ref.global_origin
        else:
            # Any unvoxelized record -> (re)voxelize everything on one lattice.
            clustercalc.global_origin = np.min([rec.bbox[0] for rec in records], axis=0)
            for rec in records:
                rec.calculator = clustercalc

            if max_proc == 1 or len(records) == 1:
                for rec in records:
                    rec._prepareVoxels()
            else:
                chunksize = max(1, len(records) // (max_proc * 4))
                with ctx.Pool(processes=max_proc) as pool:
                    records = pool.map(clustercalc._prepareVoxelsWorker, records,
                                       chunksize=chunksize)
                for rec in records:
                    if not rec.channel.getSplines()[0]:
                        rec.channel._buildSplines()
            LOGGER.info("All Records Voxelized")

    n = len(records)
    conf_ids = np.array([rec.conformation_id for rec in records])
    unique_frames = np.unique(conf_ids)
    num_frames = len(unique_frames)
    # Index arrays (not slices) so record order doesn't have to be frame-sorted
    frame_indices = [np.nonzero(conf_ids == f)[0] for f in unique_frames]

    D = np.ones((n, n), dtype=float)
    if interframe_only:
        for idx in frame_indices:
            D[np.ix_(idx, idx)] = np.nan      # same-frame = "not computed"
    np.fill_diagonal(D, 0.0)

    bboxes_mn = np.array([rec.bbox[0] for rec in records])
    bboxes_mx = np.array([rec.bbox[1] for rec in records])

    tasks = []
    for a in range(num_frames):
        for b in range(a, num_frames):
            if a == b and interframe_only:
                continue
            idx_a, idx_b = frame_indices[a], frame_indices[b]

            overlap = clustercalc._pairwiseBBoxOverlap(
                bboxes_mn[idx_a], bboxes_mx[idx_a],
                bboxes_mn[idx_b], bboxes_mx[idx_b])
            overlap = np.asarray(overlap, dtype=bool)

            if a == b:
                overlap &= np.triu(np.ones(overlap.shape, dtype=bool), k=1)

            local_i, local_j = np.nonzero(overlap)
            for i, j in zip(idx_a[local_i].tolist(), idx_b[local_j].tolist()):
                tasks.append((i, j, records[i], records[j]))
    
    counts = [clustercalc._voxelOverlapCounts(task) for task in tasks]
    LOGGER.info("All Channel Overlap Counts Computed")

    cols = np.array(counts, dtype=float).reshape(-1, 5)     # safe when there are no tasks
    ii, jj = cols[:, 0].astype(int), cols[:, 1].astype(int)
    inter, a, b = cols[:, 2], cols[:, 3], cols[:, 4]

    def _fill(alpha_v, beta_v):
        Dk = D.copy()
        score = clustercalc._voxelScoreFromCounts(inter, a, b, alpha_v, beta_v)
        Dk[ii, jj] = score
        Dk[jj, ii] = score
        return Dk

    if scan:
        grid = _scanGrid(scan, ['alpha', 'beta'], defaults={'alpha': alpha, 'beta': beta})
        return {(av, bv): _fill(av, bv) for av, bv in grid}, records

    return _fill(alpha, beta), records


def buildLiningSequenceDistanceMatrix(atoms,channels=None,trajectory=None,records=None,**kwargs):
    """
    Build a distance matrix for channel lining residue sets.

    Parameters
    ----------
    atoms : :class:`~prody.atomic.Atomic`
        Atomic structure or trajectory used to determine lining residues.
    channels : list, optional
        Channel groups to process. If provided, ``trajectory`` must also be
        supplied so that lining residues can be computed for each frame.
    trajectory : object, optional
        Trajectory associated with ``channels`` when computing lining residue
        information across conformations.
    records : list, optional
        Precomputed channel records. Pass either ``channels`` or ``records``,
        not both.
    **kwargs
        Additional keyword arguments forwarded to the lining-residue extraction
        machinery and to the clustering parameters used to build the channel
        records.

    Returns
    -------
    tuple
        A pair ``(distance_matrix, records)`` where ``distance_matrix`` is the
        Tversky distance matrix between lining residue sets and ``records`` is
        the list of channel records used to build it.

    Notes
    -----
    This function encodes each lining residue set as a packed binary vector
    over the unique residue indices in ``atoms`` and computes pairwise Tversky
    distances between all channel records.
    """

    if records is not None and channels is not None:
        raise ValueError("Pass either channels or records, not both.")
    if channels and not trajectory:
        raise ValueError("If channels is passed, trajectory must also be passed") 
    
    alpha = kwargs.pop('alpha',0.5)
    beta = kwargs.pop('beta',0.5)
    scan = kwargs.pop('scan',None)
    clustercalc = ClusterCalculator(atoms = atoms, alpha=alpha,beta=beta)
    # clustercalc.n_residues = len(np.unique(atoms.getResindices()))
    if records is None:
        if channels is None or len(channels) == 0:
            return None
        records = []
        for  conf_id, channel_group in enumerate(channels):
            for channel in channel_group:
                records.append(clustercalc._buildChannelRecord(channel=channel, conformation_id=conf_id))

        resinds = getObjectLiningInfoMultipleFrames(atoms,channels,trajectory,
                                                    callbacks={'ResIndices':getLiningResIndices},
                                                    **kwargs)
        residue_sets = [channel_group[channel]['ResIndices'] for channel_group in resinds for channel in channel_group]
        for residue_set, record in zip(residue_sets,records):
            record.lining_residues = residue_set

        packed = np.array([clustercalc._pack(residue_set) for residue_set in residue_sets])
    else:
        packed = np.array([clustercalc._pack(rec.lining_residues) for rec in records])

    if scan:
        inter, lo, hi = clustercalc._TverskyCounts(packed)
        results = {(a, b): clustercalc._TverskyFromCounts(inter, lo, hi, a, b)
                   for a, b in _scanGrid(scan, ['alpha', 'beta'],defaults={'alpha':alpha,'beta':beta})}
        return results, records
    
    return clustercalc._TverskyDistanceMatrix(packed), records


def _sklearnToLinkage(model):
    """Convert a fitted sklearn AgglomerativeClustering into a scipy-style
    linkage matrix [child_1, child_2, distance, sample_count]."""
    if getattr(model, "distances_", None) is None:
        return None
    n_samples = len(model.labels_)
    counts = np.zeros(model.children_.shape[0])
    for i, merge in enumerate(model.children_):
        c = 0
        for child in merge:
            c += 1 if child < n_samples else counts[child - n_samples]
        counts[i] = c
    return np.column_stack([model.children_, model.distances_, counts]).astype(float)

def _fitAgglomerative(distance_matrix, **kwargs):
    """Returns (model, labels, linkage_matrix). Labels are 0-based."""
    n_clusters = kwargs.pop("n_clusters", None)
    distance_threshold = kwargs.pop("distance_threshold", None)
    linkage_name = kwargs.pop("linkage", "average")

    # sklearn: n_clusters and distance_threshold are mutually exclusive
    if n_clusters is not None:
        distance_threshold = None
    elif distance_threshold is None:
        distance_threshold = 0.4  # assumes distances on a 0-1 scale

    try:
        from sklearn.cluster import AgglomerativeClustering
    except ImportError:
        AgglomerativeClustering = None

    if AgglomerativeClustering is not None:
        try:
            model = AgglomerativeClustering(n_clusters=n_clusters,
                                            distance_threshold=distance_threshold,
                                            metric="precomputed",
                                            linkage=linkage_name,
                                            compute_distances=True,
                                            **kwargs)
            model.fit(distance_matrix)
        except Exception as exc:
            raise RuntimeError(f"AgglomerativeClustering fitting failed: {exc}") from exc
        return model, model.labels_, _sklearnToLinkage(model)

    try:
        from scipy.cluster.hierarchy import linkage as sp_linkage, fcluster
        from scipy.spatial.distance import squareform
    except ImportError as exc:
        raise ImportError("AgglomerativeClustering requires scikit-learn or scipy") from exc

    Z = sp_linkage(squareform(np.asarray(distance_matrix), checks=False),
                   method=linkage_name)
    if n_clusters is not None:
        labels = fcluster(Z, t=n_clusters, criterion="maxclust") - 1
    else:
        labels = fcluster(Z, t=distance_threshold, criterion="distance") - 1
    return "scipy_hierarchical", labels, Z


def clusterChannels(distance_matrix, method, records, update_record_ids=True, **kwargs):
    """
    Cluster channels/samples from a precomputed distance matrix.

    Parameters
    ----------
    distance_matrix : array-like (n_samples, n_samples)
    method : {'HDBSCAN', 'DBSCAN', 'agglomerative'}
    records : sequence of ChannelRecord, length n_samples
    update_record_ids : bool, default True
        If True, append each label to record.cluster_ids. subcluster() passes
        False and writes the full hierarchical path itself.
    **kwargs
        sparse : bool, default False
        Everything else is forwarded to the estimator
        (HDBSCAN: min_cluster_size, min_samples, ...; DBSCAN: eps, min_samples;
        agglomerative: n_clusters, linkage, distance_threshold).

    Returns
    -------
    labels : ndarray (n_samples,)   (noise = -1 for HDBSCAN/DBSCAN)
    clustering : ClusteringRecords
        linkage_matrix is populated for agglomerative.
        Noise points (label -1) get their own Cluster with label -1.
    """
    from scipy.sparse import csr_matrix

    sparse = kwargs.pop("sparse", False)
    if not sparse:
        distance_matrix = np.array(distance_matrix, copy=True)

    shape = getattr(distance_matrix, "shape", None)
    if shape is None or len(shape) != 2:
        raise ValueError("distance_matrix must be a 2D array-like object")
    n_rows, n_cols = shape
    if n_rows != n_cols:
        raise ValueError( f"distance_matrix must be square (n_samples x n_samples); got shape {shape}")
    if len(records) != n_rows:
        raise ValueError( f"distance_matrix has {n_rows} rows but {len(records)} records were given")

    valid_methods = {"HDBSCAN", "DBSCAN", "agglomerative"}
    if method not in valid_methods:
        raise ValueError(f"Unknown method '{method}'. Expected one of {sorted(valid_methods)}")

    if sparse:
        distance_matrix = csr_matrix(distance_matrix)

    linkage_matrix = None
    if method == "HDBSCAN":
        try:
            from sklearn.cluster import HDBSCAN
        except ImportError as exc:
            raise ImportError("HDBSCAN requires scikit-learn") from exc
        kwargs.setdefault("metric", "precomputed")
        try:
            model = HDBSCAN(**kwargs)
            model.fit(distance_matrix)
        except Exception as exc:
            raise RuntimeError(f"HDBSCAN fitting failed: {exc}") from exc
        labels = model.labels_
    elif method == "DBSCAN":
        try:
            from sklearn.cluster import DBSCAN
        except ImportError as exc:
            raise ImportError("DBSCAN requires scikit-learn") from exc
        kwargs.setdefault("metric", "precomputed")
        kwargs.setdefault("eps", 0.5)
        kwargs.setdefault("min_samples", 5)
        try:
            model = DBSCAN(**kwargs)
            model.fit(distance_matrix)
        except Exception as exc:
            raise RuntimeError(f"DBSCAN fitting failed: {exc}") from exc
        labels = model.labels_
    else:
        model, labels, linkage_matrix = _fitAgglomerative(distance_matrix, **kwargs)

    labels = np.asarray(labels)
    if update_record_ids:
        for rec, label in zip(records, labels):
            rec.cluster_ids.append(int(label))

    # Stable sort so records inside a cluster keep their original order; this
    # keeps cluster.records aligned with distance_matrix[np.ix_(mask, mask)].
    order = np.argsort(labels, kind="stable")
    sorted_recs = np.empty(len(records), dtype=object)
    for i, j in enumerate(order):
        sorted_recs[i] = records[j]

    unique, counts = np.unique(labels, return_counts=True)
    clusters = np.empty(len(unique), dtype=object)
    start = 0
    for k, (label, count) in enumerate(zip(unique, counts)):
        end = start + count
        cl = Cluster(int(label), sorted_recs[start:end])
        cl.details["is_noise"] = bool(label == -1)
        cl._calcInternalDistanceMatrix(distance_matrix, labels)
        cl._calcClusterStats(labels)
        clusters[k] = cl
        start = end

    clustering = ClusteringRecords(model, records, clusters, labels, distance_matrix,
                                   linkage_matrix=linkage_matrix)
    return labels, clustering

def renameChannelFilesWithClusterLabels(channels, labels, output_path,
                                        stem='channels', ext='pqr',
                                        dry_run=False):
    """Tag calcChannelsMultipleFrames()'s per-channel files with a cluster id.

    Appends "_clid<label>" to the name of every channel file
    calcChannelsMultipleFrames() already wrote to disk, so a downstream viewer
    (see vis_channels_traj.py) can group and colour channels by cluster
    identity across frames instead of by their per-frame rank. Rank is not a
    stable identity: cavity detection runs independently per frame, so "rank 2"
    in frame 50 is not necessarily the same physical channel as "rank 2" in
    frame 0, whereas a cluster label computed across all frames at once is.

    Files must have been generated from a run with separate=True.
    mmCIF is not supported and raises immediately - every
    channel of a frame lives inside ONE shared .cif file there (one row per
    channel in the _sb_ncbr_channel_profile loop, rather than one file per
    channel), so there is no single file to rename per channel. Tagging a
    channel's cluster in that format means writing the label into the CIF's
    own fields, not renaming a file, and is a different piece of code to the
    one this function is.

    :arg channels: per-frame channel lists, in the same frame-major order
        calcChannelsMultipleFrames() returns as `channels_all` - channels[i] is
        the list of channel objects found in frame i, each carrying the
        `.origin` attribute calcChannels() sets (the search-site index behind
        the "sp<n>" tag in its filenames, or None where the run was not tagged
        by site - see below).
    :type channels: list of lists

    :arg labels: one cluster id per channel, flattened in exactly the same
        frame-major, then within-frame order `channels` itself iterates in -
        i.e. labels[n] must be the cluster id of the n-th channel walking
        `channels` as `for i in range(len(channels)): for j in
        range(len(channels[i])): ...`. Length must equal the total channel
        count across all frames; a length mismatch is the single most common
        way this silently renames the wrong files, so it is checked up front
        rather than discovered partway through a rename pass.

        A label of -1 (HDBSCAN's noise convention) is tagged "noise" rather
        than "-1", both to read clearly in a directory listing and because
        vis_channels_traj.py already special-cases the literal string "noise"
        to group it separately, in grey, switched off by default.
    :type labels: sequence of int

    :arg output_path: the SAME value passed as `output_path` to the original
        calcChannelsMultipleFrames() run, so the files this function looks for
        match the files that run actually wrote. Two forms are accepted:

        - An existing directory (e.g. "vmt_channels/"): files are assumed to
          sit directly inside it, named
          "<stem><frame>[_sp<origin>]_chl<index>.<ext>". `stem` defaults to
          "channels", the label calcChannelsMultipleFrames() hard-codes
          internally (`_frameOutputPath(output_path, j0, "channels")`);
          override `stem=` only if that call used a different label.
        - A directory-and-stem prefix with no separator between them (e.g.
          "frames/frame", giving files like "frames/frame12_chl0.pqr"),
          mirroring calcChannelsMultipleFrames() being invoked with its own
          `output_path` already carrying the intended stem. Detected
          automatically: if `output_path` is not itself an existing
          directory, its final path component becomes the stem and its
          parent becomes the directory, and `stem=` is ignored.

        Renaming happens after the files already exist, so "is it a
        directory" is an unambiguous test at this point even though it would
        not have been before the original run had written anything.
    :type output_path: str or os.PathLike

    :arg stem: filename stem used when `output_path` is a bare directory (see
        above); ignored when `output_path` already carries its own stem.
        Default "channels".
    :type stem: str

    :arg ext: extension of the files actually written - "pqr" (default) or
        "pdb", matching whatever the original run's output_path used. "cif"
        is rejected; see above.
    :type ext: str

    :arg dry_run: if True, print each (old name -> new name) pair without
        calling os.rename, so a wrong `stem`/`output_path` guess is caught
        before anything on disk is touched. Recommended on the first run
        against any given directory. Default False.
    :type dry_run: bool

    :returns: number of files renamed (or, with dry_run=True, that would be).
    :rtype: int

    :raises ValueError: if `ext` is "cif" or otherwise unsupported, or if
        `len(labels)` does not match the total channel count in `channels`.
    """
    from pathlib import Path
    ext = str(ext).lower().lstrip('.')
    if ext == 'cif':
        raise ValueError(
            "renameChannelFilesWithClusterLabels() does not support mmCIF "
            "output: every channel of a frame lives in one shared .cif file "
            "there, so there is no single per-channel file to rename. Tag "
            "clusters into the CIF's own fields instead, or re-run "
            "calcChannelsMultipleFrames() with output_format='pqr' (default) "
            "or a .pdb output_path, and separate=True.")
    if ext not in ('pqr', 'pdb'):
        raise ValueError("ext must be 'pqr' or 'pdb', got {0!r}.".format(ext))

    total_channels = sum(len(frame_channels) for frame_channels in channels)
    if len(labels) != total_channels:
        raise ValueError(
            "labels has {0} entries but channels holds {1} channels in "
            "total - they must line up 1:1, flattened frame by frame, "
            "channel index by channel index, in the order channels itself "
            "iterates.".format(len(labels), total_channels))

    output_path = Path(output_path)
    if output_path.is_dir():
        directory = output_path
        # `stem` keeps its default/passed-in value
    else:
        directory, stem = output_path.parent, output_path.name

    index = 0
    renamed = 0
    for frame_index, frame_channels in enumerate(channels):
        for channel_index, channel in enumerate(frame_channels):
            origin = getattr(channel, 'origin', None)
            # Matches _numberedPath's own rule: a single-site run tags
            # nothing with "sp<n>" at all, rather than writing "spNone".
            site_tag = '' if origin is None else '_sp{0}'.format(origin)
            name = '{0}{1}{2}_chl{3}'.format(
                stem, frame_index, site_tag, channel_index)

            label = int(labels[index])
            tag = 'noise' if label == -1 else str(label)
            index += 1

            old_path = directory / '{0}.{1}'.format(name, ext)
            new_path = directory / '{0}_clid{1}.{2}'.format(name, tag, ext)

            if not old_path.exists():
                # A filtered-out channel, a different run's naming, or a
                # mismatched output_path/stem guess -- reported and skipped
                # rather than raised, so one bad assumption doesn't abort a
                # rename pass that is otherwise correct for every other file.
                _warn("Expected channel file not found, skipped: {0}"
                     .format(old_path))
                continue

            if dry_run:
                print("{0} -> {1}".format(old_path.name, new_path.name))
            else:
                os.rename(str(old_path), str(new_path))
            renamed += 1

    LOGGER.info("{0}{1} channel file(s) with cluster labels{2}.".format(
        'Would rename ' if dry_run else 'Renamed ', renamed,
        ' (dry run)' if dry_run else ''))

    return renamed

def computeBootstrapStability(D, method, records, sub_size=0.80, n_bootstrap=100,
                                reference="full",random_state=123,**kwargs):
    """
    Evaluate clustering stability using bootstrap/subsampled datasets.
    Parameters
    ----------
    D : ndarray, shape (n_samples, n_samples)
        Precomputed distance matrix.
    method : str
        Clustering method passed to clusterChannels().
    sub_size : float, default=0.80
        Fraction of samples to retain in each replicate.
    n_bootstrap : int, default=5
        Number of subsampling replicates.
    reference : {"full", "first"}
        "full": compare each subsample to clustering of all samples.
        "first": compare each subsample to the first subsample.
    random_state : int, default=123
        Random seed.
    **kwargs
        Additional arguments passed to clusterChannels().
    Returns
    -------
    stability : dict
        Contains ARI scores, mean, standard deviation, and
        the full-data reference labels if applicable.
    """
    if not checkAndImport('sklearn'):
        raise ImportError('To run computeBootstrapStability, please install scikit-learn')
    
    from sklearn.metrics import adjusted_rand_score

    D = np.asarray(D)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square distance matrix.")

    n = len(D)
    if not 0 < sub_size <= 1:
        raise ValueError("sub_size must be in (0, 1].")

    if reference not in {"full", "first"}:
        raise ValueError("reference must be 'full' or 'first'.")

    rng = np.random.default_rng(random_state)
    n_keep = max(2, int(round(sub_size * n)))
    # Reference clustering
    if reference == "full":
        labels_ref = np.asarray(
            clusterChannels(D, method, **kwargs))
        if len(labels_ref) != n:
            raise ValueError("clusterChannels returned invalid labels.")
    else:
        labels_ref = None

    scores = []
    bootstrap_labels = []
    for i in range(n_bootstrap):
        # Random subsample without replacement
        keep = rng.choice(n, size=n_keep, replace=False)
        keep.sort()
        D_sub = D[np.ix_(keep, keep)]
        recs_sub = np.asarray(records)[keep]
        # Recluster
        update_record_ids = kwargs.pop('update_record_ids',False)
        labels_sub = np.asarray(clusterChannels(D_sub, method, recs_sub,
                                    update_record_ids=update_record_ids,**kwargs)[0])
        labels_full = np.full(n, -2, dtype=int)
        labels_full[keep] = labels_sub
        bootstrap_labels.append(labels_full)

        if len(labels_sub) != len(keep):
            raise ValueError("clusterChannels returned invalid labels.")

        # First replicate becomes reference
        if reference == "first" and labels_ref is None:
            labels_ref = labels_sub.copy()
            keep_ref = keep.copy()
            continue

        if reference == "full":
            labels_compare = labels_ref[keep]

        else:
            # Compare only channels present in both subsets
            common, idx_sub, idx_ref = np.intersect1d(
                keep, keep_ref, return_indices=True)

            if len(common) < 2:
                continue

            labels_compare = labels_ref[idx_ref]
            labels_sub = labels_sub[idx_sub]

        # ARI ignores the numerical names of clusters
        ari = adjusted_rand_score(labels_compare,labels_sub,)
        scores.append(ari)
    scores = np.asarray(scores)
    C = coassociationMatrix(bootstrap_labels, exclude_noise=True)
    return {"scores": scores,
            "mean": np.mean(scores),
            "std": np.std(scores),
            "median": np.median(scores),
            "reference_labels": labels_ref,
            "coassociation_matrix": C}


def compareClusterings(labels1,labels2):
    """
    Compare two clustering label assignments.

    Parameters
    ----------
    labels1, labels2 : array-like
        Cluster labels for the same set of samples. The two label arrays must
        have the same length and represent equivalent assignments to compare.

    Returns
    -------
    dict
        Dictionary with the adjusted Rand index (``ARI``) and adjusted mutual
        information (``AMI``) between the two clusterings.

    Notes
    -----
    This function requires scikit-learn to be installed.
    """
    
    if not checkAndImport('sklearn'):
        raise ImportError('To run computeBootstrapStability, please install scikit-learn')
    if not checkAndImport('sklearn'):
        errorMsg = 'To run getClusterStats, please install scikit-learn'
        raise ImportError(errorMsg)
    from sklearn.metrics import adjusted_rand_score,adjusted_mutual_info_score
    
    return {"ARI": adjusted_rand_score(labels1, labels2),
            "AMI": adjusted_mutual_info_score(labels1, labels2)}


def coassociationMatrix(label_list, exclude_noise=True):
    """
    Calculate a co-association matrix from multiple clustering results.

    Parameters
    ----------
    label_list : list of array-like
        Each element contains cluster labels for the same samples.
        Noise should be labeled as -1.

    exclude_noise : bool, default=True
        If True, pairs are normalized only over runs where both samples
        are assigned to a non-noise cluster.

        If False, all runs are included in the denominator, so a run in
        which either sample is noise contributes zero.

    Returns
    -------
    C : ndarray
        Symmetric co-association matrix.

        With exclude_noise=True:
            C[i, j] = fraction of informative runs in which i and j
                      were assigned to the same cluster.

        With exclude_noise=False:
            C[i, j] = fraction of all runs in which i and j were assigned
                      to the same cluster.
    """
    labels0 = np.asarray(label_list[0])
    n = len(labels0)
    C = np.zeros((n, n), dtype=float)
    counts = np.zeros((n, n), dtype=float)

    for labels in label_list:
        labels = np.asarray(labels)
        if len(labels) != n:
            raise ValueError("All label arrays must have the same length.")
        
        present = labels != -2
        if exclude_noise:
            valid = present & (labels >= 0)
        else:
            valid = present

        valid_pair = valid[:, None] & valid[None, :]
        same = ((labels[:, None] == labels[None, :]) & (labels[:, None] >= 0))

        C += same & valid_pair
        counts += valid_pair

    C = np.divide(C,counts,out=np.zeros_like(C), where=counts > 0)
    return C

def clusterConsensusStrength(C, cluster_indices):
    """
    Calculate consensus strength for a cluster.
    Parameters
    ----------
    C : ndarray, shape (n_samples, n_samples)
        Co-association matrix.
    cluster_indices : array-like
        Indices of samples belonging to the cluster.
    Returns
    -------
    float
        Mean pairwise co-association among cluster members.
        1.0 = every pair was always co-clustered
        0.0 = no pair was ever co-clustered
    """
    cluster_indices = np.asarray(cluster_indices)
    if len(cluster_indices) < 2:
        return np.nan
    sub = C[np.ix_(cluster_indices, cluster_indices)]
    i, j = np.triu_indices(len(cluster_indices), k=1)
    return np.mean(sub[i, j])


def writeClusterLabels(filename, labels, names=None):
    """
    Write one or more clustering label sets to a tab-separated text file.

    Parameters
    ----------
    filename : str or os.PathLike
        Output filename. A '.tsv' extension is recommended.

    labels : array-like or sequence of array-like
        Cluster labels to write.

        For a single clustering, this may be a 1-D array of shape
        ``(n_channels,)``.

        For multiple clusterings, this may be either:

        - a sequence of 1-D arrays, one array per clustering, or
        - a 2-D NumPy array of shape ``(n_channels, n_clusterings)``.

        All clustering label sets must contain the same number of channels.

    names : sequence of str, optional
        Column names for the clustering label sets. The number of names must
        equal the number of clustering columns. If omitted, names are generated
        as ``labels_1``, ``labels_2``, etc.

    Returns
    -------
    str
        The output filename.

    Notes
    -----
    Label -1 is conventionally used for noise/unclustered channels.

    Rows must be in the same channel order expected by the PyMOL
    visualization script: ascending frame followed by ascending channel rank
    within each frame.
    """

    if isinstance(labels, np.ndarray) and labels.ndim == 1:
        arr = labels[:, None]

    elif isinstance(labels, np.ndarray) and labels.ndim == 2:
        arr = labels

    else:
        label_sets = [np.asarray(x) for x in labels]

        if not label_sets:
            raise ValueError("At least one label set must be supplied.")

        if any(x.ndim != 1 for x in label_sets):
            raise ValueError("Each label set must be one-dimensional.")

        lengths = {len(x) for x in label_sets}
        if len(lengths) != 1:
            raise ValueError("All label sets must contain the same number of channels.")

        arr = np.column_stack(label_sets)

    if arr.ndim != 2:
        raise ValueError("Labels must form a 2-D label table.")

    n_channels, n_clusterings = arr.shape

    if names is None:
        names = ["labels_%d" % (i + 1) for i in range(n_clusterings)]
    else:
        names = list(names)

    if len(names) != n_clusterings:
        raise ValueError("Received %d column names for %d clustering label sets." % (len(names), n_clusterings))

    try:
        arr = arr.astype(int)
    except (TypeError, ValueError):
        raise ValueError("All cluster labels must be integers.")

    with open(filename, "w") as out:
        out.write("# " + "\t".join(map(str, names)) + "\n")
        for row in arr:
            out.write("\t".join(map(str, row)) + "\n")
    return str(filename)


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def _cached(obj, key, build):
    cache = obj.__dict__.setdefault("_dist_cache", {})
    if key not in cache:
        cache[key] = build()
    return cache[key]


def _arc(obj, n):
    """(t, p, r, s, kd-tree) arc-length samples, cached."""
    ch = getattr(obj, "channel", obj)
    t, p, r, s = _cached(obj, ("arc", n), lambda: _sampleSplineByArcLength(ch, n))
    p = np.asarray(p)
    tree = _cached(obj, ("arc_tree", n), lambda: _kdTree(p)) if len(p) else None
    return t, p, r, s, tree


def _dense(obj, n_dense):
    """(spline, t, points, unit tangents, kd-tree) dense samples of the spline, cached."""
    spline = getattr(obj, "channel", obj).getSplines()[0]

    def build():
        td = np.linspace(spline.x[0], spline.x[-1], n_dense)
        pd = np.asarray(spline(td))
        return td, pd, _unit(spline(td, 1)), _kdTree(pd)

    return (spline,) + _cached(obj, ("dense", n_dense), build)
  
def _frechetTable(D, full=False):
    """
    Discrete Frechet distance from a point-distance matrix D (na x nb); with
    full=True returns the whole (na, nb) coupling table (needed for the path).
 
    Same recurrence as the double loop, swept one anti-diagonal at a time on a
    skewed copy of D (S[k, i] = D[i, k-i]) so every step is a few contiguous-slice
    numpy ops with no fancy indexing: na+nb steps instead of na*nb iterations.
    """
    na, nb = D.shape
    I, K = np.arange(na)[None, :], np.arange(na + nb - 1)[:, None]
    J = K - I
    S = np.where((J >= 0) & (J < nb), D[I, np.clip(J, 0, nb - 1)], np.inf)
    p2 = np.full(na + 1, np.inf)            # diagonal k-2; index i+1 holds cell i, index 0 stays inf
    p1 = np.full(na + 1, np.inf)
    p1[1] = S[0, 0]                         # cell (0, 0)
    rows = [p1[1:].copy()] if full else None
    tmp = np.empty(na)
    for k in range(1, na + nb - 1):
        cur = np.full(na + 1, np.inf)
        np.minimum(p1[1:], p1[:-1], out=tmp)        # (i, j-1) and (i-1, j)
        np.minimum(tmp, p2[:-1], out=tmp)           # (i-1, j-1)
        np.maximum(S[k], tmp, out=cur[1:])
        if full:
            rows.append(cur[1:].copy())
        p2, p1 = p1, cur
    if not full:
        return float(p1[na])                # cell (na-1, nb-1) is on the last diagonal
    ii, jj = np.indices((na, nb))
    return np.array(rows)[ii + jj, ii]

class ClusterCalculator:
    
    def __init__(self, atoms=None, voxel_size=0.25,pad_factor=1.02,h=10,angle_deg=5.0,
                q=1.0,n_samples=5, global_origin=None,alpha=1.0,beta=2.0):
        
        self.atoms = atoms
        self.n_residues = len(np.unique(atoms.getResindices())) if atoms else None
        self.voxel_size=voxel_size
        self.pad_factor=pad_factor
        self.h=h
        self.angle_deg=angle_deg
        self.q=q
        self.n_samples=n_samples
        self.global_origin = global_origin
        self.alpha = alpha
        self.beta = beta
        POPCOUNT_TABLE = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)
        self.POPCOUNT_TABLE = POPCOUNT_TABLE

        self.scan = None

    class ChannelRecord:
        """Cached geometry for one channel"""

        def __init__(self,calculator, channel, conformation_id, lining_residues = None, rep_points=None,
                    rep_radii=None, centerline_points=None,centerline_radii=None,
                    tree=None,bbox=None,mask=None,offset=None,dims=None,volume=None,
                    endpoint=None,r_T=None,kept=True,cluster_ids=None,shell_weights=None):

            if not channel.getSplines()[0]:
                channel._buildSplines()
            self.channel = channel
            self.conformation_id = conformation_id
            self.calculator = calculator

            #calculated later
            self.centerline_points = centerline_points
            self.centerline_radii = centerline_radii
            self.dims=dims
            self.bbox = bbox
            self.tree = tree
            self.mask = mask
            self.offset = offset
            self.volume = volume
            self.lining_residues = lining_residues
            self.packed_lining = self.calculator._pack(lining_residues) if lining_residues is not None else None

            #calculated during center sampling 
            self.rep_points = rep_points
            self.rep_radii = rep_radii
            self.endpoint = endpoint
            self.r_T = r_T
            self.kept = kept
            self.cluster_ids = [] if cluster_ids is None else cluster_ids
            self.shell_weights = shell_weights
            self.voxelized = False
            self.__postInit__()


        def __postInit__(self):
            """Prepare geometry that does not depend on the Metric Calcs"""
            self.centerline_points, self.centerline_radii = _sampleObjectSpheres(self.channel,
                                                 self.calculator.n_samples)
            self.tree = _kdTree(self.centerline_points)
            self.bbox = self.calculator.channelBoundBox(self.channel)                       

        def _prepareVoxels(self):
            """Prepare Voxel Representation Geometry"""
            global_origin = self.calculator.global_origin
            voxel_size= self.calculator.voxel_size
            mn, mx = self.bbox
            i0 = np.floor((mn - global_origin) / voxel_size).astype(int)
            i1 = np.ceil((mx - global_origin) / voxel_size).astype(int)
            self.offset = i0
            self.dims = np.maximum(i1 - i0, 1)
            grid_mn = global_origin + i0 * voxel_size
            self.mask = self.calculator.voxelizeChannelAligned(self.centerline_points,
                                                            self.centerline_radii,
                                                            self.tree,
                                                            self.dims,
                                                            grid_mn)
        
            self.volume = (np.count_nonzero(self.mask) * self.calculator.voxel_size ** 3)
            self.voxelized = True
                
        @property
        def bbox_min(self):
            return self.bbox[0]

        @property
        def bbox_max(self):
            return self.bbox[1]

        @property
        def n_shells(self):
            if self.rep_points is None:
                return 0
            return len(self.rep_points)

        @property
        def n_voxels(self):
            if self.mask is None:
                return 0
            return int(np.count_nonzero(self.mask))

    def _buildChannelRecord(self,channel, conformation_id, voxelize=False):
        """Top-level (picklable) worker for parallel channel voxelization."""
        if voxelize:
            record = self.ChannelRecord(self, channel, conformation_id)
            record._prepareVoxels()
            return record
        return self.ChannelRecord(self, channel, conformation_id)

    def _prepareVoxelsWorker(self,record):
        record._prepareVoxels()
        return record

    def _buildChannelRecordWorker(self,args):
        channel, conformation_id, voxelize = args
        return self._buildChannelRecord(channel, conformation_id, voxelize=voxelize)
    
    class DendrogramStability:
        """
        Precomputes everything needed to score any dendrogram node in O(1).

        Parameters
        ----------
        Z : (n-1, 4) linkage matrix
        leaf_weight : optional (n,) array; leaves with weight <= 0 are ignored
            (noise points). Default: all leaves count.
        condense_size : int
            A sub-branch must hold at least this many counted points to count as
            a genuine split; smaller branches are treated as points falling out
            (HDBSCAN's condensing step).

        Attributes (all indexed by node id, leaves are 0..n-1)
        ----------
        raw    : (n-1,) merge heights made monotone (index = node - n)
        dist   : raw floored to a small positive value (safe for logs / 1/d)
        ncl    : counted leaves below each node
        parent : parent node id or -1
        """

        def __init__(self, Z, leaf_weight=None, condense_size=3):
            Z = np.asarray(Z, dtype=float)
            n = Z.shape[0] + 1
            self.n = n
            self.left = Z[:, 0].astype(int)
            self.right = Z[:, 1].astype(int)
            raw = Z[:, 2].copy()
            # Enforce monotone heights (centroid/median linkage can invert).
            for i in range(n - 1):
                for c in (self.left[i], self.right[i]):
                    if c >= n:
                        raw[i] = max(raw[i], raw[c - n])
            self.raw = raw
            pos = raw[raw > 0]
            self.floor = 0.5 * pos.min() if pos.size else 1e-12
            self.dist = np.maximum(raw, self.floor)
            self.condense_size = max(1, int(condense_size))

            total = 2 * n - 1
            self.parent = np.full(total, -1, dtype=int)
            self.ncl = np.zeros(total, dtype=int)
            if leaf_weight is None:
                self.ncl[:n] = 1
            else:
                self.ncl[:n] = (np.asarray(leaf_weight) > 0).astype(int)

            W = np.zeros(total)   # points counted
            A = np.zeros(total)   # sum w * log(d_shed)
            B = np.zeros(total)   # sum w / d_shed
            cond = self.condense_size
            for i in range(n - 1):
                v, a, b, d = n + i, self.left[i], self.right[i], self.dist[i]
                self.parent[a] = self.parent[b] = v
                self.ncl[v] = self.ncl[a] + self.ncl[b]
                big_a, big_b = self.ncl[a] >= cond, self.ncl[b] >= cond
                if big_a and big_b:               # genuine split: cluster ends
                    w = self.ncl[v]
                    W[v], A[v], B[v] = w, w * np.log(d), w / d
                    continue
                for ch, big in ((a, big_a), (b, big_b)):
                    if big:                        # keep following the big branch
                        W[v] += W[ch]; A[v] += A[ch]; B[v] += B[ch]
                    elif self.ncl[ch] > 0:         # small branch falls out here
                        w = self.ncl[ch]
                        W[v] += w; A[v] += w * np.log(d); B[v] += w / d
            self.W, self.A, self.B = W, A, B

        def height(self, node):
            return self.raw[node - self.n] if node >= self.n else 0.0

        def stability(self, node, ext):
            """(log-stability, 1/d-stability) of `node` born at height `ext`."""
            e = max(float(ext), self.floor)
            return (self.W[node] * np.log(e) - self.A[node],
                    self.B[node] - self.W[node] / e)
        
    @staticmethod
    def stabilityNoiseDecision(size, stability, ratio, min_cluster_size,
                                min_persistence_ratio, singleton_noise=True,
                                min_ratio_floor=1.0):
        """
        Shared noise rule. A cluster is kept if

            stability >= min_cluster_size * log(min_persistence_ratio)

        i.e. it persists at least as much as `min_cluster_size` points each
        surviving at the minimum ratio. Bigger clusters therefore need a smaller
        per-point ratio, smaller ones a larger one.
        """
        if singleton_noise and size == 1:
            return True
        if stability is None or np.isnan(stability):      # nothing to separate from
            return size < min_cluster_size
        threshold = min_cluster_size * np.log(min_persistence_ratio)
        return bool(stability < threshold) or (
            np.isfinite(ratio) and ratio < min_ratio_floor)

##### Bounding Box Functions #####
    def channelBoundBox(self,channel):
        """
        Axis-aligned bbox padded by the channel's maximum radius.
        """
        r = channel.radii.max() * self.pad_factor
        mn = channel.centers.min(axis=0) - r
        mx = channel.centers.max(axis=0) + r
        return mn, mx

    def _pairwiseBBoxOverlap(self,mn_a, mx_a, mn_b, mx_b):
        """
        Vectorized bbox-overlap test between two GROUPS of bboxes.

        mn_a, mx_a : (m, 3) arrays - bbox corners for group A
        mn_b, mx_b : (k, 3) arrays - bbox corners for group B

        Returns an (m, k) boolean matrix, computed with a single pair of
        broadcasted min/max calls instead of a python-level double loop.
        Kept as its own small block (m, k, 3) rather than doing this across
        *all* n channels at once, so memory stays O(frame_size^2) instead
        of O(n^2).
        """
        lo = np.maximum(mn_a[:, None, :], mn_b[None, :, :])   # (m, k, 3)
        hi = np.minimum(mx_a[:, None, :], mx_b[None, :, :])   # (m, k, 3)

        return np.all(hi > lo, axis=-1)                       # (m, k)
##### Bounding Box Functions #####


##### Functions Related to the Channel Lining Residue Pipeline #####
    def _pack(self,residue_set):
        arr = np.zeros(self.n_residues, dtype=np.uint8)
        arr[list(residue_set)] = 1
        return np.packbits(arr)
    
    def _TverskyCounts(self, packed):
        n = packed.shape[0]
        inter = np.zeros((n, n), dtype=np.int64)
        lo = np.zeros((n, n), dtype=np.int64)
        hi = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            inter[i] = self.POPCOUNT_TABLE[packed[i] & packed].sum(axis=1)
            x_minus_y = self.POPCOUNT_TABLE[packed[i] & ~packed].sum(axis=1)
            y_minus_x = self.POPCOUNT_TABLE[~packed[i] & packed].sum(axis=1)
            lo[i] = np.minimum(x_minus_y, y_minus_x)
            hi[i] = np.maximum(x_minus_y, y_minus_x)
        return inter, lo, hi

    @staticmethod
    def _TverskyFromCounts(inter, lo, hi, alpha, beta):
        denominator = inter + beta * (alpha * lo + (1 - alpha) * hi)
        with np.errstate(divide='ignore', invalid='ignore'):
            dist = 1.0 - inter / denominator
        dist[denominator == 0] = 0.0      # both sets empty -> identical
        return dist

    def _TverskyDistanceMatrix(self, packed):
        inter, lo, hi = self._TverskyCounts(packed)
        return self._TverskyFromCounts(inter, lo, hi, self.alpha, self.beta)
##### Functions Related to the Channel Lining Residue Pipeline #####


##### Functions Related to the Voxelization Pipeline #####
    def _segmentProject(self,p, a, b):
        """
        Project points p onto segments a -> b.

        Returns
        -------
        dist2 : ndarray
            Squared distance from each point to the segment.
        t : ndarray
            Segment parameter in [0, 1].
        """
        ab = b - a

        denom = np.einsum('ij,ij->i', ab, ab)
        denom = np.maximum(denom, 1e-12)

        t = np.einsum('ij,ij->i', p - a, ab) / denom
        t = np.clip(t, 0.0, 1.0)
        closest = a + t[:, None] * ab
        diff = p - closest
        dist2 = np.einsum('ij,ij->i', diff, diff)
        return dist2, t


    def _pointsInsideMask(self,query_points, points, radii, tree, refine=True):
        """
        Vectorized point-in-channel test, given an already-built tree.

        Single tree.query() call over all query_points at once (uses
        scipy's multi-threaded query), followed by the same vectorized
        segment-refinement used previously - no per-point or per-segment
        Python loop.
        """
        d_nn, idx = tree.query(query_points, k=1, workers=-1)

        if not refine or len(points) < 3:
            return d_nn <= radii[idx]

        n = len(points)

        idx_prev = np.clip(idx - 1, 0, n - 1)
        idx_next = np.clip(idx + 1, 0, n - 1)

        dist2_a, t_a = self._segmentProject(query_points,
                                        points[idx_prev],
                                        points[idx])

        dist2_b, t_b =  self._segmentProject(query_points,
                                        points[idx],
                                        points[idx_next])

        use_a = dist2_a < dist2_b
        r_a = (radii[idx_prev] + t_a * (radii[idx] - radii[idx_prev]))
        r_b = (radii[idx] + t_b * (radii[idx_next] - radii[idx]))
        radius = np.where(use_a, r_a, r_b)
        dist2 = np.where(use_a, dist2_a, dist2_b)

        return dist2 <= radius * radius


    def _insideChannelMask(self,query_points, prepared, refine=True):
        """
        Determine which query points are inside a channel, via KD-tree.

        Kept for callers that need point-in-channel tests for arbitrary,
        unstructured query points (e.g. lining queries). Not used by the
        volume/intersection pipeline below.
        """
        query_points = np.asarray(query_points, dtype=float)

        return self._pointsInsideMask(query_points,
                                    prepared.centerline_points,
                                    prepared.centerline_radii,
                                    prepared.tree,
                                    refine=refine)


    def _gridPoints(self,dims, grid_mn):
        """
        World-space coordinates of every voxel center in a `dims`-shaped
        grid whose first voxel center sits at `grid_mn`, built as a single
        vectorized meshgrid (no per-voxel or per-segment Python loop).
        """
        voxel_size=self.voxel_size

        xs = grid_mn[0] + (np.arange(dims[0]) + 0.5) * voxel_size
        ys = grid_mn[1] + (np.arange(dims[1]) + 0.5) * voxel_size
        zs = grid_mn[2] + (np.arange(dims[2]) + 0.5) * voxel_size

        X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

        return np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)

    def voxelizeChannelAligned(self,centers, radii, tree, dims, grid_mn):
        """
        Rasterize a channel onto a lattice-aligned boolean voxel grid.

        Parameters
        ----------
        centers : array-like
            Centerline coordinates for the channel, shaped ``(n, 3)``.
        radii : array-like
            Radius for each centerline point, shaped ``(n,)``.
        tree : scipy.spatial.cKDTree or similar
            KD-tree built from ``centers`` for efficient point-in-channel tests.
        dims : tuple[int, int, int]
            Grid dimensions in voxels along each axis.
        grid_mn : array-like
            World-space coordinates of the lower corner of the grid origin.

        Returns
        -------
        numpy.ndarray
            Boolean mask with shape ``dims`` marking voxels whose centers lie
            inside the channel volume.

        Notes
        -----
        This method first finds a compact set of candidate voxels from the
        union of padded segment bounding boxes, then evaluates only those
        candidates against the channel geometry using the point-inside-mask
        predicate. This avoids a full voxel-by-voxel scan while preserving the
        exact channel volume on the lattice.
        """
        
        voxel_size=self.voxel_size
        if len(centers) < 2:
            return np.zeros(dims, dtype=bool)

        # Vectorized per-segment padded index boxes
        a, b = centers[:-1], centers[1:]
        r = np.maximum(radii[:-1], radii[1:])

        seg_mn = np.minimum(a, b) - r[:, None]
        seg_mx = np.maximum(a, b) + r[:, None]

        i0 = np.clip(np.floor((seg_mn - grid_mn) / voxel_size).astype(int), 0, dims)
        i1 = np.clip(np.ceil((seg_mx - grid_mn) / voxel_size).astype(int), 0, dims)

        # Union of segment windows -> candidate voxels only
        active = np.zeros(dims, dtype=bool)
        for k in range(len(a)):
            if np.any(i1[k] <= i0[k]):
                continue
            active[i0[k, 0]:i1[k, 0], i0[k, 1]:i1[k, 1], i0[k, 2]:i1[k, 2]] = True

        idx = np.argwhere(active)
        if len(idx) == 0:
            return np.zeros(dims, dtype=bool)

        pts = grid_mn + (idx + 0.5) * voxel_size

        # One batched query, but only over the reduced candidate set.
        inside = self._pointsInsideMask(pts, centers, radii, tree, refine=True)

        grid = np.zeros(dims, dtype=bool)
        grid[idx[:, 0], idx[:, 1], idx[:, 2]] = inside
        return grid   

    def _maskVoxelOverlapVolume(self,prepared_a, prepared_b):
        """
        Overlap volume of two already-voxelized, lattice-aligned channels.

        Pure integer index arithmetic to find the overlapping voxel-index
        range, then a slice + bitwise AND + count_nonzero. No rasterization,
        no floating point, no search of any kind happens here.
        """
        
        start = np.maximum(prepared_a.offset, prepared_b.offset)
        end = np.minimum(prepared_a.offset + prepared_a.dims,
                            prepared_b.offset + prepared_b.dims)

        if np.any(end <= start):
            return 0.0

        sa0 = start - prepared_a.offset
        sa1 = end - prepared_a.offset
        sb0 = start - prepared_b.offset
        sb1 = end - prepared_b.offset

        sub_a = prepared_a.mask[sa0[0]:sa1[0], sa0[1]:sa1[1], sa0[2]:sa1[2]]
        sub_b = prepared_b.mask[sb0[0]:sb1[0], sb0[1]:sb1[1], sb0[2]:sb1[2]]

        intersection = np.count_nonzero(sub_a & sub_b)
        a_minus_b =  np.count_nonzero(sub_a & ~sub_b)
        b_minus_a =  np.count_nonzero(~sub_a & sub_b)
        a = np.minimum(a_minus_b, b_minus_a)
        b = np.maximum(a_minus_b, b_minus_a)
        u3 = self.voxel_size ** 3
        return intersection*u3, a*u3,b*u3

    def _voxelOverlapCounts(self, task):
        i, j, rec_a, rec_b = task
        inter, a, b = self._maskVoxelOverlapVolume(rec_a, rec_b)
        return i, j, inter, a, b

    @staticmethod
    def _voxelScoreFromCounts(inter, a, b, alpha, beta):
        """Works on scalars or arrays. No overlap volume at all -> 1.0."""
        inter = np.asarray(inter, dtype=float)
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        denominator = inter + beta * (alpha * a + (1 - alpha) * b)
        with np.errstate(divide='ignore', invalid='ignore'):
            score = 1.0 - inter / denominator
        return np.where(denominator > 0, score, 1.0)

    def _scoreVoxelOverlapPrepared(self, task):
        i, j, inter, a, b = self._voxelOverlapCounts(task)
        return i, j, float(self._voxelScoreFromCounts(inter, a, b, self.alpha, self.beta))

##### Functions Related to the Voxelization Pipeline #####


##### Functions Related to the Center Sampling Pipeline #####


    def computeAverageStart(self,records):
        """S_average = centroid of starting vertices from all conformations."""
        starts = np.array([np.asarray(r.channel.centers[0], dtype=float) for r in records])
        return starts.mean(axis=0)


    def computeRepresentativePoints(self,records,S,h=10):
        """
        Compute centerline representative points for all ChannelRecords.
        For each channel and each radial shell, multiple continuous centerline
        segments are allowed. Results are stored directly on each ChannelRecord:

            record.rep_points[i]  -> (n_segments, 3)
            record.rep_radii[i]   -> (n_segments,)

        The representative point of each segment is an arc-length-weighted
        centroid, and the representative radius is an arc-length-weighted
        mean radius.

        Parameters
        ----------
        records : list[ChannelRecord]
            Tunnel records.

        S : array-like
            Common start point shared by every channel/conformation. The radial
            shell coordinate is calculated relative to THIS point (not each
            channel's own mouth) -- that's what makes "shell i" comparable
            across different channels in Eq. 2.

        h : int
            Number of radial shells.

        angle_threshold_deg, agg_amax_deg : float
            Retained for compatibility with the original CAVER API.

        n_samples : int
            Number of points sampled along each centerline spline.

        Returns
        -------
        records
            The same records with representative data populated.
        """
        S = np.asarray(S, dtype=np.float64)

        for record in records:
            channel = record.channel
            centerline_spline = channel.centerline_spline
            radius_spline = channel.radius_spline
            n_centers = len(channel.centers)

            if n_centers < 2:
                record.rep_points = [np.asarray(channel.centers, dtype=np.float32) for _ in range(h)]
                record.rep_radii = [np.asarray(channel.radii, dtype=np.float32) for _ in range(h)]
                record.shell_weights = np.ones(h, dtype=np.float32)
                continue

            t = np.linspace(0.0, float(n_centers - 1), self.n_samples, dtype=np.float64)

            points = np.asarray(centerline_spline(t), dtype=np.float64)
            if points.ndim == 2 and points.shape == (3, self.n_samples):
                points = points.T
            points = np.asarray(points, dtype=np.float64)

            radii = np.asarray(radius_spline(t), dtype=np.float64).reshape(-1)
            radii = np.maximum(radii, 0.0)

            record.centerline_points = points.astype(np.float32, copy=False)
            record.centerline_radii = radii.astype(np.float32, copy=False)

            dxyz = np.diff(points, axis=0)
            ds = np.linalg.norm(dxyz, axis=1)

            # Radial coordinate relative to the SHARED start point S -- shell i
            # must mean the same physical distance for every channel, or the
            # per-shell comparison in Eq. 2 is comparing unrelated things.
            radial_distance = np.linalg.norm(points - S[None, :], axis=1)

            rmin = radial_distance.min()
            rmax = radial_distance.max()

            if rmax <= rmin:
                shell_edges = np.full(h + 1, rmin, dtype=np.float64)
            else:
                shell_edges = np.linspace(rmin, rmax, h + 1, dtype=np.float64)

            rep_points = []
            rep_radii = []

            for shell_i in range(h):
                lo = shell_edges[shell_i]
                hi = shell_edges[shell_i + 1]

                if shell_i == h - 1:
                    in_shell = (radial_distance >= lo) & (radial_distance <= hi)
                else:
                    in_shell = (radial_distance >= lo) & (radial_distance < hi)

                if not np.any(in_shell):
                    rep_points.append(np.empty((0, 3), dtype=np.float32))
                    rep_radii.append(np.empty((0,), dtype=np.float32))
                    continue

                change = np.flatnonzero(in_shell[1:] != in_shell[:-1]) + 1
                starts = np.concatenate((np.array([0]), change))
                stops = np.concatenate((change, np.array([len(in_shell)])))

                shell_rep_points = []
                shell_rep_radii = []

                for start, stop in zip(starts, stops):
                    if not in_shell[start] or stop <= start:
                        continue

                    p = points[start:stop]
                    r = radii[start:stop]

                    if len(p) == 1:
                        shell_rep_points.append(p[0])
                        shell_rep_radii.append(r[0])
                        continue

                    segment_ds = ds[start:stop - 1]
                    weights = np.empty(len(p), dtype=np.float64)
                    weights[0] = segment_ds[0] * 0.5
                    weights[-1] = segment_ds[-1] * 0.5
                    if len(p) > 2:
                        weights[1:-1] = (segment_ds[:-1] + segment_ds[1:]) * 0.5

                    total_weight = weights.sum()
                    if total_weight <= 0.0:
                        centroid = p.mean(axis=0)
                        mean_radius = r.mean()
                    else:
                        centroid = np.sum(p * weights[:, None], axis=0) / total_weight
                        mean_radius = np.sum(r * weights) / total_weight

                    shell_rep_points.append(centroid)
                    shell_rep_radii.append(mean_radius)

                if shell_rep_points:
                    rep_points.append(np.asarray(shell_rep_points, dtype=np.float32))
                    rep_radii.append(np.asarray(shell_rep_radii, dtype=np.float32))
                else:
                    rep_points.append(np.empty((0, 3), dtype=np.float32))
                    rep_radii.append(np.empty((0,), dtype=np.float32))

            record.rep_points = rep_points
            record.rep_radii = rep_radii
            record.shell_weights = np.ones(h, dtype=np.float32)

        return records

    def computeDistanceWeights(self, h, q= 1.0):
        """w(x) = k1*x + k2, chosen so w(0.5) = 1 and w(1)/w(0) = q."""
        if h == 1:
            return np.array([1.0])
        denom = (q - 1) * 0.5 + 1.0
        k2 = 1.0 / denom if abs(denom) > 1e-12 else 1.0
        k1 = k2 * (q - 1)
        x = np.arange(h) / (h - 1)
        return k1 * x + k2


    # module-level worker state for optional multiprocessing 
    _worker_shell_data = None
    _worker_weights = None
    _worker_n = None
    _worker_normalize_radii = True
    _worker_clip_distance = True
    _worker_centroids=None
    _worker_radii_means = None

    def _initDistanceWorker(self, shell_data, weights, n, centroids, radii_means,
                            normalize_radii=True,clip_distance=True):
        """Populate worker-local state used by _distanceBlock(). Only the
        small per-shell flattened representative arrays are sent to each
        worker (via ProcessPoolExecutor's initializer/initargs), not the full
        Channel objects."""

        global _worker_shell_data,_worker_weights, _worker_n
        global _worker_normalize_radii, _worker_clip_distance
        global _worker_centroids, _worker_radii_means
       
        _worker_shell_data = shell_data
        _worker_weights = np.asarray(weights, dtype=np.float32)
        _worker_n = n
        _worker_normalize_radii = normalize_radii
        _worker_clip_distance = clip_distance
        _worker_centroids = np.asarray(centroids, dtype=np.float32)
        _worker_radii_means = np.asarray(radii_means,dtype=np.float32)

    def _prepareAllShells(self, records, h):
        """
        Flatten every shell's representatives across all channels, ONCE, and
        also precompute the TARGET-side grouping (group_starts/group_ids)
        ONCE per shell -- these don't depend on which row-block is being
        processed, so recomputing them per block (as a naive port would) is
        pure waste. Only the SOURCE side is grouped per block, inside
        _distance_block, since that genuinely differs per block.

        Returns a list of length h, each element:
            (points, radii, channel_ids, group_starts, group_ids)

        channel_ids is sorted (built via np.repeat over np.arange(n)), and a
        channel with an empty shell simply doesn't appear in it.
        """
        n = len(records)
        per_shell = []

        for s in range(h):
            pts_list = [r.rep_points[s] for r in records]
            lens = np.fromiter((len(p) for p in pts_list), dtype=np.int64, count=n)

            if lens.sum() == 0:
                per_shell.append((np.empty((0, 3), dtype=np.float32),
                                    np.empty((0,), dtype=np.float32),
                                    np.empty((0,), dtype=np.int32),
                                    np.empty((0,), dtype=np.int64),
                                    np.empty((0,), dtype=np.int32)))
                continue

            points = np.concatenate(pts_list, axis=0)
            radii = np.concatenate([r.rep_radii[s] for r in records], axis=0)
            channel_ids = np.repeat(np.arange(n, dtype=np.int32), lens)

            group_starts = np.r_[0, np.flatnonzero(np.diff(channel_ids)) + 1]
            group_ids = channel_ids[group_starts]

            per_shell.append((points, radii, channel_ids, group_starts, group_ids))

        return per_shell

    def _prepareChannelCentroids(self, records):
        """One 3D centroid + mean radius per channel, from ALL its representative
        points/radii across every shell (falls back to the raw centerline data
        for channels with zero populated shells)."""
        n = len(records)
        centroids = np.empty((n, 3), dtype=np.float32)
        radii_means = np.empty(n, dtype=np.float32)          # <-- (n,), not (n,3)
        for i, r in enumerate(records):
            pts = [p for p in r.rep_points if len(p) > 0]
            rds = [rad for rad in r.rep_radii if len(rad) > 0]
            if pts:
                centroids[i] = np.concatenate(pts, axis=0).mean(axis=0)
                radii_means[i] = np.concatenate(rds, axis=0).mean()
            else:
                centroids[i] = np.asarray(r.channel.centers, dtype=np.float32).mean(axis=0)
                radii_means[i] = np.asarray(r.channel.radii, dtype=np.float32).mean()
        return centroids, radii_means

    def _distanceBlock(self, block):
        """
        Calculate one row-block of the channel-distance matrix.

        SOURCE-side grouping (this block's channels) is computed here since it
        depends on `block`. TARGET-side grouping (group_starts/group_ids) is
        read from the shared, precomputed shell data -- computed once in
        _prepare_all_shells, not once per (shell, block) pair.
        """
        from scipy.spatial.distance import cdist
        start, end = block
        
        shell_data = _worker_shell_data
        weights = _worker_weights
        n = _worker_n
        normalize_radii = _worker_normalize_radii
        clip_distance = _worker_clip_distance
        centroids = _worker_centroids
        radii_means = _worker_radii_means   
        h = len(weights)
        block_n = end - start

        Db = np.zeros((block_n, n), dtype=np.float32)
        Wb = np.zeros((block_n, n), dtype=np.float32)

        for shell in range(h):
            w = weights[shell]
            if w == 0:
                continue

            points, radii, channel_ids, group_starts, group_ids = shell_data[shell]
            if len(points) == 0:
                continue

            # channel_ids is globally sorted, so searchsorted is cheaper than a
            # boolean mask over the full array.
            first = np.searchsorted(channel_ids, start, side="left")
            last = np.searchsorted(channel_ids, end, side="left")
            if first >= last:
                continue

            A = points[first:last]
            rA = radii[first:last]
            idA = channel_ids[first:last] - start   

            B = points
            rB = radii

            spatial = cdist(A, B)
            if normalize_radii:
                spatial = spatial / np.maximum(rA[:, None] + rB[None, :], 1e-8)
                if clip_distance:
                    np.minimum(spatial, 1.0, out=spatial)

            # Reduce SOURCE representatives -> one row per source channel.
            source_change = np.flatnonzero(np.diff(idA)) + 1
            source_starts = np.r_[0, source_change]
            source_ids = idA[source_starts]
            source_reduced = np.minimum.reduceat(spatial, source_starts, axis=0)

            # Reduce TARGET representatives using the CACHED grouping.
            target_reduced = np.minimum.reduceat(source_reduced, group_starts, axis=1)

            Db[np.ix_(source_ids, group_ids)] += w * target_reduced
            Wb[np.ix_(source_ids, group_ids)] += w

        valid_pairs = Wb > 0
        Db[valid_pairs] /= Wb[valid_pairs]
        #no shell populated by both = maximally different
        #radii normalization and centroid sampling may be inconsistent from rest of
        #algorithm, but this is only a fallback case use to estimate the distance 
        #between channels that do not occupy any of the same shells
        #and therefore should have a large distance, approximation of channels
        #that should never be clustered
        invalid = ~valid_pairs
        if np.any(invalid):
            fallback = cdist(centroids[start:end], centroids)
            if normalize_radii:
                rA_c = radii_means[start:end]
                rB_c = radii_means
                fallback = fallback / np.maximum(rA_c[:, None] + rB_c[None, :], 1e-8)
                if clip_distance:
                    np.minimum(fallback, 1.0, out=fallback)   
            Db[invalid] = fallback[invalid]

        return start, end, Db

    def _distanceRect(self, src_shells, tgt_shells, n_src, n_tgt, weights,
                  src_centroids, src_radii_means, tgt_centroids, tgt_radii_means,
                  normalize_radii=True, clip_distance=True):
        """Distances (n_src, n_tgt) between two independently prepared record sets.
        Same math as _distanceBlock, but source and target are separate."""
        from scipy.spatial.distance import cdist

        Db = np.zeros((n_src, n_tgt), dtype=np.float32)
        Wb = np.zeros((n_src, n_tgt), dtype=np.float32)

        for shell, w in enumerate(weights):
            if w == 0:
                continue
            A, rA, idA, _, _ = src_shells[shell]
            B, rB, _, t_starts, t_ids = tgt_shells[shell]
            if len(A) == 0 or len(B) == 0:
                continue

            spatial = cdist(A, B)
            if normalize_radii:
                spatial /= np.maximum(rA[:, None] + rB[None, :], 1e-8)
                if clip_distance:
                    np.minimum(spatial, 1.0, out=spatial)

            s_starts = np.r_[0, np.flatnonzero(np.diff(idA)) + 1]
            s_ids = idA[s_starts]

            red = np.minimum.reduceat(spatial, s_starts, axis=0)
            red = np.minimum.reduceat(red, t_starts, axis=1)

            Db[np.ix_(s_ids, t_ids)] += w * red
            Wb[np.ix_(s_ids, t_ids)] += w

        valid = Wb > 0
        Db[valid] /= Wb[valid]
        invalid = ~valid
        if invalid.any():
            fb = cdist(src_centroids, tgt_centroids)
            if normalize_radii:
                fb /= np.maximum(src_radii_means[:, None] + tgt_radii_means[None, :], 1e-8)
                if clip_distance:
                    np.minimum(fb, 1.0, out=fb)
            Db[invalid] = fb[invalid]
        return Db

    def computeThroughput(self, records):
        """throughput = exp(-cost). Stored as `channel.throughput`"""
        costs = np.array([r.channel.cost if r.channel.cost is not None else np.inf for r in records])
        throughputs = np.exp(-costs)
        for r, t in zip(records, throughputs):
            r.channel.throughput = float(t)

##### centerline sampling functions######

####### POST CLUSTERING METRICS #######


    @staticmethod
    def liningSequenceDistance(rec_a, rec_b, alpha=0.5, beta=2.0):
        """
        Tversky distance between the lining-residue sets of two ChannelRecords,
        using their cached packed_lining bit vectors.

        alpha=0.5, beta=2.0 -> Jaccard distance
        alpha=0.5, beta=1.0 -> Dice distance
        """
        pa = rec_a.packed_lining
        pb = rec_b.packed_lining
        if pa is None or pb is None:
            raise ValueError("records must have packed_lining; prepare them with "
                             "metric='lining' (or pass atoms and trajectory)")

        table = ClusterCalculator.POPCOUNT_TABLE
        intersection = int(table[pa & pb].sum())
        a_only = int(table[pa & ~pb].sum())
        b_only = int(table[~pa & pb].sum())

        denominator = intersection + beta * (alpha * min(a_only, b_only)
                                             + (1 - alpha) * max(a_only, b_only))
        if denominator == 0:          
            return 0.0
        return float(1.0 - intersection / denominator)
    
    @staticmethod
    def voxelOverlapDistance(rec_a, rec_b):
        """
        1 - weighted voxel overlap score between two voxelized ChannelRecords.
        Non-overlapping channels give 1.0. alpha/beta/voxel_size/global_origin
        come from the records' calculator (set in prepareChannelRecords).
        """
        if not (rec_a.voxelized and rec_b.voxelized):
            raise ValueError("records must be voxelized; prepare them with "
                             "metric='voxel' or voxelize=True")
        calc = rec_a.calculator
        if rec_b.calculator is not calc:
            if (not np.isclose(calc.voxel_size, rec_b.calculator.voxel_size) or
                    not np.allclose(calc.global_origin, rec_b.calculator.global_origin)):
                raise ValueError("records were voxelized on different lattices")
        _, _, score = calc._scoreVoxelOverlapPrepared((0, 1, rec_a, rec_b))
        return float(score)

    def prepareCenterlineRecords(self, records, h=None, q=None):
        """
        Populate everything centerlineDistance needs on each record:
        rep_points / rep_radii (shared start point S from the full set),
        shell_weights (distance weights), centroid, radius_mean.
        Distances depend on S, so re-run this if the record set changes.
        """
        h = self.h if h is None else h
        q = self.q if q is None else q

        S = self.computeAverageStart(records)
        self.computeRepresentativePoints(records, S, h=h)

        weights = self.computeDistanceWeights(h, q).astype(np.float32)
        centroids, radii_means = self._prepareChannelCentroids(records)

        for rec, c, rm in zip(records, centroids, radii_means):
            rec.shell_weights = weights      # overwrites the placeholder ones
            rec.rep_points = c
            rec.rep_radii = float(rm)
        return records
    
    def centerlineDistanceToReps(self, new_records, rep_records, h=None, q=1.0,
                             rep_cache=None, normalize_radii=True, clip_distance=True,
                             batch_size=2000):
        """Vectorized (len(new_records), len(rep_records)) centerline distance."""
        if h is None:
            h = min(len(r.rep_points) for r in rep_records)
        weights = self.computeDistanceWeights(h, q).astype(np.float32)

        # Target side: build once and reuse across calls
        if rep_cache is None:
            rep_cache = (self._prepareAllShells(rep_records, h),
                        *self._prepareChannelCentroids(rep_records))
        tgt_shells, tgt_c, tgt_r = rep_cache

        out = np.empty((len(new_records), len(rep_records)), dtype=np.float32)
        for start in range(0, len(new_records), batch_size):
            batch = new_records[start:start + batch_size]
            src_shells = self._prepareAllShells(batch, h)
            src_c, src_r = self._prepareChannelCentroids(batch)
            out[start:start + len(batch)] = self._distanceRect(
                src_shells, tgt_shells, len(batch), len(rep_records), weights,
                src_c, src_r, tgt_c, tgt_r, normalize_radii, clip_distance)
        return out

    @staticmethod
    def centerlineDistance(rec_a, rec_b, normalize_radii=True, clip_distance=True):
        """
        Weighted shell-wise centerline distance between two ChannelRecords
        prepared with prepareCenterlineRecords. Matches _distanceBlock.
        """
        from scipy.spatial.distance import cdist

        if rec_a.shell_weights is None or rec_b.shell_weights is None:
            raise ValueError("records must be prepared with prepareCenterlineRecords "
                             "(or metric='centerline' in prepareChannelRecords)")

        weights = rec_a.shell_weights
        h = min(len(weights), len(rec_a.rep_points), len(rec_b.rep_points))

        d_sum = 0.0
        w_sum = 0.0
        for s in range(h):
            w = float(weights[s])
            if w == 0:
                continue
            A, B = rec_a.rep_points[s], rec_b.rep_points[s]
            if len(A) == 0 or len(B) == 0:
                continue

            spatial = cdist(A, B)
            if normalize_radii:
                rsum = rec_a.rep_radii[s][:, None] + rec_b.rep_radii[s][None, :]
                spatial = spatial / np.maximum(rsum, 1e-8)
                if clip_distance:
                    np.minimum(spatial, 1.0, out=spatial)

            # your block code does min over source reps, then min over target reps
            d_sum += w * spatial.min()
            w_sum += w

        if w_sum > 0:
            return float(d_sum / w_sum)

        # no shell populated by both channels -> centroid fallback
        d = float(np.linalg.norm(rec_a.centroid - rec_b.centroid))
        if normalize_radii:
            d /= max(rec_a.radius_mean + rec_b.radius_mean, 1e-8)
            if clip_distance:
                d = min(d, 1.0)
        return d


    @staticmethod
    def symmetricChamferDistance(channel_a, channel_b, n_samples=300, return_details=False):
        """
        Symmetric mean closest-centerline distance between two channels
        (Channel or ChannelRecord). Samples are uniform in arc length, so the
        plain mean approximates the length-weighted mean.
        """
        _, pa, _, _, tree_a = _arc(channel_a, n_samples)
        _, pb, _, _, tree_b = _arc(channel_b, n_samples)

        if len(pa) == 0 or len(pb) == 0:
            return (np.inf, {}) if return_details else np.inf

        d_ab, idx_ab = tree_b.query(pa, k=1)
        d_ba, idx_ba = tree_a.query(pb, k=1)

        mean_ab = float(np.mean(d_ab))
        mean_ba = float(np.mean(d_ba))
        distance = 0.5 * (mean_ab + mean_ba)

        if not return_details:
            return float(distance)

        return float(distance), {
            "distance": float(distance),
            "a_to_b_mean": mean_ab, "b_to_a_mean": mean_ba,
            "a_to_b_median": float(np.median(d_ab)), "b_to_a_median": float(np.median(d_ba)),
            "a_to_b_q90": float(np.percentile(d_ab, 90)), "b_to_a_q90": float(np.percentile(d_ba, 90)),
            "a_to_b_max": float(np.max(d_ab)), "b_to_a_max": float(np.max(d_ba)),
            "a_to_b_distances": d_ab, "b_to_a_distances": d_ba,
            "a_to_b_indices": idx_ab, "b_to_a_indices": idx_ba,
        }
    
    @staticmethod
    def splineProjectionDistance(channel_a, channel_b, n_samples=300,angle_threshold_deg=60.0,
                                 n_dense=1000, n_candidates=8, return_details=False):
        """
        Mean distance from arc-length samples of the shorter channel to the
        spline of the longer channel.

        A dense KD-tree lookup identifies candidate points on the longer channel.
        For each sample on the shorter channel, the nearest dense point is used
        as an initial match, then the sample is refined by projecting onto the
        two neighboring dense spline segments and evaluating the spline at the
        projected position. The computation is fully vectorized; n_dense controls
        the sampling resolution (chord error is about spacing^2 * curvature / 8,
        about 1e-4 A at the default setting).

        angle_threshold_deg : among the `n_candidates` nearest dense points, use
            the nearest one whose tangent direction is within this angle of the
            reference tangent; otherwise fall back to the plain nearest neighbor.
            None disables the tangent-direction constraint.
        """
        la = getattr(channel_a, "channel", channel_a).length
        lb = getattr(channel_b, "channel", channel_b).length
        short, long_ = (channel_a, channel_b) if la <= lb else (channel_b, channel_a)
        short_ch = getattr(short, "channel", short)
 
        ts, ps, rs, ss, _ = _arc(short, n_samples)
        n = len(ps)
        if n == 0:
            return (np.inf, {}) if return_details else np.inf
 
        spline, td, pd, tan_dense, tree = _dense(long_, n_dense)
        tan_a = _cached(short, ("arc_tan", n_samples),
                        lambda: _unit(short_ch.centerline_spline(ts, 1)))
 
        # initial correspondence (optionally tangent-constrained)
        use_angle = angle_threshold_deg is not None and n_candidates > 1
        _, idx = tree.query(ps, k=min(n_candidates, len(td)) if use_angle else 1)
        idx = np.asarray(idx).reshape(n, -1)
        n_fallback = 0
        if use_angle:
            cos = np.abs(np.einsum("nd,nkd->nk", tan_a, tan_dense[idx]))
            ok = cos >= np.cos(np.radians(angle_threshold_deg))
            has = ok.any(axis=1)
            sel = idx[np.arange(n), np.where(has, ok.argmax(axis=1), 0)]
            n_fallback = int((~has).sum())
        else:
            sel = idx[:, 0]
 
        # refine: closest point on the two dense segments next to the start
        # (chord error ~ h^2*curvature/8, negligible for n_dense ~ 1000), then
        # ONE spline evaluation at the resulting parameter
        nd = len(td)
        left_seg = np.clip(sel - 1, 0, nd - 2)
        best_d2 = np.full(n, np.inf)
        t_best = td[sel].astype(float)
        for seg in (left_seg, np.clip(sel, 0, nd - 2)):
            a, ab = pd[seg], pd[seg + 1] - pd[seg]
            u = np.clip(np.einsum("ij,ij->i", ps - a, ab)
                        / np.maximum(np.einsum("ij,ij->i", ab, ab), 1e-24), 0.0, 1.0)
            delta = a + u[:, None] * ab - ps
            d2 = np.einsum("ij,ij->i", delta, delta)
            better = d2 < best_d2
            best_d2 = np.where(better, d2, best_d2)
            t_best = np.where(better, td[seg] + u * (td[seg + 1] - td[seg]), t_best)
        matched_points = np.asarray(spline(t_best))
        distances = np.linalg.norm(matched_points - ps, axis=1)
        distance = float(np.mean(distances))
 
        if not return_details:
            return distance
 
        # split displacement along / across the reference tangent
        tangential = np.abs(np.einsum("ij,ij->i", matched_points - ps, tan_a))
        perpendicular = np.sqrt(np.maximum(distances ** 2 - tangential ** 2, 0.0))
 
        return distance, {
            "distance": distance,
            "reference_is_a": short is channel_a,
            "reference_parameters": ts, "reference_points": ps,
            "reference_radii": rs, "reference_arc_length": ss,
            "matched_parameters": t_best, "matched_points": matched_points,
            "distances": distances,
            "mean": distance,
            "median": float(np.median(distances)),
            "q90": float(np.percentile(distances, 90)),
            "q95": float(np.percentile(distances, 95)),
            "max": float(np.max(distances)),
            "long_channel_length": float(max(la, lb)),
            "short_channel_length": float(min(la, lb)),
            "tangential_component": tangential,     # (these two names were swapped before)
            "normal_component": perpendicular,
            "n_angle_fallback": n_fallback,
        }

    @staticmethod
    def frechetDistance(channel_a, channel_b, n_samples=300, return_details=False):
        """Discrete Frechet distance between arc-length-sampled centerlines."""
        from scipy.spatial.distance import cdist
        
        _, pa, ra, sa, _ = _arc(channel_a, n_samples)
        _, pb, rb, sb, _ = _arc(channel_b, n_samples)
        na, nb = len(pa), len(pb)
        if na == 0 or nb == 0:
            return (np.inf, {}) if return_details else np.inf

        D = cdist(pa, pb)
        if not return_details:
            return _frechetTable(D)
        ca = _frechetTable(D, full=True)
        distance = float(ca[-1, -1])

        i, j = na - 1, nb - 1
        path = [(i, j)]
        while i > 0 or j > 0:
            if i == 0:
                j -= 1
            elif j == 0:
                i -= 1
            else:
                _, i, j = min([(ca[i - 1, j], i - 1, j),
                               (ca[i, j - 1], i, j - 1),
                               (ca[i - 1, j - 1], i - 1, j - 1)], key=lambda x: x[0])
            path.append((i, j))

        return distance, {"distance": distance,
                        "a_points": pa, "b_points": pb,
                        "a_radii": ra, "b_radii": rb,
                        "a_arc_length": sa, "b_arc_length": sb,
                        "coupling_path": np.asarray(path[::-1], dtype=int),
                        "n_a": na, "n_b": nb}

####### POST CLUSTERING METRICS #######

#### post-clustering #####
    _PAIRWISE_ITEMS_A = None
    _PAIRWISE_ITEMS_B = None
    _PAIRWISE_FN = None
    _PAIRWISE_KWARGS = None


    def _initPairwiseWorker(self, items_a, items_b, fn, kwargs):
        global _PAIRWISE_ITEMS_A, _PAIRWISE_ITEMS_B
        global _PAIRWISE_FN, _PAIRWISE_KWARGS

        _PAIRWISE_ITEMS_A = items_a
        _PAIRWISE_ITEMS_B = items_b
        _PAIRWISE_FN = fn
        _PAIRWISE_KWARGS = kwargs


    def _pairwiseDistanceChunk(self, pair_chunk):
        ii, jj = pair_chunk
        distances = np.empty(len(ii), dtype=float)
        for k, (i, j) in enumerate(zip(ii, jj)):
            distances[k] = _PAIRWISE_FN(_PAIRWISE_ITEMS_A[i],
                                        _PAIRWISE_ITEMS_B[j],
                                        **_PAIRWISE_KWARGS)

        return ii, jj, distances


    def pairwiseDistance(self,items_a, items_b=None, metric='chamfer', max_proc=None, mp_context=None, dtype=float, **kwargs):
        """
        Calculate pairwise distances between ChannelRecords or Channels.

        Three modes are supported.

        Symmetric
        ---------
        pairwiseDistance(items_a)

            Computes only the upper triangle and returns an (n, n)
            symmetric distance matrix.

        Rectangular
        -----------
        pairwiseDistance(items_a, items_b)

            Computes all distances between the two collections and returns
            an (len(items_a), len(items_b)) matrix.

        One-to-many
        -----------
        pairwiseDistance(item, items_b)

            Computes the distance from one item to every item in items_b
            and returns a 1-D array of shape (len(items_b),).

        Parameters
        ----------
        items_a : ChannelRecord, Channel, or iterable
            First item or collection.

        items_b : iterable, optional
            Second collection. If None, symmetric pairwise distances are
            calculated within items_a.

        metric : str or callable
            One of:
                'chamfer'
                'frechet'
                'projection'

            or a callable:
                metric(a, b, **kwargs) -> float

        max_proc : int or None
            Number of worker processes.
            1 performs serial calculation.
            None uses approximately half available CPUs.

        mp_context : str or None
            Multiprocessing context, e.g. 'fork', 'spawn'.

        dtype : numpy dtype
            Output dtype.

        **kwargs
            Forwarded to the distance function.

        Returns
        -------
        ndarray
            Symmetric mode:
                (n, n)

            Rectangular mode:
                (len(items_a), len(items_b))

            One-to-many mode:
                (len(items_b),)
        """
        import multiprocessing as mp

        metrics = {'chamfer': self.symmetricChamferDistance,
                    'frechet': self.frechetDistance,
                    'projection': self.splineProjectionDistance}

        if callable(metric):
            fn = metric
        elif metric in metrics:
            fn = metrics[metric]
        else:
            raise ValueError(f"unknown metric {metric!r}; "
                            f"use one of {list(metrics)} or a callable")

        kwargs.pop("return_details", None)

        symmetric = items_b is None

        if symmetric:
            items_a = list(items_a)
            items_b = items_a
            one_to_many = False

        else:
            try:
                iter(items_a)
                items_a = list(items_a)
                one_to_many = False
            except TypeError:
                items_a = [items_a]
                one_to_many = True

            items_b = list(items_b)

        n_a = len(items_a)
        n_b = len(items_b)

        if symmetric:
            if n_a == 0:
                return np.empty((0, 0), dtype=dtype)

        elif one_to_many:
            if n_b == 0:
                return np.empty(0, dtype=dtype)

        elif n_a == 0 or n_b == 0:
            return np.empty((n_a, n_b), dtype=dtype)

        if symmetric:
            ii, jj = np.triu_indices(n_a,k=1)
            D = np.zeros( (n_a, n_a), dtype=dtype)
        else:
            ii = np.repeat(np.arange(n_a),n_b)
            jj = np.tile( np.arange(n_b),n_a)
            D = np.empty( (n_a, n_b),dtype=dtype)
        n_pairs = len(ii)
        if n_pairs == 0:
            if one_to_many:
                return D[0]
            return D

        if max_proc == 1:
            for i, j in zip(ii, jj):
                D[i, j] = fn(items_a[i],items_b[j],**kwargs)

        else:
            if mp_context is None:
                ctx = mp.get_context()
            else:
                ctx = mp.get_context(mp_context)

            if max_proc is None:
                max_proc = max(1, ctx.cpu_count() // 2)

            chunksize = max(1,n_pairs // (max_proc * 4))
            pair_chunks = ((ii[start:start + chunksize], jj[start:start + chunksize])
                                    for start in range(0, n_pairs,chunksize))

            with ctx.Pool(processes=max_proc,initializer=self._initPairwiseWorker,
                            initargs=(items_a,items_b,fn,kwargs)) as pool:
                for ci, cj, distances in pool.imap_unordered(self._pairwiseDistanceChunk,pair_chunks,chunksize=1):
                    D[ci, cj] = distances
        if symmetric:
            D += D.T
        if one_to_many:
            return D[0]
        return D

    def pairwiseDistanceMatrix(self,items,metric='chamfer',max_proc=None,mp_context=None,**kwargs):
        return self.pairwiseDistance(items, metric=metric, max_proc=max_proc, mp_context=mp_context,**kwargs)


    def clusterDistanceMatrix(self, cluster, metric='chamfer', store=False,D=None, **kwargs):
        metrics = {'chamfer': self.symmetricChamferDistance,
                'frechet': self.frechetDistance,
                'projection': self.splineProjectionDistance}
        if callable(metric):
            fn = metric
        elif metric in metrics:
            fn = metrics[metric]
        else:
            raise ValueError(f"unknown metric {metric!r}; use one of {list(metrics)} or a callable")

        if D is None:
            D = self.pairwiseDistanceMatrix(cluster.records, fn, **kwargs)
        if store:
            n = len(D)
            nn = (D + np.diag(np.full(n, np.inf))).min(axis=1) if n > 1 else np.zeros(n)
            cluster.internal_distance_matrix = D
            cluster.mean_pairwise_distance = float(D[np.triu_indices(n, 1)].mean()) if n > 1 else 0.0
            cluster.max_nearest_neighbor_distance = float(nn.max()) if n else 0.0
            cluster.details["medoid_local"] = int(np.argmin(D.sum(axis=1))) if n else None
            cluster.details["distance_metric"] = getattr(fn, "__name__", "custom")
        return D

    def clusteringDistanceMatrices(self, clustering, metric=None,
                                   reuse_clustering_matrix=False, **kwargs):
        """
        clusterDistanceMatrix for every cluster in a ClusteringRecords.
        Returns {cluster.label: matrix} and sets cluster.mediod to the GLOBAL
        index into clustering.records when all its records are in there.

        reuse_clustering_matrix=True slices the square clustering.distance_matrix
        instead of recomputing -- only valid if it used the same metric/params.
        """
        clusters = clustering.clusters
        clusters = clusters.values() if isinstance(clusters, dict) else clusters
        gidx = {id(r): i for i, r in enumerate(clustering.records)}
        full = np.asarray(clustering.distance_matrix) if reuse_clustering_matrix else None
        out = {}
        for cl in clusters:
            ids = [gidx.get(id(r)) for r in cl.records]
            known = None not in ids
            D = full[np.ix_(ids, ids)] if full is not None and known else None
            out[cl.label] = self.clusterDistanceMatrix(cl, metric, D=D, **kwargs)
        return out

class Cluster:
    """ Stores information of individual clusters """

    def __init__(self, label, records, details=None, parent=None, depth=0):
        self.label = label
        self.records = records
        self.details = details if details is not None else {}   # was a shared mutable default
        self.conformation_ids = [rec.conformation_id for rec in records]
        self.approximated = None
        self.soft_assigned = None
        self.approximated_indices = []
        self.soft_assigned_indices = []

        self.parent = parent
        self.children = []
        self.local_label = None

        self.size = len(records)
        self.persistence = None
        self.mediod = None              # GLOBAL index into the fitted subsample
        self.mean_pairwise_distance = None
        self.density = None
        self.max_nearest_neighbor_distance = None
        self.internal_distance_matrix = None
        self.mean_probability = None    # only if hdbscan
        self.voxelized_graph = None
        self.graph_occupancy = None
        #ClusterCalculator.ClusteringRecords()
        # set by ClusteringRecords.fit(); refer to the ORIGINAL subsample,
        # so they stay valid after extendCluster()
        self.rep_idx = None
        self.assignment_rep_idx = None
        self.rep_threshold = None     # acceptance radius computed on the reps

        self.depth = depth
        self.distance_matrix = None #if subcluster
        self.linkage_matrix = None #if subcluster
        self.model = None #if subcluster


    def _threshold(self):
        """Representative acceptance radius or fallback cluster scale."""
        if self.rep_threshold is not None:
            return self.rep_threshold
        vals = [x for x in (self.mean_pairwise_distance, self.max_nearest_neighbor_distance) if x is not None]
        return max(vals) if vals else 0.0

    def _resetStats(self):
        """Clear summaries invalidated when cluster membership changes."""
        self.persistence = None
        self.mean_pairwise_distance = None
        self.density = None
        self.max_nearest_neighbor_distance = None
        self.internal_distance_matrix = None
        self.mean_probability = None
        self.voxelized_graph = None
        self.graph_occupancy = None
    
    def extendCluster(self, records, indices=None,approximated=False, soft_assigned=False):
        """Append records and optionally record their indices in the owning ClusteringRecords.records collection. """
        records = list(records)
        self.records = list(self.records) + records
        self.conformation_ids.extend(rec.conformation_id for rec in records)
        if indices is not None:
            indices = [int(i) for i in indices]
            if approximated:
                self.approximated_indices.extend(indices)
            if soft_assigned:
                self.soft_assigned_indices.extend(indices)
        if approximated:
            self.approximated = True
        if soft_assigned:
            self.soft_assigned = True
        self.size = len(self.records)
        self._resetStats()

    def removeRecords(self, records, indices=None):
        """Remove records (matched by identity) from this cluster. Returns count removed.
        `indices` are the records' global indices, used to clean the index lists."""
        drop = {id(r) for r in records}
        keep = [r for r in self.records if id(r) not in drop]
        n_removed = len(self.records) - len(keep)
        self.records = keep
        self.conformation_ids = [r.conformation_id for r in keep]
        self.size = len(keep)

        if indices is not None:
            gone = {int(i) for i in indices}
            self.approximated_indices = [i for i in self.approximated_indices if i not in gone]
            self.soft_assigned_indices = [i for i in self.soft_assigned_indices if i not in gone]
            if self.mediod is not None and int(self.mediod) in gone:
                LOGGER.warn('Medoid of cluster {0} was reassigned; medoid is now None'.format(self.label))
                self.mediod = None
            if self.rep_idx is not None and gone & {int(i) for i in self.rep_idx}:
                LOGGER.warn('Representative(s) of cluster {0} were reassigned'.format(self.label))

        self._resetStats()
        return n_removed

    def _calcInternalDistanceMatrix(self,distance_matrix,labels):
        """Store the pairwise distances for members belonging to this cluster."""
        D = np.asarray(distance_matrix)
        labels = np.asarray(labels)
        idx = np.flatnonzero(labels == self.label)
        if len(idx) == 0:
            raise ValueError("Cluster {0} does not exist.".format(self.label))
        Dc = D[np.ix_(idx, idx)]
        self.internal_distance_matrix = Dc
    
    def _calcPersistance(self,num_conformations):
        """Compute the cluster persistence fraction relative to the full sample."""
        self.persistance = self.size/num_conformations

    def _calcClusterStats(self,labels):
        """Compute medoid, mean pairwise distance, max nearest-neighbor distance, and density."""
        labels = np.asarray(labels)
        idx = np.flatnonzero(labels == self.label)
        n = len(idx)            # not self.size: size grows after extendCluster()
        Dc = self.internal_distance_matrix
        self.mediod = int(idx[np.argmin(Dc.sum(axis=1))])
        if n > 1:
            upper = Dc[np.triu_indices(n, k=1)]

            mean_pairwise_distance = upper.mean()
        else:
            mean_pairwise_distance = 0.0
        self.mean_pairwise_distance = mean_pairwise_distance
        # For every point, find its closest OTHER point in the cluster.
        # Then take the largest of those nearest-neighbor distances.
        if n > 1:
            Dc_no_diag = Dc.copy()
            np.fill_diagonal(Dc_no_diag, np.inf)
            nearest_neighbor_distances = Dc_no_diag.min(axis=1)
            self.max_nearest_neighbor_distance = nearest_neighbor_distances.max()
        else:
            self.max_nearest_neighbor_distance = 0.0
        # Here density is defined as inverse mean pairwise distance.
        # Higher = more compact.
        if n > 1 and mean_pairwise_distance > 0:
            self.density = 1.0 / mean_pairwise_distance
        else:
            self.density = np.inf

    def buildVoxelizedGraph(self, min_occupancy=0.0, return_coordinates=False):
        """Build a cluster occupancy grid by aggregating the voxel masks of all records.

        Parameters
        ----------
        min_occupancy : float, optional
            Minimum occupancy value required to keep a voxel in the graph.
        return_coordinates : bool, optional
            If True, return the occupancy array together with the global origin and
            voxel size.

        Returns
        -------
        occupancy : numpy.ndarray
            Occupancy map for the cluster, optionally accompanied by origin and voxel
            size metadata.
        """
        if not np.any(self.records):
            self.graph_occupancy = None
            self.voxelized_graph = None
            if return_coordinates:
                return (np.zeros((0, 0, 0), dtype=np.float32),None,None)
            
        calculator = self.records[0].calculator
        voxel_size = calculator.voxel_size
        global_origin = np.asarray(calculator.global_origin, dtype=float)
        for record in self.records:
            if not record.voxelized:
                record._prepareVoxels()

        global_min = np.min(np.vstack([record.offset 
                            for record in self.records]),
                            axis=0)

        global_max = np.max(np.vstack([ record.offset + record.dims
                            for record in self.records
                            ]),axis=0)

        global_dims = global_max - global_min
        occupancy_count = np.zeros(global_dims,dtype=np.uint16)
        for record in self.records:
            start = record.offset - global_min
            end = start + record.dims
            sub = occupancy_count[start[0]:end[0],
                                start[1]:end[1],
                                start[2]:end[2]]

            sub += record.mask.astype(np.uint16)
        occupancy = occupancy_count.astype(np.float32) / self.size
        if min_occupancy > 0:
            occupancy[occupancy < min_occupancy] = 0.0
        occupancy_origin = (global_origin + global_min * voxel_size)
        self.graph_occupancy = occupancy
        self.voxelized_graph = {
            "occupancy": occupancy,
            "origin": occupancy_origin,
            "voxel_size": voxel_size,
            "dims": global_dims,
            "global_offset": global_min,
            "n_channels": self.size}

        if return_coordinates:
            return occupancy, occupancy_origin, voxel_size

        return occupancy
    
    def writeVoxelOccupancyPQR(self,output_file, min_occupancy=0.0, occupancy_as_charge=True):
        """
        Write the cluster voxel occupancy as a PQR-like FIL file.

        Each occupied voxel is written as one FIL pseudoatom. The occupancy
        value is stored in the PQR charge field.

        Parameters
        ----------
        output_file : str
            Output PQR filename.

        min_occupancy : float, optional
            Only voxels with occupancy >= this value are written.
            Default is 0.0.

        occupancy_as_charge : bool, optional
            Store occupancy in the PQR charge field. Default True.

        Returns
        -------
        output_file : str
        """
        if self.graph_occupancy is None:
            self.buildVoxelizedGraph(min_occupancy=min_occupancy)

        occupancy = self.graph_occupancy
        info = self.voxelized_graph
        if occupancy is None:
            return None
        origin = np.asarray(info["origin"], dtype=float)
        voxel_size = float(info["voxel_size"])
        indices = np.argwhere(occupancy >= min_occupancy)
        if len(indices) == 0:
            with open(output_file, "w"):
                pass
            return output_file
        values = occupancy[indices[:, 0],
                            indices[:, 1],
                            indices[:, 2]]
        coordinates = (origin + (indices + 0.5) * voxel_size)

        with open(output_file, "w") as out:
            for atom_id, (xyz, occ) in enumerate(zip(coordinates, values),start=1):
                x, y, z = xyz
                charge = float(occ) if occupancy_as_charge else 0.0
                # Use voxel radius as the PQR radius.
                radius = voxel_size * 0.5
                label = self.label if self.label is not None else -1
                out.write(
                    f"ATOM  {atom_id:5d}  H   FIL T{label:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{charge:8.4f}{radius:7.3f}\n"
                )
                
        return output_file

class ClusteringRecords:
    """ Stores the results of a clustering run, including the fitted model, the
    records, the clusters, the labels, and the distance matrix. Also stores
    the linkage matrix if hierarchical clustering was used. Provides methods
    for reassignment of clusters, subclustering, noise assignment
    and other post-clustering analysis. """

    def __init__(self,model, records, clusters, labels, distance_matrix, linkage_matrix=None, levels=None):
        self.model = model
        self.records = records
        self.linkage_matrix = linkage_matrix
        self.clusters = clusters
        self.labels = labels
        self.distance_matrix = distance_matrix
        self.num_clusters = len(clusters)
        self.levels = [0] if levels is None else levels
        self.between_distances = None
        self.approximated = False
        self.approximated_distance_matrices = []
        self._by_label = {}
        self._n_sub = len(self.labels)
        self._fitted = False
        self.sub_clusterings = None
        self._rep_idx = None
        self._rep_records = None
        self._rep_lookup = None

    def getClusterByLabel(self, label):
        if self._by_label:
            cluster = self._by_label.get(int(label))
            if cluster is not None:
                return cluster
        for cluster in self.clusters:
            if cluster.label == label:
                return cluster
        LOGGER.warn('No Cluster With Label {0} found'.format(label))

        return None

    @staticmethod
    def _clusterPath(cluster):
        """Labels from the top-level cluster down to `cluster`."""
        path, node = [], cluster
        while node is not None:
            path.append(node.label)
            node = node.parent
        path.reverse()
        return path

    def clusteringDistanceStats(self):
        """
        Summarize pairwise distances within and between clusters.

        Parameters
        ----------
        D : array-like
            Square distance matrix of shape ``(n_channels, n_channels)``.
        labels : array-like
            Cluster assignments for each channel. Entries with negative labels are
            ignored.

        Returns
        -------
        dict
            Dictionary containing the mean and median distance values for
            within-cluster and between-cluster pairs, along with the separation
            ratio defined as the mean between-cluster distance divided by the mean
            within-cluster distance.
        """
        D = self.distance_matrix
        labels = np.asarray(self.labels)
        valid = labels >= 0
        labels = labels[valid]
        D = D[np.ix_(valid, valid)]
        clusters = np.unique(labels)
        within = []
        between = []
        for i, c1 in enumerate(clusters):
            idx1 = labels == c1
            # within-cluster
            if idx1.sum() > 1:
                vals = D[np.ix_(idx1, idx1)]
                within.append(vals[np.triu_indices_from(vals, k=1)])
            # between-cluster
            for c2 in clusters[i+1:]:
                idx2 = labels == c2
                between.append(D[np.ix_(idx1, idx2)].ravel())
        within = np.concatenate(within) if within else np.array([])
        between = np.concatenate(between) if between else np.array([])
        return {"mean_within": np.mean(within),
                "median_within": np.median(within),
                "mean_between": np.mean(between),
                "median_between": np.median(between),
                "separation_ratio":np.mean(between) / np.mean(within)}


    def silhouetteScore(self,ignore_noise=True):
        """
        Calculate the silhouette coefficient for each sample using a
        precomputed pairwise distance matrix.
        Parameters
        ----------
        distance_matrix : (N, N) array_like
            Symmetric pairwise distance matrix. Diagonal should be zero.
        labels : (N,) array_like
            Cluster label for each sample.
        Returns
        -------
        silhouettes : (N,) ndarray
            Silhouette coefficient for each sample.
        """
        def silhouetteSamples(distance_matrix, labels):

            D = np.asarray(distance_matrix, dtype=float)
            labels = np.asarray(labels)
            n = len(labels)
            if D.shape != (n, n):
                raise ValueError("distance_matrix must have shape (N, N)")

            unique_labels = np.unique(labels)
            if len(unique_labels) < 2:
                raise ValueError("Silhouette score requires at least 2 clusters")

            silhouettes = np.zeros(n, dtype=float)
            for label in unique_labels:
                idx = np.flatnonzero(labels == label)
                n_cluster = len(idx)
                #singleton clusters get silhouette = 0
                if n_cluster == 1:
                    silhouettes[idx] = 0.0
                    continue
                # Mean intra-cluster distance.
                # Diagonal is zero, so divide by n_cluster - 1.
                a = D[np.ix_(idx, idx)].sum(axis=1) / (n_cluster - 1)
                # Mean distance to the nearest other cluster
                b = np.full(n_cluster, np.inf)
                for other_label in unique_labels:
                    if other_label == label:
                        continue

                    other_idx = np.flatnonzero(labels == other_label)
                    mean_dist = D[np.ix_(idx, other_idx)].mean(axis=1)
                    b = np.minimum(b, mean_dist)

                denom = np.maximum(a, b)
                silhouettes[idx] = np.divide(b - a, denom, out=np.zeros_like(a), where=denom != 0)
            return silhouettes
        
        D = np.asarray(self.distance_matrix)
        labels = np.asarray(self.labels)
        if ignore_noise:
            keep = labels != -1
            D = D[np.ix_(keep, keep)]
            labels = labels[keep]

        return silhouetteSamples(D, labels).mean()

    def frameClusterDistributions(self):
        labels = np.asarray(self.labels)
        frames = np.asarray([record.conformation_id for record in self.records])

        unique_frames = np.unique(frames)
        clusters = np.unique(labels[labels >= 0])

        distributions = np.zeros((len(unique_frames), len(clusters)))

        for i, f in enumerate(unique_frames):
            lf = labels[frames == f]

            # Only clustered records contribute to distribution
            lf = lf[lf >= 0]

            if len(lf):
                counts = np.array([np.sum(lf == c) for c in clusters])
                distributions[i] = counts / counts.sum()

        return distributions, clusters

    def calcTemporalCoherence(self,by_cluster=False):
        P, clusters = self.frameClusterDistributions()
        if by_cluster:
            if len(P) < 2:
                return {c: np.nan for c in clusters}
            delta = np.abs(np.diff(P, axis=0))
            coherence = 1.0 - np.mean(delta, axis=0)
            return {cluster: score for cluster, score in zip(clusters, coherence)}
        from scipy.spatial.distance import jensenshannon
        js = [jensenshannon(P[i], P[i+1]) for i in range(len(P)-1)]
        return 1 - np.mean(js)
    

    def fit(self, max_reps_per_cluster=500, random_state=123, keep_internal=True):
        """Choose representatives and an acceptance threshold per cluster,
        using only the (n_sub, n_sub) distance matrix.

        Cluster-level stats (medoid etc.) are computed here if the Cluster
        doesn't have them yet. Representatives are a random subset that
        always includes the medoid.
        """
        D = np.asarray(self.distance_matrix, dtype=float)
        if D.ndim != 2 or D.shape[0] != D.shape[1]:
            raise ValueError("distance_matrix must be square (n_sub, n_sub)")
        if D.shape[0] != len(self.labels):
            raise ValueError("distance_matrix size must match len(labels)")

        rng = np.random.default_rng(random_state)
        self._n_sub = D.shape[0]
        self._by_label = {}

        for cluster in self.clusters:
            member_idx = np.flatnonzero(self.labels == cluster.label)
            if len(member_idx) == 0:
                LOGGER.warn('Cluster {0} has no members in labels'.format(cluster.label))
                continue

            if cluster.internal_distance_matrix is None or cluster.mediod is None:
                cluster._calcInternalDistanceMatrix(D, self.labels)
                cluster._calcClusterStats(self.labels)

            n = len(member_idx)
            # cluster.mediod is a global index -> position within this cluster
            medoid_local = int(np.flatnonzero(member_idx == cluster.mediod)[0])

            if max_reps_per_cluster is not None and n > max_reps_per_cluster:
                others = np.delete(np.arange(n), medoid_local)
                pick = rng.choice(others, size=max_reps_per_cluster - 1, replace=False)
                rep_local = np.concatenate([[medoid_local], pick])
            else:
                rep_local = np.arange(n)
            rep_idx = member_idx[rep_local]

            if len(rep_idx) > 1:
                D_rep = D[np.ix_(rep_idx, rep_idx)]          # fancy indexing copies
                mean_pw = float(D_rep[np.triu_indices(len(rep_idx), k=1)].mean())
                D_nn = D_rep.copy()
                np.fill_diagonal(D_rep, np.inf)
                max_nn = float(D_rep.min(axis=1).max())
            else:
                mean_pw, max_nn = 0.0, 0.0

            cluster.rep_idx = rep_idx
            cluster.rep_threshold = max(mean_pw, max_nn)
            cluster.approximated = len(rep_idx) < n
            if not keep_internal:
                cluster.internal_distance_matrix = None

            self._by_label[int(cluster.label)] = cluster

        if self._by_label:
            self._rep_idx = np.unique(np.concatenate([cluster.rep_idx 
                                for cluster in self._by_label.values()])).astype(int)
        else:
            self._rep_idx = np.empty(0,dtype=int)
        self._rep_records = [self.records[i] for i in self._rep_idx]
        self._rep_lookup = {int(global_idx): local_idx for local_idx, global_idx in enumerate(self._rep_idx)}

        for cluster in self._by_label.values():
            cluster.assignment_rep_idx = np.array([self._rep_lookup[int(i)]
                                            for i in cluster.rep_idx ], dtype=int)

        self.num_clusters = len(self._by_label)
        self.approximated = True
        self._fitted = True
        return self

    def _checkFitted(self):
        if not self._fitted or not self._by_label:
            raise RuntimeError("call fit() first")
        
    def _getCalculator(self, calculator=None):
        """
        Resolve the ClusterCalculator used to calculate new distances.
        Explicit calculator takes priority. Otherwise use the calculator
        cached by the first representative ChannelRecord.
        """
        if calculator is not None:
            return calculator
        if not self._rep_records:
            raise RuntimeError("No representative records available")
        calculator = getattr(self._rep_records[0], "calculator", None)
        if calculator is None:
            raise RuntimeError("Representative ChannelRecords do not have a calculator")
        return calculator

    def _distanceBatches(self, new_records, calculator, metric, batch_size, max_proc, mp_context, **metric_kwargs):
        """
        Yield distances from batches of new records to the compressed
        representative set.

        Yields
        ------
        start, end, D

        where:
            D.shape == (end-start, n_representatives)
        """
        n = len(new_records)
        #####
        if metric == 'centerline':
            # cache the target side once per fit
            if getattr(self, '_rep_cache', None) is None:
                h = min(len(r.rep_points) for r in self._rep_records)
                self._rep_cache = (calculator._prepareAllShells(self._rep_records, h),
                                *calculator._prepareChannelCentroids(self._rep_records))
            for start in range(0, n, batch_size):
                end = min(start + batch_size, n)
                D = calculator.centerlineDistanceToReps(
                    new_records[start:end], self._rep_records,
                    q=metric_kwargs.get('q', calculator.q),
                    rep_cache=self._rep_cache, batch_size=end - start,
                    normalize_radii=metric_kwargs.get('normalize_radii', True),
                    clip_distance=metric_kwargs.get('clip_distance', True))
                yield start, end, D.astype(float)
            return
    # original generic path for other metrics
        #####
        for start in range(0,n,batch_size):
            end = min(start + batch_size, n)
            batch = new_records[start:end]
            D = calculator.pairwiseDistance(batch, self._rep_records, metric=metric, max_proc=max_proc,
                                            mp_context=mp_context, **metric_kwargs)

            yield (start, end, np.asarray(D, dtype=float))

    def _distToClusters(self, block, labels_list):
        """(b, n_sub) block -> (b, k): min distance to each cluster's reps."""
        n_reps = len(self._rep_records)
        if block.shape[1] != self._n_sub:
            raise ValueError(f"expected {self._n_sub} columns (one per subsample "
                                f"point), got {block.shape[1]}")
        return np.stack([block[:, self._by_label[l].rep_idx].min(axis=1)
                            for l in labels_list], axis=1)

    def assign(self, new_records, metric = 'centerline', calculator=None, n_new=None, batch_size=2000, strict=True,
                slack=1.0, atol=1e-9, max_proc=None, mp_context=None, **metric_kwargs):
        """Hard-assign new points; -1 = unassigned. See earlier docstring."""
        self._checkFitted()
        new_records = list(new_records)
        n = len(new_records)
        if n == 0: 
            return np.empty(0,dtype=int)
        calculator = self._getCalculator(calculator)
        labels_list = list(self._by_label)
        lbl_arr = np.array(labels_list)
        thr = np.array([self._by_label[l]._threshold() for l in labels_list]) * slack
        thr = np.where(thr == 0.0, atol, thr)
        out = np.full(n, -1, dtype=int)

        for start, end, block in self._distanceBatches( new_records, calculator, metric,
                                        batch_size, max_proc, mp_context, **metric_kwargs):
            d = self._distToClusters(block, labels_list)
            best_j = d.argmin(axis=1)
            best_d = d[np.arange(len(d)), best_j]
            accept = np.ones(len(d), bool) if not strict else best_d <= thr[best_j]
            out[start:end] = np.where(accept, lbl_arr[best_j], -1)
        return out

    def softAssign(self, new_records, metric ='centerline', calculator = None, 
                max_distance=None, bandwidth=None, top_k=None, batch_size=2000,
                    max_proc=None, mp_context=None, **metric_kwargs):
        """Soft-assign borderline points. Returns (labels_list, probs)."""
        self._checkFitted()
        new_records = list(new_records)
        n = len(new_records)
        if n == 0: 
            return np.empty(0,dtype=int)
        calculator = self._getCalculator(calculator)
        labels_list = list(self._by_label)
        k = len(labels_list)
        thr = np.array([self._by_label[l]._threshold() for l in labels_list])
        thr = thr[thr > 0]
        if max_distance is None:
            max_distance = float(2.0 * thr.max()) if thr.size else 1.0
        if bandwidth is None:
            bandwidth = max(float(np.median(thr)) if thr.size else 1.0, 1e-6)

        probs = np.zeros((n, k))
        for start, end, block in self._distanceBatches(new_records, calculator, metric,
                                        batch_size, max_proc, mp_context, **metric_kwargs):
            d = self._distToClusters(block, labels_list)
            w = np.where(d <= max_distance, np.exp(-d / bandwidth), 0.0)

            if top_k is not None and top_k < k:
                keep = np.argpartition(-w, top_k - 1, axis=1)[:, :top_k]
                mask = np.zeros_like(w, dtype=bool)
                np.put_along_axis(mask, keep, True, axis=1)
                w = np.where(mask, w, 0.0)

            s = w.sum(axis=1, keepdims=True)
            probs[start:end] = np.divide(w, s, out=np.zeros_like(w), where=s > 0)
        return labels_list, probs

    def extend(self, new_records, metric='chamfer',calculator = None,
                            batch_size=2000, strict=True, slack=1.0, 
                            max_proc = None, mp_context = None, reassign=False, **metric_kwargs):
        """Hard-assign records and add them to their Clusters.

        reassign=False : new_records are new; appended to self.records.
        reassign=True  : new_records are already in self.records; those that get a
                        different, valid label are moved out of their old cluster.
        reassign=<lbl> : as True, but only records currently in cluster <lbl>
                        (e.g. -1 for noise) are considered.
        Returns the final label per input record (-1 = unassigned / unchanged noise).
        """

        new_records = list(new_records)
        if not new_records:
            return np.empty(0,dtype=int)
        
        labels = self.assign( new_records, metric=metric, calculator=calculator, 
                                batch_size=batch_size, strict=strict, slack=slack,
                                max_proc=max_proc, mp_context=mp_context, **metric_kwargs)
       
        if not reassign:
            offset = len(self.records)
            self.records.extend(new_records)
            new_global_idx = offset + np.arange(len(new_records))
            for lbl in np.unique(labels):
                local_idx = np.flatnonzero(labels == lbl)
                self._by_label[lbl].extendCluster( [new_records[i] for i in local_idx],
                                    indices=new_global_idx[local_idx], approximated=True)
            return labels

        pos = {id(r): i for i, r in enumerate(self.records)}
        try:
            global_idx = np.array([pos[id(r)] for r in new_records])
        except KeyError:
            raise ValueError("reassign requires new_records to already be in self.records")

        # current label of every record, built from cluster membership
        current = np.full(len(self.records), -1, dtype=int)
        for lbl, cl in self._by_label.items():
            current[np.asarray(cl.indices, dtype=int)] = lbl
        cur = current[global_idx]

        result = cur.copy()                      # default: everything stays put

        # restrict to one source cluster if a label was given
        if reassign is True:
            cand = np.arange(len(new_records))
        else:
            cand = np.flatnonzero(cur == reassign)
        if cand.size == 0:
            return result

        cand_records = [new_records[i] for i in cand]
        new_lbl = self.assign(cand_records, **metric_kwargs)

        # only move records that got a real, different label
        move = (new_lbl != -1) & (new_lbl != cur[cand])
        if not move.any():
            return result

        m_cand = cand[move]                      # positions in new_records
        m_global = global_idx[m_cand]
        m_old = cur[m_cand]
        m_new = new_lbl[move]

        # 1) remove from old clusters
        for old in np.unique(m_old):
            sel = m_old == old
            self._by_label[old].removeFromCluster(m_global[sel])   # see note below

        # 2) add to new clusters
        for lbl in np.unique(m_new):
            sel = np.flatnonzero(m_new == lbl)
            self._by_label[lbl].extendCluster(
                [new_records[m_cand[j]] for j in sel],
                indices=m_global[sel], approximated=True)

        result[m_cand] = m_new
        return result

            # for lbl in np.unique(labels):
            #     if lbl == -1:
            #         continue
            #     local_idx = np.flatnonzero(labels == lbl)
            #     members = [new_records[i] for i in local_idx]
            #     global_idx = np.arange(len(self.records))
            #     indices = global_idx[local_idx]
            #     self._by_label[lbl].extendCluster(members, indices=indices,approximated=True)
                
        # offset= len(self.records)
        # self.records.extend(new_records)
        # new_global_idx = (offset + np.arange(len(new_records)))
        # for lbl in np.unique(labels):
        #     local_idx = np.flatnonzero(labels == lbl)
        #     members = [new_records[i] for i in local_idx]
        #     indices = new_global_idx[local_idx]
        #     self._by_label[lbl].extendCluster(members, indices=indices,approximated=True)
        # return labels

    
    def _bestLabels(self, probs, labels_list, min_probability):
        """Highest-probability label per row if >= min_probability, else -1."""
        lbl_arr = np.asarray(labels_list, dtype=int)
        best_j = probs.argmax(axis=1)
        best_p = probs[np.arange(len(probs)), best_j]
        labels = np.where(best_p >= min_probability, lbl_arr[best_j], -1)
        labels[probs.sum(axis=1) == 0] = -1
        return labels


    def softExtend(self, new_records, metric='chamfer', calculator=None,
                min_probability=0.5, max_distance=None, bandwidth=None, top_k=None,
                batch_size=2000, max_proc=None, mp_context=None, reassign=False,
                **metric_kwargs):
        """
        Soft-assign records and add each to its highest-probability cluster when
        that probability is >= min_probability.

        reassign=False : new_records are new; appended to self.records. Records
                        below the threshold stay in self.records, in no Cluster.
        reassign=True  : new_records are already in self.records; those whose best
                        accepted label differs from their current cluster are moved.
        reassign=<lbl> : as True, but only records currently in cluster <lbl>
                        (e.g. -1 for noise) are considered.

        Returns
        -------
        labels : ndarray   Final label per input record (-1 = unassigned).
        probabilities : ndarray   Soft-assignment matrix, one row per input record.
                        In reassign mode, rows for records that were not candidates
                        are NaN (they were not scored).
        labels_list : list   Column -> cluster label mapping for probabilities.
        """
        new_records = list(new_records)
        if not new_records:
            return (np.empty(0, dtype=int), np.empty((0, len(self._by_label))),
                    list(self._by_label))

        kw = dict(metric=metric, calculator=calculator, max_distance=max_distance,
                bandwidth=bandwidth, top_k=top_k, batch_size=batch_size,
                max_proc=max_proc, mp_context=mp_context, **metric_kwargs)

        # ---------- brand-new records ----------
        if reassign is False or reassign is None:
            labels_list, probs = self.softAssign(new_records, **kw)
            labels = self._bestLabels(probs, labels_list, min_probability)

            offset = len(self.records)
            self.records.extend(new_records)
            gidx = offset + np.arange(len(new_records))
            for lbl in np.unique(labels):
                if lbl == -1:
                    continue
                sel = np.flatnonzero(labels == lbl)
                self._by_label[int(lbl)].extendCluster(
                    [new_records[i] for i in sel], indices=gidx[sel], soft_assigned=True)
            return labels, probs, labels_list

        # ---------- reassignment of existing records ----------
        pos = {id(r): i for i, r in enumerate(self.records)}
        try:
            gidx = np.array([pos[id(r)] for r in new_records])
        except KeyError:
            raise ValueError("reassign requires new_records to already be in self.records")

        # current label per input record from cluster membership (default -1)
        owner = {}
        for cl in self.clusters:
            for r in cl.records:
                owner[id(r)] = cl.label
        cur = np.array([owner.get(id(r), -1) for r in new_records], dtype=int)

        labels = cur.copy()
        labels_list = list(self._by_label)
        probs = np.full((len(new_records), len(labels_list)), np.nan)

        cand = np.arange(len(new_records)) if reassign is True else np.flatnonzero(cur == reassign)
        if cand.size == 0:
            return labels, probs, labels_list

        labels_list, cand_probs = self.softAssign([new_records[i] for i in cand], **kw)
        if cand_probs.shape[1] != probs.shape[1]:
            probs = np.full((len(new_records), cand_probs.shape[1]), np.nan)
        probs[cand] = cand_probs

        new_lbl = self._bestLabels(cand_probs, labels_list, min_probability)
        move = (new_lbl != -1) & (new_lbl != cur[cand])
        if not move.any():
            return labels, probs, labels_list

        m_pos = cand[move]          # positions in new_records
        m_gidx = gidx[m_pos]
        m_old = cur[m_pos]
        m_new = new_lbl[move]

        # 1) remove from old clusters (noise may have no Cluster object)
        for old in np.unique(m_old):
            cl = self.getClusterByLabel(old) if (old != -1 or self._hasLabel(-1)) else None
            if cl is None:
                continue
            sel = np.flatnonzero(m_old == old)
            cl.removeRecords([new_records[m_pos[j]] for j in sel], indices=m_gidx[sel])

        # 2) add to new clusters
        for lbl in np.unique(m_new):
            sel = np.flatnonzero(m_new == lbl)
            self.getClusterByLabel(int(lbl)).extendCluster(
                [new_records[m_pos[j]] for j in sel], indices=m_gidx[sel], soft_assigned=True)

        # 3) keep fitted-subsample labels in sync
        in_fit = m_gidx < self._n_sub
        if in_fit.any():
            self.labels = np.asarray(self.labels).copy()
            self.labels[m_gidx[in_fit]] = m_new[in_fit]

        labels[m_pos] = m_new
        return labels, probs, labels_list
        


    def updateRecordClusterIds(self):
        """
        Rebuild every record's `cluster_ids` from `self.clusters`.

        Each record gets the labels of all clusters that contain it, ordered
        coarse -> fine (by cluster depth). Records in no cluster get []. Because it
        is rebuilt from the clusters, it is safe to call after any change to
        `self.clusters` (new hierarchy, noise removal, ...).
        """
        ids = {id(r): [] for r in self.records}
        ordered = sorted(self.clusters, key=lambda c: (getattr(c, "depth", 0), c.label))
        for cl in ordered:
            for r in cl.records:
                ids[id(r)].append(cl.label)
        for r in self.records:
            r.cluster_ids = ids[id(r)]


    def removeClusters(self, noise_labels, update_labels=True):
        """
        Declare clusters noise: drop them from the clustering and update everything
        that refers to them.

        Parameters
        ----------
        noise_labels : iterable of int
            Labels (Cluster.label) of the clusters to remove.
        update_labels : bool, default=True
            Also write -1 into `self.labels` for the removed clusters' points.
            This assumes Cluster.label values are the values used in self.labels
            (true for flat clusterings). Pass False for hierarchy clusters, whose
            labels live in `self.hierarchy_labels` instead.

        What is updated
        ---------------
        * self.clusters            removed clusters dropped
        * parent / children        surviving children of a removed cluster are
                                re-attached to its nearest surviving ancestor
                                (or become roots); local_label renumbered
        * hierarchy_roots / hierarchy_levels / hierarchy_labels / hierarchy_noise
                                updated if a hierarchy exists
        * self.labels              -1 for removed points (new array, not in place)
        * record.cluster_ids       rebuilt via updateRecordClusterIds()

        Returns
        -------
        removed : list of Cluster
        """
        if np.isscalar(noise_labels):
            noise_labels = [noise_labels]
        noise = {int(x) for x in noise_labels}
        clusters = getattr(self, "clusters", None)
        clusters = [] if clusters is None else list(clusters)
        removed = [c for c in clusters if c.label in noise]
        if not removed:
            return []
        removed_ids = {id(c) for c in removed}
        survivors = [c for c in clusters if id(c) not in removed_ids]

        # Re-attach survivors to their nearest surviving ancestor.
        touched = {}
        for cl in survivors:
            p = cl.parent
            if p is None or id(p) not in removed_ids:
                continue
            while p is not None and id(p) in removed_ids:
                p = p.parent
            cl.parent = p
            if p is not None:
                p.children.append(cl)
                touched[id(p)] = p
        # Parents that lost a child.
        for cl in survivors:
            if any(id(ch) in removed_ids for ch in cl.children):
                cl.children = [ch for ch in cl.children if id(ch) not in removed_ids]
                touched[id(cl)] = cl
        for p in touched.values():
            p.children = sorted(p.children, key=lambda c: c.label)
            for k, ch in enumerate(p.children):
                ch.local_label = k
        rec_noise = []
        for cl in removed:                       # detach the dropped objects
            cl.children = []
            cl.parent = None
            rec_noise.extend(cl.records)
        self.clusters = survivors
        self.clusters.append(Cluster(-1,rec_noise))

        # Hierarchy bookkeeping, if present.
        if getattr(self, "hierarchy_roots", None) is not None:
            self.hierarchy_roots = [c for c in survivors if c.parent is None]
        if getattr(self, "hierarchy_levels", None) is not None:
            self.hierarchy_levels = [[c for c in lvl if id(c) not in removed_ids]
                                    for lvl in self.hierarchy_levels]
        hl = getattr(self, "hierarchy_labels", None)
        if hl is not None and getattr(hl, "size", 0):
            mask = np.isin(hl, list(noise))
            hl = hl.copy()
            hl[mask] = -1
            self.hierarchy_labels = hl
            hn = getattr(self, "hierarchy_noise", None)
            if hn is not None and hn.shape == mask.shape:
                self.hierarchy_noise = hn | mask

        if update_labels and getattr(self, "labels", None) is not None:
            labels = np.asarray(self.labels)
            self.labels = np.where(np.isin(labels, list(noise)), -1, labels)

        self.updateRecordClusterIds()
        return removed


    def agglomerativeClusterQuality(self,labels=None,linkage_matrix=None,min_cluster_size=3,
                min_persistence_gap = None, min_persistence_ratio=2.0, singleton_noise=True,
                condense_size=None, external="cluster", linkage_name='average', min_ratio_floor=1.0):
        """
        Parameters
        ----------
        labels : array (n_samples,)   flat labels, negative = noise.
        linkage_matrix : (n_samples - 1, 4) scipy/sklearn linkage matrix.
        min_cluster_size, min_persistence_ratio : see stability_noise_decision.
        singleton_noise : singletons are always noise.
        condense_size : min points for a sub-branch to count as a real split
            (default min_cluster_size).
        external : "cluster" -> birth = first ancestor merge joining another
            non-noise cluster (noise-only merges ignored); "parent" -> birth =
            merge with the immediate parent node.
        min_ratio_floor : extra guard, noise if ext/int is below this
            (1.0 = disabled; try 1.2).

        Returns DataFrame: cluster, size, node_size, coherence,
            internal_merge_distance, external_merge_distance, persistence_ratio,
            persistence_gap, stability, mean_persistence, stability_lambda,
            stability_ratio, recommended_noise
        """
        import pandas as pd
        _NONE, _MIXED = -1, -2

        if linkage_matrix is None:
            if (linkage_matrix := self.linkage_matrix) is None:
                from scipy.spatial.distance import squareform
                from scipy.cluster.hierarchy import linkage as sp_linkage
                linkage_matrix =  sp_linkage(squareform(np.asarray(self.distance_matrix), checks=False), method=linkage_name)
        if labels is None:
            labels = self.labels

        n = len(labels)
        Z = np.asarray(linkage_matrix, dtype=float)
        if Z.shape != (n - 1, 4):
            raise ValueError("linkage_matrix must have shape (n_samples - 1, 4).")
        if external not in ("cluster", "parent"):
            raise ValueError("external must be 'cluster' or 'parent'.")
        cond = int(condense_size or min_cluster_size)

        D = ClusterCalculator.DendrogramStability(Z, leaf_weight=labels >= 0, condense_size=cond)
        parent, ncl, raw, floor = D.parent, D.ncl, D.raw, D.floor
        left, right = D.left, D.right

        total = 2 * n - 1
        state = np.full(total, _NONE, dtype=int)   # pure label / NONE / MIXED
        for i in range(n):
            if labels[i] >= 0:
                state[i] = labels[i]
        for i in range(n - 1):
            v, a, b = n + i, left[i], right[i]
            sa, sb = state[a], state[b]
            if sa == _NONE:
                state[v] = sb
            elif sb == _NONE:
                state[v] = sa
            elif sa == sb and sa >= 0:
                state[v] = sa
            else:
                state[v] = _MIXED

        # Highest pure node per cluster (the one holding most of its members).
        best = {}
        for v in range(total):
            c = state[v]
            if c < 0:
                continue
            p = parent[v]
            if p != -1 and state[p] == c:
                continue
            if c not in best or ncl[v] > ncl[best[c]]:
                best[c] = v

        threshold = min_cluster_size * np.log(min_persistence_ratio)
        rows = []
        for c in np.unique(labels[labels >= 0]):
            size = int((labels == c).sum())
            v = best.get(c)
            node_size = int(Z[v - n, 3]) if v >= n else 1
            coherence = ncl[v] / size  # <1 => cluster is split across the tree

            internal = D.height(v)
            ext = np.nan
            p = parent[v]
            while p != -1:
                if external == "parent" or state[p] == _MIXED:
                    ext = raw[p - n]
                    break
                p = parent[p]

            if np.isfinite(ext):
                if internal > 0:
                    ratio, gap = ext / internal, ext - internal
                else:
                    ratio, gap = np.inf, ext
                if size == 1:
                    s_log = s_lam = np.inf
                else:
                    s_log, s_lam = D.stability(v, ext)
            else:  # only cluster / nothing to separate from
                ratio = gap = s_log = s_lam = np.nan

            noise = ClusterCalculator.stabilityNoiseDecision(size, s_log, ratio, min_cluster_size,
                                            min_persistence_ratio,
                                            singleton_noise, min_ratio_floor)
            
            if (not noise and min_persistence_gap is not None and np.isfinite(gap) and gap < min_persistence_gap):
                noise = True
            
            rows.append({"cluster": c,
                        "node":int(v),
                        "size": size,
                        "node_size": node_size,
                        "coherence": coherence,
                        "internal_merge_distance": internal,
                        "external_merge_distance": ext,
                        "persistence_ratio": ratio,
                        "persistence_gap": gap,
                        "stability": s_log,
                        "mean_persistence": s_log / size if np.isfinite(s_log) else s_log,
                        "stability_lambda": s_lam,
                        "stability_ratio": s_log / threshold if threshold > 0 else np.nan,
                        "recommended_noise": noise})

        return pd.DataFrame(rows).sort_values("cluster").reset_index(drop=True)


    def assignAggClusteringNoise(self, linkage_matrix=None, min_cluster_size=3,
                min_persistence_ratio=2.0, singleton_noise=True, condense_size=None,
                external="cluster",linkage_name='average', min_ratio_floor=1.0):
        """
        Assign noise points based on agglomerative clustering results.

        Parameters
        ----------
        linkage_matrix : ndarray, optional
            The linkage matrix from agglomerative clustering. If None, it will be computed.
        min_cluster_size : int, optional
            Minimum size of clusters to be considered non-noise.
        min_persistence_ratio : float, optional
            Minimum persistence ratio for small clusters to be retained.
        singleton_noise : bool, optional
            Whether singletons should be considered noise.

        Returns
        -------
        labels : ndarray
            Array of cluster labels with noise points marked as -1.
        """
        if linkage_matrix is None:
            if (linkage_matrix := self.linkage_matrix) is None:
                from scipy.spatial.distance import squareform
                from scipy.cluster.hierarchy import linkage as sp_linkage
                linkage_matrix =  sp_linkage(squareform(np.asarray(self.distance_matrix), checks=False), method=linkage_name)

        quality_df = self.agglomerativeClusterQuality( labels=self.labels, 
                    linkage_matrix=linkage_matrix, min_cluster_size=min_cluster_size,
                    min_persistence_ratio=min_persistence_ratio,
                    singleton_noise=singleton_noise, condense_size=condense_size,
                    external=external, min_ratio_floor=min_ratio_floor)
        
        noise_clusters = quality_df.loc[quality_df["recommended_noise"], "cluster"].values
        labels = self.labels
        labels[np.isin(labels, noise_clusters)] = -1
        if len(noise_clusters):
            self.removeClusters(noise_clusters)  
        
        return labels

    def buildAgglomerativeHierarchy(self, linkage_matrix=None, cut_heights=None, n_levels=4,
                                    min_cluster_size=3, singleton_noise=True,
                                    min_persistence_ratio=2.0, min_persistence_gap=None,
                                    condense_size=None, min_ratio_floor=1.0,
                                    external="cluster", linkage_name='average',
                                    identify_stable_branches=False,
                                    stable_min_size=None, stable_min_ratio=2.0,
                                    stable_min_gap=None, update_record_ids=False):
        """
        Build a hierarchy of FLAT CUTS through an agglomerative dendrogram.

        A cluster at level k is attached to the non-noise cluster at level k-1
        containing its members. Noise is evaluated independently at every cut by
        agglomerativeClusterQuality (the same decision used by
        assignAggClusteringNoise), so a group may be noise at a fine cut but
        belong to a valid cluster at a coarser cut.

        Parameters
        ----------
        linkage_matrix : ndarray (n-1, 4), optional
        cut_heights : sequence of float, optional
            Coarse -> fine order is enforced. If None, `n_levels` heights are
            taken from merge-distance quantiles.
        n_levels : int, default=4
        min_cluster_size, min_persistence_ratio, singleton_noise, condense_size,
        external, min_ratio_floor, min_persistence_gap :
            Passed straight to agglomerativeClusterQuality (see there).
            min_persistence_ratio must not be None.
        identify_stable_branches : bool
            Independently find dendrogram branches with
                size >= stable_min_size and
                stability >= stable_min_size * log(stable_min_ratio)
            (and gap >= stable_min_gap if given) -> self.agglomerative_stable_branches.
        stable_min_size : int or None  (default min_cluster_size)
        stable_min_ratio : float, default=2.0
        stable_min_gap : float or None
        update_record_ids : bool
            Replace record.cluster_ids with the non-noise cluster path.

        Returns
        -------
        clusters : list of Cluster (non-noise, coarse -> fine)

        Each Cluster carries in `.details`: stability, mean_persistence,
        stability_lambda, persistence ratio/gap, dendrogram_node, etc.
        `.persistence` is the ratio; `.stability` mirrors details["stability"].
        """
        from scipy.cluster.hierarchy import fcluster

        if min_persistence_ratio is None:
            raise ValueError("min_persistence_ratio is required.")

        if linkage_matrix is None:
            if (linkage_matrix := self.linkage_matrix) is None:
                from scipy.spatial.distance import squareform
                from scipy.cluster.hierarchy import linkage as sp_linkage
                linkage_matrix =  sp_linkage(squareform(np.asarray(self.distance_matrix), checks=False), method=linkage_name)

        Z = np.asarray(linkage_matrix, dtype=float)
        n = len(self.records)

        if n == 0:
            self.clusters = []
            self.hierarchy_roots = []
            self.hierarchy_levels = []
            self.hierarchy_labels = np.empty((0, 0), dtype=int)
            self.hierarchy_noise = np.empty((0, 0), dtype=bool)
            self.hierarchy_cut_heights = np.array([])
            self.agglomerative_stable_branches = []
            return []

        if n == 1:
            if Z.size != 0:
                raise ValueError("linkage_matrix must be empty for one sample.")
        elif Z.shape != (n - 1, 4):
            raise ValueError(f"linkage_matrix must have shape ({n - 1}, 4), got {Z.shape}")

        cond = int(condense_size or min_cluster_size)
        D = ClusterCalculator.DendrogramStability(Z, condense_size=cond) if n > 1 else None

        def nodeDistance(node):
            return 0.0 if node < n else float(D.raw[node - n])

        def nodeSize(node):
            return 1 if node < n else int(Z[node - n, 3])

        def parentOf(node):
            p = D.parent[node] if D is not None else -1
            return None if p == -1 else int(p)

        def leavesOf(node):
            out, stack = [], [node]
            while stack:
                u = stack.pop()
                if u < n:
                    out.append(u)
                else:
                    stack.append(int(D.left[u - n]))
                    stack.append(int(D.right[u - n]))
            return np.sort(np.asarray(out, dtype=int))

        def branchStats(node):
            """A branch exists for internal <= h < external."""
            internal = nodeDistance(node)
            p = parentOf(node)
            ext = np.nan if p is None else nodeDistance(p)
            if np.isfinite(ext):
                gap = ext - internal
                if internal > 0:
                    ratio = ext / internal
                elif ext > 0:
                    ratio = np.inf
                else:
                    ratio = np.nan
            else:
                gap = ratio = np.nan
            size = nodeSize(node)
            if np.isfinite(ext) and size > 1:
                s_log, s_lam = D.stability(node, ext)
            elif np.isfinite(ext):
                s_log = s_lam = np.inf
            else:
                s_log = s_lam = np.nan
            return {"node": int(node), "size": int(size),
                    "internal_distance": float(internal),
                    "external_distance": float(ext),
                    "persistence_gap": float(gap),
                    "persistence_ratio": float(ratio),
                    "stability": float(s_log),
                    "mean_persistence": float(s_log / size) if np.isfinite(s_log) else float(s_log),
                    "stability_lambda": float(s_lam)}

        if n == 1:
            heights = np.array([0.0])
        elif cut_heights is not None:
            heights = np.asarray(cut_heights, dtype=float)
            if heights.ndim != 1:
                raise ValueError("cut_heights must be one-dimensional.")
            if len(heights) == 0:
                raise ValueError("cut_heights cannot be empty.")
            if not np.all(np.isfinite(heights)):
                raise ValueError("cut_heights contains non-finite values.")
            if np.any(heights < 0):
                raise ValueError("cut_heights must be >= 0.")
            heights = np.unique(heights)[::-1]                     # coarse -> fine
        else:
            if n_levels < 1:
                raise ValueError("n_levels must be >= 1.")
            merge_distances = np.unique(Z[:, 2])
            if len(merge_distances) <= n_levels:
                heights = merge_distances[::-1]
            else:
                qs = np.linspace(0.85, 0.15, n_levels)
                heights = np.unique(np.quantile(merge_distances, qs))[::-1]

        all_clusters = []
        hierarchy_levels = []
        hierarchy_roots = []
        hierarchy_labels = np.full((len(heights), n), -1, dtype=int)
        hierarchy_noise = np.zeros((len(heights), n), dtype=bool)

        global_counter = 0
        previous_by_label = {}

        for depth, h in enumerate(heights):
            flat = np.ones(1, dtype=int) if n == 1 else fcluster(Z, t=float(h), criterion="distance")

            # Same convention as assignAggClusteringNoise: noise points are -1
            cut_labels = flat.copy()
            if singleton_noise:
                cut_labels[np.bincount(flat)[flat] == 1] = -1

            if n > 1 and (cut_labels >= 0).any():
                q = self.agglomerativeClusterQuality(
                        labels=cut_labels, linkage_matrix=Z,
                        min_cluster_size=min_cluster_size,
                        min_persistence_ratio=min_persistence_ratio,
                        singleton_noise=singleton_noise, condense_size=condense_size,
                        external=external, min_ratio_floor=min_ratio_floor,
                        min_persistence_gap=min_persistence_gap).set_index("cluster")
            else:
                q = None

            order = np.argsort(flat, kind="stable")
            uniq, starts = np.unique(flat[order], return_index=True)
            groups = np.split(order, starts[1:])

            current_level, current_by_label = [], {}
            for local_label, idx in zip(uniq, groups):
                if q is not None:
                    if local_label not in q.index or q.at[local_label, "recommended_noise"]:
                        hierarchy_noise[depth, idx] = True
                        continue
                    row = q.loc[local_label]
                    stats = {"node": int(row["node"]), "size": int(row["size"]),
                             "internal_distance": float(row["internal_merge_distance"]),
                             "external_distance": float(row["external_merge_distance"]),
                             "persistence_gap": float(row["persistence_gap"]),
                             "persistence_ratio": float(row["persistence_ratio"]),
                             "stability": float(row["stability"]),
                             "mean_persistence": float(row["mean_persistence"]),
                             "stability_lambda": float(row["stability_lambda"])}
                elif n == 1 and not singleton_noise:
                    # trivial single-point dendrogram: nothing to separate from
                    nan = float("nan")
                    stats = {"node": 0, "size": 1, "internal_distance": 0.0,
                             "external_distance": nan, "persistence_gap": nan,
                             "persistence_ratio": nan, "stability": nan,
                             "mean_persistence": nan, "stability_lambda": nan}
                else:
                    hierarchy_noise[depth, idx] = True
                    continue

                # Parent = cluster of the PREVIOUS (coarser) cut holding these points.
                parent = None
                if depth > 0:
                    prev_label = hierarchy_labels[depth - 1, idx[0]]
                    if prev_label >= 0:
                        parent = previous_by_label.get(int(prev_label))

                cl = Cluster(label=global_counter,
                             records=[self.records[i] for i in idx],
                             parent=parent,
                             depth=depth)
                global_counter += 1
                cl.local_label = len(current_level)
                cl.details["cut_height"] = float(h)
                cl.details["cut_level"] = depth
                cl.details["flat_cut_label"] = int(local_label)
                cl.details["dendrogram_node"] = stats["node"]
                cl.details["point_indices"] = idx
                cl.details.update(stats)
                cl.persistence = stats["persistence_ratio"]
                cl.stability = stats["stability"]

                if parent is None:
                    hierarchy_roots.append(cl)
                else:
                    cl.local_label = len(parent.children)
                    parent.children.append(cl)

                current_level.append(cl)
                current_by_label[cl.label] = cl
                all_clusters.append(cl)
                hierarchy_labels[depth, idx] = cl.label

            hierarchy_levels.append(current_level)
            previous_by_label = current_by_label

        stable_branches = []
        if identify_stable_branches and n > 1:
            if stable_min_size is None:
                stable_min_size = min_cluster_size
            stable_threshold = stable_min_size * np.log(stable_min_ratio)
            for node in range(n, 2 * n - 1):
                stats = branchStats(node)
                if stats["size"] < stable_min_size:
                    continue
                strong = (np.isfinite(stats["stability"])
                          and stats["stability"] >= stable_threshold)
                gap_ok = (stable_min_gap is None
                          or (np.isfinite(stats["persistence_gap"])
                              and stats["persistence_gap"] >= stable_min_gap))
                if not (strong and gap_ok):
                    continue
                idx = leavesOf(node)
                stable_branches.append({**stats,
                                        "point_indices": idx,
                                        "records": [self.records[i] for i in idx],
                                        "height_min": stats["internal_distance"],
                                        "height_max": stats["external_distance"]})
            stable_branches.sort(key=lambda x: (x["stability"], x["persistence_gap"], x["size"]),
                                 reverse=True)

        if update_record_ids:
            for i, record in enumerate(self.records):
                col = hierarchy_labels[:, i]
                record.cluster_ids = [int(c) for c in col[col >= 0]]

        self.hierarchy_roots = hierarchy_roots
        self.hierarchy_levels = hierarchy_levels
        self.hierarchy_labels = hierarchy_labels
        self.hierarchy_noise = hierarchy_noise
        self.hierarchy_cut_heights = heights
        self.agglomerative_stable_branches = stable_branches
        self.levels = list(heights)

        return all_clusters
        
    def subcluster(self, cluster, distance_matrix, method="agglomerative",
                   linkage_method="average", n_clusters=None,
                   distance_threshold=None, update_record_ids=True, **kwargs):
        """
        Subdivide `cluster` using a new distance matrix, via clusterChannels.

        Parameters
        ----------
        cluster : Cluster
            Cluster (top level or already a child) whose records belong to self.
        distance_matrix : ndarray (cluster.size, cluster.size)
            Ordered like cluster.records.
        method : {'agglomerative', 'HDBSCAN', 'DBSCAN'}
        linkage_method : str
            Agglomerative only; used unless `linkage` is passed in kwargs.
        n_clusters, distance_threshold :
            Agglomerative only, mutually exclusive, exactly one required.
        update_record_ids : bool
            Set record.cluster_ids to [*ancestor labels, cluster.label, child.label].
            Noise points get -1 as their last entry.
        **kwargs
            Forwarded to the estimator (eps, min_samples, min_cluster_size, ...).

        Returns
        -------
        children : list of Cluster (includes a label -1 noise cluster if any)
        sub_clustering : ClusteringRecords for this cluster's records
            (linkage_matrix is set for agglomerative, else None)
        labels : ndarray, local labels aligned with cluster.records

        Side effects
        ------------
        cluster.children / .distance_matrix / .linkage_matrix / .model / .details
        child.parent / .depth / .local_label / .details
        self.sub_clusterings[tuple(path_to_cluster)] = sub_clustering
        self.levels gains cluster.depth + 1
        record.cluster_ids (if update_record_ids)
        """
        if not isinstance(cluster, Cluster):
            raise TypeError("cluster must be a Cluster instance.")
        if method not in {"agglomerative", "HDBSCAN", "DBSCAN"}:
            raise ValueError(f"Unknown method '{method}'.")

        records = cluster.records
        n = len(records)
        if n < 2:
            raise ValueError("A cluster must contain at least two records to be subclustered.")

        owned = {id(r) for r in self.records}
        if not all(id(r) in owned for r in records):
            raise ValueError("cluster does not belong to this ClusteringRecords.")

        distance_matrix = np.asarray(distance_matrix, dtype=float)
        if distance_matrix.shape != (n, n):
            raise ValueError(f"distance_matrix must have shape ({n}, {n}) for this "
                             f"cluster, got {distance_matrix.shape}.")
        if not np.all(np.isfinite(distance_matrix)):
            raise ValueError("distance_matrix contains non-finite values.")
        if not np.allclose(distance_matrix, distance_matrix.T, atol=1e-8):
            raise ValueError("distance_matrix must be symmetric.")
        if np.any(np.diag(distance_matrix) != 0):
            raise ValueError("distance_matrix must have a zero diagonal.")

        fit_kwargs = dict(kwargs)
        if method == "agglomerative":
            if n_clusters is not None and distance_threshold is not None:
                raise ValueError("Specify either n_clusters or distance_threshold, not both.")
            if n_clusters is None and distance_threshold is None:
                raise ValueError("Specify either n_clusters or distance_threshold.")
            fit_kwargs.setdefault("linkage", linkage_method)
            if n_clusters is not None:
                fit_kwargs["n_clusters"] = n_clusters
                fit_kwargs["distance_threshold"] = distance_threshold
        elif n_clusters is not None or distance_threshold is not None:
            raise ValueError("n_clusters / distance_threshold only apply to 'agglomerative'; "
                             "use eps / min_samples / min_cluster_size for DBSCAN/HDBSCAN.")

        # Nothing on the records or cluster is touched until the fit succeeds.
        labels, sub = clusterChannels(distance_matrix, method, records,
                                      update_record_ids=False, **fit_kwargs)

        path = self._clusterPath(cluster)
        children = list(sub.clusters)
        for child in children:
            child.parent = cluster
            child.depth = cluster.depth + 1
            child.local_label = int(child.label)
            child.details["parent_label"] = cluster.label
            child.details["local_label"] = int(child.label)
            child.details["member_indices"] = np.flatnonzero(labels == child.label)

        n_noise = int(np.sum(labels == -1))
        cluster.children = children
        cluster.distance_matrix = distance_matrix
        cluster.linkage_matrix = sub.linkage_matrix
        cluster.model = sub.model
        cluster.details["subclustering_method"] = method
        cluster.details["subclustering_params"] = dict(fit_kwargs)
        cluster.details["subclustering_n_clusters"] = len(children) - (1 if n_noise else 0)
        cluster.details["subclustering_n_noise"] = n_noise
        cluster.details["subclustering_distance_threshold"] = distance_threshold

        sub.levels = [cluster.depth + 1]
        sub.parent_cluster = cluster
        if self.sub_clusterings is None:
            self.sub_clusterings = {}
        self.sub_clusterings[tuple(path)] = sub
        self.levels = sorted(set(self.levels) | {cluster.depth + 1})

        if update_record_ids:
            for child in children:
                child_path = path + [int(child.label)]
                for record in child.records:
                    record.cluster_ids = list(child_path)

        return children, sub, labels
    
    
#: Colour of each rank - the 0-based position of an object among those a viewer
#: script draws. One table for every PyMOL, VMD and ChimeraX script this module
#: writes, so that a channel is the same colour in all of them; it is
#: CaviTracerMD's viewer palette, unchanged.
#:
#: Ranks 0-5 are CAVER 3's first six colours, from its out/pymol/modules/rgb.py
#: in the order its view.py hands them to tunnel clusters. They are what makes a
#: CAVER figure recognisable, so they are kept verbatim.
_CAVER_PRIMARIES = [(0.0, 0.0, 1.0),    # blue
                    (0.0, 1.0, 0.0),    # green
                    (1.0, 0.0, 0.0),    # red
                    (0.0, 1.0, 1.0),    # cyan
                    (1.0, 1.0, 0.0),    # yellow
                    (1.0, 0.0, 1.0)]    # magenta

# Ranks 6-199 are a table, not generated: the rule below searches many thousands
# of candidates, which a viewer should not be doing before it draws.
#
# Taken farthest-first: each rank is the candidate farthest from every colour
# before it and from the reserved greys, so the closest pair among the first n
# colours falls as slowly as the candidates allow and the worst pairs come last.
# Distance is OKLab x100 between versions of each colour dimmed in linear light
# to ten levels from 0.55 to 1, the closest pair of versions counting: a sphere
# runs from lit to shadowed, so two colours are distinct only if no shade of one
# matches a shade of the other. Candidates are HSV at every degree of hue, with
# saturation 0.45-1 and value 0.60-1 in five steps each, less those within 7 of
# either reserved grey - 0.45, and the grey80 a protein is drawn in. The value
# floor keeps out colours that turn near-black in shadow, which the distance
# alone would rank highly for being far from everything light. Regenerate on
# those terms or not at all.
_PALETTE = _CAVER_PRIMARIES + [
    (0.430, 0.000, 0.600),  # 6
    (0.000, 0.533, 1.000),  # 7
    (0.782, 0.550, 1.000),  # 8
    (1.000, 0.580, 0.100),  # 9
    (0.000, 0.600, 0.260),  # 10
    (0.700, 0.000, 0.432),  # 11
    (1.000, 0.550, 0.670),  # 12
    (0.150, 0.255, 0.600),  # 13
    (0.600, 0.300, 0.000),  # 14
    (0.550, 1.000, 0.617),  # 15
    (0.400, 0.100, 1.000),  # 16
    (0.600, 0.560, 0.000),  # 17
    (0.550, 0.715, 1.000),  # 18
    (0.688, 0.250, 1.000),  # 19
    (0.479, 0.330, 0.600),  # 20
    (1.000, 0.400, 0.850),  # 21
    (0.000, 0.642, 0.700),  # 22
    (0.600, 0.000, 0.100),  # 23
    (1.000, 0.880, 0.550),  # 24
    (0.470, 0.400, 1.000),  # 25
    (0.230, 0.000, 0.600),  # 26
    (0.600, 0.330, 0.406),  # 27
    (1.000, 0.250, 0.462),  # 28
    (0.240, 0.414, 0.600),  # 29
    (0.683, 1.000, 0.000),  # 30
    (1.000, 0.400, 0.250),  # 31
    (1.000, 0.670, 0.550),  # 32
    (0.000, 0.333, 1.000),  # 33
    (0.000, 0.000, 0.600),  # 34
    (0.330, 0.600, 0.456),  # 35
    (0.870, 0.400, 1.000),  # 36
    (0.600, 0.150, 0.585),  # 37
    (0.550, 0.887, 1.000),  # 38
    (1.000, 0.100, 0.760),  # 39
    (0.000, 0.700, 1.000),  # 40
    (1.000, 0.767, 0.000),  # 41
    (0.550, 0.565, 1.000),  # 42
    (0.785, 0.900, 0.495),  # 43
    (0.000, 0.900, 0.450),  # 44
    (0.550, 1.000, 0.843),  # 45
    (0.600, 0.456, 0.240),  # 46
    (0.567, 0.000, 1.000),  # 47
    (1.000, 0.550, 0.865),  # 48
    (1.000, 0.400, 0.660),  # 49
    (0.490, 0.700, 0.175),  # 50
    (0.600, 0.240, 0.474),  # 51
    (0.336, 0.240, 0.600),  # 52
    (0.817, 0.000, 1.000),  # 53
    (0.000, 1.000, 0.767),  # 54
    (0.650, 0.400, 1.000),  # 55
    (0.250, 0.100, 1.000),  # 56
    (0.600, 0.150, 0.292),  # 57
    (0.250, 0.425, 1.000),  # 58
    (0.700, 0.448, 0.385),  # 59
    (0.700, 0.245, 0.175),  # 60
    (1.000, 0.400, 0.460),  # 61
    (1.000, 0.720, 0.400),  # 62
    (0.581, 0.280, 0.700),  # 63
    (0.400, 0.200, 0.800),  # 64
    (1.000, 0.250, 0.662),  # 65
    (0.505, 0.600, 0.330),  # 66
    (0.330, 0.339, 0.600),  # 67
    (0.955, 0.550, 1.000),  # 68
    (0.000, 0.170, 0.600),  # 69
    (0.760, 1.000, 0.400),  # 70
    (0.400, 0.620, 1.000),  # 71
    (0.320, 0.632, 0.800),  # 72
    (1.000, 0.940, 0.400),  # 73
    (1.000, 0.250, 0.912),  # 74
    (0.400, 0.500, 1.000),  # 75
    (1.000, 0.000, 0.517),  # 76
    (1.000, 0.000, 0.300),  # 77
    (0.000, 0.700, 0.478),  # 78
    (0.385, 0.700, 0.674),  # 79
    (0.000, 0.340, 0.600),  # 80
    (0.250, 0.200, 0.800),  # 81
    (1.000, 0.520, 0.400),  # 82
    (0.887, 0.250, 1.000),  # 83
    (0.550, 1.000, 0.985),  # 84
    (0.000, 0.700, 0.618),  # 85
    (0.390, 0.700, 0.385),  # 86
    (0.394, 0.900, 0.225),  # 87
    (1.000, 0.467, 0.000),  # 88
    (0.700, 0.385, 0.679),  # 89
    (0.330, 0.000, 0.600),  # 90
    (0.250, 0.887, 1.000),  # 91
    (0.683, 0.000, 1.000),  # 92
    (1.000, 0.250, 0.275),  # 93
    (0.700, 0.000, 0.548),  # 94
    (0.446, 0.175, 0.700),  # 95
    (0.000, 0.700, 0.000),  # 96
    (0.700, 0.280, 0.329),  # 97
    (0.700, 0.350, 0.280),  # 98
    (0.569, 0.175, 0.700),  # 99
    (1.000, 0.400, 0.980),  # 100
    (0.580, 1.000, 0.400),  # 101
    (0.000, 0.195, 0.900),  # 102
    (1.000, 0.267, 0.000),  # 103
    (0.163, 0.000, 0.700),  # 104
    (0.850, 1.000, 0.000),  # 105
    (0.700, 0.537, 0.000),  # 106
    (1.000, 0.550, 0.550),  # 107
    (0.250, 1.000, 0.875),  # 108
    (0.250, 1.000, 0.637),  # 109
    (0.550, 0.000, 0.600),  # 110
    (0.700, 0.280, 0.665),  # 111
    (0.700, 0.647, 0.385),  # 112
    (0.000, 0.617, 1.000),  # 113
    (0.000, 0.460, 0.600),  # 114
    (0.700, 0.280, 0.441),  # 115
    (1.000, 0.667, 0.000),  # 116
    (0.294, 0.600, 0.060),  # 117
    (0.000, 0.800, 1.000),  # 118
    (0.700, 0.443, 0.000),  # 119
    (1.000, 0.000, 0.633),  # 120
    (1.000, 0.883, 0.000),  # 121
    (0.700, 0.175, 0.525),  # 122
    (0.662, 0.550, 1.000),  # 123
    (0.550, 0.640, 1.000),  # 124
    (1.000, 0.610, 0.400),  # 125
    (0.700, 0.385, 0.574),  # 126
    (0.483, 0.000, 1.000),  # 127
    (0.280, 0.350, 0.700),  # 128
    (0.700, 0.000, 0.257),  # 129
    (1.000, 0.000, 0.883),  # 130
    (1.000, 0.752, 0.550),  # 131
    (0.700, 0.175, 0.236),  # 132
    (0.630, 0.700, 0.000),  # 133
    (1.000, 0.000, 0.417),  # 134
    (0.940, 1.000, 0.400),  # 135
    (1.000, 0.250, 0.562),  # 136
    (0.250, 0.337, 1.000),  # 137
    (0.294, 0.280, 0.700),  # 138
    (0.700, 1.000, 0.550),  # 139
    (0.495, 0.900, 0.657),  # 140
    (1.000, 0.550, 0.767),  # 141
    (0.400, 0.690, 1.000),  # 142
    (1.000, 0.800, 0.400),  # 143
    (0.700, 0.140, 0.000),  # 144
    (0.788, 0.250, 1.000),  # 145
    (0.000, 1.000, 0.350),  # 146
    (0.700, 0.455, 0.280),  # 147
    (1.000, 0.400, 0.560),  # 148
    (0.280, 0.413, 0.700),  # 149
    (0.700, 0.175, 0.438),  # 150
    (0.413, 0.250, 1.000),  # 151
    (0.483, 0.280, 0.700),  # 152
    (0.600, 0.250, 1.000),  # 153
    (0.750, 0.400, 1.000),  # 154
    (1.000, 0.400, 0.750),  # 155
    (0.000, 0.467, 1.000),  # 156
    (0.600, 0.230, 0.000),  # 157
    (1.000, 0.970, 0.550),  # 158
    (0.700, 0.665, 0.280),  # 159
    (1.000, 0.250, 0.750),  # 160
    (0.988, 0.250, 1.000),  # 161
    (0.070, 0.469, 0.700),  # 162
    (0.700, 0.000, 0.350),  # 163
    (0.900, 0.484, 0.225),  # 164
    (0.385, 0.458, 0.700),  # 165
    (0.400, 1.000, 0.530),  # 166
    (0.175, 0.236, 0.700),  # 167
    (1.000, 0.440, 0.400),  # 168
    (0.550, 0.782, 1.000),  # 169
    (0.700, 0.000, 0.665),  # 170
    (0.350, 0.700, 0.280),  # 171
    (0.000, 0.360, 0.900),  # 172
    (0.333, 0.000, 1.000),  # 173
    (0.560, 0.400, 1.000),  # 174
    (0.400, 1.000, 0.940),  # 175
    (0.560, 0.700, 0.280),  # 176
    (0.895, 0.100, 1.000),  # 177
    (0.700, 0.385, 0.385),  # 178
    (0.248, 0.150, 0.600),  # 179
    (0.225, 0.439, 0.900),  # 180
    (0.280, 0.616, 0.700),  # 181
    (0.400, 0.560, 1.000),  # 182
    (0.400, 1.000, 0.790),  # 183
    (0.250, 0.575, 1.000),  # 184
    (0.469, 0.385, 0.700),  # 185
    (0.400, 1.000, 0.700),  # 186
    (1.000, 0.325, 0.250),  # 187
    (0.400, 0.440, 1.000),  # 188
    (0.000, 0.230, 0.600),  # 189
    (1.000, 0.367, 0.000),  # 190
    (0.145, 0.100, 1.000),  # 191
    (1.000, 0.250, 0.375),  # 192
    (0.400, 0.960, 1.000),  # 193
    (1.000, 0.738, 0.250),  # 194
    (1.000, 0.610, 0.550),  # 195
    (1.000, 0.650, 0.250),  # 196
    (0.060, 0.096, 0.600),  # 197
    (1.000, 0.250, 0.825),  # 198
    (0.800, 0.370, 0.200),  # 199
]


def _paletteIndex(rank):
    """Row of :data:`_PALETTE` holding the colour of the 0-based *rank*.

    Past rank 199 the ranks cycle through 6-199 and the primaries never repeat.
    That far out no palette keeps every pair apart, so colour stops identifying
    an object there and its name has to."""

    first = len(_CAVER_PRIMARIES)
    if rank < len(_PALETTE):
        return rank
    return first + (rank - first) % (len(_PALETTE) - first)


def _hexColour(rank):
    """*rank*'s colour as ``#rrggbb``, the form a ChimeraX script gives it in."""

    return '#%02x%02x%02x' % tuple(int(round(255 * value))
                                   for value in _PALETTE[_paletteIndex(rank)])


#: The VMD colour IDs ranks 0-29 take, in order: every fixed ID from 0 to 32 but
#: 2 (gray), 8 (white) and 16 (black), which the protein, the background and the
#: labels keep. Their RGB is redefined to the palette's, and the first six were
#: already VMD's nearest to the CAVER primaries. Later ranks take the IDs from 33
#: on, VMD's colour-scale slots: the whole palette fits there, but VMD regenerates
#: them whenever its colour scale is changed, which it never does to a fixed ID.
_VMD_COLOR_IDS = (0, 7, 1, 10, 4, 11, 3, 9, 12, 13, 14, 15, 5, 6, 17, 18, 19,
                  20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)


def _vmdColorID(rank):
    """The VMD colour ID carrying *rank*'s colour (see :data:`_VMD_COLOR_IDS`)."""

    index = _paletteIndex(rank)
    if index < len(_VMD_COLOR_IDS):
        return _VMD_COLOR_IDS[index]
    return 33 + index - len(_VMD_COLOR_IDS)


def _vmdPalette(count):
    """Tcl giving ranks ``0 .. count-1`` their palette colours, and listing the
    IDs that carry them, rank by rank, as ``cavitracer_colors``."""

    lines = ['# The CaviTracer palette, the same as in the PyMOL and ChimeraX',
             '# scripts. Ranks 0-29 redefine fixed colour IDs; later ranks use the',
             '# colour-scale IDs from 33 on, which VMD regenerates if the colour',
             '# scale is changed.']
    ids, defined = [], set()
    for rank in range(count):
        color_id = _vmdColorID(rank)
        if color_id not in defined:
            defined.add(color_id)
            lines.append('color change rgb %d %.3f %.3f %.3f'
                         % ((color_id,) + tuple(_PALETTE[_paletteIndex(rank)])))
        ids.append(str(color_id))
    lines.append('set cavitracer_colors {%s}' % ' '.join(ids))
    return '\n'.join(lines)


#: The palette as a PyMOL script carries it: the table, and ``caverColour`` to
#: register a rank's colour on first use.
_VIS_PALETTE = r'''# --- Palette ---
# The CaviTracer palette, the same as in the VMD and ChimeraX scripts, so that a
# channel is the same colour in all three. Ranks 0-5 are CAVER 3's first six
# colours; ranks 6-199 were taken farthest-first for how they read on shaded
# spheres. Past rank 199 the ranks cycle through 6-199 and the primaries never
# repeat.
PALETTE = [
__ROWS__]
CAVER_PRIMARIES = PALETTE[:6]

def caverColour(rank):
    """Name of the colour for a 0-based rank, registered on first use.

    A lookup and nothing else: a rank is the same colour in every structure,
    every run and every viewer, whatever was loaded beside it."""
    first = len(CAVER_PRIMARIES)
    if rank >= len(PALETTE):
        rank = first + (rank - first) % (len(PALETTE) - first)
    name = "caver%d" % (rank + 1) if rank < first else "gen%d" % rank
    cmd.set_color(name, list(PALETTE[rank]))
    return name
# --- End of palette ---
'''.replace('__ROWS__', ''.join('    (%.3f, %.3f, %.3f),  # %d\n' % (rgb + (rank,))
                                for rank, rgb in enumerate(_PALETTE)))


#: Source of the PyMOL viewer that :func:`_writeVisScript` leaves beside what
#: :func:`writePyMolCaviTracerScript` writes. Held inline so that this module
#: carries everything it writes, and raw so the rank patterns keep their
#: backslashes.
#: Source of the PyMOL viewer that :func:`_writeVisScript` leaves beside what
#: :func:`writePyMolCaviTracerScript` writes. Held inline so that this module
#: carries everything it writes, and raw so the rank patterns keep their
#: backslashes.
_VIS_CHANNELS_SCRIPT = r'''import glob
import os
import re
import shlex
import sys

# --- Parse command-line args ---
# Invoke as:  pymol vis_channels.py -- protein.pdb "por*chl*.pqr"
#         or: pymol vis_channels.py -- protein.pdb channels.cif
#         or: pymol vis_channels.py -- protein.pdb channels.pqr  (every channel)
#         or: pymol vis_channels.py -- protein.pdb frames.pqr    (a state per frame)
#         or: pymol vis_channels.py -- channels.cif        (structure inside)
# The regex MUST be quoted so the shell doesn't glob-expand it before PyMOL sees it.
#
# One script for both formats rather than one each. The two write the same
# spheres and are coloured off the same 0-based rank, so a channel is the same
# colour whichever way the run was written, and a directory holding both opens
# with either in view.
def holdsChannels(path):
    """Whether *path* is a PQR of channels rather than a structure: every
    sphere this module writes is a FIL residue, which no structure holds."""
    with open(path) as handle:
        return any(line.startswith(("ATOM", "HETATM")) and line[17:20] == "FIL"
                   for line in handle)

protein_file = None
cif_file = None
pqr_files = []
channel_regex = None
outfile=None
for arg in sys.argv[1:]:
    # Without a "--" PyMOL leaves its own flags and this script in argv, and the
    # script is an existing file, so unguarded it gets loaded as the protein.
    # That only bites when nothing else is passed, which is exactly the
    # "just find the mmCIF yourself" invocation.
    if arg == "--" or arg.startswith("-") or arg.lower().endswith(".py"):
        continue
    if os.path.isfile(arg):
        if arg.lower().endswith(".cif"):
            if cif_file is None:
                cif_file = arg
        # A PQR named outright is a file of channels, split below, unless it
        # holds none: a structure can be a PQR as well.
        elif arg.lower().endswith(".pqr") and holdsChannels(arg):
            pqr_files.append(arg)
        elif protein_file is None:
            protein_file = argi
    elif arg.lower().endswith(".pse"):
        outfile=arg
    elif channel_regex is None:
        channel_regex = arg

if channel_regex is None:
    channel_regex = "*chl*.pqr"   # fallback default

# Nothing named and no PQRs about: an mmCIF run leaves a single file, so look
# for one before giving up.
if cif_file is None and not pqr_files and not glob.glob(channel_regex):
    found = sorted(glob.glob("*.cif"))
    if found:
        cif_file = found[0]

print(f"Using channel regex: {channel_regex}")
if cif_file:
    print(f"Using mmCIF: {cif_file}")
for path in pqr_files:
    print(f"Using PQR: {path}")

''' + _VIS_PALETTE + r'''
def cifLoops(path):
    """category -> (columns, rows-as-token-lists). Enough CIF for what we write.

    Deliberately not a general CIF parser: this reads back files this module
    wrote, whose loops are one row per line with no multi-line values. shlex
    does the splitting so a quoted value ('ProDy 2.6.1') stays one token.
    """
    loops = {}
    with open(path) as handle:
        lines = handle.read().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        i += 1
        columns, category = [], None
        while i < len(lines) and lines[i].startswith("_"):
            category, _, column = lines[i].strip()[1:].partition(".")
            columns.append(column)
            i += 1
        rows = []
        while i < len(lines) and lines[i] and not lines[i][0] in "#_" \
                and not lines[i].startswith(("loop_", "data_")):
            try:
                tokens = shlex.split(lines[i])
            except ValueError:
                tokens = lines[i].split()
            if len(tokens) == len(columns):
                rows.append(tokens)
            i += 1
        loops[category] = (columns, rows)
    return loops

# The protein: a file named on the command line, or the _atom_site the mmCIF
# carries when it was written with one. Loading the .cif for its structure is
# safe either way -- PyMOL simply finds no atoms when there is none -- but an
# empty object beside the channels is noise, so it is checked for first.
structure_source = protein_file
if structure_source is None and cif_file:
    with open(cif_file) as handle:
        if "_atom_site." in handle.read(200000):
            structure_source = cif_file

if structure_source:
    protein_name = os.path.splitext(os.path.basename(structure_source))[0]
    cmd.load(structure_source, protein_name)
    cmd.hide("everything", protein_name)
    cmd.show("cartoon", protein_name)
    cmd.show("surface", protein_name)
    cmd.color("grey80", protein_name)
    cmd.set("transparency", 0.5, protein_name)
    print(f"Loaded protein: {protein_name}")
else:
    print("No protein file found in args (use: pymol vis_channels.py -- your.pdb)")

# --- Load channels/tunnels ---
def natural_sort_key(s):
    parts = re.split(r'(\d+\.\d+|\d+)', s)
    key = []
    for text in parts:
        try:
            key.append(float(text))
        except ValueError:
            key.append(text.lower())
    return key

def loadSpheres(filename, colour):
    """Show a PQR as its real probe spheres, in one colour.

    The radius is read from the file text and assigned with `alter`. Do NOT use
    `vdw=b`: PyMOL parses a .pqr extension with its PQR reader, which does not
    put this column in the B-factor, so `vdw=b` sets every radius to zero and
    the spheres vanish while the object is still loaded.
    """
    obj_name = os.path.splitext(os.path.basename(filename))[0]
    cmd.load(filename, obj_name)
    radii_list = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                radii_list.append(float(line.split()[-1]))
    if radii_list:
        cmd.alter(obj_name, "vdw = radii_list.pop(0)",
                  space={'radii_list': radii_list})
    cmd.hide("everything", obj_name)
    cmd.show("spheres", obj_name)
    cmd.color(colour, obj_name)
    return obj_name

# Which number in the name is the rank, tried in the order the two producers
# write them. chl0.pqr counts from 0; CAVER writes tun_cl_001_1.pdb and colours
# by the cluster 001, not by the trailing tunnel index within it, and numbers
# its clusters from 1 -- so that one is shifted down to share the 0-based scale.
RANK_PATTERNS = ((r'chl(\d+)', 0),
                 (r'cl_(\d+)', 1),
                 (r'(\d+)', 0))      # anything else: the first number in the name

UNRANKED = "grey60"   # for a file whose name carries no number

def channelRank(filename):
    """Rank deciding the colour: the number the producer put in the file name.

    Read from the name rather than from the position in the loaded list, so a
    channel keeps its colour whether or not its lower-numbered siblings were
    written and whatever else the glob picked up alongside it. None when the
    name holds no number at all.
    """
    base = os.path.basename(filename)
    for pattern, first in RANK_PATTERNS:
        match = re.search(pattern, base)
        if match:
            return max(0, int(match.group(1)) - first)
    return None

def freeName(name):
    """`name`, or the first numbered variant of it no loaded object has taken.

    The all-channels dump written beside the per-channel files loads as an
    object named after the set, and PyMOL refuses to make a group over a name
    an ordinary object already holds.
    """
    taken = set(cmd.get_names("objects"))
    suffix = ""
    while name + suffix in taken:
        suffix = str(int(suffix or 1) + 1)
    return name + suffix

def loadCifChannels(path):
    """Draw every object in a tunnels-schema mmCIF, one PyMOL object each.

    The same spheres the PQR holds, from the same profile: x/y/z and radius per
    sample. They are built as a PDB string rather than one pseudoatom per
    sphere, which for a thousand samples is the difference between instant and
    a visible wait, and the radius is then assigned exactly as loadSpheres does
    for a PQR.
    """
    loops = cifLoops(path)
    if "sb_ncbr_channel_profile" not in loops:
        print(f"  {path}: no _sb_ncbr_channel_profile loop")
        return {}

    columns, rows = loops["sb_ncbr_channel_profile"]
    index = {name: columns.index(name) for name in columns}
    needed = ("channel_id", "x", "y", "z", "radius")
    if any(name not in index for name in needed):
        print(f"  {path}: profile loop is missing one of {needed}")
        return {}

    samples = {}
    for row in rows:
        samples.setdefault(row[index["channel_id"]], []).append(
            (float(row[index["x"]]), float(row[index["y"]]),
             float(row[index["z"]]), float(row[index["radius"]])))

    # The type each id was written under, so channels, pores and links can be
    # grouped apart the way the PQR output keeps them in separate files.
    #
    # The group names are this module's, not the schema's: a Tunnel goes into
    # chnl_grp and a Path into link_grp, so that a directory opens with the same
    # group names whichever format it was written in. Naming them after the
    # schema instead would give one run chnl_grp and the other tunnel_grp for
    # the very same channels.
    kinds = {}
    if "sb_ncbr_channel" in loops:
        cols, crows = loops["sb_ncbr_channel"]
        if "id" in cols and "type" in cols:
            for row in crows:
                kinds[row[cols.index("id")]] = row[cols.index("type")]

    groups = {}
    for channel_id in sorted(samples, key=natural_sort_key):
        spheres = samples[channel_id]
        rank = channelRank(channel_id)
        colour = caverColour(rank if rank is not None else 0)

        text = "".join(
            "ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
            % (i + 1, i + 1, x, y, z, 1.00, radius)
            for i, (x, y, z, radius) in enumerate(spheres))

        obj = channel_id
        cmd.read_pdbstr(text, obj)
        radii_list = [radius for _, _, _, radius in spheres]
        cmd.alter(obj, "vdw = radii_list.pop(0)",
                  space={'radii_list': radii_list})
        cmd.hide("everything", obj)
        cmd.show("spheres", obj)
        cmd.color(colour, obj)

        kind = kinds.get(channel_id, "Tunnel")
        groups.setdefault(kind, []).append(obj)
        print(f"  {obj:<20s} {colour}  ({kind}, {len(spheres)} spheres)")

    return groups

# The label each object's REMARK carries in a PQR, in the schema's words where
# it has them, so that a PQR's objects are grouped exactly as an mmCIF's are.
PQR_KINDS = {"channel": "Tunnel", "pore": "Pore", "link": "Path",
             "cavity": "Cavity"}

def loadPqrChannels(path):
    """Draw every object in a PQR holding several, one PyMOL object each.

    The file keeps them apart by residue number alone, so PyMOL loads it whole
    as a single object whose channels cannot be hidden or coloured apart. Split
    here the way loadCifChannels splits an mmCIF, with the same names, colours
    and groups, so a run opens the same way whichever format it was written in.

    A multi-model PQR, a MODEL per frame as mergeFramesPQR writes it, gives
    every object a state per frame, titled with the frame's number: the object
    of rank n holds the n-th channel of each frame, and stepping through the
    states steps through the frames. In a frame with fewer channels than that
    the object has nothing to show.

    Surface cavities are drawn as the surface around their markers, as the VMD
    and ChimeraX scripts draw them. A file of connected cavities and channels
    holds the cavities on chain C and the channels on chain H, each numbered
    from 1, so the chain is part of what tells one object from another.
    """
    # (MODEL number, {(chain, index): spheres}) per frame; a file without MODEL
    # records is a single frame, numbered None. placeholders: state -> the NIL
    # atom mergeFramesPQR puts in a frame that found nothing.
    models, kinds, label, placeholders = [], {}, "channel", {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                models.append((line[5:].strip(), {}))
            elif line.startswith("REMARK"):
                words = line.split()
                if len(words) > 2 and words[1] in PQR_KINDS:
                    label = words[1]
            elif line.startswith(("ATOM", "HETATM")) and line[17:20] == "NIL":
                if not models:
                    models.append((None, {}))
                placeholders[len(models)] = line
            elif line.startswith(("ATOM", "HETATM")) and line[17:20] == "FIL":
                if not models:
                    models.append((None, {}))
                # Fixed columns for the coordinates: three %8.3f values run
                # together without a space once one of them reaches -100.
                key = (line[21], int(line[22:26]) - 1)
                kinds.setdefault(key, "cavity" if key[0] == "C" else label)
                models[-1][1].setdefault(key, []).append(
                    (float(line[30:38]), float(line[38:46]),
                     float(line[46:54]), float(line.split()[-1])))

    # Colours by rank, as in the VMD and ChimeraX scripts: a file holding both
    # cavities and channels colours the cavities in turn and the channels after
    # them, so that the two are not given the same colours.
    cavities = sorted(key for key in kinds if kinds[key] == "cavity")
    others = sorted(key for key in kinds if kinds[key] != "cavity")
    if cavities and others:
        colours = dict((key, caverColour(position))
                       for position, key in enumerate(cavities + others))
    else:
        colours = dict((key, caverColour(key[1])) for key in kinds)

    groups = {}
    for key in cavities + others:
        index = key[1]
        colour = colours[key]
        obj = freeName(f"{kinds[key]}{index}")
        count = frames = 0

        for state, (number, samples) in enumerate(models, start=1):
            spheres = samples.get(key)
            if not spheres:
                continue
            text = "".join(
                "ATOM  %5d  H   FIL T%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                % (i + 1, i + 1, x, y, z, 1.00, radius)
                for i, (x, y, z, radius) in enumerate(spheres))
            # Discrete, so that each state keeps atoms of its own: the channel
            # of one frame is as long as it is, not as the first frame's.
            cmd.read_pdbstr(text, obj, state=state, discrete=1)
            if number is not None:
                cmd.set_title(obj, state, f"frame {number}")
            count += len(spheres)
            frames += 1

        # Read as PDB text, which keeps the radius column as the B-factor in
        # every state (unlike PyMOL's PQR reader, see loadSpheres).
        cmd.alter(obj, "vdw = b")
        cmd.hide("everything", obj)
        if kinds[key] == "cavity":
            # A cavity is a cloud of markers of one size, not a route of probe
            # spheres, so it is drawn as the surface around them. Every marker
            # is a hydrogen, which a surface leaves out unless told otherwise.
            cmd.set("surface_mode", 1, obj)
            cmd.show("surface", obj)
        else:
            cmd.show("spheres", obj)
        cmd.color(colour, obj)

        kind = PQR_KINDS[kinds[key]]
        groups.setdefault(kind, []).append(obj)
        note = f" in {frames} of {len(models)} frames" if len(models) > 1 else ""
        print(f"  {obj:<20s} {colour}  ({kind}, {count} spheres{note})")

    # A frame that found nothing holds only its placeholder, which goes into an
    # object of its own, drawn as nothing. The frame then keeps its state even
    # where no object has anything in it - an empty last frame would otherwise
    # leave PyMOL counting one frame fewer than the file holds.
    if placeholders:
        empty = freeName("empty_frames")
        for state, line in sorted(placeholders.items()):
            cmd.read_pdbstr(line, empty, state=state, discrete=1)
            number = models[state - 1][0]
            if number is not None:
                cmd.set_title(empty, state, f"frame {number}, nothing found")
        cmd.hide("everything", empty)
        print(f"  {empty:<20s} {len(placeholders)} frame(s) that found nothing")

    if len(models) > 1:
        print(f"  {len(models)} frames, a state each: step through them with the "
              f"arrow keys or the movie controls.")

    return groups

# both sets read their rank off the same 0-based scale, so the first channel of
# either program is blue and the two stay comparable side by side
sets = [("chnl_grp", sorted(glob.glob(channel_regex), key=natural_sort_key), False),
        ("tun_grp", sorted(glob.glob("tun_*"), key=natural_sort_key), True)]

CIF_GROUPS = {"Tunnel": "chnl_grp", "Pore": "pore_grp", "Path": "link_grp",
              "Cavity": "cav_grp"}

if cif_file or pqr_files:
    groups = loadCifChannels(cif_file) if cif_file else {}
    for path in pqr_files:
        for kind, objects in loadPqrChannels(path).items():
            groups.setdefault(kind, []).extend(objects)
    for kind, objects in sorted(groups.items()):
        group = CIF_GROUPS.get(kind, kind.lower() + "_grp")
        cmd.group(freeName(group), " ".join(objects))
    cmd.rebuild()
    cmd.set("sphere_scale", 1.0)
    cmd.set("sphere_quality", 2)
    cmd.bg_color("white")
    cmd.zoom()
    if outfile:
        cmd.save(outfile)
        print(f"Saved to {outfile}")
    print(f"Success: {sum(len(o) for o in groups.values())} object(s) "
          f"loaded from {', '.join(([cif_file] if cif_file else []) + pqr_files)}.")
elif not any(files for _, files, _ in sets):
    print("Error: No channel files found. Check your working directory (pwd).")
else:
    for group, files, start_off in sets:
        objects, unranked = [], []
        ranks = [channelRank(f) for f in files]
        # channels.pqr, the dump of every channel at once, sits beside the
        # per-channel files and carries no number, so there is no rank to
        # colour it by -- and being a copy of all of them it would hide them.
        # Grey it and start it switched off. Unless nothing in the set has a
        # rank, in which case the set is dumps only and position will do.
        grey_the_unranked = any(rank is not None for rank in ranks)
        for i, (filename, rank) in enumerate(zip(files, ranks)):
            grey = rank is None and grey_the_unranked
            colour = UNRANKED if grey else caverColour(i if rank is None else rank)
            obj = loadSpheres(filename, colour)
            objects.append(obj)
            print(f"  {obj:<20s} {colour}{'  (no rank in name, off)' if grey else ''}")
            if grey:
                unranked.append(obj)
        if objects:
            name = freeName(group)
            cmd.group(name, " ".join(objects))
            if start_off:
                cmd.disable(name)
        for obj in unranked:          # after the group, which enables its members
            cmd.disable(obj)
    cmd.rebuild()
    cmd.set("sphere_scale", 1.0)
    cmd.set("sphere_quality", 2)
    cmd.bg_color("white")
    cmd.zoom()
    if outfile:
        cmd.save(outfile)
        print(f"Saved to {outfile}")
    print("Success: Files loaded in perfect sequential order with custom radii!")
'''

_VIS_CHANNELS_FLAT_TRAJ_SCRIPT = r"""

import re
import glob
import sys
import colorsys
import os
import argparse
import time


parser = argparse.ArgumentParser(description=("Load trajectory channel PQR files into PyMOL, color the channels "
        "according to externally generated cluster labels, and save one PyMOL "
        "session (.pse) for each supplied clustering."),formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("file_pattern",help=("Glob pattern selecting the trajectory channel PQR files. "
        "Quote the pattern so that it is expanded by this script rather than "
        "by the shell, e.g. 'channels/frame*_channel*.pqr'."))

parser.add_argument("--protein-file", dest="protein_file", default=None,
    help=("Optional protein structure to load together with the channels. "
        "The protein is displayed as a cartoon and included in each saved "
        "PyMOL session."))

parser.add_argument("--labels",nargs="+",default=None,metavar="LABEL_FILE",
    help=("One or more text files containing cluster labels. Each file may "
        "contain either one clustering or multiple clustering columns. "
        "Rows correspond to channels in ascending frame number, then "
        "ascending channel rank. Label -1 is noise. Label -2 excludes the "
        "channel completely from that visualization."))

parser.add_argument("--out-file",dest="out_file", default="traj_out.pse",
    help=("Output PyMOL session filename. For a single clustering this is the "
        "requested filename. If multiple clusterings are supplied, the "
        "clustering name is appended to the stem. Existing files are not "
        "overwritten."))

args = parser.parse_args()

file_pattern = args.file_pattern
protein_file = args.protein_file
label_files = args.labels
out_file = args.out_file

files = sorted(glob.glob(file_pattern))
if not files:
    parser.error( "No trajectory PQR files matched pattern: %s" % file_pattern)

if protein_file:
    protein_name = os.path.splitext( os.path.basename(protein_file))[0]
    cmd.load(protein_file, protein_name)
    cmd.hide("everything", protein_name)
    cmd.show("cartoon", protein_name)
    cmd.set("all_states", 1, protein_name)

CAVER_PRIMARIES = [(0.0, 0.0, 1.0),
                    (0.0, 1.0, 0.0),
                    (1.0, 0.0, 0.0),
                    (0.0, 1.0, 1.0),
                    (1.0, 1.0, 0.0),
                    (1.0, 0.0, 1.0),]

GOLDEN_ANGLE = (3.0 - 5.0 ** 0.5) / 2.0
HUE_OFFSET = 0.796

SATURATION_VALUE = ((0.95, 0.55),
                    (0.65, 0.65),
                    (0.65, 0.95))


def caverColour(rank):
    rank = max(0, int(rank))
    if rank < len(CAVER_PRIMARIES):
        name = "caver%d" % (rank + 1)
        rgb = CAVER_PRIMARIES[rank]
    else:
        step = rank - len(CAVER_PRIMARIES)
        saturation, value = SATURATION_VALUE[step % len(SATURATION_VALUE)]
        name = "gen%d" % rank
        rgb = colorsys.hsv_to_rgb((HUE_OFFSET+ (step + 1) * GOLDEN_ANGLE) % 1.0,
                                    saturation,
                                    value)
    cmd.set_color(name, list(rgb))

    return name

UNCLUSTERED = "grey60"

def clusterColour(cluster):
    if str(cluster).lower() in ("noise", "-1"):
        return UNCLUSTERED
    try:
        return caverColour(int(cluster))
    except (TypeError, ValueError):
        raise ValueError("Cluster labels must be integers or 'noise'; got %r" % cluster)

FRAME_PREFIX = re.compile(r'^[A-Za-z]+(\d+)')
CHANNEL_RANK = re.compile(r'chl(\d+)')
LINK_RANK = re.compile(r'lnk(\d+)')
CLUSTER_ID = re.compile(r'clid(\w+)')

def channelRank(rest):
    match = CHANNEL_RANK.search(rest)
    return int(match.group(1)) if match else None

def linkRank(rest):
    match = LINK_RANK.search(rest)
    return int(match.group(1)) if match else None

def clusterId(rest):
    match = CLUSTER_ID.search(rest)
    return match.group(1) if match else None

def splitFrame(filename):
    base = os.path.basename(filename)
    match = FRAME_PREFIX.match(base)
    if not match:
        return None, base

    return int(match.group(1)), base[match.end():]

by_chl_frame = {}
frame_ranks = {}

for fname in files:
    frame, rest = splitFrame(fname)
    if frame is None:
        continue
    chl_rank = channelRank(rest)
    if chl_rank is not None:
        frame_ranks.setdefault(frame, []).append(chl_rank)
        by_chl_frame.setdefault(frame, {})[chl_rank] = fname

if not by_chl_frame:
    parser.error( "Files matched the glob, but no frame/channel information could "
        "be parsed from their filenames.")

channel_order = []
for frame in sorted(by_chl_frame):
    frame_channels = by_chl_frame[frame]
    for rank in sorted(frame_channels):
        channel_order.append(frame_channels[rank])

n_channels = len(channel_order)
print("Found %d channels across %d frames." % (n_channels, len(by_chl_frame)))

def readLabelFile(path):
    import csv
    if not os.path.isfile(path):
        raise ValueError("Label file does not exist: %s" % path)

    raw_rows = []
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle,delimiter="\t")
        for row in reader:
            row = [field.strip() for field in row]
            if (not row or all(field == "" for field in row)):
                continue

            if row[0].startswith("#"):
                continue

            raw_rows.append(row)

    if not raw_rows:
        raise ValueError( "No labels found in %s" % path)

    try:
        [int(value) for value in raw_rows[0]]
        header = None
        data_rows = raw_rows

    except ValueError:
        header = raw_rows[0]
        data_rows = raw_rows[1:]

    if not data_rows:
        raise ValueError("Label file %s contains a header " "but no label rows." % path)

    ncols = len(data_rows[0])
    if ncols == 0:
        raise ValueError("No label columns found in %s"% path)

    for row_num, row in enumerate(data_rows,start=2 if header is not None else 1):
        if len(row) != ncols:
            raise ValueError("Inconsistent number of columns in %s ""at row %d: expected %d, found %d."% (
                                path, row_num, ncols ,len(row)))
                                
    rows = []

    for row_num, row in enumerate(data_rows,start=2 if header is not None else 1):
        try:
            rows.append([int(value) for value in row])

        except ValueError:
            raise ValueError("Non-integer cluster label in %s " 
                                "at row %d:\n%s" % (path,
                                                    row_num,
                                                    "\t".join(row)))

    if len(rows) != n_channels:
        raise ValueError("%s contains %d label rows, "
                         "but %d channels were found."
                        % (path, len(rows), n_channels))

    if header is not None:
        if len(header) != ncols:
            raise ValueError("Header in %s contains %d names "
                            "but the data contain %d columns."
                            % (path, len(header),ncols))

        if any(not name for name in header):
            raise ValueError("Header in %s contains an empty "
                            "column name." % path)

        if len(set(header)) != len(header):
            raise ValueError("Header in %s contains duplicate "
                                "clustering names." % path)

    base = os.path.splitext(os.path.basename(path))[0]

    clusterings = []
    for col in range(ncols):
        if header is not None:
            name = header[col]

        elif ncols == 1:
            name = base

        else:
            name = "%s_%d" % (base,col + 1)

        labels = [row[col] for row in rows]
        clusterings.append((name, labels))

    return clusterings


def freePath(path):
    if not os.path.exists(path):
        return path

    stem, ext = os.path.splitext(path)
    n = 1
    while os.path.exists("%s_%d%s" % (stem, n, ext)):
        n += 1

    return "%s_%d%s" % (stem, n, ext)


def safeName(name):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name))
    return name.strip("_.") or "clustering"

def makeOutputPath(base_path, clustering_name=None, multiple=False,):
    stem, ext = os.path.splitext(base_path)
    if ext.lower() != ".pse":
        if ext:
            stem = base_path
        ext = ".pse"

    if (multipleand clustering_name is not None):
        path = "%s_%s%s" % (stem, safeName(clustering_name), ext)
    else:
        path = stem + ext

    return freePath(path)


def applyClustering(channel_objects,labels,):
    objects_by_label = {}
    for name, label in zip(channel_objects,labels,):
        label = int(label)
        # -2 means the channel is excluded from this visualization entirely.
        if label == -2:
            continue
        objects_by_label.setdefault(label,[]).append(name)
    group_names = []
    for label, names in objects_by_label.items():
        if label == -1:
            group_name = "noise"
        else:
            group_name = "cluster_%d" % label

        selection = " or ".join(names)

        cmd.group(group_name," ".join(names))
        cmd.color( clusterColour(label),selection)
        cmd.alter(selection,"segi='%s'" % label)
        group_names.append(group_name)

    return group_names

clusterings = []

if label_files:
    for label_file in label_files:
        clusterings.extend( readLabelFile(label_file))

if clusterings:
    load_mask = [False] * n_channels
    for clustering_name, labels in clusterings:
        for i, label in enumerate(labels):
            if label != -2:
                load_mask[i] = True
else:
    load_mask = [True] * n_channels

n_load = sum(load_mask)
print("Loading %d of %d channels." % ( n_load, n_channels))

cmd.set("auto_zoom", 0)
cmd.set("defer_builds_mode", 3)
channel_objects = []
channel_indices = []
t0 = time.perf_counter()

for global_index, fname in enumerate(channel_order):
    if not load_mask[global_index]:
        continue

    basename = os.path.splitext(os.path.basename(fname))[0]
    name = "chn_%06d_%s" % (global_index,basename)
    cmd.load(fname,name)
    channel_objects.append(name)
    channel_indices.append(global_index)

print("Time to load %d files: %.2f s" % (len(channel_objects), time.perf_counter() - t0))


if not clusterings:
    for name in channel_objects:
        cmd.color(UNCLUSTERED, name)
        cmd.alter(name,"segi='ZZZZ'")
    cmd.bg_color("white")
    output_path = makeOutputPath(out_file)
    cmd.save(output_path)
    print("Saved session: %s" % output_path)

else:
    multiple = len(clusterings) > 1
    for clustering_name, full_labels in clusterings:
        print("Applying clustering '%s' (%d global labels)." % (clustering_name,len(full_labels)))
        labels = [full_labels[i] for i in channel_indices]
        excluded_objects = [name for name, label in zip(channel_objects,labels) if label == -2]
        included_objects = [name for name, label in zip(channel_objects,labels) if label != -2]

        if excluded_objects:
            cmd.disable(" ".join(excluded_objects))

        if included_objects:
            cmd.enable(" ".join(included_objects))

        t0 = time.perf_counter()

        group_names = applyClustering(channel_objects,labels)
        t1 = time.perf_counter()
        print("Clustering applied in %.2f s" % (t1 - t0))

        output_path = makeOutputPath(out_file, clustering_name=clustering_name, multiple=multiple,)
        cmd.bg_color("white")
        cmd.save(output_path)
        print( "Saved clustering '%s': %s" % (clustering_name, output_path))

        for group_name in group_names:
            cmd.ungroup(group_name)
            cmd.delete(group_name)

        if excluded_objects:
            cmd.enable(" ".join(excluded_objects))

"""

_VIS_CHANNELS_TRAJ_SCRIPT = r"""
import glob
import os
import re
import sys
import colorsys

# Unified state-based viewer
# This viewer accepts BOTH trajectory representations used by ProDy/CaviTracer:

#   1. One multi-MODEL PQR:
#        channels_frames.pqr
#
#   2. One PQR per frame:
#        channels100_sp0_chl0.pqr
#        channels100_sp0_chl1.pqr
#        channels101_sp0_chl0.pqr
#        ...
#
# The two paths end in the same PyMOL representation: one object per channel
# rank (or cluster), one state per frame.

# For per-frame cluster files, several files may legitimately have the same
# cluster id AND frame number.  They are accumulated in a list and all of their
# geometry is placed into that cluster's state.  Nothing is discarded.

TRAJ_EXTENSIONS = ('.dcd', '.xtc', '.trr', '.nc', '.dtr')

args = [a for a in sys.argv[1:] if a != "--"]
file_pattern = None
protein_file = None
traj_file = None

# Existing files are inputs; the first non-file argument is the glob pattern.
for arg in args:
    lower = arg.lower()
    if os.path.isfile(arg) and lower.endswith(TRAJ_EXTENSIONS):
        if traj_file is None:
            traj_file = arg
    elif os.path.isfile(arg) and lower.endswith(".pse"):
        out_file = arg
    elif os.path.isfile(arg):
        if protein_file is None:
            protein_file = arg
    elif file_pattern is None:
        file_pattern = arg
if file_pattern is None:
    file_pattern = "*.pqr"

files = sorted(glob.glob(file_pattern))
if not files:
    print("No files matched %r" % file_pattern)
    raise SystemExit

# If the command line supplied a glob plus a protein, the protein is the input
# file that is not itself a channel PQR.  Otherwise use the first non-PQR file.

if protein_file:
    protein_name = os.path.splitext(os.path.basename(protein_file))[0]
    cmd.load(protein_file, protein_name)
    cmd.hide("everything", protein_name)
    cmd.show("cartoon", protein_name)
    cmd.set("all_states", 1, protein_name)


CAVER_PRIMARIES = [(0.0, 0.0, 1.0), (0.0, 1.0, 0.0), (1.0, 0.0, 0.0),
                   (0.0, 1.0, 1.0), (1.0, 1.0, 0.0), (1.0, 0.0, 1.0)]
GOLDEN_ANGLE = (3.0 - 5.0 ** 0.5) / 2.0
HUE_OFFSET = 0.796
SATURATION_VALUE = ((0.95, 0.55), (0.65, 0.65), (0.65, 0.95))


def caverColour(rank):
    rank = max(0, int(rank))
    if rank < len(CAVER_PRIMARIES):
        name, rgb = "caver%d" % (rank + 1), CAVER_PRIMARIES[rank]
    else:
        step = rank - len(CAVER_PRIMARIES)
        saturation, value = SATURATION_VALUE[step % len(SATURATION_VALUE)]
        name = "gen%d" % rank
        rgb = colorsys.hsv_to_rgb(
            (HUE_OFFSET + (step + 1) * GOLDEN_ANGLE) % 1.0,
            saturation, value)
    cmd.set_color(name, list(rgb))
    return name


UNCLUSTERED = "grey60"


def clusterColour(cluster):
    if str(cluster).lower() in ("noise", "-1"):
        cmd.set_color("noise", [0.6, 0.6, 0.6])
        return UNCLUSTERED
    try:
        return caverColour(int(cluster))
    except (TypeError, ValueError):
        return caverColour(abs(hash(str(cluster))) % 1000)

FRAME_PREFIX = re.compile(r'^[A-Za-z]+(\d+)')
CHANNEL_RANK = re.compile(r'chl(\d+)')
LINK_RANK = re.compile(r'lnk(\d+)')
CLUSTER_ID = re.compile(r'clid(\w+)')


def splitFrame(filename):
    base = os.path.basename(filename)
    match = FRAME_PREFIX.match(base)
    if not match:
        return None, base
    return int(match.group(1)), base[match.end():]


def channelRank(rest):
    match = CHANNEL_RANK.search(rest)
    return int(match.group(1)) if match else None


def linkRank(rest):
    match = LINK_RANK.search(rest)
    return int(match.group(1)) if match else None


def clusterId(rest):
    match = CLUSTER_ID.search(rest)
    return match.group(1) if match else None

# gather, per rank (and per cluster id, if present), the file per frame 
by_chl_rank = {}      # rank -> {frame: filename}, used when no file carries a clid tag
by_cluster = {}       # cluster id (str) -> {frame: filename}, used once any file does
by_lnk_rank = {}      # rank -> {frame: filename}
skipped_dumps = 0
collisions = []          # (cluster_id, frame, kept_fname, dropped_fname)
unclustered_chl = []     # (rank, frame, fname) for chl files with no clid tag
any_clustered = False



PQR_KINDS = {"channel": "Tunnel", "pore": "Pore", "link": "Path",
             "cavity": "Cavity"}


def isMultiModelPQR(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                return True
    return False


def pqrToPdbText(lines):
    '''Convert PQR atom records to PDB-like records with radius in B-factor.

    read_pdbstr() then gives us exact per-atom vdw values via ``vdw = b``.
    '''
    out = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20] not in ("FIL", "NIL"):
            continue
        try:
            radius = float(line.split()[-1])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            serial = int(line[6:11]) if line[6:11].strip() else len(out) + 1
            resid = int(line[22:26]) if line[22:26].strip() else len(out) + 1
            chain = line[21:22] or " "
        except (ValueError, IndexError):
            continue
        out.append(
            "ATOM  %5d  H   FIL %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
            % (serial, chain, resid, x, y, z, 1.00, radius))
    return "".join(out)


def loadPqrObject(fname, obj):
    '''Load one per-frame PQR into a temporary PyMOL object and set radii.'''
    cmd.load(fname, obj)
    radii = []
    with open(fname) as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    radii.append(float(line.split()[-1]))
                except ValueError:
                    pass
    if radii:
        cmd.alter(obj, "vdw = radii.pop(0)", space={'radii': radii})
    return obj


def addFilesToState(fnames, obj, state):
    '''Put all PQRs in *fnames* into one state of *obj*.

    Multiple files are deliberately combined rather than treated as a
    collision. This is the important case for clustered trajectories: several
    channels from one frame can belong to the same cluster.
    '''
    temps = []
    try:
        for i, fname in enumerate(fnames):
            tmp = "__vis_tmp_%d_%d" % (state, i)
            # Avoid stale names from a partially interrupted PyMOL run.
            if tmp in cmd.get_names("objects"):
                cmd.delete(tmp)
            loadPqrObject(fname, tmp)
            temps.append(tmp)
        if temps:
            selection = " or ".join("%s" % name for name in temps)
            cmd.create(obj, selection, 1, state, 1)
    finally:
        for tmp in temps:
            if tmp in cmd.get_names("objects"):
                cmd.delete(tmp)


def buildMultistate(file_map, obj, colour, frame_to_state, n_states):
    '''file_map: frame -> filename or frame -> list of filenames.'''
    n_written = 0
    for frame in sorted(file_map):
        value = file_map[frame]
        fnames = value if isinstance(value, (list, tuple)) else [value]
        fnames = [f for f in fnames if f]
        if not fnames:
            continue
        addFilesToState(fnames, obj, frame_to_state[frame])
        n_written += 1

    cmd.hide("everything", obj)
    cmd.show("spheres", obj)
    cmd.color(colour, obj)
    print("  %-20s %s  present in %d/%d frame(s)" %
          (obj, colour, n_written, n_states))
    return n_written

def loadMultiModelPQR(path):
    '''Build one stateful object per channel/link/cavity from a multi-MODEL PQR.'''
    models = []
    current = None
    current_kind = "channel"
    kind_by_key = {}

    with open(path) as handle:
        for line in handle:
            if line.startswith("MODEL"):
                number_text = line[5:].strip()
                try:
                    number = int(number_text)
                except ValueError:
                    number = len(models) + 1
                current = {"number": number, "records": {}}
                models.append(current)
                continue

            if line.startswith("REMARK"):
                words = line.split()
                if len(words) > 2 and words[1] in PQR_KINDS:
                    current_kind = words[1]
                continue

            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[17:20] != "FIL":
                continue

            if current is None:
                current = {"number": 1, "records": {}}
                models.append(current)

            try:
                chain = line[21:22] or " "
                resid = int(line[22:26])
                radius = float(line.split()[-1])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue

            key = (chain, resid)
            kind = "cavity" if chain == "C" else current_kind
            kind_by_key.setdefault(key, kind)
            current["records"].setdefault(key, []).append((x, y, z, radius))

    if not models:
        print("No FIL records found in %s" % path)
        return [], [], 0

    frame_numbers = [m["number"] for m in models]
    frame_to_state = {frame: i + 1 for i, frame in enumerate(frame_numbers)}

    # Stable object order: cavities first, then everything else by residue id.
    keys = sorted(kind_by_key, key=lambda k: (kind_by_key[k] != "cavity", k[0], k[1]))
    built = []
    groups = {}

    for key in keys:
        kind = kind_by_key[key]
        rank = max(0, key[1] - 1)
        obj_kind = PQR_KINDS.get(kind, "Tunnel")
        obj = "%s%d" % (kind, rank)
        if obj in cmd.get_names("objects"):
            suffix = 1
            base = obj
            while obj in cmd.get_names("objects"):
                obj = base + str(suffix)
                suffix += 1

        colour = caverColour(rank)
        written = 0
        for state, model in enumerate(models, start=1):
            samples = model["records"].get(key)
            if not samples:
                continue
            lines = []
            for i, (x, y, z, radius) in enumerate(samples, start=1):
                lines.append(
                    "ATOM  %5d  H   FIL %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n"
                    % (i, key[0], i, x, y, z, 1.00, radius))
            cmd.read_pdbstr("".join(lines), obj, state=state, discrete=1)
            cmd.alter(obj, "vdw = b")
            cmd.set_title(obj, state, "frame %s" % model["number"])
            written += 1

        cmd.hide("everything", obj)
        if kind == "cavity":
            cmd.set("surface_mode", 1, obj)
            cmd.show("surface", obj)
        else:
            cmd.show("spheres", obj)
        cmd.color(colour, obj)
        groups.setdefault(obj_kind, []).append(obj)
        built.append(obj)
        print("  %-20s %s  (%s, %d/%d frame(s))" %
              (obj, colour, obj_kind, written, len(models)))

    return built, groups, len(models)


def loadPerFramePQRs(files):
    by_chl_rank = {}
    by_cluster = {}
    by_lnk_rank = {}
    unclustered_chl = []
    skipped_dumps = 0
    any_clustered = False

    for fname in files:
        frame, rest = splitFrame(fname)
        if frame is None:
            continue

        chl_rank = channelRank(rest)
        if chl_rank is not None:
            by_chl_rank.setdefault(chl_rank, {}).setdefault(frame, []).append(fname)
            cluster = clusterId(rest)
            if cluster is not None:
                any_clustered = True
                by_cluster.setdefault(cluster, {}).setdefault(frame, []).append(fname)
            else:
                unclustered_chl.append((chl_rank, frame, fname))
            continue

        lnk_rank = linkRank(rest)
        if lnk_rank is not None:
            by_lnk_rank.setdefault(lnk_rank, {}).setdefault(frame, []).append(fname)
            continue

        skipped_dumps += 1

    all_frames = sorted({frame for maps in (by_chl_rank, by_lnk_rank, by_cluster)
                         for rank_map in maps.values() for frame in rank_map})
    if not all_frames:
        print("No ranked channel/link files found among %d file(s); "
              "%d unranked dump(s) skipped." % (len(files), skipped_dumps))
        return [], [], 0

    frame_to_state = {frame: i + 1 for i, frame in enumerate(all_frames)}
    n_states = len(all_frames)

    print("Found %d frame(s), %d channel rank(s), %d cluster(s), %d link rank(s) "
          "(%d unranked dump(s) skipped)." %
          (n_states, len(by_chl_rank), len(by_cluster), len(by_lnk_rank),
           skipped_dumps))

    if any_clustered:
        # Preserve partially tagged runs. Untagged channels get a unique
        # synthetic cluster key based on their rank, and values are lists just
        # like the real clusters.
        for rank, frame, fname in unclustered_chl:
            key = "rank%d_untagged" % rank
            by_cluster.setdefault(key, {}).setdefault(frame, []).append(fname)

    built_chl, built_lnk = [], []

    if any_clustered:
        for cluster in sorted(by_cluster, key=lambda c: (str(c).lower() == "noise", str(c))):
            obj = "cl%s" % cluster
            colour = clusterColour(cluster)
            if buildMultistate(by_cluster[cluster], obj, colour,
                               frame_to_state, n_states):
                built_chl.append(obj)
                if str(cluster).lower() in ("noise", "-1"):
                    cmd.disable(obj)
    else:
        for rank in sorted(by_chl_rank):
            obj = "chl%d" % rank
            if buildMultistate(by_chl_rank[rank], obj, caverColour(rank),
                               frame_to_state, n_states):
                built_chl.append(obj)

    for rank in sorted(by_lnk_rank):
        obj = "lnk%d" % rank
        if buildMultistate(by_lnk_rank[rank], obj, caverColour(rank),
                           frame_to_state, n_states):
            built_lnk.append(obj)

    if built_chl:
        cmd.group("chnl_grp", " ".join(built_chl))
    if built_lnk:
        cmd.group("link_grp", " ".join(built_lnk))
        cmd.disable("link_grp")

    cmd.mset("1 -%d" % n_states)
    return built_chl, built_lnk, n_states

if len(files) == 1 and files[0].lower().endswith(".pqr") and isMultiModelPQR(files[0]):
    print("Detected multi-MODEL PQR: %s" % files[0])
    built, groups, n_states = loadMultiModelPQR(files[0])
    for kind, objects in sorted(groups.items()):
        group_name = {"Tunnel": "chnl_grp", "Pore": "pore_grp",
                      "Path": "link_grp", "Cavity": "cav_grp"}.get(
                          kind, kind.lower() + "_grp")
        if objects:
            cmd.group(group_name, " ".join(objects))
    cmd.mset("1 -%d" % n_states)
    print("Success: %d object(s) across %d state(s)." % (len(built), n_states))
else:
    # If a single ordinary PQR was supplied, it is treated as one frame.
    built_chl, built_lnk, n_states = loadPerFramePQRs(files)
    print("Success: %d channel/cluster object(s), %d link object(s) across %d state(s)." %
          (len(built_chl), len(built_lnk), n_states))

cmd.rebuild()
cmd.set("sphere_scale", 1.0)
cmd.set("sphere_quality", 2)
cmd.bg_color("white")
cmd.zoom()

def freePath(path):
    if not os.path.exists(path):
        return path
    stem, ext = os.path.splitext(path)
    n = 1
    while os.path.exists("%s_%d%s" % (stem, n, ext)):
        n += 1
    return "%s_%d%s" % (stem, n, ext)

out_pse = freePath('traj_out.pse')
cmd.save(out_file if 'out_file' in locals() else out_pse)
print("Saved session: %s" % out_pse)
"""


def _writeVisTrajScript(directory, pattern='channels*_sp*_chl*.pqr', flat=False):
    """Write the unified state-based PyMOL trajectory viewer.

    The viewer accepts either a multi-MODEL PQR or a glob of per-frame PQRs.
    This function is retained as the trajectory-specific writer API for
    backwards compatibility; the generated viewer is now the same state
    pathway used for both representations.
    """
    import os

    path = os.path.join(str(directory), 'vis_channels_traj.py')
    usage = ('`pymol vis_channels_traj.py -- "{0}" [protein.pdb]`'.format(pattern))
    script = None
    if flat:
        script = _VIS_CHANNELS_FLAT_TRAJ_SCRIPT
    else:
        script = _VIS_CHANNELS_TRAJ_SCRIPT

    if os.path.exists(path):
        LOGGER.info('View the trajectory with the PyMOL viewer already in {0}: '
                    '{1}.'.format(os.path.dirname(path), usage))
        return
    try:
        with open(path, 'w') as script_file:
            script_file.write(script)
    except (IOError, OSError) as err:
        _warn("Could not write the PyMOL trajectory viewer {0}: {1}".format(path, err))
    else:
        LOGGER.info('Wrote the PyMOL trajectory viewer {0}. View the output with '
                    '{1}.'.format(path, usage))


def _writeVisScript(directory, pattern='chl*.pqr'):
    """Leave ``vis_channels.py`` in ``directory`` unless it is already there.

    :func:`writePyMolCaviTracerScript` leaves the viewer beside what it writes, as
    CAVER leaves ``view.py`` beside its clusters; a channel run itself writes no
    script, that being what the viewer-script functions are for. An existing file
    is never overwritten: edits made to one copy survive a rerun, and so does a
    newer script left by an earlier one.

    One script serves both output formats. *pattern* is what the log tells the
    reader to pass, a glob for a PQR run and the file itself for an mmCIF one;
    the script reads whichever it is handed, and finds the other on its own when
    handed nothing, so a directory holding both opens either way.
    """
    import os

    path = os.path.join(str(directory), 'vis_channels.py')
    usage = '`pymol vis_channels.py -- <protein>.pdb "{0}"`'.format(pattern)

    if os.path.exists(path):
        # Left as it is, but still announced. The hint is the useful half of this
        # function, and a rerun into the same directory - the usual way to try a
        # different setting - used to get a script and no word on how to use it,
        # the one run that said so having scrolled away.
        LOGGER.info('View the output with the PyMOL viewer already in {0}: '
                    '{1}.'.format(os.path.dirname(path), usage))
        return
    try:
        with open(path, 'w') as script_file:
            script_file.write(_VIS_CHANNELS_SCRIPT)
    except (IOError, OSError) as err:
        # a viewer that cannot be written is no reason to lose the run
        _warn("Could not write the PyMOL viewer {0}: {1}".format(path, err))
    else:
        LOGGER.info('Wrote the PyMOL viewer {0}. View the output with '
                    '{1}.'.format(path, usage))


def _writeViewerPQR(result_file, objects, object_type, surface, num_samples):
    """Write *objects* to the PQR a viewer script draws: channels, pores or links
    as their spheres, surface cavities as their markers, and connected results as
    both, the cavities on chain C and the channels on chain H.

    One writer for the VMD, PyMOL and ChimeraX functions, which all draw these
    same records. *surface* is what :func:`calcSurfaceCavities` returned, needed
    only for the two cavity types."""

    if object_type in ('channels', 'pores', 'links'):

        labels = {'channels': 'channel', 'pores': 'pore', 'links': 'link'}
        label = labels[object_type]

        with open(str(result_file), 'w') as handle:
            atom_index = 1

            for object_index, obj in enumerate(objects):
                lines, written = ChannelCalculator._channelRecords(
                    object_index, obj, atom_index, num_samples,
                    label=label, name_sites=True)

                handle.writelines(lines)
                handle.write("\n")
                atom_index += written

    elif object_type == 'surface_cavities':

        vertices = np.asarray(surface[4])
        atom_index = 1
        drawn = []

        for cavity in objects:
            tetrahedra = cavity.tetrahedra

            if tetrahedra is None or len(tetrahedra) == 0:
                continue

            drawn.append((len(drawn), cavity,
                          vertices[np.asarray(tetrahedra, dtype=np.intp)]))

        if not drawn:
            raise ValueError("No surface cavity contains points that can be visualized")

        with open(str(result_file), 'w') as handle:
            for cavity_index, cavity, centers in drawn:
                lines, written = ChannelCalculator._cavityRecords(
                    cavity_index, cavity, centers, atom_index)

                handle.writelines(lines)
                atom_index += written

    else:
        _saveConnectedCavityChannels(objects, surface, result_file, separate=False,
                                    num_samples=num_samples)


def writeVmdMultiModelScript(pqr_file, atoms, output_path='.',
                             object_type='channels'):
    """Write a VMD script for a PQR already on disk, such as the multi-model
    one :func:`mergeFramesPQR` writes.

    The script reads the PQR itself, one molecule per frame, each named for its
    frame and coloured by rank as :func:`writeVmdCaviTracerScript` colours one
    structure, and VMD's animation slider shows one frame at a time. A PQR of
    one structure opens just as that function's own would.

    Rank is all that ties an object to its namesakes in other frames: channel 3
    of one frame and channel 3 of the next are the fourth cheapest of each, not
    one channel followed through the trajectory.

    :arg pqr_file: the PQR to draw, one structure or several frames.
    :type pqr_file: str

    :arg atoms: structure to draw the objects against. Its active coordinate
        set is written, a single backdrop for every frame.
    :type atoms: :class:`.Atomic`

    :arg output_path: directory the structure and the script are written to.
        The script finds the PQR where it is. Default is the current directory.
    :type output_path: str

    :arg object_type: what the PQR holds: ``'channels'`` (the default),
        ``'pores'``, ``'links'``, ``'surface_cavities'`` or
        ``'connected_cavities_channels'``, as in
        :func:`writeVmdCaviTracerScript`.
    :type object_type: str

    :returns: paths to the PQR, the structure PDB and the VMD script
    :rtype: tuple

    Usage:
    mergeFramesPQR('frames', 'channels_frames.pqr')
    writeVmdMultiModelScript('channels_frames.pqr', atoms)

    Next (bash console):
    $ vmd -e vis_channels_frames.tcl"""

    from pathlib import Path

    _requireCoords(atoms)

    object_types = ('channels', 'pores', 'links', 'surface_cavities',
                    'connected_cavities_channels')
    if object_type not in object_types:
        raise ValueError('object_type must be one of {0}'.format(
            ', '.join(repr(name) for name in object_types)))

    pqr_file = Path(pqr_file)
    if not pqr_file.is_file():
        raise ValueError('no PQR at {0}'.format(pqr_file))

    output_path = Path(output_path)
    if not output_path.exists():
        output_path.mkdir(parents=True)
    if not output_path.is_dir():
        raise ValueError("output_path must be a directory")

    number = 0
    while True:
        suffix = '' if number == 0 else '-{0}'.format(number)
        protein_file = output_path / ('protein' + suffix + '.pdb')
        script_file = output_path / ('vis_' + pqr_file.stem + suffix + '.tcl')
        if not (protein_file.exists() or script_file.exists()):
            break
        number += 1

    # The active frame only: written whole, a trajectory's structure would come
    # out as a model per frame, all of them drawn at once behind the channels.
    writePDB(str(protein_file), atoms, csets=atoms.getACSIndex())
    _writeVmdScript(script_file, protein_file, pqr_file, object_type)

    LOGGER.info("CaviTracer VMD files written:")
    LOGGER.info("    protein: {0}".format(protein_file))
    LOGGER.info("    results: {0}".format(pqr_file))
    LOGGER.info("    VMD script: {0}".format(script_file))
    LOGGER.info("View the result with: vmd -e {0}".format(script_file.name))

    return str(pqr_file), str(protein_file), str(script_file)
    
    
def writePyMolCaviTracerScript(objects, atoms, object_type='channels',
                               output_path='.', num_samples=5,
                               output_format=None, surface=None):
    """Prepare CaviTracer results for visualization in PyMOL.

    This function writes CaviTracer objects in the tunnels-schema mmCIF format
    or as a PQR, saves the supplied molecular structure as a PDB file, and
    creates the CaviTracer PyMOL visualization script, which reads either.

    The mmCIF is the one :func:`writeChannelsCIF` writes, and the PQR holds the
    same records as the one :func:`writeVmdCaviTracerScript` and
    :func:`writeChimeraXCaviTracerScript` write: every object in one file, told
    apart by its residue number, which the viewer splits it by. Surface cavities
    are drawn as the surface around their markers, as in those two.

    :arg objects: CaviTracer objects to visualize: channels, pores, links,
        surface cavities, or connected cavity-channel results.
    :type objects: list or CaviTracer object

    :arg atoms: Molecular structure to display together with the CaviTracer
        results.
    :type atoms: :class:`.Atomic`

    :arg object_type: Type of CaviTracer result. Accepted values are
        ``'channels'``, ``'pores'``, ``'links'``,
        ``'surface_cavities'`` and ``'connected_cavities_channels'``.
        Default is ``'channels'``.
    :type object_type: str

    :arg output_path: Directory in which the mmCIF or PQR, the PDB and the
        PyMOL viewer are saved. Default is the current directory.
    :type output_path: str or pathlib.Path

    :arg num_samples: Number of samples per tetrahedron used for the channel,
        pore, or link profile. Default is 5.
    :type num_samples: int

    :arg output_format: ``'mmcif'`` or ``'pqr'``, what the objects are written
        as. Default is mmCIF for channels, pores and links, and PQR for surface
        cavities and connected results, which have no mmCIF form.
    :type output_format: str

    :arg surface: Surface information returned by :func:`calcSurfaceCavities`.
        Required for ``'surface_cavities'`` and ``'connected_cavities_channels'``.
        Ignored for channels, pores and links.
    :type surface: list or None

    :returns: Paths to the mmCIF or PQR file, protein PDB file and PyMOL script.
    :rtype: tuple

    Usage:
    atoms = parsePDB('1tqn').select("protein")
    channels, surface = calcChannels(atoms)

    For channels:
    writePyMolCaviTracerScript(channels, protein)

    For pores:
    writePyMolCaviTracerScript(pores, protein, object_type='pores')

    As a PQR:
    writePyMolCaviTracerScript(channels, protein, output_format='pqr')

    Surface cavities:
    writePyMolCaviTracerScript(cavities, protein, object_type='surface_cavities',
                               surface=cavity_surface)

    Next (bash console):
    $ pymol vis_channels.py -- protein.pdb channels.cif
    $ pymol vis_channels.py -- protein.pdb channels.pqr   """

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    _requireCoords(atoms)

    aliases = {
        'channel': 'channels',
        'channels': 'channels',
        'pore': 'pores',
        'pores': 'pores',
        'link': 'links',
        'links': 'links',
        'cavity': 'surface_cavities',
        'cavities': 'surface_cavities',
        'surface_cavity': 'surface_cavities',
        'surface_cavities': 'surface_cavities',
        'connected': 'connected_cavities_channels',
        'connected_cavity_channel': 'connected_cavities_channels',
        'connected_cavities_channels': 'connected_cavities_channels'}

    object_type = aliases.get(str(object_type).lower())

    if object_type is None:
        raise ValueError("object_type must be 'channels', 'pores', 'links', "
                         "'surface_cavities' or 'connected_cavities_channels'")

    cavity_types = ('surface_cavities', 'connected_cavities_channels')
    if object_type in cavity_types:
        if surface is None or len(surface) < 5:
            raise ValueError("surface returned by calcSurfaceCavities() must be "
                             "provided for object_type='{0}'.".format(object_type))

    # Raises on anything but mmCIF or PQR, before a file is written.
    if output_format is None:
        output_format = 'pqr' if object_type in cavity_types else 'mmcif'
    mmcif = _isMmcifFormat(output_format)
    if mmcif and object_type in cavity_types:
        raise ValueError("object_type='{0}' has no mmCIF form; leave output_format "
                         "unset, or pass 'pqr'.".format(object_type))

    if objects is None:
        raise ValueError("objects cannot be None")

    if isinstance(objects, list):
        if not objects:
            raise ValueError("objects contains no CaviTracer objects to visualize")
    else:
        objects = [objects]

    output_path = Path(output_path)

    if not output_path.exists():
        output_path.mkdir(parents=True)

    if not output_path.is_dir():
        raise ValueError("output_path must be a directory")

    result_stems = {'channels': 'channels', 'pores': 'pores', 'links': 'links',
                    'surface_cavities': 'cavities',
                    'connected_cavities_channels': 'connected_cavities_channels'}
    result_stem = result_stems[object_type]
    number = 0

    while True:
        suffix = '' if number == 0 else '-{0}'.format(number)

        result_file = output_path / (
            result_stem + suffix + ('.cif' if mmcif else '.pqr'))
        protein_file = output_path / (
            'protein' + suffix + '.pdb')

        if not result_file.exists() and not protein_file.exists():
            break

        number += 1

    writePDB(str(protein_file), atoms)

    if not mmcif:
        # The records the VMD and ChimeraX writers write: every object in one
        # file, told apart by its residue number, which the viewer splits by.
        _writeViewerPQR(result_file, objects, object_type, surface, num_samples)
        written = result_file

    elif object_type == 'channels':
        written = writeChannelsCIF(result_file, objects,
                            atoms=atoms, num_samples=num_samples)

    elif object_type == 'pores':
        written = writeChannelsCIF(result_file, objects, atoms=atoms,
            object_type='pore', num_samples=num_samples)

    else:
        written = writeChannelsCIF(result_file, [], 
                                    atoms=atoms, 
                                    links=objects,
                                    num_samples=num_samples)

    if written is None:
        raise ValueError("No CaviTracer objects were written to the mmCIF file")

    written = Path(written)
    _writeVisScript(written.parent, written.name)
    script_file = written.parent / 'vis_channels.py'

    LOGGER.info("CaviTracer PyMOL files written:")
    LOGGER.info("    protein: {0}".format(protein_file))
    LOGGER.info("    results: {0}".format(written))
    LOGGER.info("    PyMOL script: {0}".format(script_file))
    LOGGER.info(
        "View the result with: pymol {0} -- {1} {2}".format(
            script_file.name, protein_file.name, written.name))

    return str(written), str(protein_file), str(script_file)
    
    
def writePyMolMultiModelScript(atoms, output_path='.', object_type='channels', flat=False):
    """Write a PyMOL trajectory visualization script.

    atoms : Atomic
        Molecular structure to display with the trajectory.

    output_path : str or pathlib.Path, optional
        Directory in which the protein PDB and PyMOL visualization script
        are written. Default is the current directory.

    object_type : str, optional
        Type of CaviTracer object stored in the PQR. Accepted values are
        'channels', 'pores', 'links', 'surface_cavities', and
        'connected_cavities_channels'. Default is 'channels'.

    flat : bool, optional
        Select the trajectory visualization layout.

        If False (default), generate the state-based viewer. Channels of the
        same rank or cluster are combined into PyMOL objects whose states
        correspond to trajectory frames.

        If True, generate the flat viewer. Each channel instance remains a
        separate PyMOL object, allowing externally generated clustering labels
        to be applied to individual channels. The generated script accepts
        one or more TSV label files and can save separate PSE sessions for
        multiple clustering solutions.

    Returns
    -------
    tuple
        protein PDB file, and PyMOL script.
    """

    import os
    from pathlib import Path

    _requireCoords(atoms)

    object_types = ('channels', 'pores', 'links', 'surface_cavities',
                    'connected_cavities_channels')
    if object_type not in object_types:
        raise ValueError('object_type must be one of {0}'.format(
            ', '.join(repr(name) for name in object_types)))

    output_path = Path(output_path)
    if not output_path.exists():
        output_path.mkdir(parents=True)
    if not output_path.is_dir():
        raise ValueError("output_path must be a directory")

    # Use the same state-based viewer as per-frame trajectory files.  The
    # generated script understands both a single multi-MODEL PQR and a glob of
    # per-frame PQRs, so there is no longer a separate visualization pathway.
    number = 0
    while True:
        suffix = '' if number == 0 else '-{0}'.format(number)
        protein_file = output_path / ('protein' + suffix + '.pdb')
        script_file = output_path / ('vis_channels_traj' + suffix + '.py')
        if flat:
            script_file = output_path / ('vis_channels_traj_flat' + suffix + '.py')
        if not (protein_file.exists() or script_file.exists()):
            break
        number += 1

    writePDB(str(protein_file), atoms, csets=atoms.getACSIndex())

    if flat:
        script = _VIS_CHANNELS_FLAT_TRAJ_SCRIPT
    else:
        script = _VIS_CHANNELS_TRAJ_SCRIPT


    with open(str(script_file), 'w') as handle:
        handle.write(script)

    LOGGER.info("CaviTracer PyMOL files written:")
    LOGGER.info("    protein: {0}".format(protein_file))
    LOGGER.info("    PyMOL script: {0}".format(script_file))
    LOGGER.info("View the result from {0} with: pymol {1} -- {2}".format(
        output_path, script_file.name, protein_file.name))
    
    return str(protein_file), str(script_file) 

def writeChimeraXCaviTracerScript(objects, atoms, object_type='channels',
                                  surface=None, output_path='.', num_samples=5):
    """Write CaviTracer results and a ChimeraX command script.

    The function saves the supplied CaviTracer objects to a combined PQR file,
    saves the molecular structure to a PDB file, and creates a ChimeraX command
    file (.cxc) that loads and visualizes both files.

    Channels, pores and chamber links are displayed as spheres using the
    radii stored in the CaviTracer PQR file. Surface cavities are displayed
    as molecular surfaces generated around the cavity marker atoms.
    Connected surface-cavity/channel results display the cavity regions as
    surfaces and the connected channels as spheres.

    The protein is displayed as a cartoon together with a transparent
    molecular surface.

    :arg objects: CaviTracer objects to visualize.
    :type objects: list or CaviTracer object

    :arg atoms: Molecular structure displayed with the CaviTracer results.
    :type atoms: :class:`.Atomic`

    :arg object_type: Type of CaviTracer result. Accepted values are
        ``'channels'``, ``'pores'``, ``'links'``,
        ``'surface_cavities'`` and ``'connected_cavities_channels'``.
        Default is ``'channels'``.
    :type object_type: str

    :arg surface: Surface information returned by
        :func:`calcSurfaceCavities`. Required for ``'surface_cavities'`` and
        ``'connected_cavities_channels'``.
    :type surface: list or None

    :arg output_path: Directory in which the PQR, PDB and ChimeraX command
        files are saved. Default is the current directory.
    :type output_path: str or pathlib.Path

    :arg num_samples: Number of samples per tetrahedron used to represent
        channels, pores and links. Default is 5.
    :type num_samples: int

    :returns: Paths to the result PQR file, protein PDB file and ChimeraX
        command file.
    :rtype: tuple  
    
    Usage:
    p = parsePDB('1tqn')
    atoms = p.select("protein")
    channels, surface = calcChannels(atoms)
    writeChimeraXCaviTracerScript(channels, atoms)
    
    Other CaviTracer objects:
    writeChimeraXCaviTracerScript(pores, atoms, object_type='pores')

    writeChimeraXCaviTracerScript(cavities, atoms, object_type='surface_cavities', 
                                    surface=surface)

    writeChimeraXCaviTracerScript(connected, atoms,
                    object_type='connected_cavities_channels', surface=cavity_surface) 
                    
    Next (bash console): 
    $ chimerax vis_channels.cxc  """

    if PY3K:
        from pathlib import Path
    else:
        from pathlib2 import Path

    _requireCoords(atoms)

    aliases = {
        'channel': 'channels',
        'channels': 'channels',
        'pore': 'pores',
        'pores': 'pores',
        'link': 'links',
        'links': 'links',
        'cavity': 'surface_cavities',
        'cavities': 'surface_cavities',
        'surface_cavity': 'surface_cavities',
        'surface_cavities': 'surface_cavities',
        'connected': 'connected_cavities_channels',
        'connected_cavity_channel': 'connected_cavities_channels',
        'connected_cavities_channels': 'connected_cavities_channels'}

    object_type = aliases.get(str(object_type).lower())

    if object_type is None:
        raise ValueError("object_type must be 'channels', 'pores', 'links', "
                        "'surface_cavities' or 'connected_cavities_channels'")

    if object_type in ('surface_cavities', 'connected_cavities_channels'):
        if surface is None or len(surface) < 5:
            raise ValueError("surface returned by calcSurfaceCavities() must be provided "
                            "for object_type='{0}'.".format(object_type))

    if objects is None:
        raise ValueError("objects cannot be None")

    if isinstance(objects, list):
        if not objects:
            raise ValueError("objects contains no CaviTracer objects to visualize")
    else:
        objects = [objects]

    output_path = Path(output_path)

    if not output_path.exists():
        output_path.mkdir(parents=True)

    if not output_path.is_dir():
        raise ValueError("output_path must be a directory")

    names = {'channels': ('channels', 'vis_channels'),
            'pores': ('pores', 'vis_pores'),
            'links': ('links', 'vis_links'),
            'surface_cavities': ('cavities', 'vis_surface_cavities'),
            'connected_cavities_channels':
                ('connected_cavities_channels',
                'vis_connected_cavities_channels')}

    result_stem, script_stem = names[object_type]
    number = 0

    while True:
        suffix = '' if number == 0 else '-{0}'.format(number)

        result_file = output_path / (result_stem + suffix + '.pqr')
        protein_file = output_path / ('protein' + suffix + '.pdb')
        script_file = output_path / (script_stem + suffix + '.cxc')

        if not (result_file.exists() or
                protein_file.exists() or
                script_file.exists()):
            break

        number += 1

    writePDB(str(protein_file), atoms)
    _writeViewerPQR(result_file, objects, object_type, surface, num_samples)

    _writeChimeraXScript(script_file, protein_file, result_file, object_type)

    LOGGER.info("CaviTracer ChimeraX files written:")
    LOGGER.info("    protein: {0}".format(protein_file))
    LOGGER.info("    results: {0}".format(result_file))
    LOGGER.info("    ChimeraX script: {0}".format(script_file))
    LOGGER.info("View the result with: chimerax {0}".format(script_file.name))

    return str(result_file), str(protein_file), str(script_file)


def _writeChimeraXScript(script_file, protein_file, result_file, object_type):
    """Write the ChimeraX command file drawing the PQR *result_file* over the
    structure *protein_file*.

    What it colours is read from the PQR itself, the residues on each of its
    chains, so a file written by :func:`writeChimeraXCaviTracerScript` and one
    already on disk open alike. So do a single structure and a multi-model file
    from :func:`mergeFramesPQR`: ChimeraX opens the latter as one submodel per
    ``MODEL``, the colours reach every one of them, and the script names each
    for its frame and shows the first, with a slider to step through the rest."""

    from pathlib import Path

    residues, frames = {}, []
    with open(str(result_file)) as handle:
        for line in handle:
            if line.startswith('MODEL'):
                frames.append(line[5:].strip())
            elif line.startswith(('ATOM', 'HETATM')):
                residues.setdefault(line[21], set()).add(int(line[22:26]))

    protein_path = Path(protein_file).resolve().as_posix().replace('"', '\\"')
    result_path = Path(result_file).resolve().as_posix().replace('"', '\\"')

    lines = [
        '# CaviTracer visualization for ChimeraX',
        '#',
        '# Open this file in ChimeraX or start ChimeraX with this file.',
        '',
        '# Start from a clean visualization.',
        'close',
        '',
        '# Display settings.',
        'graphics bgColor white',
        'camera ortho',
        '',
        # coordsets: a structure of several models, such as a trajectory's
        # frames, is read as one structure with a coordinate set per model.
        # ChimeraX otherwise opens a submodel per model and draws them all at
        # once, one on top of another.
        '# Protein.',
        'open "{0}" id #1 name Protein autoStyle false coordsets true'.format(
            protein_path),
        'hide #1 atoms,bonds',
        'show #1 cartoons',
        'color #1 lightgray target c',
        'surface #1',
        'color #1 lightgray target s',
        'transparency #1 65 target s',
        '',
        # combineSymAtoms: by default ChimeraX merges atoms of one name that sit
        # at nearly one position, as on a symmetry axis. Every sphere here is a
        # hydrogen of its channel's residue, so it would merge spheres that
        # merely lie close, and the copies of a cavity written once per channel.
        '# CaviTracer result.',
        'open "{0}" id #2 name "CaviTracer {1}" autoStyle false atomic false '
        'combineSymAtoms false'.format(
            result_path, object_type.replace('_', ' '))]

    # ChimeraX has no PQR reader of its own: the PDB reader takes the file, and
    # the radius column is where a PDB keeps the B-factor. One command maps that
    # onto the atom radius, as the line through 0.01:0.01 and 1000:1000 is the
    # identity (a waypoint must be above zero). One size command per sphere did
    # the same a command at a time, which took over a minute once there were
    # thousands of spheres.
    lines.extend(['', '# Set the CaviTracer radii, read into the B-factor.',
                  'size byattribute bfactor #2 0.01:0.01 1000:1000', ''])

    if object_type in ('channels', 'pores', 'links'):
        lines.extend([
            'hide #2 atoms,bonds',
            'show #2 atoms',
            'style #2 sphere',
            ''])

        # The palette of the PyMOL and VMD scripts, so that a channel is the
        # same colour whichever of them shows it. An object's residue number is
        # its rank plus one, in every frame.
        for resid in sorted(residues.get('T', ())):
            lines.append('color #2/T:{0} {1} target a'.format(
                resid, _hexColour(resid - 1)))

    elif object_type == 'surface_cavities':

        lines.extend(['hide #2 atoms,bonds', ''])

        for resid in sorted(residues.get('T', ())):
            lines.append('surface #2/T:{0}'.format(resid))
            lines.append('color #2/T:{0} {1} target s'.format(
                resid, _hexColour(resid - 1)))

    else:
        lines.extend(['hide #2 atoms,bonds', 'show #2/H atoms', 'style #2/H sphere', ''])

        # Cavities on chain C, channels on chain H, numbered by their own
        # indices; coloured in turn, the channels after the cavities.
        cavity_resids = sorted(residues.get('C', ()))
        channel_resids = sorted(residues.get('H', ()))

        lines.append('# Connected surface cavities.')

        for i, resid in enumerate(cavity_resids):
            lines.append('surface #2/C:{0}'.format(resid))
            lines.append('color #2/C:{0} {1} target s'.format(resid, _hexColour(i)))

        lines.append('')
        lines.append('# Connected channels.')

        offset = len(cavity_resids)

        for i, resid in enumerate(channel_resids):
            lines.append('color #2/H:{0} {1} target a'.format(
                resid, _hexColour(i + offset)))

    if len(frames) > 1:
        lines.extend(['', '# One submodel per frame: each named for its frame, the '
                      'first shown and a slider to step through the rest.'])
        for number, frame in enumerate(frames, start=1):
            lines.append('rename #2.{0} "frame {1}"'.format(number, frame))
        lines.extend(['hide #2.2-{0} models'.format(len(frames)), 'mseries slider #2'])

    lines.extend(['', '# Fit the complete system to the window.',
        'view', '', 'log text CaviTracer visualization ready.'])

    with open(str(script_file), 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def writeChimeraXMultiModelScript(pqr_file, atoms, output_path='.',
                                  object_type='channels'):
    """Write a ChimeraX command file for a PQR already on disk, such as the
    multi-model one :func:`mergeFramesPQR` writes.

    ChimeraX opens a multi-model PQR as one submodel per frame. The script names
    each for its frame, shows the first and adds a slider to step through the
    rest, and colours every object by its rank in its own frame, as
    :func:`writeChimeraXCaviTracerScript` does for one structure. A PQR of one
    structure opens just as that function would have drawn it.

    Rank is all that ties an object to its namesakes in other frames: channel 3
    of one frame and channel 3 of the next are the fourth cheapest of each, not
    one channel followed through the trajectory.

    :arg pqr_file: the PQR to draw, one structure or several frames.
    :type pqr_file: str

    :arg atoms: structure to draw the objects against. Its active coordinate
        set is written, a single backdrop for every frame.
    :type atoms: :class:`.Atomic`

    :arg output_path: directory the structure and the script are written to.
        The script opens the PQR where it is. Default is the current directory.
    :type output_path: str

    :arg object_type: what the PQR holds: ``'channels'`` (the default),
        ``'pores'``, ``'links'``, ``'surface_cavities'`` or
        ``'connected_cavities_channels'``, as in
        :func:`writeChimeraXCaviTracerScript`.
    :type object_type: str

    :returns: paths to the PQR, the structure PDB and the ChimeraX script
    :rtype: tuple

    Usage:
    mergeFramesPQR('frames', 'channels_frames.pqr')
    writeChimeraXMultiModelScript('channels_frames.pqr', atoms)

    Next (bash console):
    $ chimerax vis_channels_frames.cxc"""

    from pathlib import Path

    _requireCoords(atoms)

    object_types = ('channels', 'pores', 'links', 'surface_cavities',
                    'connected_cavities_channels')
    if object_type not in object_types:
        raise ValueError('object_type must be one of {0}'.format(
            ', '.join(repr(name) for name in object_types)))

    pqr_file = Path(pqr_file)
    if not pqr_file.is_file():
        raise ValueError('no PQR at {0}'.format(pqr_file))

    output_path = Path(output_path)
    if not output_path.exists():
        output_path.mkdir(parents=True)
    if not output_path.is_dir():
        raise ValueError("output_path must be a directory")

    number = 0
    while True:
        suffix = '' if number == 0 else '-{0}'.format(number)
        protein_file = output_path / ('protein' + suffix + '.pdb')
        script_file = output_path / ('vis_' + pqr_file.stem + suffix + '.cxc')
        if not (protein_file.exists() or script_file.exists()):
            break
        number += 1

    # The active frame only: written whole, a trajectory's structure would come
    # out as a model per frame, all of them drawn at once behind the channels.
    writePDB(str(protein_file), atoms, csets=atoms.getACSIndex())
    _writeChimeraXScript(script_file, protein_file, pqr_file, object_type)

    LOGGER.info("CaviTracer ChimeraX files written:")
    LOGGER.info("    protein: {0}".format(protein_file))
    LOGGER.info("    results: {0}".format(pqr_file))
    LOGGER.info("    ChimeraX script: {0}".format(script_file))
    LOGGER.info("View the result with: chimerax {0}".format(script_file.name))

    return str(pqr_file), str(protein_file), str(script_file)