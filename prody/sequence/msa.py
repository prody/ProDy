# -*- coding: utf-8 -*-
"""This module defines MSA analysis functions."""

from numbers import Integral

from numpy import all, zeros, dtype, array, char, cumsum, ceil, reshape
from numpy import where, sort, concatenate, vstack, isscalar, chararray

from prody import LOGGER, PY3K
from prody.atomic import Atomic
from prody.utilities import toChararray, pystr, splitSeqLabel
from .sequence import Sequence

import sys

__all__ = ['MSA', 'refineMSA', 'mergeMSA', 'specMergeMSA',]

class MSA(object):

    """Store and manipulate multiple sequence alignments."""

    def __init__(self, msa, title='Unknown', labels=None, **kwargs):
        """*msa* must be a 2D Numpy character array. *labels* is a list of
        sequence labels (or titles).  *mapping* should map label or part of
        label to sequence index in *msa* array. If *mapping* is not given,
        one will be build from *labels*."""

        self._aligned = aligned = kwargs.get('aligned', True)
        msa = toChararray(msa, aligned)
        numseq = msa.shape[0]

        if labels and len(labels) != numseq:
            raise ValueError('len(labels) must be equal to number of '
                             'sequences')
        
        if labels is None:
            labels = [str(i+1) for i in range(numseq)]

        if PY3K:
            for i, label in enumerate(labels):
                if not isinstance(label, str):
                    labels[i] = label.decode()
                    
        self._labels = labels
        
        mapping = kwargs.get('mapping')
        self._map(mapping)
        self._msa = msa
        self._title = str(title) or 'Unknown'
        self._split = bool(kwargs.get('split', True))

    def _map(self, mapping=None):

        labels = self._labels
        if mapping is not None:
            try:
                mapping['isdict']
            except KeyError:
                pass
            except Exception:
                raise TypeError('mapping must be a dictionary')
            
            for key, value in mapping.items():
                values = [value] if isscalar(value) else value
                for i in values:
                    if not key in labels[i]:
                        labels[i] = key
        
        self._mapping = mapping = {}
        for index, label in enumerate(labels):
            label = splitSeqLabel(label)[0]
            try:
                value = mapping[label]
            except KeyError:
                mapping[label] = index
            else:
                try:
                    value.append(index)
                except AttributeError:
                    mapping[label] = [value, index]
        return mapping

    def __str__(self):

        return 'MSA ' + self._title

    def __repr__(self):

        if self._aligned:
            return ('<MSA: {0} ({1} sequences, {2} residues)>'
                    ).format(self._title, self.numSequences(),
                             self.numResidues())
        else:
            return ('<MSA: {0} ({1} sequences, not aligned)>'
                    ).format(self._title, self.numSequences())

    def __len__(self):

        return len(self._msa)

    def __add__(self, other):
        """Concatenate two MSAs."""

        A = self.getArray()
        B = other.getArray()
        try:
            AB = vstack((A, B))
        except:
            raise ValueError('failed to add {0} and {1} together: shapes do not match'
                             .format(repr(self), repr(other)))
        aligned = self._aligned or other._aligned

        labels = list(self._labels)
        labels.extend(other._labels)

        split = self._split or other._split

        msa = MSA(AB, title='%s + %s'%(self.getTitle(), other.getTitle()), 
                  aligned=aligned, labels=labels, split=split)
        return msa

    def __getitem__(self, index):

        if isinstance(index, Integral):
            return Sequence(self, index)

        if isinstance(index, str):
            try:
                rows = self._mapping[index]
            except KeyError:
                raise KeyError('label {0} is not mapped to a sequence'
                               .format(index))
            else:
                msa = self._msa[rows]
                if msa.ndim == 1:
                    return Sequence(self, rows)
                else:
                    if msa.base is not None:
                        msa = msa.copy()
                    labels = self._labels
                    return MSA(msa, title='{0}[{1}]'.format(self._title,
                               index), labels=[labels[i] for i in rows])
        elif isinstance(index, tuple):
            if len(index) == 1:
                return self[index[0]]
            elif len(index) == 2:
                rows, cols = index
            else:
                raise IndexError('invalid index: ' + str(index))
        else:
            rows, cols = index, None

        # handle list of labels
        if isinstance(rows, list):
            rows = self.getIndex(rows) or rows
        elif isinstance(rows, int):
            return Sequence(self._msa[rows, cols].tobytes(),
                            self._labels[rows])
        elif isinstance(rows, str):
            try:
                rows = self._mapping[rows]
            except KeyError:
                raise KeyError('label {0} is not mapped to a sequence'
                               .format(index))
            else:
                if isinstance(rows, int):
                    return Sequence(self._msa[rows, cols].tobytes(),
                                    self._labels[rows])

        if cols is None:
            msa = self._msa[rows]
        else:
            try:
                msa = self._msa[rows, cols]
            except TypeError:
                raise IndexError('invalid index: ' + str(index))

        try:
            lbls = list(array(self._labels)[rows])
        except TypeError:
            labels = self._labels
            lbls = [labels[i] for i in rows]
        finally:
            if not isinstance(lbls, list):
                lbls = [lbls]

        if msa.ndim == 0:
            msa = msa.reshape((1, 1))
        elif msa.ndim == 1:
            msa = msa.reshape((1, len(msa)))
        if msa.base is not None:
            msa = msa.copy()

        return MSA(msa=msa, title=self._title, labels=lbls,
                   aligned=self._aligned)

    def __iter__(self):

        for i in range(len(self._msa)):
            yield Sequence(self, i)

    def __contains__(self, key):

        try:
            return key in self._mapping
        except Exception:
            pass
        return False

    def __eq__(self, other):

        try:
            other = other._getArray()
        except AttributeError:
            return False

        try:
            return all(other == self._msa)
        except Exception:
            pass
        return False

    def _getSplit(self):

        return self._split

    def _setSplit(self, split):

        self._split = bool(split)

    split = property(_getSplit, _setSplit,
                     doc='Return split label when iterating or indexing.')

    def extend(self, other):
        """Adds *other* to this MSA."""

        A = self.getArray()
        if isinstance(other, MSA):
            B = other.getArray()
            otherlabels = other._labels
        elif isinstance(other, Sequence):
            B = other.getArray()
            B = reshape(B, (1, B.shape[0]))
            otherlabels = other.getLabel(full=True)
        else:
            try:
                B = toChararray(other, aligned=False)
                numA, numB = A.shape[0], B.shape[0]
                otherlabels = [str(i+1) for i in range(numA, numA+numB)]
            except:
                raise ValueError('failed to add {1} to {0}'
                             .format(repr(self), repr(other)))
        try:
            AB = vstack((A, B))
        except:
            raise ValueError('failed to add {1} to {0}: shapes do not match'
                             .format(repr(self), repr(other)))

        labels = list(self._labels)
        labels.extend(otherlabels)

        self._msa = AB
        self._labels = labels
        self._map()

    def isAligned(self):
        """Returns **True** if MSA is aligned."""

        return self._aligned

    def numSequences(self):
        """Returns number of sequences."""

        return self._msa.shape[0]

    def numResidues(self):
        """Returns number of residues (or columns in the MSA), if MSA is
        aligned."""

        if self._aligned:
            return self._msa.shape[1]

    def numIndexed(self):
        """Returns number of sequences that are indexed using the identifier
        part or all of their labels.  The return value should be equal to
        number of sequences."""

        count = len(self._mapping)
        if len(self._msa) == count:
            return count
        else:
            count = len(self._mapping)
            for val in self._mapping.values():
                try:
                    count += len(val) - 1
                except TypeError:
                    pass
            return count

    def getTitle(self):
        """Returns title of the instance."""

        return self._title

    def setTitle(self, title):
        """Set title of the instance."""

        self._title = str(title)

    def getLabel(self, index, full=False):
        """Returns label of the sequence at given *index*.  Residue numbers will
        be removed from the sequence label, unless *full* is **True**."""

        index = self._mapping.get(index, index)
        if full:
            return self._labels[index]
        else:
            return splitSeqLabel(self._labels[index])[0]

    def getResnums(self, index):
        """Returns starting and ending residue numbers (:term:`resnum`) for the
        sequence at given *index*."""

        index = self._mapping.get(index, index)
        return splitSeqLabel(self._labels[index])[1:]

    def getArray(self):
        """Returns a copy of the MSA character array."""

        return self._msa.copy()

    def _getArray(self):
        """Returns MSA character array."""

        return self._msa

    def getIndex(self, label):
        """Returns index of the sequence that *label* maps onto.  If *label*
        maps onto multiple sequences or *label* is a list of labels, a list
        of indices is returned.  If an index for a label is not found,
        return **None**."""

        try:
            index = self._mapping[label]
        except KeyError:
            try:
                return [v for k, v in self._mapping.items() if label in k][0]
            except:
                return None
        except TypeError:
            mapping = self._mapping
            indices = []
            append, extend = indices.append, indices.extend
            for key in label:
                try:
                    index = mapping[key]
                except KeyError:
                    return None
                try:
                    extend(index)
                except TypeError:
                    append(index)
            return indices
        else:
            try:
                return list(index)
            except TypeError:
                return index

    def iterLabels(self, full=False):
        """Yield sequence labels.  By default the part of the label used for
        indexing sequences is yielded."""

        if full:
            for label in self._labels:
                yield label
        else:
            for label in self._labels:
                yield splitSeqLabel(label)[0]

    def getLabels(self, full=False):
        """Returns all labels"""

        labels = [label for label in self.iterLabels(full=full)]
        return labels

    def countLabel(self, label):
        """Returns the number of sequences that *label* maps onto."""

        try:
            return len(self._mapping[label])
        except KeyError:
            return 0
        except TypeError:
            return 1

def refineMSA(msa, index=None, label=None, rowocc=None, seqid=None, colocc=None, **kwargs):
    """Refine *msa* by removing sequences (rows) and residues (columns) that
    contain gaps.

    :arg msa: multiple sequence alignment
    :type msa: :class:`.MSA`

    :arg index: remove columns that are gaps in the sequence with that index
    :type index: int

    :arg label: remove columns that are gaps in the sequence matching label,
        ``msa.getIndex(label)`` must return a sequence index, a PDB identifier
        is also acceptable
    :type label: str

    :arg rowocc: row occupancy, sequences with less occupancy will be
        removed after *label* refinement is applied
    :type rowocc: float

    :arg seqid: keep unique sequences at specified sequence identity level,
        unique sequences are identified using :func:`.uniqueSequences`
    :type seqid: float

    :arg colocc: column occupancy, residue positions with less occupancy
        will be removed after other refinements are applied
    :type colocc: float

    :arg keep: keep columns corresponding to residues not resolved in the PDB
        structure, default is **False**, applies when *label* is a PDB
        identifier
    :arg type: bool

    For Pfam MSA data, *label* is UniProt entry name for the protein.  You may
    also use PDB structure and chain identifiers, e.g. ``'1p38'`` or
    ``'1p38A'``, for *label* argument and UniProt entry names will be parsed
    using :func:`.parsePDBHeader` function (see also :class:`.Polymer` and
    :class:`.DBRef`).

    The order of refinements are applied in the order of arguments.  If *label*
    and *unique* is specified, sequence matching *label* will
    be kept in the refined :class:`.MSA` although it may be similar to some
    other sequence."""

    # if msa is a char array, it will be refined but label won't work
    try:
        ndim, dtype_ = msa.ndim, msa.dtype
    except AttributeError:
        try:
            arr = msa._getArray()
        except AttributeError:
            raise TypeError('msa must be a character array or an MSA instance')
        ndim, dtype_ = arr.ndim, arr.dtype
    else:
        arr, msa = msa, None

    if dtype('|S1') != dtype_:
        raise ValueError('msa must be a character array or an MSA instance')
    if ndim != 2:
        raise ValueError('msa must be a 2D array or an MSA instance')

    title = []
    cols = None

    if index is not None:
        before = arr.shape[1]
        LOGGER.timeit('_refine')
        cols = char.isalpha(arr[index]).nonzero()[0]
        arr = arr.take(cols, 1)
        title.append('index=' + str(index))
        LOGGER.report('Index refinement reduced number of columns from {0} to '
                      '{1} in %.2fs.'.format(before, arr.shape[1]), '_refine')

    if label is not None:
        if index is not None:
            LOGGER.info('An index was provided so the label will be ignored.')

        else:
            before = arr.shape[1]
            LOGGER.timeit('_refine')
            try:
                upper, lower = label.upper(), label.lower()
            except AttributeError:
                raise TypeError('label must be a string')

            if msa is None:
                raise TypeError('msa must be an MSA instance, '
                                'label cannot be used')

            index = msa.getIndex(label)
            if index is None:
                index = msa.getIndex(upper)
            if index is None:
                index = msa.getIndex(lower)

            chain = None
            if index is None and (len(label) == 4 or len(label) == 5):
                from prody import parsePDB
                try:
                    structure, header = parsePDB(label[:4], header=True)
                except Exception as err:
                    raise IOError('failed to parse header for {0} ({1})'
                                  .format(label[:4], str(err)))

                chid = label[4:].upper()
                for poly in header['polymers']:
                    if chid and poly.chid != chid:
                        continue
                    for dbref in poly.dbrefs:
                        if index is None:
                            index = msa.getIndex(dbref.idcode)
                            if index is not None:
                                LOGGER.info('{0} idcode {1} for {2}{3} '
                                            'is found in chain {4}.'.format(
                                            dbref.database, dbref.idcode,
                                            label[:4], poly.chid, str(msa)))
                                break
                        if index is None:
                            index = msa.getIndex(dbref.accession)
                            if index is not None:
                                LOGGER.info('{0} accession {1} for {2}{3} '
                                            'is found in chain {4}.'.format(
                                            dbref.database, dbref.accession,
                                            label[:4], poly.chid, str(msa)))
                                break
                if index is not None:
                    chain = structure[poly.chid]
                    resnums = chain.ca.getResnums()
    
            if index is None:
                raise ValueError('label is not in msa, or msa is not indexed')
            try:
                len(index)
            except TypeError:
                pass
            else:
                raise ValueError('label {0} maps onto multiple sequences, '
                                 'so cannot be used for refinement'.format(label))

            title.append('label=' + label)
            cols = char.isalpha(arr[index]).nonzero()[0]
            arr = arr.take(cols, 1)
            LOGGER.report('Label refinement reduced number of columns from {0} to '
                          '{1} in %.2fs.'.format(before, arr.shape[1]), '_refine')

            if chain is not None and not kwargs.get('keep', False):
                before = arr.shape[1]
                LOGGER.timeit('_refine')
                
                from prody.utilities import MATCH_SCORE, MISMATCH_SCORE, alignBioPairwise
                from prody.utilities import GAP_PENALTY, GAP_EXT_PENALTY, ALIGNMENT_METHOD

                chseq = chain.getSequence()
                algn = alignBioPairwise(pystr(arr[index].tobytes().upper()), pystr(chseq),
                                        "local",
                                        MATCH_SCORE, MISMATCH_SCORE,
                                        GAP_PENALTY, GAP_EXT_PENALTY,
                                        max_alignments=1)
                torf = []
                for s, c in zip(*algn[0][:2]):
                    if s == '-':
                        continue
                    elif c != '-':
                        torf.append(True)
                    else:
                        torf.append(False)
                torf = array(torf)
                tsum = torf.sum()
                assert tsum <= before, 'problem in mapping sequence to structure'
                if tsum < before:
                    arr = arr.take(torf.nonzero()[0], 1)
                    resnums = resnums.take(torf.nonzero()[0]-torf.nonzero()[0][0]+1)
                    LOGGER.report('Structure refinement reduced number of '
                                  'columns from {0} to {1} in %.2fs.'
                                  .format(before, arr.shape[1]), '_refine')
                else:
                    LOGGER.debug('All residues in the sequence are contained in '
                                 'PDB structure {0}.'.format(label))

                labels = msa._labels
                labels[index] = splitSeqLabel(labels[index])[0] + '/' + str(resnums[0]) + '-' + str(resnums[-1])

    from .analysis import calcMSAOccupancy, uniqueSequences

    rows = None
    if rowocc is not None:
        before = arr.shape[0]
        LOGGER.timeit('_refine')
        try:
            rowocc = float(rowocc)
        except Exception as err:
            raise TypeError('rowocc must be a float ({0})'.format(str(err)))
        assert 0. <= rowocc <= 1., 'rowocc must be between 0 and 1'

        rows = calcMSAOccupancy(arr, 'row') >= rowocc
        if index is not None:
            index = rows[:index].sum()
        rows = (rows).nonzero()[0]
        arr = arr[rows]
        title.append('rowocc>=' + str(rowocc))
        LOGGER.report('Row occupancy refinement reduced number of rows from '
                      '{0} to {1} in %.2fs.'.format(before, arr.shape[0]),
                      '_refine')

    if seqid is not None:
        before = arr.shape[0]
        LOGGER.timeit('_refine')
        unique = uniqueSequences(arr, seqid)
        if index is not None:
            unique[index] = True
        unique = unique.nonzero()[0]
        arr = arr[unique]
        title.append('seqid>=' + str(seqid))
        if rows is not None:
            rows = rows[unique]
        else:
            rows = unique
        LOGGER.report('Sequence identity refinement reduced number of rows '
                      'from {0} to {1} in %.2fs.'.format(before, arr.shape[0]),
                      '_refine')

    if colocc is not None:
        before = arr.shape[1]
        LOGGER.timeit('_refine')
        try:
            colocc = float(colocc)
        except Exception as err:
            raise TypeError('colocc must be a float ({0})'.format(str(err)))
        assert 0. <= colocc <= 1., 'colocc must be between 0 and 1'

        cols = (calcMSAOccupancy(arr, 'col') >= colocc).nonzero()[0]
        arr = arr.take(cols, 1)
        title.append('colocc>=' + str(colocc))
        LOGGER.report('Column occupancy refinement reduced number of columns '
                      'from {0} to {1} in %.2fs.'.format(before, arr.shape[1]),
                      '_refine')

    if not title:
        raise ValueError('label, index, seqid, rowocc, colocc all cannot be None')

    # depending on slicing of rows, arr may not have it's own memory
    if arr.base is not None:
        arr = arr.copy()

    if msa is None:
        return arr
    else:
        if rows is None:
            from copy import copy
            labels = copy(msa._labels)
        else:
            labels = msa._labels
            labels = [labels[i] for i in rows]
        return MSA(arr, title=msa.getTitle() + ' refined ({0})'
                   .format(', '.join(title)), labels=labels)


def mergeMSA(*msa, **kwargs):
    """Returns an :class:`.MSA` obtained from merging parts of the sequences
    of proteins present in multiple *msa* instances.  Sequences are matched
    based on protein identifiers found in the sequence labels.  Order of
    sequences in the merged MSA will follow the order of sequences in the
    first *msa* instance.  Note that protein identifiers that map to multiple
    sequences will be excluded.
    
    MSAs with different identifiers can be merged with the *ignore_ids* kwarg.
    This only works when all MSAs have the same number of sequences.
    
    :arg msa: a set of :class:`.MSA` objects to be analysed
    :type msa: list, tuple, :class:`~numpy.ndarray`
    
    :arg ignore_ids: where to ignore identifiers instead of matching them
        Default is False
    :type ignore_ids: bool
    """
    ignore_ids = kwargs.pop('ignore_ids', False)

    if len(msa) <= 1:
        raise ValueError('more than one msa instances are needed')

    try:
        arrs = [m._getArray() for m in msa]
        sets = []
        for m in msa:
            aset = set([])
            count = m.countLabel
            for label in m.iterLabels():
                if count(label) == 1:
                    aset.add(label)
            sets.append(aset)
    except AttributeError:
        raise TypeError('all msa arguments must be MSA instances')

    num_seqs = [len(aset) for aset in sets]
    sets = iter(sets)
    common = next(sets)
    if not all(array(num_seqs) == num_seqs[0]):
        LOGGER.warn("Can only ignore identifiers if all MSAs have the same number of sequences")
        return None
    
    for aset in sets:
        if not ignore_ids:
            common = common.intersection(aset)
        else:
            [common.add(item) for item in aset]
    if not common:
        LOGGER.warn("No sequences with common identifiers to merge")
        return None

    lens = [m.numResidues() for m in msa]
    rngs = [0]
    rngs.extend(cumsum(lens))
    rngs = [(start, end) for start, end in zip(rngs[:-1], rngs[1:])]

    if not ignore_ids:
        idx_arr_rng = list(zip([m.getIndex for m in msa], arrs, rngs))
        merger = zeros((len(common), sum(lens)), '|S1')
        labels = []
        for index, label in enumerate(common):
            if label not in common:
                continue
            for idx, arr, (start, end) in idx_arr_rng:
                merger[index, start:end] = arr[idx(label)]
            labels.append(label)
    else:
        idx_arr_rng = list(zip([i for i in range(len(msa))], arrs, rngs))
        merger = zeros((len(num_seqs), sum(lens)), '|S1')
        labels = msa[0].getLabels()
        for index, label in enumerate(msa[0].getLabels()):
            for i, arr, (start, end) in idx_arr_rng:
                    merger[index, start:end] = arr[index]
            labels[index] += '_' + msa[i].getLabels()[index]
            
    merger = MSA(merger, labels=labels,
                 title=' + '.join([m.getTitle() for m in msa]))
    return merger

def specMergeMSA(*msa, **kwargs):
    """Returns an :class:`.MSA` obtained from merging parts of the sequences
    of proteins present in multiple *msa* instances. Sequences are matched
    based on species section of protein identifiers found in the sequence labels.  
    Order of sequences in the merged MSA will follow the order of sequences in the
    first *msa* instance. Note that protein identifiers that map to multiple
    sequences will be excluded."""

    if len(msa) <= 1:
        raise ValueError('more than one msa instances are needed')
    lbl={}
    try:
        arrs = [m._getArray() for m in msa]
        sets = []
        labells = []
        for m in msa:
            aset = set([])
            labell = {}
            count = m.countLabel
            for label in m.iterLabels():
                lbl[label]=label.rsplit('_')[1]
                if count(label) == 1 and lbl[label] not in aset:
                    aset.add(lbl[label])
                    labell[lbl[label]]=label
            sets.append(aset)
            labells.append(labell)
    except AttributeError:
        raise TypeError('all msa arguments must be MSA instances')  
    sets = iter(sets)
    common = next(sets)
    for aset in sets:
        common = common.intersection(aset)
    if not common:
        return None
    lens = [m.numResidues() for m in msa]
    rngs = [0]
    rngs.extend(cumsum(lens))
    rngs = [(start, end) for start, end in zip(rngs[:-1], rngs[1:])]
    
    idx_arr_rng = list(zip([m.getIndex for m in msa], arrs, rngs))

    merger = zeros((len(common), sum(lens)), '|S1')
    labels = []
    for index, lbl in enumerate(common):
        for idx, arr, (start, end) in idx_arr_rng:
            merger[index, 0:start]=list(str(msa[0][msa[0].getIndex(labells[0][lbl])]))
            merger[index, start:end]=list(str(msa[1][msa[1].getIndex(labells[1][lbl])]))
        label = labells[0][lbl]

        labels.append(label)
    merger = MSA(merger, labels=labels,
                 title=' + '.join([m.getTitle() for m in msa]))
    return merger
