"""This module defines functions and constants related to sequences."""

import re

__all__ = ['MATCH_SCORE', 'MISMATCH_SCORE', 'GAP_PENALTY',
           'GAP_EXT_PENALTY', 'ALIGNMENT_METHOD', 'splitSeqLabel',
           'alignBioPairwise']

MATCH_SCORE = 1.0
MISMATCH_SCORE = 0.0
GAP_PENALTY = -1.
GAP_EXT_PENALTY = -0.1
ALIGNMENT_METHOD = 'local'

SPLITLABEL = re.compile(r'[/-]+').split


def splitSeqLabel(label):
    """Returns label, starting residue number, and ending residue number parsed
    from sequence label."""

    try:
        if label.strip() == '':
            raise Exception
        idcode, start, end = SPLITLABEL(label)
    except Exception:
        return label, None, None
    else:
        try:
            return idcode, int(start), int(end)
        except Exception:
            return label, None, None


def alignBioPairwise(a_sequence, b_sequence,
                     ALIGNMENT_METHOD=ALIGNMENT_METHOD,
                     MATCH_SCORE=MATCH_SCORE, MISMATCH_SCORE=MISMATCH_SCORE,
                     GAP_PENALTY=GAP_PENALTY, GAP_EXT_PENALTY=GAP_EXT_PENALTY,
                     max_alignments=1):
    """
    Wrapper function to align two sequences using Biopython to support deprecation
    and associated warnings.

    It first attempts to use :mod:`Bio.Align.PairwiseAligner` and formats the output
    to match :mod:`Bio.pairwise2` methods and if it doesn't find it then it defaults
    back to :mod:`Bio.pairwise2` methods.

    :arg a_sequence: first sequence to align
    :type a_sequence: str

    :arg b_sequence: second sequence to align
    :type b_sequence: str

    :arg ALIGNMENT_METHOD: method for pairwise2 alignment
        Possible values are ``"local"`` and ``"global"``
    :type ALIGNMENT_METHOD: str

    :arg MATCH_SCORE: a positive integer, used to reward finding a match
    :type MATCH_SCORE: int

    :arg MISMATCH_SCORE: a negative integer, used to penalise finding a mismatch
    :type MISMATCH_SCORE: int

    :arg GAP_PENALTY: a negative integer, used to penalise opening a gap
    :type GAP_PENALTY: int

    :arg GAP_EXT_PENALTY: a negative integer, used to penalise extending a gap
    :type GAP_EXT_PENALTY: int

    :arg max_alignments: maximum number of alignments to extract
    :type max_alignments: int
    """
    try:
        from Bio.Align import PairwiseAligner
        
        try:
            aligner = PairwiseAligner()
            aligner.mode = ALIGNMENT_METHOD
            aligner.match_score = MATCH_SCORE
            aligner.mismatch_score = MISMATCH_SCORE
            aligner.open_internal_gap_score = GAP_PENALTY
            aligner.extend_internal_gap_score = GAP_EXT_PENALTY
            alns = aligner.align(a_sequence, b_sequence)
        except AttributeError:
            try:
                aligner = PairwiseAligner()
                aligner.mode = ALIGNMENT_METHOD
                aligner.match_score = MATCH_SCORE
                aligner.mismatch_score = MISMATCH_SCORE
                aligner.open_gap_score = GAP_PENALTY
                aligner.extend_gap_score = GAP_EXT_PENALTY
                alns = aligner.align(a_sequence, b_sequence)
            except AttributeError:
                aligner = PairwiseAligner()
                aligner.mode = ALIGNMENT_METHOD
                aligner.match_score = MATCH_SCORE
                aligner.mismatch_score = MISMATCH_SCORE
                aligner.open_internal_insertion_score = GAP_PENALTY
                aligner.extend_internal_insertion_score = GAP_EXT_PENALTY
                alns = aligner.align(a_sequence, b_sequence)

        results = []

        for i, aln in enumerate(alns):
            if i == max_alignments:
                break

            split_aln = aln.format().split('\n')

            if split_aln[0].startswith('target'):
                # new biopython
                import re

                m1 = re.search('[0-9]+', split_aln[0])
                begin = m1.group()

                m2 = re.search('[0-9]+', split_aln[-4])
                m3 = re.search('[0-9]+', split_aln[-4][m2.end():])
                end = m3.group()

                row_1 = ''.join([line[m1.end()+1:] for line in split_aln[0:-4:4]])
                row_1 += split_aln[-4][m2.end()+1:m2.end()+m3.start()-1]

                row_2 = ''.join([line[m1.end()+1:] for line in split_aln[2:-4:4]])
                row_2 += split_aln[-2][m2.end()+1:m2.end()+m3.start()-1]
            else:
                begin = split_aln[1].find('|')
                end = len(split_aln[1])-1

                row_1 = split_aln[0].replace(" ", "-")
                row_2 = split_aln[2].replace(" ", "-")

                if len(row_1) < len(row_2):
                    row_1 += "-"*(len(row_2)-len(row_1))
                elif len(row_2) < len(row_1):
                    row_2 += "-"*(len(row_1)-len(row_2))

            results.append((row_1, row_2, aln.score, begin, end))

        return results
        
    except (ImportError, AttributeError):
        from Bio import pairwise2

        if ALIGNMENT_METHOD == "local":
            return pairwise2.align.localms(a_sequence, b_sequence,
                                            MATCH_SCORE, MISMATCH_SCORE,
                                            GAP_PENALTY, GAP_EXT_PENALTY,
                                            one_alignment_only=1)
        elif ALIGNMENT_METHOD == "global":
            return pairwise2.align.globalms(a_sequence, b_sequence,
                                            MATCH_SCORE, MISMATCH_SCORE,
                                            GAP_PENALTY, GAP_EXT_PENALTY,
                                            one_alignment_only=1)
        else:
            raise ValueError("method should be local or global")
