"""This module defines functions and constants related to sequences."""

import re

__all__ = ['MATCH_SCORE', 'MISMATCH_SCORE', 'GAP_PENALTY',
           'GAP_EXT_PENALTY', 'ALIGNMENT_METHOD', 'splitSeqLabel']

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