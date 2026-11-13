"""Reading back alignments written in each supported MSA file format, and
reading SELEX, Stockholm, CLUSTAL and PIR files written by other programs."""

import gzip
import os
import re
from io import StringIO

import numpy as np
import pytest

from prody import MSA, MSAFile, parseMSA, writeMSA
from prody.tests.datafiles import pathDatafile

ROUNDTRIP_LABELS = [
    'sp|P01.2|ABC_HUMAN/1-12',
    'short',
    '/query',
    '/lead/2-9',
    'sp/inner/4-11',
    '.dot',
    ':colon',
    '*star',
    '|pipe',
    '%pct',
    'x/1/2',
    'foo/bar-1-2',
    'L' + 'abcdefghij' * 50 + '/7-19',
    'chain:A*model.2',
    'exactly_16_chars',
    'seq31_abcdefghijklmnopqrstuvwxy',
    'a_very_long_label_exceeding_the_selex_label_field/5-16',
    'Q9XYZ1',
]

ROUNDTRIP_RESIDUES = 'ACDEFGHIKLMNPQRSTVWY-.'

ROUNDTRIP_EXTENSIONS = ['.slx', '.sth', '.aln', '.ali',
                        '.slx.gz', '.sth.gz', '.aln.gz', '.ali.gz']

ROUNDTRIP_FORMATS = {'.slx': 'selex', '.sth': 'stockholm',
                     '.aln': 'clustal', '.ali': 'pir'}

EXPECTED_FIXTURES = {
    'blocks.sth': [
        ('P12345.1/3-24', 'ACDEFGHIKLMNPQRSTVWY'),
        ('short', 'AC-EFG.IKLMN-QRSTV-Y'),
        ('seq:2*x|y', 'ACDEFGHIK-mnpqrsTVWY'),
    ],
    'blocks.slx': [
        ('Q9XYZ1', 'ACDEFGHIKLMNPQRS'),
        ('long_label_with.dots_and:colons*star', 'AC-EFGHIK-MNPQ-S'),
        ('x/1-12', 'acdefghiklmnpqrs'),
    ],
    'clustalw2.aln': [
        ('sp|P01.2|ABC_HUMAN',
         'ACDEFGHIKLMNPQRSTVWY' * 6 + 'MNPQ'),
        ('tr:Q8.1*',
         'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKL-NPQRSTVWY'
         '--DEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY'
         'MN--'),
        ('seq.3',
         'ACDEFGHIKLMNPQRSTVWYACDEF-HIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY'
         + 'ACDEFGHIKLMNPQRSTVWY' * 3 + 'MNAC'),
    ],
    'clustalw1.aln': [
        ('1abc.A', 'ACDEFGHIKLMNPQ'),
        ('model:2', 'AC-EFGHIKLMN-Q'),
        ('target*', 'ACDEF-HIKLMNPQ'),
    ],
    'two_alignments.sth': [
        ('/query', 'ACDEFGHIKLMNPQR'),
        ('P99999.2/1-15', 'AC-EFGHIK-mn-qr'),
    ],
    'modeller.ali': [
        ('1abc', 'ACDEFGHIKLMNPQRSTVWY' * 3 + 'ACDEFGHIKL'),
        ('target.2', 'AC-EFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEF'
                     'SHIKLMNPQRSTVWYACDEFGHIK-'),
        ('model_B', 'SCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM'
                    'NPQRSTVWYAC-EFGHIKL'),
        ('custom.4', 'sDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM'),
    ],
}


def roundtrip_sequences(ncols):
    """Returns one aligned sequence of *ncols* columns per label."""

    seqs = []
    for i in range(len(ROUNDTRIP_LABELS)):
        seq = ''.join(ROUNDTRIP_RESIDUES[(5 * i + 3 * j) %
                                         len(ROUNDTRIP_RESIDUES)]
                      for j in range(ncols))
        seqs.append(seq.lower() if i == 1 else seq)
    return seqs


def roundtrip_msa(ncols):

    seqs = roundtrip_sequences(ncols)
    return MSA(np.array([list(seq) for seq in seqs], '|S1'),
               title='roundtrip', labels=list(ROUNDTRIP_LABELS))


def roundtrip_expected(ncols):

    return list(zip(ROUNDTRIP_LABELS, roundtrip_sequences(ncols)))


def msa_pairs(msa):
    """Returns full label and sequence pairs of a parsed MSA."""

    return [(msa.getLabel(i, full=True), str(msa[i]))
            for i in range(msa.numSequences())]


def msafile_pairs(msafile):
    """Returns full label and sequence pairs read from an MSAFile."""

    return [(seq.getLabel(True), str(seq)) for seq in msafile]


def formats_fixture(name, tmp_path, compressed=False):
    """Returns the path of a fixture file, or of a gzipped copy of it."""

    path = pathDatafile('msa_formats_' + name)
    if not compressed:
        return path
    copy = str(tmp_path / (name + '.gz'))
    with open(path, 'rb') as inp, gzip.open(copy, 'wb') as out:
        out.write(inp.read())
    return copy


def read_text(path):

    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as inp:
            return inp.read()
    with open(path) as inp:
        return inp.read()


def join_labelled_lines(lines):
    """Returns label and sequence pairs from lines with a label, whitespace
    and a sequence piece, joining pieces by label in order of appearance."""

    labels = []
    pieces = {}
    for line in lines:
        items = line.split()
        assert len(items) == 2, 'expected label and sequence: ' + repr(line)
        label, piece = items
        if label not in pieces:
            labels.append(label)
            pieces[label] = []
        pieces[label].append(piece)
    return [(label, ''.join(pieces[label])) for label in labels]


@pytest.mark.parametrize('ncols', [1, 12, 60, 150, 12000])
@pytest.mark.parametrize('ext', ROUNDTRIP_EXTENSIONS)
def test_writemsa_then_parsemsa(tmp_path, ext, ncols):

    filename = writeMSA(str(tmp_path / ('roundtrip' + ext)),
                        roundtrip_msa(ncols))
    assert msa_pairs(parseMSA(filename)) == roundtrip_expected(ncols)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('ext', ['.slx', '.sth', '.aln', '.ali'])
def test_writemsa_then_msafile(tmp_path, ext, ncols):

    filename = writeMSA(str(tmp_path / ('roundtrip' + ext)),
                        roundtrip_msa(ncols))
    assert msafile_pairs(MSAFile(filename)) == roundtrip_expected(ncols)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('ext', ROUNDTRIP_EXTENSIONS)
def test_msafile_writer_then_parsemsa(tmp_path, ext, ncols):

    filename = str(tmp_path / ('roundtrip' + ext))
    with MSAFile(filename, 'w') as out:
        for seq in roundtrip_msa(ncols):
            out.write(seq)
    assert msa_pairs(parseMSA(filename)) == roundtrip_expected(ncols)
    assert msafile_pairs(MSAFile(filename)) == roundtrip_expected(ncols)
    if ROUNDTRIP_FORMATS[ext.replace('.gz', '')] == 'stockholm':
        lines = [line for line in read_text(filename).splitlines()
                 if line.strip()]
        assert lines[-1] == '//'


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('format', ['selex', 'clustal', 'pir'])
def test_msafile_stream_roundtrip(format, ncols):

    stream = StringIO()
    with MSAFile(stream, 'w', format=format) as out:
        for seq in roundtrip_msa(ncols):
            out.write(seq)
    stream.seek(0)
    pairs = msafile_pairs(MSAFile(stream, format=format))
    assert pairs == roundtrip_expected(ncols)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
def test_msafile_stockholm_stream_ends_alignment(ncols):

    stream = StringIO()
    with MSAFile(stream, 'w', format='stockholm') as out:
        for seq in roundtrip_msa(ncols):
            out.write(seq)
    text = stream.getvalue()
    lines = [line for line in text.splitlines() if line.strip()]
    assert lines[-1] == '//'
    pairs = msafile_pairs(MSAFile(StringIO(text), format='stockholm'))
    assert pairs == roundtrip_expected(ncols)


@pytest.mark.parametrize('ext', ['.slx', '.sth'])
def test_selex_writer_separates_labels(tmp_path, ext):

    filename = writeMSA(str(tmp_path / ('layout' + ext)), roundtrip_msa(70))
    lines = [line for line in read_text(filename).splitlines()
             if line.strip() and not line.startswith(('#', '//'))]
    assert join_labelled_lines(lines) == roundtrip_expected(70)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('ext', ['.aln', '.aln.gz'])
def test_clustal_writer_layout(tmp_path, ext, ncols):

    filename = writeMSA(str(tmp_path / ('layout' + ext)), roundtrip_msa(ncols))
    lines = read_text(filename).splitlines()
    assert lines[0].startswith('CLUSTAL')
    lines = [line for line in lines[1:]
             if line.strip() and not line[0].isspace()]
    assert join_labelled_lines(lines) == roundtrip_expected(ncols)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('ext', ['.ali', '.ali.gz'])
def test_pir_writer_layout(tmp_path, ext, ncols):

    filename = writeMSA(str(tmp_path / ('layout' + ext)), roundtrip_msa(ncols))
    entries = read_text(filename).split('>P1;')
    assert entries[0].strip() == ''
    pairs = []
    for entry in entries[1:]:
        lines = entry.splitlines()
        sequence = ''.join(line.strip() for line in lines[2:])
        assert sequence.endswith('*'), repr(entry)
        pairs.append((lines[0].strip(), sequence[:-1]))
    assert pairs == roundtrip_expected(ncols)


@pytest.mark.parametrize('ncols', [1, 12, 70, 150])
@pytest.mark.parametrize('ext', ['.slx', '.sth', '.slx.gz', '.sth.gz'])
def test_label_lookup_after_selex_roundtrip(tmp_path, ext, ncols):

    filename = writeMSA(str(tmp_path / ('lookup' + ext)), roundtrip_msa(ncols))
    msa = parseMSA(filename)
    for label, seq in roundtrip_expected(ncols):
        name = re.sub(r'/\d+-\d+$', '', label)
        assert str(msa[name]) == seq, name


def test_label_lookup_in_stockholm_file():

    msa = parseMSA(pathDatafile('msa_formats_blocks.sth'))
    assert str(msa['short']) == 'AC-EFG.IKLMN-QRSTV-Y'
    assert str(msa['seq:2*x|y']) == 'ACDEFGHIK-mnpqrsTVWY'
    assert str(msa['P12345.1']) == 'ACDEFGHIKLMNPQRSTVWY'


@pytest.mark.parametrize('ext', ['.slx', '.sth'])
def test_repository_selex_fixtures_labels(ext):

    fasta = parseMSA(pathDatafile('msa_Cys_knot.fasta'))
    other = parseMSA(pathDatafile('msa_Cys_knot' + ext))
    assert other.getLabels(full=True) == fasta.getLabels(full=True)
    assert msa_pairs(other) == msa_pairs(fasta)


@pytest.mark.parametrize('compressed', [False, True], ids=['plain', 'gz'])
@pytest.mark.parametrize('name', ['blocks.sth', 'blocks.slx',
                                  'two_alignments.sth',
                                  'clustalw2.aln', 'clustalw1.aln',
                                  'modeller.ali'])
def test_parsemsa_reads_fixture(tmp_path, name, compressed):

    filename = formats_fixture(name, tmp_path, compressed)
    assert msa_pairs(parseMSA(filename)) == EXPECTED_FIXTURES[name]


@pytest.mark.parametrize('name', ['blocks.sth', 'blocks.slx',
                                  'two_alignments.sth',
                                  'clustalw2.aln', 'clustalw1.aln',
                                  'modeller.ali'])
def test_msafile_reads_fixture(name):

    msafile = MSAFile(pathDatafile('msa_formats_' + name))
    assert msafile_pairs(msafile) == EXPECTED_FIXTURES[name]


@pytest.mark.parametrize('compressed', [False, True], ids=['plain', 'gz'])
@pytest.mark.parametrize('bad, good', [
    ('blocks_reordered.sth', 'blocks.sth'),
    ('blocks_extra_row.sth', 'blocks.sth'),
    ('blocks_label_only.sth', 'blocks.sth'),
    ('blocks_missing_row.slx', 'blocks.slx'),
    ('blocks_renamed.slx', 'blocks.slx'),
])
def test_mismatched_blocks_raise_ioerror(tmp_path, bad, good, compressed):

    filename = formats_fixture(good, tmp_path, compressed)
    assert msa_pairs(parseMSA(filename)) == EXPECTED_FIXTURES[good]
    with pytest.raises(IOError):
        parseMSA(formats_fixture(bad, tmp_path, compressed))
    with pytest.raises(IOError):
        list(MSAFile(formats_fixture(bad, tmp_path, compressed)))


@pytest.mark.parametrize('ext', ['.aln', '.aln.gz'])
def test_clustal_roundtrip_labels_with_spaces(tmp_path, ext):

    labels = ['sp|P1|X_HUMAN some description', 'second 2 words']
    msa = MSA(np.array([list('ACDEFGHIKL'), list('MNPQRSTVWY')], dtype='|S1'),
              labels=labels, aligned=True)
    filename = writeMSA(str(tmp_path / ('spaces' + ext)), msa)
    assert msa_pairs(parseMSA(filename)) == [
        (labels[0], 'ACDEFGHIKL'), (labels[1], 'MNPQRSTVWY')]


def test_clustal_output_is_read_by_biopython(tmp_path):

    from Bio import AlignIO

    filename = writeMSA(str(tmp_path / 'bio.aln'), roundtrip_msa(70))
    alignment = AlignIO.read(filename, 'clustal')
    assert len(alignment) == len(ROUNDTRIP_LABELS)
    assert alignment.get_alignment_length() == 70


def test_clustal_reads_other_aligner_headers(tmp_path):

    path = tmp_path / 'kalign.aln'
    path.write_text('Kalign (2.0) alignment in ClustalW format\n\n'
                    'seq1  ACDE-F\nseq2  AC-EGF\n')
    assert msa_pairs(parseMSA(str(path))) == [('seq1', 'ACDE-F'),
                                              ('seq2', 'AC-EGF')]


@pytest.mark.parametrize('ext', ['.slx', '.slx.gz'])
def test_selex_non_ascii_labels_across_blocks(tmp_path, ext):

    text = 'café ACD\nb EFG\n\ncafé HI\nb KL\n'
    path = str(tmp_path / ('utf8' + ext))
    opener = gzip.open if ext.endswith('.gz') else open
    with opener(path, 'wt', encoding='utf-8') as out:
        out.write(text)
    expected = [('café', 'ACDHI'), ('b', 'EFGKL')]
    assert msa_pairs(parseMSA(path)) == expected
    assert msafile_pairs(MSAFile(path)) == expected


def test_failed_msafile_init_does_not_raise_on_cleanup(monkeypatch):

    import gc
    import sys

    unraisable = []
    monkeypatch.setattr(sys, 'unraisablehook', unraisable.append)
    stream = StringIO()
    stream.close()
    with pytest.raises(ValueError):
        MSAFile(stream, format='fasta')
    gc.collect()
    assert unraisable == []


def test_fasta_label_ends_at_tab(tmp_path):

    path = tmp_path / 'tabs.fasta'
    path.write_text('>id1\tsome desc\nACDEF\n>id2/1-5\tx\nAC-EF\n')
    msa = parseMSA(str(path))
    assert [msa.getLabel(i, full=True) for i in range(2)] == ['id1', 'id2/1-5']
    assert msa.getResnums(1) == (1, 5)
    assert str(msa['id1']) == 'ACDEF'
