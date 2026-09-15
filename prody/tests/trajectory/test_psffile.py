"""This module contains unit tests for :mod:`~prody.trajectory.psffile`."""

from os.path import join

from numpy.testing import assert_equal

from prody import parsePSF, writePSF
from prody.tests import TEMPDIR, TestCase
from prody.tests.datafiles import DATA_FILES, pathDatafile


PSF = pathDatafile('topology_psf')
COUNTS = DATA_FILES['topology_psf']

# the topology of the fixture, in the order it appears in the file and using
# 1-based indices as the file does
BONDS = [[2, 1], [3, 1], [4, 3], [5, 3], [6, 5]]
ANGLES = [[2, 1, 3], [1, 3, 5], [4, 3, 5], [3, 5, 6]]
DIHEDRALS = [[2, 1, 3, 5], [1, 3, 5, 6], [4, 3, 5, 6]]
IMPROPERS = [[5, 3, 7, 6]]
DONORS = [[1, 2]]
ACCEPTORS = [[6, 5]]
CROSSTERMS = [[5, 7, 9, 10, 1, 3, 5, 7]]

SECTIONS = ('bonds', 'angles', 'dihedrals', 'impropers',
            'donors', 'acceptors', 'crossterms')


def topology(ag):
    """Returns the topology of *ag* section by section, as 1-based lists."""

    return {name: (getattr(ag, '_' + name) + 1).tolist() for name in SECTIONS}


class TestPSFFile(TestCase):

    def setUp(self):

        self.ag = parsePSF(PSF)
        self.top = topology(self.ag)

    def testCounts(self):
        """Every section is read, and read whole."""

        self.assertEqual(self.ag.numAtoms(), COUNTS['atoms'])
        for name in SECTIONS:
            self.assertEqual(len(self.top[name]), COUNTS[name],
                             'wrong number of {0}'.format(name))

    def testCrosstermWidth(self):
        """A CMAP cross-term is the eight atoms of two coupled dihedrals."""

        assert_equal(self.ag._crossterms.shape, (COUNTS['crossterms'], 8))

    def testAtomOrderWithinTerms(self):
        """The atom order within a term identifies it and must be preserved.

        An angle's middle index is its vertex, a torsion's four indices are a
        sequence, and a cross-term's eight are two dihedrals; sorting them
        turns each into a different term."""

        assert_equal(self.top['angles'], ANGLES)
        assert_equal(self.top['dihedrals'], DIHEDRALS)
        assert_equal(self.top['impropers'], IMPROPERS)
        assert_equal(self.top['crossterms'], CROSSTERMS)
        # donors and acceptors are ordered pairs too
        assert_equal(self.top['donors'], DONORS)
        assert_equal(self.top['acceptors'], ACCEPTORS)

    def testBondsAreCanonicalised(self):
        """A bond is symmetric, so its pair is sorted and deduplicated."""

        assert_equal(self.top['bonds'], sorted(sorted(b) for b in BONDS))

    def testAngleVertex(self):
        """The vertex an Angle reports is the one the file gives it."""

        vertices = [angle.getAtoms()[1].getIndex() + 1
                    for angle in self.ag.iterAngles()]
        assert_equal(vertices, [a[1] for a in ANGLES])

    def testRoundTrip(self):
        """Writing a parsed topology and reading it back changes nothing."""

        out = writePSF(join(TEMPDIR, 'topology_roundtrip.psf'), self.ag)
        again = topology(parsePSF(out))
        for name in SECTIONS:
            assert_equal(again[name], self.top[name],
                         '{0} did not round-trip'.format(name))

    def testRoundTripIsIdempotent(self):
        """A second write/read pass is a no-op, not a further drift."""

        first = parsePSF(writePSF(join(TEMPDIR, 'topology_rt1.psf'), self.ag))
        second = parsePSF(writePSF(join(TEMPDIR, 'topology_rt2.psf'), first))
        assert_equal(topology(second), topology(first))


class TestWritePSFSubset(TestCase):
    """A selection is what cropping a solvated system produces, and it has to come
    out as a PSF that stands on its own."""

    def setUp(self):

        self.ag = parsePSF(PSF)

    def crop(self, selstr, name):

        sel = self.ag.select(selstr)
        return sel, parsePSF(writePSF(join(TEMPDIR, name), sel))

    def testIndicesAreRenumbered(self):
        """Terms are written as positions in the new file, not in the parent group."""

        sel, back = self.crop('index 2 to 5', 'crop_bonds.psf')
        self.assertEqual(back.numAtoms(), sel.numAtoms())
        assert_equal((back._bonds + 1).tolist(), [[1, 2], [1, 3], [3, 4]])
        # parent angles 4-3-5 and 3-5-6 are the two wholly inside; atom 3 is local 1,
        # so the vertex stays in the middle
        assert_equal((back._angles + 1).tolist(), [[2, 1, 3], [1, 3, 4]])
        self.assertLessEqual(back._bonds.max() + 1, back.numAtoms(),
                             'a term points past the end of the file')

    def testCrosstermIsRenumberedInOrder(self):
        """The parent term 5 7 9 10 1 3 5 7 becomes local 3 4 5 6 1 2 3 4."""

        _, back = self.crop('index 0 2 4 6 8 9', 'crop_cmap.psf')
        assert_equal((back._crossterms + 1).tolist(), [[3, 4, 5, 6, 1, 2, 3, 4]])

    def testOrderSurvivesSelection(self):
        """Selecting must not reorder a term.  A neighbour map records which atoms
        share a term, not where in it they sit, so rebuilding a term from one put the
        pivot atom first and renamed an angle's vertex."""

        _, back = self.crop('index 0 to 5', 'crop_order.psf')
        assert_equal((back._angles + 1).tolist(), ANGLES)
        assert_equal((back._acceptors + 1).tolist(), ACCEPTORS)

    def testPartialTermsAreDropped(self):
        """A term reaching outside the selection cannot be written at all."""

        _, back = self.crop('index 0 1', 'crop_partial.psf')
        assert_equal((back._bonds + 1).tolist(), [[1, 2]])
        for name in ('_angles', '_dihedrals', '_crossterms'):
            self.assertIsNone(getattr(back, name),
                              '{0} kept a term that is not wholly selected'.format(name))


class TestTopologyUnderSelection(TestCase):
    """A selection has to carry the topology two ways: by delegation, where terms are
    read through the parent group and keep its numbering, and through copy(), which
    makes an independent AtomGroup numbered from one."""

    def setUp(self):

        self.ag = parsePSF(PSF)
        self.sel = self.ag.select('index 0 to 5')

    def testNumbondsDataCoversEveryAtom(self):
        """The num* arrays are per-atom, so they are as long as the group even when
        the last atoms take part in no term -- `numbonds 0` selects ions by relying
        on that."""

        for label in ('numbonds', 'numangles', 'numdihedrals', 'numimpropers',
                      'numdonors', 'numacceptors', 'numcrossterms'):
            data = self.ag.getData(label)
            self.assertEqual(len(data), self.ag.numAtoms(),
                             '{0} is not per-atom'.format(label))

    def testSelectionDelegatesTerms(self):
        """Iterating a selection yields the terms wholly inside it, in parent
        numbering and in the order they were set."""

        assert_equal([list(t) for t in self.sel._iterAngles()], self.ag._angles.tolist())
        assert_equal([list(t) for t in self.sel._iterAcceptors()],
                     self.ag._acceptors.tolist())

    def testSelectionNumBonds(self):
        """A Bond must be built on the parent group; on the pointer it has no
        _bondOrders and getBonds raised AttributeError."""

        self.assertEqual(self.sel.numBonds(), 5)

    def testCopyCarriesEveryTopologySection(self):
        """copy() used to carry the bonds and silently drop everything else."""

        c = self.sel.copy()
        self.assertEqual(c.numAtoms(), 6)
        assert_equal((c._bonds + 1).tolist(), sorted(sorted(b) for b in BONDS))
        assert_equal((c._angles + 1).tolist(), ANGLES)
        assert_equal((c._dihedrals + 1).tolist(), DIHEDRALS)
        assert_equal((c._donors + 1).tolist(), DONORS)
        assert_equal((c._acceptors + 1).tolist(), ACCEPTORS)
        # the improper and the cross-term reach outside the selection
        self.assertIsNone(c._impropers)
        self.assertIsNone(c._crossterms)

    def testWholeGroupCopyIsUnchanged(self):

        c = self.ag.copy()
        for name in SECTIONS:
            assert_equal(getattr(c, '_' + name), getattr(self.ag, '_' + name))

    def testCopyThenWriteMatchesWritingTheSelection(self):
        """The two routes to a cropped PSF must agree."""

        viaCopy = parsePSF(writePSF(join(TEMPDIR, 'via_copy.psf'), self.sel.copy()))
        viaSel = parsePSF(writePSF(join(TEMPDIR, 'via_sel.psf'), self.sel))
        for name in SECTIONS:
            assert_equal(getattr(viaCopy, '_' + name), getattr(viaSel, '_' + name),
                         '{0} differs between the two routes'.format(name))
