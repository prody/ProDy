from numpy import array
from numpy.testing import assert_allclose

from prody.tests import TestCase

from prody.utilities import rangeString
from prody.utilities.misctools import getMasses


class TestRangeString(TestCase):

    def testContinuous(self):

        self.assertEqual(rangeString(list(range(10))), '0 to 9')

    def testNegative(self):

        self.assertEqual(rangeString(list(range(-5, 10)), pos=False),
                         '-5 to 9')

    def testGapped(self):

        self.assertEqual(rangeString(list(range(-5, 10)) +
                                     list(range(15, 20)) +
                                     list(range(25, 30)), pos=False),
                                     '-5 to 9 15 to 19 25 to 29')

    def testRepeated(self):

        self.assertEqual(rangeString(list(range(10, 20)) +
                                     list(range(15, 20)) +
                                     list(range(30))), '0 to 29')


class TestGetMasses(TestCase):

    def testKnownElements(self):

        assert_allclose(getMasses(['C', 'N', 'O', 'S']),
                        [12.0107, 14.0067, 15.9994, 32.065])

    def testCaseInsensitive(self):
        """Element symbols are matched regardless of case."""

        assert_allclose(getMasses(['c', 'n', 'FE']),
                        getMasses(['C', 'N', 'Fe']))

    def testUnknownElementIsZero(self):
        """An unrecognised symbol contributes zero mass, it does not raise."""

        assert_allclose(getMasses(['C', 'Xx', 'O']), [12.0107, 0., 15.9994])

    def testRepeatedElements(self):
        """Repeated symbols must all map back to their own mass.

        The lookup is done once per distinct symbol and mapped back onto the
        atoms, so a wrong mapping would show up here as masses landing on the
        wrong atoms.
        """

        elements = ['O', 'C', 'C', 'N', 'O', 'S', 'C', 'N']
        expected = [15.9994, 12.0107, 12.0107, 14.0067,
                    15.9994, 32.065, 12.0107, 14.0067]
        assert_allclose(getMasses(elements), expected)

    def testEmpty(self):

        self.assertEqual(len(getMasses([])), 0)

    def testString(self):
        """A single symbol still returns a scalar."""

        self.assertAlmostEqual(getMasses('c'), 12.0107)

    def testArrayInput(self):

        assert_allclose(getMasses(array(['C', 'O'])), [12.0107, 15.9994])
