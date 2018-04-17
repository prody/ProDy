from prody.tests import TestCase

from prody.utilities import rangeString


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
