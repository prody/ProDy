"""This module contains unit tests for general functions in
:mod:`~prody.dynamics.analysis`, i.e. ones that are not tied to a
particular model or to a more specific test module."""

import numpy as np

from prody import calcCollectivity, LOGGER
from prody.tests import unittest

LOGGER.verbosity = 'none'


class TestCalcCollectivity(unittest.TestCase):
    """The masses path must divide by ``sum(u2in)``, not ``sqrt(sum(u2in))``.

    Dividing by the square root breaks the probability-distribution
    requirement of the collectivity definition in equation 5 of
    Bruschweiler 1995.
    """

    def test_masses_path_matches_bruschweiler_definition(self):
        v = np.array([0.6, 0.8])
        masses = np.array([1.0, 4.0])

        u2 = v ** 2 / masses
        p = u2 / u2.sum()
        self.assertAlmostEqual(p.sum(), 1.0, places=12)
        expected = np.exp(-(p * np.log(p)).sum()) / len(v)

        actual = calcCollectivity(v, masses=masses, is3d=False)
        self.assertAlmostEqual(actual, expected, places=8)

    def test_default_path_without_masses_is_unaffected(self):
        # Unit-norm vector: sum(u2in) == sqrt(sum(u2in)) == 1, so this path
        # is identical before and after the fix. Guards against regressions.
        v = np.array([0.6, 0.8])
        u2 = v ** 2
        p = u2 / u2.sum()
        expected = np.exp(-(p * np.log(p)).sum()) / len(v)

        actual = calcCollectivity(v, is3d=False)
        self.assertAlmostEqual(actual, expected, places=8)


if __name__ == '__main__':
    unittest.main(verbosity=2)
