import sys
import unittest
import test_proteins
import test_select


def test(verbosity=2, descriptions=True, stream=sys.stderr):
    testrunner = unittest.TextTestRunner(stream, descriptions, verbosity)
    for module in [test_proteins, test_select]:
        testrunner.run(unittest.defaultTestLoader.loadTestsFromModule(module))
    
