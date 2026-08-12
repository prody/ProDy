"""This module contains some unit tests for :mod:`prody.utilities.catchall` module,
starting with tree-related ones."""

import os
import tempfile
import numpy as np

from prody.tests import unittest
from prody.utilities import calcTree, findSubgroups, writeTree, parseTree
from prody.tests.datafiles import parseDatafile, pathDatafile

class TestCalcTree(unittest.TestCase):

    def testCalcTreeUPGMA(self):
        """Test calcTree with UPGMA method."""
        names = ['A', 'B', 'C', 'D']
        distance_matrix = np.array([[0,   1,   2, 1],
                                    [1,   0, 1.5, 2],
                                    [2, 1.5,   0, 2],
                                    [1,   2,   2, 0]])
        tree = calcTree(names, distance_matrix, method='upgma')
        self.assertIsNotNone(tree)
        # Check that tree has 4 leaves and they include the names
        leaves = tree.get_terminals()
        self.assertEqual(len(leaves), 4)
        self.assertEqual(set([leaf.name for leaf in leaves]), set(names))
        # Check that the tree has split evenly as expected for UPGMA
        self.assertEqual(len(tree.root.clades), 2)

    def testCalcTreeNJ(self):
        """Test calcTree with NJ method."""
        names = ['A', 'B', 'C', 'D']
        distance_matrix = np.array([[0,   1,   2, 1],
                                    [1,   0, 1.5, 2],
                                    [2, 1.5,   0, 2],
                                    [1,   2,   2, 0]])
        tree = calcTree(names, distance_matrix, method='nj')
        self.assertIsNotNone(tree)
        leaves = tree.get_terminals()
        # Check that tree has 4 leaves and they include the names
        self.assertEqual(len(leaves), 4)
        self.assertEqual(set([leaf.name for leaf in leaves]), set(names))
        # Check that the tree has split unevenly as expected for NJ
        self.assertEqual(len(tree.root.clades), 3)

    def testCalcTreeMismatchSize(self):
        """Test calcTree with mismatched names and matrix sizes."""
        names = ['A', 'B']
        distance_matrix = np.array([[0, 1, 2],
                                    [1, 0, 1.5],
                                    [2, 1.5, 0]])
        with self.assertRaises(ValueError):
            calcTree(names, distance_matrix)


class TestFindSubgroups(unittest.TestCase):

    def setUp(self):
        """Set up a test tree for findSubgroups tests."""
        # Create a simple distance matrix with clear clustering
        # Points A,B are close (distance 0.5), C,D are close (distance 0.5)
        # But A,B are far from C,D (distance 5)
        self.names = ['A', 'B', 'C', 'D']
        self.distance_matrix = np.array([[0.0, 0.5, 5.0, 5.0],
                                         [0.5, 0.0, 5.0, 5.0],
                                         [5.0, 5.0, 0.0, 0.5],
                                         [5.0, 5.0, 0.5, 0.0]])
        self.tree = calcTree(self.names, self.distance_matrix, method='upgma')

    def testFindSubgroupsNaiveMethod(self):
        """Test findSubgroups with naive method."""
        # Using cutoff 2.0 should separate into 2 subgroups
        subgroups = findSubgroups(self.tree, 2.0, method='naive')
        self.assertIsNotNone(subgroups)
        self.assertEqual(len(subgroups), 2)
        # Check that subgroups contain the expected names
        all_names = [name for subgroup in subgroups for name in subgroup]
        self.assertEqual(set(all_names), set(self.names))

    def testFindSubgroupsNaiveLargeCutoff(self):
        """Test findSubgroups with naive method and large cutoff."""
        # Using cutoff 10.0 should keep everything in one subgroup
        subgroups = findSubgroups(self.tree, 10.0, method='naive')
        self.assertEqual(len(subgroups), 1)
        self.assertEqual(set(subgroups[0]), set(self.names))

    def testFindSubgroupsNaiveTinyCutoff(self):
        """Test findSubgroups with naive method and tiny cutoff."""
        # Using cutoff 0.1 should separate all into individual subgroups
        subgroups = findSubgroups(self.tree, 0.1, method='naive')
        self.assertEqual(len(subgroups), 4)
        # Each subgroup should have one member
        for subgroup in subgroups:
            self.assertEqual(len(subgroup), 1)

    def testFindSubgroupsReturnsListOfLists(self):
        """Test that findSubgroups returns a list of lists."""
        subgroups = findSubgroups(self.tree, 2.0, method='naive')
        self.assertIsInstance(subgroups, list)
        for subgroup in subgroups:
            self.assertIsInstance(subgroup, list)


class TestParseTree(unittest.TestCase):

    def testParseUPGMATree(self):
        """Test parsing an UPGMA tree from a file."""
        tree_fn = pathDatafile('upgma_tree')
        tree = parseTree(tree_fn)
        self.assertIsNotNone(tree)
        # Check that tree has expected number of leaves
        leaves = tree.get_terminals()
        self.assertEqual(len(leaves), 4)
        # Check that tree has expected number of top-level clades
        self.assertEqual(len(tree.root.clades), 2)

    def testParseNJTree(self):
        """Test parsing a neighbor-joining tree from a file."""
        tree_fn = pathDatafile('nj_tree')
        tree = parseTree(tree_fn)
        self.assertIsNotNone(tree)
        # Check that tree has expected number of leaves
        leaves = tree.get_terminals()
        self.assertEqual(len(leaves), 4)
        # Check that tree has expected number of top-level clades
        self.assertEqual(len(tree.root.clades), 3)

    def testParseTreeTreeType(self):
        """Test that parseTree returns a Biopython Tree object."""
        try:
            from Bio import Phylo
            tree = parseDatafile('upgma_tree')
            self.assertIsInstance(tree, Phylo.BaseTree.Tree)
        except ImportError:
            self.skipTest("Biopython not available")

    def testParseTreeWrongFilepath(self):
        """Test parseTree with non-existent file."""
        with self.assertRaises((AssertionError, FileNotFoundError)):
            parseTree('/nonexistent/path/to/tree.nwk')

    def testParseTreeWrongFileType(self):
        """Test parseTree with invalid filename argument."""
        with self.assertRaises(TypeError):
            parseTree(123)


class TestWriteTree(unittest.TestCase):

    def setUp(self):
        """Set up test trees for writing."""
        self.upgma_tree = parseDatafile('upgma_tree')
        self.nj_tree = parseDatafile('nj_tree')
        # Create a temporary directory for test files
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary test files."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def testWriteUPGMATree(self):
        """Test writing an UPGMA tree to a file."""
        output_file = os.path.join(self.temp_dir, 'test_upgma.nwk')
        try:
            writeTree(output_file, self.upgma_tree)
            self.assertTrue(os.path.exists(output_file))
            # Verify the file is not empty
            with open(output_file, 'r') as f:
                content = f.read()
                self.assertTrue(len(content) > 0)
                # Check for Newick format markers
                self.assertIn(';', content)
        except ImportError:
            self.skipTest("Biopython not available")

    def testWriteNJTree(self):
        """Test writing a neighbor-joining tree to a file."""
        output_file = os.path.join(self.temp_dir, 'test_nj.nwk')
        try:
            writeTree(output_file, self.nj_tree)
            self.assertTrue(os.path.exists(output_file))
            # Verify the file is not empty
            with open(output_file, 'r') as f:
                content = f.read()
                self.assertTrue(len(content) > 0)
                # Check for Newick format markers
                self.assertIn(';', content)
        except ImportError:
            self.skipTest("Biopython not available")

    def testWriteTreeWrongFilename(self):
        """Test writeTree with invalid filename argument."""
        with self.assertRaises(TypeError):
            writeTree(123, self.upgma_tree)

    def testWriteTreeWrongTreeType(self):
        """Test writeTree with invalid tree argument."""
        output_file = os.path.join(self.temp_dir, 'test.nwk')
        with self.assertRaises(TypeError):
            writeTree(output_file, "not a tree")

    def testWriteTreeWrongFormat(self):
        """Test writeTree with invalid format argument."""
        output_file = os.path.join(self.temp_dir, 'test.nwk')
        with self.assertRaises(TypeError):
            writeTree(output_file, self.upgma_tree, format_str=123)

    def testWriteAndParseRoundtrip(self):
        """Test writing a tree and then parsing it back."""
        output_file = os.path.join(self.temp_dir, 'roundtrip.nwk')
        try:
            # Write the tree
            writeTree(output_file, self.upgma_tree)
            # Parse it back
            parsed_tree = parseTree(output_file)
            # Verify the parsed tree is valid
            self.assertIsNotNone(parsed_tree)
            leaves = parsed_tree.get_terminals()
            self.assertEqual(len(leaves), 4)
            self.assertEqual(len(parsed_tree.root.clades), 2)
        except ImportError:
            self.skipTest("Biopython not available")

if __name__ == '__main__':
    unittest.main()
