# Copyright (C) 2015 by Per Unneberg
import unittest
import random
from gaqtk.mpld3.cumulative_contig_plot import plot_cumulative_contigs

class TestCumulativeContigs(unittest.TestCase):
    def test_plot_cumulative_contig(self):
        """Test plotting cumulative contig function"""
        tmp = plot_cumulative_contigs([("c1", 10), ("c2", 5)])
        print (tmp)
