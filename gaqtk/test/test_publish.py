# Copyright (C) 2015 by Per Unneberg
import unittest
import random
from gaqtk.bokeh.publish import static_html
from gaqtk.bokeh.scatter import scatter
from bokeh.plotting import figure

class TestBokehPublish(unittest.TestCase):
    def test_static_html(self):
        """Test static html function"""
        tmp = static_html(**{'p1': figure(title="p1"), 'p2' : figure(title="p2")})
        print (tmp)

