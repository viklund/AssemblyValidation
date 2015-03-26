# Copyright (C) 2015 by Per Unneberg
import unittest
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpld3 import plugins
import mpld3

css= """
table
{
  border-collapse: collapse;
}
th
{
  color: #ffffff;
  background-color: #000000;
}
td
{
  background-color: #cccccc;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 1px solid black;
  text-align: right;
}
"""

class TestCumulativeContigs(unittest.TestCase):
    def setUp(self):
        random.seed(12345)
        self.ir = range(0, 15)
        self.contigsizes = sorted([random.randint(1, 10**6) for i in self.ir], reverse=True)
        self.cumsum = np.cumsum(self.contigsizes)
        self.contignames = ["contig_{}".format(i) for i in self.ir]

    def test_plot(self):
        fig, ax = plt.subplots()
        ax.grid(True, alpha=0.3)
        # fig, ax = plt.subplots(subplot_kw=dict(axisbg='#EEEEEE'))
        # ax.grid(True, alpha=0.3)
        scatter = ax.plot(range(0, len(self.cumsum)),
                             self.cumsum,
                            alpha=0.3)
#cmap=plt.cm.jet)
        ax.grid(color='white', linestyle='solid')

        ax.set_title("Scatter Plot (with tooltips!)", size=20)

        labels = []
        for i in range(0, 15):
            label = df.ix[[i], :].T
            label.columns = ['Row {0}'.format(i)]
            # .to_html() is unicode; so make leading 'u' go away with str()
            labels.append(str(label.to_html()))

            #labels = self.contignames
        tooltip = mpld3.plugins.PointHTMLTooltip(scatter[0], labels=labels, voffset=10, hoffset=10, css=css)
        mpld3.plugins.connect(fig, tooltip)

        mpld3.show()


    def test_html(self):
        fig, ax = plt.subplots()
        ax.grid(True, alpha=0.3)

        labels = []
        for i in range(0, 15):
            label = [["Number of contigs", i], ["Length (bp)", self.cumsum[i]], ["Contig name", self.contignames[i]]]
            # .to_html() is unicode; so make leading 'u' go away with str()
            #labels.append(str(label.to_html()))
            labels.append(str(label))

        points = ax.plot(self.ir, self.cumsum, 'o', color='b',
                    mec='k', ms=15, mew=1, alpha=.6)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('HTML tooltips', size=20)

        tooltip = plugins.PointHTMLTooltip(points[0], labels,
                                        voffset=10, hoffset=10, css=css)
        plugins.connect(fig, tooltip)
        mpld3.show()
