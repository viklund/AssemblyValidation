# Copyright (C) 2015 by Per Unneberg
import mpld3
import matplotlib.pyplot as plt
import numpy as np
import random
import unittest
from mpld3 import plugins

def plot_cumulative_contigs(contigs, title="Cumulative contigs", *args, **kwargs):
    """Generate plot of cumulative contigs

    Args:
      contigs: list of tuples, where tuples are (contigname, contigsize)

    Returns:
      pltobj: a plotting object of class ...
    """
    N = len(contigs)
    ir = range(1, N+1)
    # Unzip the tuple pairs to two tuple lists
    (contignames, contigsizes) = zip(*contigs)
    cumsum = np.cumsum(contigsizes)
    return cumsum



    # # FIXME: customize from **kwargs
    # fig, ax = plt.subplots()
    # ax.grid(True, alpha=0.3)

    # labels = []
    # for i in range(1, N+1):
    #     # FIXME: add generic function for converting list of lists to a html table. texttable?
    #     label = "<table>" + "\n".join(["<tr><td>Number of contigs: {}</td></tr>".format(i),
    #                                     "<tr><td>Length (bp): {}</td></tr>".format(cumsum[i-1]),
    #                                     "<tr><td>Contig name: {}</td></tr>".format(contignames[i-1])]) + "</table>"
    #     labels.append(label)

# points = ax.plot(ir, cumsum, 'o', color='b', mec='k', ms=15, mew=1, alpha=.6)
    # plt.xlim(0, max(self.ir)+1)

    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_title(title, size=20)

    #tooltip = plugins.PointHTMLTooltip(points[0], labels, voffset=10, hoffset=10, css=css)
    #plugins.connect(fig, tooltip)
#    return points
    
    


# Module-level tests - just for testing the plotting
# This is probably not the way to do it
with open("../static/basic.css") as fh:
    css = "\n".join(fh.readlines())
    
class TestCumulativeContigs(unittest.TestCase):
    def setUp(self):
        random.seed(12345)
        self.ir = range(1, 15)
        self.contigsizes = sorted([random.randint(1, 10**6) for i in self.ir], reverse=True)
        self.cumsum = np.cumsum(self.contigsizes)
        self.contignames = ["contig_{}".format(i) for i in self.ir]

    def test_html(self):
        fig, ax = plt.subplots()
        ax.grid(True, alpha=0.3)

        labels = []
        for i in self.ir:
            label = "<table>" + "\n".join(["<tr><td>Number of contigs: {}</td></tr>".format(i),
                                            "<tr><td>Length (bp): {}</td></tr>".format(self.cumsum[i-1]),
                                            "<tr><td>Contig name: {}</td></tr>".format(self.contignames[i-1])]) + "</table>"
            labels.append(str(label))

        points = ax.plot(self.ir, self.cumsum, 'o', color='b',
                         mec='k', ms=15, mew=1, alpha=.6)
        plt.xlim(0, max(self.ir)+1)

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('HTML tooltips', size=20)

        tooltip = plugins.PointHTMLTooltip(points[0], labels,
                                        voffset=10, hoffset=10, css=css)
        plugins.connect(fig, tooltip)
        mpld3.show()
