# Copyright (C) 2015 by Per Unneberg
import mpld3
import matplotlib.pyplot as plt
import numpy as np
import random



# Unittests
class TestCumulativeContigs(unittest.TestCase):
    def setUp(self):
        random.seed(12345)
        ir = range(0, 15)
        self.contigsizes = [random.randint(1, 10**6) for i in ir]
        self.contignames = ["contig_{}".format(i) for i in ir]

    def test_plot(self):
        fig, ax = plt.subplots(subplot_kw=dict(axisbg='#EEEEEE'))
        N = 100

        scatter = ax.scatter(np.random.normal(size=N),
                            np.random.normal(size=N),
                            c=np.random.random(size=N),
                            s=1000 * np.random.random(size=N),
                            alpha=0.3,
                            cmap=plt.cm.jet)
        ax.grid(color='white', linestyle='solid')

        ax.set_title("Scatter Plot (with tooltips!)", size=20)

        labels = ['point {0}'.format(i + 1) for i in range(N)]
        tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
        mpld3.plugins.connect(fig, tooltip)

        mpld3.show()
