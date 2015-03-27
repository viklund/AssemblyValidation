# Copyright (C) 2015 by Per Unneberg
from __future__ import division
import numpy as np
from bokeh.plotting import figure, HBox, output_file, show, VBox
from bokeh.models import Range1d

DEFAULT_TOOLS = "pan,wheel_zoom,box_zoom,reset,save"

def scatter(x, y, tools=DEFAULT_TOOLS, width=300, height=300, color="red", alpha=0.5, size=12, title="Scatter plot"):
    xr = Range1d(start=min(x), end=max(x))
    yr = Range1d(start=min(y), end=max(y))
    p = figure(x_range=xr, y_range=yr, tools=tools, plot_width=width, plot_height=height, title=title)
    p.segment(x[0:len(x)-1], y[0:len(y)-1], x[1:len(x)], y[1:len(y)], line_color="black", line_width=1)
    p.scatter(x, y, size=size, color=color, alpha=alpha)

    return p

def cumulative_contigs(contigs, tools=DEFAULT_TOOLS, width=300, height=300, color="red", alpha=0.5, size=8, title="Cumulative contig plot"):
    contigs = sorted(contigs, reverse=True, key=lambda tup: tup[1])
    (contignames, contigsizes) = zip(*contigs)
    contignames = [""] + list(contignames)
    contigsizes = [0] + list(contigsizes)
    cumsum = np.cumsum(contigsizes)
    contignames = contignames
    total_size = sum(contigsizes)
    n50 = n90 = None
    for i in range(len(cumsum)):
        cs = float(cumsum[i]) / total_size
        if float(cs) >= 0.5 and n50 is None:
            n50 = (i + 1, cumsum[i])
        if float(cs) >= 0.9 and n90 is None:
            n90 = (i + 1, cumsum[i])
            
    cumsum_frac = [float(x) / sum(contigsizes) for x in cumsum]
    xr = Range1d(start=0, end=len(contigsizes) + 1)
    yr = Range1d(start=-1, end=1.1 * float(max(cumsum)))
    x = list(range(1, len(contigsizes) + 1))
    p = figure(x_range=xr, y_range=yr, tools=tools, plot_width=width, plot_height=height, title=title)
    p.segment(x[0:len(x)-1], cumsum[0:len(cumsum)-1], x[1:len(x)], cumsum[1:len(cumsum)], line_color="black", line_width=1)
    p.segment([n50[0], n90[0]], [0, 0], [n50[0], n90[0]], [n50[1], n90[1]], line_width=1.5, alpha=1, line_color="black")
    p.scatter(x, cumsum, size=size, color=color, alpha=alpha )
    return p
