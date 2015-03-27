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
    p.scatter(x, y, size=size, color=color, alpha=alpha)
    return p

