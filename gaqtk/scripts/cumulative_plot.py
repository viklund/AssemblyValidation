#!/usr/bin/env python
# Copyright (C) 2015 by Per Unneberg
import sys
import os
import io
import jinja2
import textwrap
from gaqtk.bokeh.scatter import scatter, cumulative_contigs
from gaqtk.bokeh.publish import static_html

def main():
    x = list(range(1, 5))
    y = list(range(1, 5))
    z = ["c1", "c2", "c3", "c4", "c5"]
    p1 = scatter(x=x, y=y, title="Plot 1")
    p2 = cumulative_contigs(zip(z, y))
    kw = {'p1' : p1, 'p2' : p2}
    html = static_html(**kw)
    with open("scatter.html", "w") as fh:
        fh.write(html)

if __name__ == "__main__":
    main()

