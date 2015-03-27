#!/usr/bin/env python
# Copyright (C) 2015 by Per Unneberg
import sys
import os
import io
import jinja2
import textwrap
from gaqtk.bokeh.scatter import scatter
from gaqtk.bokeh.publish import static_html

def main():
    x = list(range(5))
    y = list(range(5))
    p1 = scatter(x=x, y=y, title="Plot 1")
    p2 = scatter(x=x, y=y, title="Plot 2")
    kw = {'p1' : p1, 'p2' : p2}
    html = static_html(**kw)
    with open("scatter.html", "w") as fh:
        fh.write(html)

if __name__ == "__main__":
    main()

