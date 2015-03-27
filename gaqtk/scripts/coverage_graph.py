#!/usr/bin/env python2.7
"""
Script that plots the coverage of a sam/bam-file.
"""

import pysam
import numpy

from bokeh.plotting import figure, output_file, save


class CoveragePlot(object):
    
    def __init__(self, assembly):
        self.assembly = assembly
        output_file("scatter.html")
        self.p1 = figure(plot_width=700, plot_height=300)
        self.sections = []
        
    def plot(self, name, length, cov, offset):
        print "%s %s %s" % (name, length, len(cov))

        local_bin_width = self.bin_width
        nbins = length / self.bin_width + 1
        if nbins < 10:
            nbins = 10
            local_bin_width = length / (nbins - 1)

        cov_bins = {
            "min":  [float('Inf')] * nbins,
            "max":  [0.0] * nbins,
            "mean": [0.0] * nbins,
        }
        for i, x in enumerate(cov):
            bini = i / local_bin_width
            if cov_bins['min'][bini] > x:
                cov_bins['min'][bini] = x
            elif cov_bins['max'][bini] < x:
                cov_bins['max'][bini] = x
            cov_bins['mean'][bini] += x

        l = local_bin_width * 1.0
        for i, x in enumerate(cov_bins['mean']):
            #print cov_bins[i][2]
            cov_bins['mean'][i] = x*1.0/l
        cov_bins['x_vals'] = [x*local_bin_width for x in range(offset/local_bin_width, offset/local_bin_width+len(cov_bins['mean']))]
        cov_bins['length'] = length
        self.sections.append(cov_bins)

    def save(self):
        mx = max( [ max(x['max']) for x in self.sections ] )
        colors = ['#D8A031', '#D37600', '#C16A00', '#A84600']
        for i, cov_bins in enumerate(self.sections):
            x_vals = cov_bins['x_vals']
            self.p1.rect(x_vals[0] + cov_bins['length']/2, 0+mx/2, cov_bins['length'], mx, color=colors[i%4], alpha=1)

            self.p1.line(x_vals, cov_bins['min'], size=12, color="green", alpha=1)
            self.p1.line(x_vals, cov_bins['max'], size=12, color="red", alpha=1)
            self.p1.line(x_vals, cov_bins['mean'], size=12, color="blue", alpha=1)

        save(self.p1)
    
    def run(self):
        ext = self.assembly.split('.')[-1]
        if ext == 'bam':
            samfile = pysam.Samfile(self.assembly, 'rb')
        else:
            samfile = pysam.Samfile(self.assembly, 'r')

        tot_length = 0
        header = []
        for i in samfile.header['SQ']:
            header.append( (i['SN'], i['LN']) )
            tot_length += i['LN']
        self.bin_width = tot_length/1000

        length_sofar = 0
        
        for name, length in header:
            cov = [0] * length
            for pileupcolumn in samfile.pileup( name, 0, length ):
                cov[pileupcolumn.reference_pos] = pileupcolumn.n
            self.plot(name, length, cov, length_sofar)
            length_sofar += length
    

if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser( description = __doc__ )
    parser.add_argument("assembly", help="Assembly file in SAM/BAM format.")
    args = parser.parse_args()
    
    plotter = CoveragePlot(args.assembly)
    plotter.run()

    plotter.save()
