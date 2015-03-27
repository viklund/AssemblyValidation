#!/usr/bin/env python2.7
"""
Script that plots the coverage of a sam/bam-file.
"""

import pysam
import numpy
#from matplotlib import pyplot

import matplotlib.pyplot as plt
import mpld3

class CoveragePlot(object):
    
    def __init__(self, assembly):
        self.assembly = assembly
        fig, ax = plt.subplots()
        self.fig = fig
        self.ax = ax
        
    def plot(self, name, length, cov, offset):
        print "%s %s %s" % (name, length, len(cov))

        local_bin_width = self.bin_width
        nbins = length / self.bin_width + 1
        if nbins < 10:
            nbins = 10
            local_bin_width = length / (nbins - 1)

        cov_bins = [0] * nbins
        for i, x in enumerate(cov):
            cov_bins[ i/local_bin_width ] += x

        cov_bins = [ x/local_bin_width for x in cov_bins ]
        
        self.ax.plot([x*local_bin_width+offset for x in range(len(cov_bins))], cov_bins)
        #mpld3.save_html(self.fig, "%s.html" % name)

    def save(self):
        mpld3.save_html(self.fig, "snubbe.html")
    
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
        self.bin_width = tot_length/10000

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
