"""
September 2014: I'm writing this script to calculate the size
of chromosomes in the C. albicans genome, in order to run JAMM
for peak-finding. JAMM requires a text file listing the sizes
of each chromosome in the analysis.
"""

import os, sys

fastapath = open(sys.argv[1])
fin = open(fastapath, "r")
chrom_size = {}
for l in fin.xreadlines():
    