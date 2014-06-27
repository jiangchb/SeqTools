#
# This script maps ChIP-Seq peaks onto regulatory sequences
#

import re, sys, os
from argParser import ArgParser
ap = ArgParser(sys.argv)

# peakpath is the filepath to a *_peaks.encodePeak
peakpath = ap.getArg("--peaks")

# the file containing the intergenic region mappings
gff = ap.getArg("--gff")

# build a library of genomic regions and 3' and 5' gene

chr_sites_genes = {} # key = chromosome name, value = hash; key = start site of gene, value = (end site, gene name)

fin = open(gff, "r")
for l in fin.xreadlines():
    if l.__len__() > 0 and False == l.startswith("#"):
        tokens = l.split()
        if tokens[2] != "gene": # We're only interested in genes
            continue
        chr = tokens[0]
        start = int( tokens[3] )
        stop = int( tokens[4] )
        gene = tokens[8].split(";")[0].split("=")[1]
        print chr, start, stop, gene

fin.close() 