#
# This script maps ChIP-Seq peaks onto regulatory sequences
#

import re, sys, os
from argParser import ArgParser
ap = ArgParser(sys.argv)

# peakpath is the filepath to a *_summits.bed
summitpath = ap.getArg("--summits")

# the file containing features of the genome
gff = ap.getArg("--gff")

outname = ap.getArg("--name") # A unique name for the output


def splash():
    print "============================================"
    print "  peaks2genes.py"
    print "  Victor Hanson-Smith, 2014"
    print "  victor.hanson-smith@ucsf.edu"
    print "============================================"


splash()

#
# Build a library of genes
#
chr_gene_sites = {} # key = chromosome name, value = hash; key = gene name, value = (Start, stop)
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
        if chr not in chr_gene_sites:
            chr_gene_sites[chr] = {}
        chr_gene_sites[chr][gene] = (start, stop)
fin.close() 

count_genes = 0
for chr in chr_gene_sites:
    count_genes += chr_gene_sites[chr].keys().__len__()
print "\n. I found", count_genes, "genes."

#
# Build a library of summit peaks
#
chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
fin = open(summitpath, "r")
for l in fin.xreadlines():
    if l.__len__() > 0 and False == l.startswith("#"):
        tokens = l.split()
        chr = tokens[0]
        site = int( tokens[1] )
        score = float( tokens[4] )
        if chr not in chr_site_score:
            chr_site_score[chr] = {}
        chr_site_score[chr][site] = score
fin.close()

count_summits = 0
for chr in chr_site_score:
    count_summits += chr_site_score[chr].keys().__len__()
print "\n. I found", count_summits, "summits."

#
# For each summit, find the closest genes
#
chr_sumsite_stats = {}

fout = open("p2g." + outname + ".txt", "w")
for chr in chr_site_score:
    if chr not in chr_sumsite_stats:
        chr_sumsite_stats[chr] = {}
    
    for sumsite in chr_site_score[chr]:
        min_up = None 
        closest_up = ""
        min_down = None
        closest_down = ""
        for gene in chr_gene_sites[chr]:
            start = chr_gene_sites[chr][gene][0]
            stop = chr_gene_sites[chr][gene][1]
            if start >= sumsite:
                d = start - sumsite
                if min_up == None:
                    min_up = d
                    closest_up = gene
                if min_up > d:
                    min_up = d
                    closest_up = gene
            elif stop <= sumsite:
                d = sumsite - stop
                if min_down == None:
                    min_down = d
                    closest_down = gene
                if min_down > d:
                    min_down = d
                    closest_down = gene
        if min_up == None:
            min_up = "n/a"
        if min_down == None:
            min_down = "n/a"
        #print chr, sumsite, chr_site_score[chr][sumsite], closest_up, min_up, closest_down, min_down
        chr_sumsite_stats[chr][sumsite] = (min_up, closest_up, min_down, closest_down)
        fout.write(chr + "\t" + sumsite.__str__() + "\t" + chr_site_score[chr][sumsite].__str__() + "\t" + closest_up + "\t" + min_up.__str__() + "\t" + closest_down + "\t -" + min_down.__str__() + "\n")
fout.close()

#
# Plot stats about the peak-gene mappings
#
from plot_histogram import *

# Distribution of distances
up_distances = []
down_distances = []
for chr in chr_sumsite_stats:
    for sumsite in chr_sumsite_stats[chr]:
        if chr_sumsite_stats[chr][sumsite][0] != "n/a":
            up_distances.append( int( chr_sumsite_stats[chr][sumsite][0]) )
        if chr_sumsite_stats[chr][sumsite][2] != "n/a":
            down_distances.append( int( chr_sumsite_stats[chr][sumsite][2]) )
plot_histogram(up_distances, outname + "up", xlab="min. distance to gene", ylab="proportion of summits", fixed_max=5000)
plot_histogram(down_distances, outname + "down", xlab="min. distaance to gene", ylab="proportion of summits", fixed_max=5000)

