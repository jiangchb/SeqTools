#
# This script maps ChIP-Seq peaks onto regulatory sequences
#

import re, sys, os
from plot_histogram import *
from plot_scatter import *
from argParser import ArgParser
ap = ArgParser(sys.argv)

#################################################################
#
# User Input:
#

# peakpath is the filepath to a *_summits.bed
summitpaths = ap.getList("--summits")

# the file containing features of the genome
gffpath = ap.getArg("--gff")

# outname will be used as a keyword in output file names.
outname = ap.getArg("--name") # A unique name for the output

#
#
#################################################################


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
def read_gff(gffpath):
    """Returns chr_gene_sites"""
    chr_gene_sites = {} # key = chromosome name, value = hash; key = gene name, value = (Start, stop)
    fin = open(gffpath, "r")
    for l in fin.xreadlines():
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            if tokens[2] != "gene": # We're only interested in genes
                continue
            chr = tokens[0]
            start = int( tokens[3] )
            stop = int( tokens[4] )
            strand = tokens[6]
            gene = tokens[8].split(";")[0].split("=")[1]
            if chr not in chr_gene_sites:
                chr_gene_sites[chr] = {}
            if strand == "+":
                chr_gene_sites[chr][gene] = (start, stop)
            else: # antisense strand
                chr_gene_sites[chr][gene] = (stop, start)
    fin.close() 
    
    count_genes = 0
    for chr in chr_gene_sites:
        count_genes += chr_gene_sites[chr].keys().__len__()
    print "\n. I found", count_genes, "genes " in gffpath + "."
    
    return chr_gene_sites
    
def build_pk2gn(summitpath, chr_gene_sites):
    """Returns chr_sumsite_stats."""
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
    print "\n. I found", count_summits, "summits in " + summitpath + "."
    
    #
    # For each summit, find the closest genes
    #
    chr_sumsite_stats = {}
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
                d = start - sumsite
                
                """Sense direction, and upstream"""
                if start < stop and start >= sumsite:
                    
                    if min_up == None:
                        min_up = d
                        closest_up = gene
                    if min_up > d:
                        min_up = d
                        closest_up = gene
                elif start > stop and start <= sumsite:
                    """Antisense and downstream"""
                    if min_down == None:
                        min_down = d
                        closest_down = gene
                    if min_down < d: # be careful here, we're looking for the largest NEGATIVR number.
                        min_down = d
                        closest_down = gene
                        
                """Note: this ignores the cases where the summit is inside a
                coding region. It also ignores the cases where the nearest gene
                is in the wrong sense orientation, such that the gene's stop site
                is placed between the start site and the peak summit."""
            
            #print chr, sumsite, chr_site_score[chr][sumsite], closest_up, min_up, closest_down, min_down
            chr_sumsite_stats[chr][sumsite] = (min_up, closest_up, min_down, closest_down, chr_site_score[chr][sumsite])    
    return chr_sumsite_stats

def build_gn2pk(chr_sumsite_stats):
    """Returns gene_peaks: gene_peaks[gene] = [ (summit site, distance to gene, summit score), (n2), (n3) ]"""
    gene_peaks = {}
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            upgene = chr_sumsite_stats[chr][sumsite][1]
            upd = chr_sumsite_stats[chr][sumsite][0]
            downgene = chr_sumsite_stats[chr][sumsite][3]
            downd = chr_sumsite_stats[chr][sumsite][2]
            if upd != None:
                if upgene not in gene_peaks:
                    gene_peaks[ upgene ] = [] # continue here
                gene_peaks[ upgene ].append(  (sumsite, -1*int(upd), chr_sumsite_stats[chr][sumsite][4])  )       
            if downd != None:
                if downgene not in gene_peaks:
                    gene_peaks[ downgene ] = [] # continue here
                gene_peaks[ downgene ].append(  (sumsite, -1*int(downd), chr_sumsite_stats[chr][sumsite][4])  )       
    
    # sort peaks by distance to gene
    for gene in gene_peaks:
        summits = gene_peaks[gene]
        hashed_summits = {}
        for s in summits:
            d = abs( int(s[1]) )
            if d in hashed_summits:
                d += 1
            hashed_summits[d] = s
        keys = hashed_summits.keys()
        keys.sort()
        gene_peaks[gene] = []
        for d in keys:
            gene_peaks[gene].append( hashed_summits[d] )    
    return gene_peaks
    
def write_pk2gn(chr_sumsite_stats, repid):
    """Writes a table with one row for each summit, and columns for the nearest 3' and 5' target genes."""
    fout = open(outname + ".rep" + repid.__str__() + ".pk2gn.txt", "w")
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            fout.write( chr + "\t" + sumsite.__str__())
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][4].__str__()            )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][0].__str__()            )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][1].__str__()            )
            fout.write( "\t -" + chr_sumsite_stats[chr][sumsite][2].__str__()          )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][3].__str__().__str__() + "\n"  )

    fout.close()

def write_gn2pk(gene_peaks, repid):
    """Writes a table with one row for each gene, and columns for the nearest 3' and 5' peaks, with scores for those peaks."""
    fout = open(outname + ".rep" + repid.__str__() + ".gn2pk.txt", "w")
    for gene in gene_peaks:
        # find the maximum peak to the this gene
        fout.write(gene + "\t")
        fout.write(gene_peaks[gene].__len__().__str__() + "\t")
        for summit in gene_peaks[gene]:
            #fout.write(summit[0].__str__() + "\t")
            fout.write(summit[1].__str__() + "\t")
            fout.write( summit[2].__str__() + "\t")
        fout.write("\n")
    fout.close()

def resolve_replicates(replicate_gn2pk):
    gn2pk = {}
    rep_genes = {}
    for rep in replicate_gn2pk:
        rep_genes[rep] = replicate_gn2pk[rep].keys()

    good_genes = [] # genes with peaks in all reps
    bad_genes = [] # genes lacking a peak in one or more reps
    for rep in rep_genes:
        for gene in rep_genes[rep]:
            if gene in good_genes:
                continue
            elif gene in bad_genes:
                continue
            else: # we haven't seen this gene yet
                missing = False
                for rep in rep_genes:
                    if gene not in rep_genes[rep]:
                        missing = True
                        bad_genes.append( gene )
                        break
                if missing == False:
                    good_genes.append( gene )
    
    gene_best_scores = {} # key = gene, value = array of scores for this rep
    for gene in (good_genes + bad_genes):
        gene_best_scores[gene] = []
        for rep in replicate_gn2pk:
            if gene in replicate_gn2pk[rep]:
                gene_best_scores[gene].append( find_max_peak(replicate_gn2pk[rep][gene]) )
        
    good_genes_scores = []
    fout = open(outname + ".genes.union.txt", "w")
    for gene in good_genes:
        fout.write( gene + "\t" + "\t".join( map(str,gene_best_scores[gene]) ) + "\n")
        good_genes_scores.append( max(gene_best_scores[gene]) )
    fout.close()
    
    bad_genes_scores = []
    fout = open(outname + ".genes.disunion.txt", "w")
    for gene in bad_genes:
        fout.write(gene + "\t" + "\t".join( map(str, gene_best_scores[gene]) ) +  "\n")
        bad_genes_scores.append( max(gene_best_scores[gene]) )
    fout.close()

    fout = open(outname + ".peak_replicates.txt", "w")
    fout.write("=============================================\n")
    fout.write("Summary for " + outname + "\n\n")
    fout.write("I found " + (good_genes.__len__() + bad_genes.__len__()).__str__() + " total ChIP-Seq peaks.\n.")
    fout.write(good_genes.__len__().__str__() + " have ChIP peaks in all replicates.\n")
    fout.write(bad_genes.__len__().__str__() + " don't have peaks in all replicates.\n")
    fout.write("Union genes summit score, mean= " + "%.3f"%mean(good_genes_scores) + "sd=" + "%.3f"%sd(good_genes_scores) + "\n")
    fout.write("Disunion genes summit score, mean= " + "%.3f"%mean(bad_genes_scores) + "sd=" + "%.3f"%sd(bad_genes_scores) + "\n")
    fout.write("\n")
    fout.close()

    print "\n." + good_genes.__len__().__str__() + " have ChIP peaks in all replicates."
    print "\n." + bad_genes.__len__().__str__() + " don't have peaks in all replicates."
    print "\n. Union genes summit score, mean= ", "%.3f"%mean(good_genes_scores), "sd=", "%.3f"%sd(good_genes_scores)
    print "\n. Disunion genes summit score, mean= ", "%.3f"%mean(bad_genes_scores), "sd=", "%.3f"%sd(bad_genes_scores)

def find_max_peak(summits):
    """Returns the summit score for the best peak in the list. THe list 'summits' is a 5-tuple, created in the method read_replicate."""
    max_score = 0
    for s in summits:
        if float(s[2]) > max_score:
            max_score = float(s[2])
    return max_score
#
# Plot stats about the peak-gene mappings
#
def make_plots(chr_sumsite_stats, repid):
    up_distances = []
    down_distances = []
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            if chr_sumsite_stats[chr][sumsite][0] != None:
                up_distances.append( int( chr_sumsite_stats[chr][sumsite][0]) )
            if chr_sumsite_stats[chr][sumsite][2] != None:
                down_distances.append( int( chr_sumsite_stats[chr][sumsite][2]) )
    
    # for debugging:
    #print up_distances
    #print down_distances
    #exit()
    
    plot_histogram(up_distances + down_distances, outname + ".rep" + repid.__str__() + ".dist", xlab="min. distance to gene", ylab="proportion of summits",)    
    
    distances = []
    scores = []
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            upd = None
            if chr_sumsite_stats[chr][sumsite][0] != None:
                upd = int( chr_sumsite_stats[chr][sumsite][0])
            dd = None
            if chr_sumsite_stats[chr][sumsite][2] != None:
                dd = int( chr_sumsite_stats[chr][sumsite][2])
            this_d = None
            if dd == None and upd != None:
                this_d = upd
            elif dd != None and upd == None:
                this_d = dd
            elif dd != None and upd != None:
                if abs(dd) < upd:
                    this_d = dd
                else:
                    this_d = upd
            distances.append( this_d )
            scores.append( float( chr_sumsite_stats[chr][sumsite][4]) )
    scatter1(distances, scores, outname + ".rep" + repid.__str__() + ".dist-v-score", xlab="min. distance to gene", ylab="summit score")

####################
#
# main
chr_gene_sites = read_gff(gffpath)

replicate_pk2gn = {}
replicate_gn2pk = {} 
for ii in range(0, summitpaths.__len__() ):
    spath = summitpaths[ii]
    
    """Parse data and build peak-2-gene library."""
    replicate_pk2gn[spath] = build_pk2gn( spath, chr_gene_sites )
    write_pk2gn( replicate_pk2gn[spath], ii )
    
    """Also build gene-2-peak library."""
    replicate_gn2pk[spath] = build_gn2pk(replicate_pk2gn[spath])
    write_gn2pk( replicate_gn2pk[spath], ii)
    
    make_plots(replicate_pk2gn[spath], ii)

resolve_replicates( replicate_gn2pk )
    


                