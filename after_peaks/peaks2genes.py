#
# This script maps ChIP-Seq summits onto regulatory sequences
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

# summitpath is the filepath to a *_summits.bed
summitpaths = ap.getList("--summits")

# peakpaths is a list of *.encodePeak files
peakpaths = ap.getOptionalList("--peaks")

# the file containing features of the genome
gffpath = ap.getArg("--gff")

# runid will be used as a keyword in output file names.
runid = ap.getArg("--runid") # A unique name for the output

#
#
#################################################################


def splash():
    print "============================================"
    print "  summits2genes.py"
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
            if strand == "+": # sense strand
                chr_gene_sites[chr][gene] = (start, stop)
            else: # antisense strand
                chr_gene_sites[chr][gene] = (stop, start)
    fin.close() 
    
    count_genes = 0
    for chr in chr_gene_sites:
        count_genes += chr_gene_sites[chr].keys().__len__()
    print "\n. I found", count_genes, "genes " in gffpath + "."
    
    return chr_gene_sites
    
def build_sm2gn(summitpath, chr_gene_sites):
    """Returns chr_sumsite_stats."""
    #
    # Build a library of summits
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
        else:
            print "91:", l
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
                
                #
                # to-do: there is a bug here. See notes from Eugenio July 2014
                #
                
                """Sense direction, and upstream"""
                if start < stop and start >= sumsite:
                    if min_up == None:
                        min_up = d
                        closest_up = gene
                    if min_up > d:
                        min_up = d
                        closest_up = gene
                        
                elif start < stop and stop < sumsite:
                    """Sense direction and downstream"""
                    if min_down == None:
                        min_down = d
                        closest_down = None
                    if min_down < d:
                        min_down = d
                        closest_down = None
                    
                elif start > stop and start <= sumsite:
                    """Antisense and downstream"""
                    if min_down == None:
                        min_down = d
                        closest_down = gene
                    if min_down < d: # be careful here, we're looking for the largest NEGATIVR number.
                        min_down = d
                        closest_down = gene
                        
                elif start > stop and stop > sumsite:
                    """Antisense and upstream"""
                    if min_up == None:
                        min_up = d
                        closest_up = None
                    if min_up > d:
                        min_up = d
                        closest_up = None
                        
                """Note: this ignores the cases where the summit is inside a
                coding region. It also ignores the cases where the nearest gene
                is in the wrong sense orientation, such that the gene's stop site
                is placed between the start site and the summit."""
            if closest_up == None:
                min_up = None
            if closest_down == None:
                min_down = None
            #print chr, sumsite, chr_site_score[chr][sumsite], closest_up, min_up, closest_down, min_down
            chr_sumsite_stats[chr][sumsite] = (min_up, closest_up, min_down, closest_down, chr_site_score[chr][sumsite])    
    return chr_sumsite_stats

def build_gn2sm(chr_sumsite_stats):
    """Returns gene_summits: gene_summits[gene] = [ (summit site, distance to gene, summit score), (n2), (n3) ]"""
    gene_summits = {}
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            upd = chr_sumsite_stats[chr][sumsite][0]
            upgene = chr_sumsite_stats[chr][sumsite][1]
            downd = chr_sumsite_stats[chr][sumsite][2]
            downgene = chr_sumsite_stats[chr][sumsite][3]
            if upgene != None and upd != None:
                if upgene not in gene_summits:
                    gene_summits[ upgene ] = [] 
                gene_summits[ upgene ].append(  (sumsite, -1*int(upd), chr_sumsite_stats[chr][sumsite][4])  )       
            if downgene != None and downd != None:
                if downgene not in gene_summits:
                    gene_summits[ downgene ] = [] 
                gene_summits[ downgene ].append(  (sumsite, -1*int(downd), chr_sumsite_stats[chr][sumsite][4])  )       
    
    # sort summits by distance to gene
    for gene in gene_summits:
        summits = gene_summits[gene]
        hashed_summits = {}
        for s in summits:
            d = abs( int(s[1]) )
            if d in hashed_summits:
                d += 1
            hashed_summits[d] = s
        keys = hashed_summits.keys()
        keys.sort()
        gene_summits[gene] = []
        for d in keys:
            gene_summits[gene].append( hashed_summits[d] )    
    return gene_summits
    
def write_sm2gn(chr_sumsite_stats, repid):
    """Writes a table with one row for each summit, and columns for the nearest 3' and 5' target genes."""
    fout = open(runid + ".rep" + repid.__str__() + ".sm2gn.txt", "w")
    for chr in chr_sumsite_stats:
        for sumsite in chr_sumsite_stats[chr]:
            fout.write( chr + "\t" + sumsite.__str__())
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][4].__str__()            )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][0].__str__()            )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][1].__str__()            )
            fout.write( "\t " + chr_sumsite_stats[chr][sumsite][2].__str__()          )
            fout.write( "\t" + chr_sumsite_stats[chr][sumsite][3].__str__().__str__() + "\n"  )

    fout.close()

def write_gn2sm(gene_summits, repid):
    """Writes a table with one row for each gene, and columns for the nearest 3' and 5' summits, with scores for those summits."""
    fout = open(runid + ".rep" + repid.__str__() + ".gn2sm.txt", "w")
    for gene in gene_summits:
        # find the maximum summit to the this gene
        fout.write(gene + "\t")
        fout.write(gene_summits[gene].__len__().__str__() + "\t")
        for summit in gene_summits[gene]:
            #fout.write(summit[0].__str__() + "\t")
            fout.write(summit[1].__str__() + "\t")
            fout.write( summit[2].__str__() + "\t")
        fout.write("\n")
    fout.close()

def resolve_replicates(replicate_gn2sm):
    gn2sm = {}
    rep_genes = {}
    for rep in replicate_gn2sm:
        rep_genes[rep] = replicate_gn2sm[rep].keys()

    good_genes = [] # genes with summits in all reps
    bad_genes = [] # genes lacking a summit in one or more reps
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
        for rep in replicate_gn2sm:
            if gene in replicate_gn2sm[rep]:
                gene_best_scores[gene].append( find_max_summit(replicate_gn2sm[rep][gene]) )
        
    good_genes_scores = []
    fout = open(runid + ".genes.union.xls", "w")
    for gene in good_genes:
        fout.write( gene + "\t" + "\t".join( map(str,gene_best_scores[gene]) ) + "\n")
        good_genes_scores.append( max(gene_best_scores[gene]) )
    fout.close()
    
    bad_genes_scores = []
    fout = open(runid + ".genes.disunion.xls", "w")
    for gene in bad_genes:
        fout.write(gene + "\t" + "\t".join( map(str, gene_best_scores[gene]) ) +  "\n")
        bad_genes_scores.append( max(gene_best_scores[gene]) )
    fout.close()

    fout = open(runid + ".summit_replicates.txt", "w")
    fout.write("=============================================\n")
    fout.write("Summary for " + runid + "\n\n")
    fout.write("I found " + (good_genes.__len__() + bad_genes.__len__()).__str__() + " total ChIP-Seq summits.\n")
    fout.write(good_genes.__len__().__str__() + " genes have nearby ChIP summits in all replicates. Mean summit score = " + "%.2f"%mean(good_genes_scores) + " s.d. = " + "%.2f"%sd(good_genes_scores) + "\n")
    fout.write(bad_genes.__len__().__str__() + " genes don't have nearby summits in all replicates. Mean summit score = " + "%.2f"%mean(bad_genes_scores) + " s.d. = " + "%.2f"%sd(bad_genes_scores) + "\n")
    fout.write("\n")
    fout.close()

    print "\n." + good_genes.__len__().__str__() + " have ChIP summits in all replicates."
    print "\n." + bad_genes.__len__().__str__() + " don't have summits in all replicates."
    print "\n. Union genes summit score, mean= ", "%.3f"%mean(good_genes_scores), "sd=", "%.3f"%sd(good_genes_scores)
    print "\n. Disunion genes summit score, mean= ", "%.3f"%mean(bad_genes_scores), "sd=", "%.3f"%sd(bad_genes_scores)
    
    return (good_genes, bad_genes)

def find_max_summit(summits):
    """Returns the summit score for the best summit in the list. THe list 'summits' is a 5-tuple, created in the method read_replicate."""
    max_score = 0
    for s in summits:
        if float(s[2]) > max_score:
            max_score = float(s[2])
    return max_score
#
# Plot stats about the summit-gene mappings
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
    
    plot_histogram(up_distances + down_distances, runid + ".rep" + repid.__str__() + ".dist", xlab="min. distance to gene", ylab="proportion of summits",)    
    
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
            
            if this_d != None:
                distances.append( this_d )
                scores.append(float( chr_sumsite_stats[chr][sumsite][4]) )
    scatter1(distances, scores, runid + ".rep" + repid.__str__() + ".dist-v-score", xlab="min. distance to gene", ylab="summit score")

def read_peaks(ppath):
    """Analyzes the peaks found in an *.encodePeak file from MACS."""
    fin = open(ppath, "r")
    chr_peaksite_stats = {} # key = chromosome, value = hashtable, key = peak start site, value = (stopsite, rank, enrichment, pval, qval)
    for l in fin.xreadlines():
        if l.__len__() > 2:
            tokens = l.split()
            chr = tokens[0]
            start = int(tokens[1])
            stop = int(tokens[2])
            enrichment = float( tokens[6])
            pval = float(tokens[7])
            qval = float(tokens[8])
            if chr not in chr_peaksite_stats:
                chr_peaksite_stats[chr] = {}
            chr_peaksite_stats[chr][start] = (stop, enrichment, pval, qval)
    fin.close()
    return chr_peaksite_stats

def make_peak_plots( chr_peaksite_stats, repid ):
    evals = []
    logpvals = []
    pvals = []
    logqvals = []
    qvals = []
    for chr in chr_peaksite_stats:
        for start in chr_peaksite_stats[chr]:
            eval = chr_peaksite_stats[chr][start][1]
            logpval = chr_peaksite_stats[chr][start][2]
            logqval = chr_peaksite_stats[chr][start][3]

            evals.append( eval )
            logpvals.append( logpval )
            logqvals.append( logqval )
            pvals.append( 10.0**(logpval/-10) )
            qvals.append( 10.0**(logqval/-10) )
    scatter1(evals, logpvals, runid + ".rep" + repid.__str__() + ".e-v-logp", xlab="mean enrichment over peak region", ylab="log(P) value of peak")
    scatter1(evals, logqvals, runid + ".rep" + repid.__str__() + ".e-v-logq", xlab="mean enrichment over peak region", ylab="log(Q) value of peak")
    scatter1(logpvals, logqvals, runid + ".rep" + repid.__str__() + ".logp-v-logq", xlab="log(P) value of the peak", ylab="log(Q) value of peak")
    scatter1(pvals, logpvals, runid + ".rep" + repid.__str__() + ".p-v-logp", xlab="P value", ylab="log(P) value of peak")
    scatter1(qvals, logqvals, runid + ".rep" + repid.__str__() + ".q-v-logq", xlab="Q value", ylab="log(Q) value of peak")

def rank_peaks( gn2sm ):
    """Returns gene_rank"""
    qvals_genes = {}
    for g in gn2sm:
        max_qval = 0.0
        for sm in gn2sm[g]:
            if float(sm[2]) > max_qval:
                max_qval = float(sm[2])
        qval = max_qval
        if qval not in qvals_genes:
            qvals_genes[qval] = [ g ]
        else:
            qvals_genes[qval].append( g )
    qvals = qvals_genes.keys()
    qvals.sort(reverse=True)
    
    gene_rank = {}
    
    rank = 0
    for qval in qvals:
        rank += 1
        for g in qvals_genes[qval]:
            gene_rank[g] = rank
    return gene_rank

def calculate_gamma_rank(replicate_gene_rank, genes):
    """See Li 2011"""
    xranks = []
    yranks = []
    for g in genes:      
        xranks.append( replicate_gene_rank[ 0 ][g] )
        yranks.append( replicate_gene_rank[ 1 ][g] )
    scatter1(xranks, yranks, runid + ".rank", xlab="Gene Rank, rep1", ylab="Gene Rank, rep2")
    
####################
#
# main
chr_gene_sites = read_gff(gffpath)
replicate_sm2gn = {}
replicate_gn2sm = {}
for ii in range(0, summitpaths.__len__() ):
    spath = summitpaths[ii]
    """Parse data and build summit-2-gene library."""
    replicate_sm2gn[spath] = build_sm2gn( spath, chr_gene_sites )
    write_sm2gn( replicate_sm2gn[spath], ii )
    """Also build gene-2-summit library."""
    replicate_gn2sm[spath] = build_gn2sm(replicate_sm2gn[spath])
    write_gn2sm( replicate_gn2sm[spath], ii)    
    make_plots(replicate_sm2gn[spath], ii)

(good_genes, bad_genes) = resolve_replicates( replicate_gn2sm )
replicate_gene_rank = {}
for ii in range(0, summitpaths.__len__() ):
    spath = summitpaths[ii]
    replicate_gene_rank[ii] = rank_peaks( replicate_gn2sm[spath] )
calculate_gamma_rank(replicate_gene_rank, good_genes)
    
if peakpaths:
    replicate_peaks = {}

    for ii in range(0, peakpaths.__len__()):
        ppath = peakpaths[ii]
        replicate_peaks[ii] = read_peaks(ppath)
        make_peak_plots( replicate_peaks[ii], ii)


                