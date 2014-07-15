#
# Create a master table of gene-peak occurances across multiple species and replicate conditions
#

import os, sys, re

from plot_venn import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

runids = ap.getList("--runids") # a list of runid names
runid = ap.getArg("--runid") # the runid for this peak comparison
pillarspath = ap.getArg("--pillars") # a file containing G
orf_names = {}

def splash():
    print "============================================"
    print "  compare_peaks.py"
    print "  Victor Hanson-Smith, 2014"
    print "  victor.hanson-smith@ucsf.edu"
    print "============================================"
    
splash()

def get_runid_reps():
    print "\n. Reading Replicates..."
    runid_reps = {}
    for t in runids:
        if t not in runid_reps:
            runid_reps[t] = []
        repid = 0
        foundit = True
        while(foundit == True):
            if os.path.exists(t + ".rep" + repid.__str__() + ".gn2pk.txt"):
                runid_reps[t].append( repid )
                repid += 1
            else:
                foundit = False
    return runid_reps

def read_pillars():
    print "\n. Reading Pillars..."
    fin = open(pillarspath, "r")
    for l in fin.xreadlines():
        if l.startswith("orf"):
            l = l.strip()
            tokens = l.split()
            orf_list = []
            for t in tokens:
                if False == t.startswith("-"):
                    orf_list.append(t)
            orf_names[ orf_list[0] ] = orf_list

def get_orf(genename):
    for g in [genename, genename.upper()]:
        for orf in orf_names:
            if orf == g:
                return orf
            if orf_names[orf].__contains__(g):
                return orf
    return genename

def make_output(runid_reps):
    print "\n. Building the output matrix..."
    """Also plots a venn diagram"""
    seen_genes = []
    runid_rep_gene_pkcount = {}
    columns = []
    for t in runid_reps:
        if t not in runid_rep_gene_pkcount:
            runid_rep_gene_pkcount[t] = {}
        for rep in runid_reps[t]:
            print "\n-->", t, rep
            columns.append(t + ":" + rep.__str__())
            
            if rep not in runid_rep_gene_pkcount[t]:
                runid_rep_gene_pkcount[t][rep] = {}
            fin = open(t + ".rep" + rep.__str__() + ".gn2pk.txt")
            lines = fin.readlines()
            for l in lines:
                if l.__len__() > 2:
                    tokens = l.split()
                    orfname = get_orf(tokens[0])
                    if orfname not in seen_genes:
                        seen_genes.append( orfname )
                    peak_count = int( tokens[1] )
                    runid_rep_gene_pkcount[t][rep][orfname] = peak_count

    outpath = "peak_count." + runid + ".xls"
    print "\n. Writing", outpath
    fout = open(outpath, "w")
    fout.write(".\t")
    fout.write( "\t".join(columns) )
    fout.write("\n")
    
    #print "87: seen_genes", seen_genes
    
    runid_genes = {} # key = runid ID, value = list of gene IDs that are in all replicates of the runid ID
    for ii in range(0, seen_genes.__len__()):
        gene = seen_genes[ii]
        fout.write(gene + "\t")
        for t in runid_reps:
            if t not in runid_genes:
                runid_genes[t] = []
            row_sum = 0.0
            rep_sum = 0.0
            for rep in runid_reps[t]:
                if gene in runid_rep_gene_pkcount[t][rep]:
                    fout.write( "%d"%runid_rep_gene_pkcount[t][rep][gene] + "\t" )
                    row_sum += runid_rep_gene_pkcount[t][rep][gene]
                    rep_sum += 1.0 / runid_reps[t].__len__()
                else:
                    fout.write( "0\t")
            if rep_sum >= 1.0:
                runid_genes[t].append( ii )
        fout.write("\n")    
    fout.close()

    plot_venn_diagram(runid_genes, runid)
    
#     in_all = [] # a list of genes with peaks in all runids:
#     for ii in range(0, seen_genes.__len__()):
#         foundit = True
#         for t in runid_genes:
#             if ii not in runid_genes[t]:
#                 foundit = False
#                 break
#         if foundit == True:
#             in_all.append( ii )
    

    #tt_overlap = compute_venn( runid_genes, in_all )
    #print "\n. 124:", tt_overlap
    #print "\n. 125:", in_all
    #plot_venn_diagram( tt_overlap, in_all.__len__(), "peak_count.venn." + runid)
    
    
#
# depricated
#
def compute_venn( runid_genes, in_all ):
    """This method prepares data for the method plot_venn_diagram.
    runid_genes[runid] = list of genes with peaks in this runid
    in_all = list of genes with peaks in all runids (i.e., the center of the Venn diagram."""
    runids = runid_genes.keys()
    runids.sort()
    
    """Compute all the 2-way overlaps"""
    runid_runid_overlap = {} # key = runid, value = hash, key = runid, value = count of overlap
    for t in runid_genes:
        runid_runid_overlap[t] = {}
        for s in runid_genes:
            count = 0
            for ii in runid_genes[t]:
                if ii in runid_genes[s]:
                    count += 1
            runid_runid_overlap[t][s] = count
    outpath = "peak_count.venn_table." + runid + ".xls"
    print "\n. Writing", outpath
    fout = open(outpath, "w")
    fout.write(".\t")
    fout.write( "\t".join( runids ) )
    fout.write("\n")
    for t in runids:
        fout.write(t + "\t")
        for s in runids:
            fout.write( runid_runid_overlap[t][s].__str__() + "\t")
        fout.write("\n")
    fout.close() 
    
    return runid_runid_overlap
    
#########################################
#
# main
#
read_pillars()
runid_reps = get_runid_reps()
print runid_reps
make_output(runid_reps)
