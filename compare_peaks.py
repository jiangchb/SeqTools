#
# Create a master table of gene-peak occurances across multiple species and replicate conditions
#

import os, sys, re

from plot_venn import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

runids = ap.getList("--runids") # a list of runid names
runid = ap.getArg("--runid") # the (output) runid for this peak comparison
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
    
    #
    # Fill these data structures...
    #
    seen_genes = [] # a list of orfnames that we've seen in any context
    runid_rep_gene_pkcount = {} # key = runid, value = hash, key = rep, value = hash, key = gene (orfname), value = count of peaks for this gene

    #
    # gn2pk.txt
    #
    for t in runid_reps:
        if t not in runid_rep_gene_pkcount:
            runid_rep_gene_pkcount[t] = {}
        for rep in runid_reps[t]:            
            if rep not in runid_rep_gene_pkcount[t]:
                runid_rep_gene_pkcount[t][rep] = {}
            gn2pk_path = t + ".rep" + rep.__str__() + ".gn2pk.txt"
            print "\n. Reading the gene-2-peak file", gn2pk_path
            fin = open(gn2pk_path)
            lines = fin.readlines()
            for l in lines:
                if l.__len__() > 2:
                    tokens = l.split()
                    orfname = get_orf(tokens[0])
                    if orfname not in seen_genes:
                        seen_genes.append( orfname )
                    peak_count = int( tokens[1] )
                    runid_rep_gene_pkcount[t][rep][orfname] = peak_count

    #
    # Write peak_count...xls 
    #
    outpath = "peak_count." + runid + ".xls"
    print "\n. Writing", outpath
    columns = []
    for t in runid_reps:
        for rep in runid_reps[t]:
            columns.append(t + ":" + rep.__str__())
    fout = open(outpath, "w")
    fout.write(".\t")
    fout.write( "\t".join(columns) )
    fout.write("\n") 
    runid_genes = {} # key = runid, value = list of gene IDs that are in all replicates of the runid ID
    runid_rep_genes = {} # key = runid, value = hash, key = rep, value = list of genes in that replicate
    for ii in range(0, seen_genes.__len__()):
        gene = seen_genes[ii]
        fout.write(gene + "\t")
        for t in runid_reps:
            if t not in runid_genes:
                runid_genes[t] = []
                runid_rep_genes[t] = {}
            row_sum = 0.0
            rep_sum = 0.0
            for rep in runid_reps[t]:
                if "rep" + rep.__str__() not in runid_rep_genes[t]:
                    runid_rep_genes[t]["rep" + rep.__str__()] = []
                if gene in runid_rep_gene_pkcount[t][rep]:
                    fout.write( "%d"%runid_rep_gene_pkcount[t][rep][gene] + "\t" )
                    row_sum += runid_rep_gene_pkcount[t][rep][gene]
                    rep_sum += 1.0 / runid_reps[t].__len__()
                    runid_rep_genes[t]["rep" + rep.__str__()].append( ii )
                else:
                    fout.write( "0\t")
            if rep_sum >= 1.0: # eliminate peaks that don't exist in all replicates
                runid_genes[t].append( ii )
        fout.write("\n")    
    fout.close()

    #
    # Write summarized, new, gn2pk.txt
    #
    gn2pk_outpath = runid + ".rep0.gn2pk.txt"
    fout = open(gn2pk_outpath, "w")
    for ii in range(0, seen_genes.__len__()):
        in_all = True
        for t in runid_genes:
            if ii not in runid_genes[t]:
                in_all = False
                break
        if in_all:
            fout.write( seen_genes[ii] + "\t1\n")
    fout.close()

    #
    # Plot Venn Diagrams
    #
    for t in runid_reps:
        plot_venn_diagram(runid_rep_genes[t], t)
    plot_venn_diagram(runid_genes, runid)

#########################################
#
# main
#
read_pillars()
runid_reps = get_runid_reps()
print runid_reps
make_output(runid_reps)
