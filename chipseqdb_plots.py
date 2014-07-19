#############################################################################
#
# A collection of functions to generate plots and assorted visualizations
# for data inside a chipseqdb SqlLite3 database.
#
#############################################################################

from chipseqdb import *
from plot_scatter import *

def correlate_two_reps(repid1, repid2, con):
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where count=" + repid1.__str__())
    data = cur.fetchone()
    species1 = data[0]
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where count=" + repid2.__str__())
    species2 = cur.fetchone()[0]    
    if species1 != species2:
        return None
    
    xvals = []
    yvals = []
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        #print "25:", chr[0]
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            x = get_max_summit_score_for_gene(gid, repid1, con)
            y = get_max_summit_score_for_gene(gid, repid2, con)
            if x != None and y != None:
                #print x, y
                xvals.append(x[0])
                yvals.append(y[0])
    scatter1(xvals,yvals,"corr.score.rep" + repid1.__str__() + ".rep" + repid2.__str__(), xlab="Replicate 1: max summit score for gene", ylab="Replicate 2: max summit score for gene")

    scoresx_genes = {}
    scoresy_genes = {}
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            x = get_max_summit_score_for_gene(gid, repid1, con)
            y = get_max_summit_score_for_gene(gid, repid2, con)
            if x != None and y != None:
                if x[0] not in scoresx_genes:
                    scoresx_genes[ x[0] ] = []
                scoresx_genes[ x[0] ].append( gid )
                
                if y[0] not in scoresy_genes:
                    scoresy_genes[ y[0] ] = []
                scoresy_genes[ y[0] ].append( gid )
    xkeys = scoresx_genes.keys()
    xkeys.sort( reverse=True )
    ykeys = scoresy_genes.keys()
    ykeys.sort( reverse=True )
    
    gene_xrank = {}
    gene_yrank = {}
    xrank = 1
    for x in xkeys:
        for gene in scoresx_genes[x]:
            gene_xrank[gene] = xrank
            xrank += 1
    yrank = 1
    for y in ykeys:
        for gene in scoresy_genes[y]:
            gene_yrank[gene] = yrank
            yrank += 1
    xvals = []
    yvals = []
    for gene in gene_xrank:
        xvals.append( gene_xrank[gene] )
        yvals.append( gene_yrank[gene] )
            
    scatter1(xvals,yvals,"corr.rank.rep" + repid1.__str__() + ".rep" + repid2.__str__(), xlab="Replicate 1: rank order of genes by max summit score", ylab="Replicate 2: rank order of genes by max summit score")
         
        
        