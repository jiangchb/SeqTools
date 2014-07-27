#############################################################################
#
# A collection of functions to generate plots and assorted visualizations
# for data inside a chipseqdb SqlLite3 database.
#
#############################################################################

from chipseqdb import *
from plot_scatter import *
from plot_venn import *

def correlate_two_reps(repid1, repid2, con):
    print "\n. Comparing replicates", repid1, "and", repid2
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where id=" + repid1.__str__())
    data = cur.fetchone()
    
    ##
#     print "19:", data
#     cur.execute("SELECT species from Replicates")
#     data = cur.fetchall()
    #print "22:", data
    ##
    
    species1 = data[0]
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where id=" + repid2.__str__())
    species2 = cur.fetchone()[0]    
    if species1 != species2:
        return None
    
    """Venn diagram of genes with/without summits in both replicates."""
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]

    venn_data = {}
    venn_data[ rep1name ] = []
    venn_data[ rep2name ] = []
    x = get_geneids(con, repid1)
    for ii in x:
        venn_data[ rep1name ].append(ii[0])
    y = get_geneids(con, repid2)
    for ii in y:
        venn_data[ rep2name ].append(ii[0])
    plot_venn_diagram( venn_data, "genes." + rep1name.__str__() + "." + rep2name.__str__())
    
    
    """ Plot score correlations."""
    xvals = []
    yvals = []
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            x = get_max_summit_score_for_gene(gid, repid1, con)
            y = get_max_summit_score_for_gene(gid, repid2, con)
            if x != None and y != None:
                xvals.append(x[0])
                yvals.append(y[0])
    scatter1(xvals,yvals,"corr.score.rep" + repid1.__str__() + ".rep" + repid2.__str__(), xlab="Replicate 1: max summit score for gene", ylab="Replicate 2: max summit score for gene")

    """Plot rank correlations."""
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
         
        
        