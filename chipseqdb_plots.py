#############################################################################
#
# A collection of functions to generate plots and assorted visualizations
# for data inside a chipseqdb SqlLite3 database.
#
#############################################################################

from chipseqdb import *
from plot_scatter import *
from plot_venn import *

def correlate_reps_in_group(rgroupid, con):
    
    repids = get_reps_in_group(rgroupid, con)
    repid1 = repids[0][0]
    repid2 = repids[1][0]
    
    """This is a long method (sorry). It correlates the summits in two replicates,
    using a variety of methods."""
    cur = con.cursor()
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]
    
    print "\n. Comparing replicates", rep1name, "and", rep2name
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where id=" + repid1.__str__())
    data = cur.fetchone()
        
    species1 = data[0]

    cur.execute("SELECT species from Replicates where id=" + repid2.__str__())
    species2 = cur.fetchone()[0]    
    if species1 != species2:
        return None
    
    """Venn diagram of genes with/without summits in both replicates."""
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
    
 
    """For the following plots, we'll use these data structures."""   
    print "\n. Building a dictionary of summits & genes for", rep1name, "and", rep2name
    rep1_gene_summitscores = {} # key = geneid, value = summit scores from rep1
    rep2_gene_summitscores = {}
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            x = get_summit_scores_for_gene(gid, repid1, con)
            if x.__len__() > 0:
                rep1_gene_summitscores[gid] = x
            y = get_summit_scores_for_gene(gid, repid2, con)
            if y.__len__() > 0:
                rep2_gene_summitscores[gid] = y
    
    
    """Write an Excel table with genes and their scores in each replicate."""
    xlpath = "genes." + rep1name.__str__() + "." + rep2name.__str__() + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    fout.write("N(" + rep1name.__str__() + ")\t")
    fout.write("mean(" + rep1name.__str__() + ")\t")
    #fout.write("mode(" + rep1name.__str__() + ")\t")
    fout.write("max(" + rep1name.__str__() + ")\t")
    fout.write("min(" + rep1name.__str__() + ")\t")
    fout.write("N(" + rep2name.__str__() + ")\t")
    fout.write("mean(" + rep2name.__str__() + ")\t")
    #fout.write("mode(" + rep2name.__str__() + ")\t")
    fout.write("max(" + rep2name.__str__() + ")\t")
    fout.write("min(" + rep2name.__str__() + ")\n")
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            """Ignore genes that aren't in at least on replicate."""
            if gid not in rep1_gene_summitscores and gid not in rep2_gene_summitscores:
                continue # skip to the next gene
            
            elif gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
                # This summit is in both replicates:
                cur.execute("SELECT COUNT(*) from RepgroupGenes where repgroupid=" + rgroupid.__str__() + " and geneid=" + gid.__str__() )
                if cur.fetchone()[0] == 0: # if this entry doesn't already exist...
                    sql = "INSERT into RepgroupGenes (repgroupid, geneid) VALUES(" + rgroupid.__str__() + "," + gid.__str__() + ")"
                    cur.execute(sql)
                    con.commit()
            
            fout.write(g[0].__str__() + "\t" + g[1].__str__() + "\t")
            if gid in rep1_gene_summitscores:
                thisn = rep1_gene_summitscores[gid].__len__()
                if thisn == 1:
                    thismode = rep1_gene_summitscores[gid][0]
                    thismean = rep1_gene_summitscores[gid][0]
                else:
                    thismean = stats.tmean( rep1_gene_summitscores[gid] ) # the [0] index is because SciPy returns an array of results
                    thismode = stats.mode( rep1_gene_summitscores[gid] )[0]
                thismax = max( rep1_gene_summitscores[gid] )
                thismin = min( rep1_gene_summitscores[gid] )
                #print "81:", thisn, thismean, thismode, thismax, thismin
                fout.write("%d"%thisn + "\t")
                fout.write("%.3f"%thismean + "\t")
                #fout.write("%.3f"%thismode + "\t")
                fout.write("%.3f"%thismax + "\t")
                fout.write("%.3f"%thismin + "\t")
            else:
                fout.write("0\t0\t0\t0\t")               
            
            if gid in rep2_gene_summitscores:
                #print "93:", gid
                #print "97:", rep2_gene_summitscores[gid]
                #print "98:", type( rep2_gene_summitscores[gid] )
                thisn = rep2_gene_summitscores[gid].__len__()
                if thisn == 1:
                    thismode = rep2_gene_summitscores[gid][0]
                    thismean = rep2_gene_summitscores[gid][0]
                else:
                    thismean = stats.tmean( rep2_gene_summitscores[gid] )
                    thismode = stats.mode( rep2_gene_summitscores[gid] )[0]
                thismax = max( rep2_gene_summitscores[gid] )
                thismin = min( rep2_gene_summitscores[gid] )
                #print "95:", thisn, thismean, thismode, thismax, thismin
                fout.write("%d"%thisn + "\t")
                fout.write("%.3f"%thismean + "\t")
                #fout.write("%.3f"%thismode + "\t")
                fout.write("%.3f"%thismax + "\t")
                fout.write("%.3f"%thismin + "\n")
            else:
                fout.write("0\t0\t0\t0\n")
    fout.close()

    
    """ Plot score correlations."""
    xvals = []
    yvals = []
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            #x = get_max_summit_score_for_gene(gid, repid1, con)
            #y = get_max_summit_score_for_gene(gid, repid2, con)
            #if x != None and y != None:
            if gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
                xvals.append( max(rep1_gene_summitscores[gid]) )
                yvals.append( max(rep2_gene_summitscores[gid]) )
    scatter1(xvals,yvals,"corr.score." + rep1name.__str__() + "." + rep2name.__str__(), xlab="Replicate 1: max summit score for gene", ylab="Replicate 2: max summit score for gene")

    """Plot rank correlations."""
    scoresx_genes = {} # bins of genes, organized by their score range
    scoresy_genes = {}
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes(con, chr[0])
        for g in genes:
            gid = g[0]
            if gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
                x = max(rep1_gene_summitscores[gid])
                y = max(rep2_gene_summitscores[gid])

                if x not in scoresx_genes:
                    scoresx_genes[ x ] = []
                scoresx_genes[ x ].append( gid )
                
                if y not in scoresy_genes:
                    scoresy_genes[ y ] = []
                scoresy_genes[ y ].append( gid )
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
    scatter1(xvals,yvals,"corr.rank." + rep1name.__str__() + "." + rep2name.__str__(), xlab="Replicate 1: rank order of genes by max summit score", ylab="Replicate 2: rank order of genes by max summit score")
         
        
        