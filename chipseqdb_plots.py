#############################################################################
#
# A collection of functions to generate plots and assorted visualizations
# for data inside a chipseqdb SqlLite3 database.
#
#############################################################################

from chipseqdb import *
from plot_scatter import *
from plot_venn import *
from plot_histogram import *

def correlate_enrichments_for_species(speciesid, con):
    cur = con.cursor()
      


def correlate_enrichments_for_union(unionid, con, keyword=None):
    
    # depricated:
    #compute_enrichments_for_union(unionid, con, keyword=keyword)
    
    cur = con.cursor()
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    rgroupids.sort()
    
    #
    # What gene IDs apply to this union?
    #
    sql = "SELECT geneid from GroupEnrichmentStats where"
    sqlbits = []
    for rgid in rgroupids:
        sqlbits.append( " rgroupid=" + rgid.__str__() )
    sql += " OR ".join( sqlbits )
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
        
    
    #
    # Write an excel table,
    # while also filling data structures for the scatterplot
    #
    rgid_maxvals = {}
    rgid_meanvals = {}
    rgid_sumvals = {}
    
    xlpath = "enrich." + unionname.__str__() + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for rgid in rgroupids:
        rgname = get_repgroup_name(rgid, con)
        fout.write("max(" + rgname + ")\t")
        fout.write("mean(" + rgname + ")\t")
        fout.write("sum(" + rgname + ")\t")

        rgid_maxvals[rgid] = []
        rgid_meanvals[rgid] = []
        rgid_sumvals[rgid] = []
    
    count = 0
    fout.write("\n")
    for geneid in genes:
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        for rgid in rgroupids:
            count += 1
            if count%50 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            
            sql = "SELECT * from GroupEnrichmentStats where rgroupid=" + rgid.__str__() + " and geneid=" + geneid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            if x != None:
                fout.write(x[2].__str__() + "\t")
                fout.write(x[3].__str__() + "\t")
                fout.write(x[4].__str__() + "\t")

                rgid_maxvals[rgid].append(x[2])
                rgid_meanvals[rgid].append(x[3])
                rgid_sumvals[rgid].append(x[4])
                    
            else:
                fout.write("0\t0\t0\t")
                rgid_maxvals[rgid].append(0)
                rgid_meanvals[rgid].append(0)
                rgid_sumvals[rgid].append(0)
        fout.write("\n")
    fout.close()
    
    #
    # Scatterplots
    #
    for ii in range(0, rgroupids.__len__() ):
        for jj in range(ii+1, rgroupids.__len__() ):
            #print "90:", ii, jj

            this_rgid = rgroupids[ii]
            that_rgid = rgroupids[jj]
            
            this_rgname = get_repgroup_name(this_rgid, con)
            that_rgname = get_repgroup_name(that_rgid, con)
            
            scatter1(rgid_maxvals[this_rgid], rgid_maxvals[that_rgid], "enrich.max." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": max enrichment", ylab=that_rgname.__str__()+": max enrichment", force_square=True)
            scatter1(rgid_meanvals[this_rgid], rgid_meanvals[that_rgid], "enrich.mean." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": mean enrichment", ylab=that_rgname.__str__()+": mean enrichment", force_square=True)
            scatter1(rgid_sumvals[this_rgid], rgid_sumvals[that_rgid], "enrich.sum." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": sum enrichment", ylab=that_rgname.__str__()+": sum enrichment", force_square=True)
            

def correlate_enrichments_for_reps_in_group(rgroupid, con):
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
    cur = con.cursor()
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]
    
    x_maxe = []
    y_maxe = []

    x_meane = []
    y_meane = []
    
    x_sume = []
    y_sume = []
    
    x = get_enrichment_stats_for_replicate( repid1, con )
    y = get_enrichment_stats_for_replicate( repid2, con )
    
    geneids_all = x.keys()
    geneids_all.sort()
    
    geneids = []
    for geneid in geneids_all:
        if geneid in y:
            geneids.append( geneid )
    
    for geneid in geneids:
        if geneid in y:
            x_maxe.append( x[geneid][0] )
            y_maxe.append( y[geneid][0] )
            x_meane.append( x[geneid][1] )
            y_meane.append( y[geneid][1] )
            x_sume.append( x[geneid][2] )
            y_sume.append( y[geneid][2] )
    
    scatter1(x_maxe,y_maxe,"enrich.max." + rep1name.__str__() + "." + rep2name.__str__(), xlab=rep1name.__str__()+": max enrichment", ylab=rep2name.__str__()+": max enrichment", force_square=True)
    scatter1(x_meane,y_meane,"enrich.mean." + rep1name.__str__() + "." + rep2name.__str__(), xlab=rep1name.__str__()+": mean enrichment", ylab=rep2name.__str__()+": mean enrichment", force_square=True)
    scatter1(x_sume,y_sume,"enrich.sum." + rep1name.__str__() + "." + rep2name.__str__(), xlab=rep1name.__str__()+": sum enrichment", ylab=rep2name.__str__()+": sum enrichment", force_square=True)

    #
    # Write an Excel Table
    #
    xlpath = "enrich." + rep1name + "." + rep2name + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for rgname in [rep1name,rep2name]:
        fout.write("max(" + rgname + ")\t")
        fout.write("mean(" + rgname + ")\t")
        fout.write("sum(" + rgname + ")\t")        
    fout.write("\n")
    for ii in range(0, geneids.__len__() ):
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        fout.write(x_maxe[ii].__str__() + "\t")
        fout.write(x_meane[ii].__str__() + "\t")
        fout.write(x_sume[ii].__str__() + "\t")
        fout.write(y_maxe[ii].__str__() + "\t")
        fout.write(y_meane[ii].__str__() + "\t")
        fout.write(y_sume[ii].__str__())
        fout.write("\n")
    fout.close()

    print "\n. Updating the table GroupEnrichmentStats"
    for ii in range(0, geneids.__len__() ):
        geneid = geneids[ii]
        sql = "SELECT count(*) from GroupEnrichmentStats where rgroupid=" + rgroupid.__str__()
        sql += " and geneid=" + geneid.__str__()
        cur.execute(sql)
        x = cur.fetchone()[0]
        if x == 0:  
            sql = "INSERT into GroupEnrichmentStats (rgroupid, geneid, maxenrich, meanenrich, sumenrich)"
            sql += " VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ","
            sql += max( [x_maxe[ii], y_maxe[ii]] ).__str__() + ","
            sql += mean( [x_meane[ii], y_meane[ii]] ).__str__() + ","
            sql += sum( [x_sume[ii],y_sume[ii]] ).__str__() + ")"
            cur.execute(sql)
            con.commit()


def correlate_summits_for_reps_in_group(rgroupid, con):
    """This is a long method (sorry). 
    It correlates the summits from all the replicates in a group,
    using a variety of methods.
    It assumes that all the replicates come from the same species."""
        
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
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
    x = get_geneids_with_summits(con, repid1)
    for ii in x:
        venn_data[ rep1name ].append(ii[0])
    y = get_geneids_with_summits(con, repid2)
    for ii in y:
        venn_data[ rep2name ].append(ii[0])
    plot_venn_diagram( venn_data, "summits." + rep1name.__str__() + "." + rep2name.__str__())
    
 
    """For the following plots, we'll use these data structures."""   
    print "\n. Building a dictionary of summits & genes for", rep1name, "and", rep2name
    rep1_gene_summitscores = {} # key = geneid, value = summit scores from rep1
    rep2_gene_summitscores = {}
    count = 0
        
    genes = get_genes_for_species(con, species1)
    for g in genes:
        count += 1
        if count%100 == 0:
            sys.stdout.write(".")
            sys.stdout.flush() 
        
        gid = g[0]
        x = get_summit_scores_for_gene(gid, repid1, con)
        if x.__len__() > 0:
            rep1_gene_summitscores[gid] = x
        y = get_summit_scores_for_gene(gid, repid2, con)
        if y.__len__() > 0:
            rep2_gene_summitscores[gid] = y
    
    
    """Write an Excel table with genes and their scores in each replicate.
    During the iteration over genes, also update the SQL tables with summary
    statistics about summits."""
    xlpath = "genes." + rep1name.__str__() + "." + rep2name.__str__() + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    fout.write("Nsites(" + rep1name.__str__() + ")\t")
    fout.write("meanE(" + rep1name.__str__() + ")\t")
    #fout.write("mode(" + rep1name.__str__() + ")\t")
    fout.write("maxE(" + rep1name.__str__() + ")\t")
    fout.write("minE(" + rep1name.__str__() + ")\t")
    fout.write("Nsites(" + rep2name.__str__() + ")\t")
    fout.write("meanE(" + rep2name.__str__() + ")\t")
    #fout.write("mode(" + rep2name.__str__() + ")\t")
    fout.write("maxE(" + rep2name.__str__() + ")\t")
    fout.write("minE(" + rep2name.__str__() + ")\n")

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
            
                sql = "INSERT INTO RepgroupSummitStats(repgroupid, geneid, maxsummit, nsummits) "
                sql += " VALUES("
                sql += rgroupid.__str__() + ","
                sql += gid.__str__() + ","
                sql += max( rep1_gene_summitscores[gid] + rep2_gene_summitscores[gid] ).__str__() + ","
                sql += ( rep1_gene_summitscores[gid] + rep2_gene_summitscores[gid] ).__len__().__str__()
                sql += ")"
                cur.execute(sql)
                con.commit()
        
        fout.write(g[0].__str__() + "\t" + g[1].__str__() + "\t")
        
        repdata = [rep1_gene_summitscores, rep2_gene_summitscores]
        for data in repdata:
            if gid in data:
                thisn = data[gid].__len__()
                if thisn == 1:
                    thismode = data[gid][0]
                    thismean = data[gid][0]
                else:
                    thismean = scipystats.tmean( data[gid] ) # the [0] index is because SciPy returns an array of results
                    thismode = scipystats.mode( data[gid] )[0]
                thismax = max( data[gid] )
                thismin = min( data[gid] )
                #print "81:", thisn, thismean, thismode, thismax, thismin
                fout.write("%d"%thisn + "\t")
                fout.write("%.3f"%thismean + "\t")
                #fout.write("%.3f"%thismode + "\t")
                fout.write("%.3f"%thismax + "\t")
                fout.write("%.3f"%thismin + "\t")
            else:
                fout.write("0\t0\t0\t0\t")               
    fout.close()

    """ Plot score correlations."""
    xvals = []
    yvals = []
    for g in genes:
        gid = g[0]
        if gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
            xvals.append( max(rep1_gene_summitscores[gid]) )
            yvals.append( max(rep2_gene_summitscores[gid]) )
    scatter1(xvals,yvals,"summits.score." + rep1name.__str__() + "." + rep2name.__str__(), xlab=rep1name+": max summit score for gene", ylab=rep2name+": max summit score for gene")

    """Plot rank correlations."""
    scoresx_genes = {} # bins of genes, organized by their score range
    scoresy_genes = {}
    chromids = get_chrom_ids(con, species1)
    for chr in chromids:
        genes = get_genes_for_chrom(con, chr[0])
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
    scatter1(xvals,yvals,"summits.rank." + rep1name.__str__() + "." + rep2name.__str__(), xlab=rep1name.__str__()+": rank order of genes by max summit score", ylab=rep2name.__str__()+": rank order of genes by max summit score")
         
    return con
         
def correlate_summits_for_union(unionid, con):
    """repgroups is a list of repgroupids. 
    Each repgroup contains multiple replicates."""

    cur = con.cursor()
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )

    seen_genes = []

    """Venn diagram of genes with/without summits in both replicates."""
    venn_data = {}
    rgroupids
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )
        genes = get_geneids_from_repgroup(con, rgroupid)
        seen_genes += genes
        venn_data[ rgroupname ] = genes
    plot_venn_diagram( venn_data, "summits." + unionname.__str__() )

    """Are all the replicates in this union from the same species?"""
    possible_species = []
    for rgroupid in rgroupids:
        cur.execute("SELECT species from Replicates where id in (SELECT replicate from GroupReplicate where rgroup=" + rgroupid.__str__() + ")")
        x = cur.fetchall()
        for ii in x:
            if ii[0] not in possible_species:
                possible_species.append( ii[0] )
    if possible_species.__len__() > 1:
        """Then the replicates in the union come from different species."""
        print "The replicates in the union", unionname, "come from multiple species."
        exit()    
        
    """
    1. Update the table UnionGenes with genes that have a peak in all repgroups in this union.
    2. Write the Excel table listing the genes and their scores in the repgroups."""
    xlpath = "genes." + unionname + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )        
        fout.write("max(" + rgroupname.__str__() + ")\t")
        fout.write("N(" + rgroupname.__str__() + ")\t")
    
    count = 0
    for geneid in seen_genes:
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        for rgroupid in rgroupids:
            count += 1
            if count%50 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            
            sql = "SELECT * from RepgroupSummitStats where repgroupid=" + rgroupid.__str__()
            sql += " and geneid=" + geneid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            if x != None:
                fout.write(x[2].__str__() + "\t")
                fout.write(x[3].__str__() + "\t")
        fout.write("\n")
    fout.close()
            
    count = 0
    genes = get_genes_for_species(con, possible_species[0])
    for g in genes:
        gid = g[0]
        
        count += 1
        if count%100 == 0:
            sys.stdout.write(".")
            sys.stdout.flush()
        
        # 1.
        # Does this gene have a peak in all replicategroups in this species?
        found_in_all = True
        for rgroupid in rgroupids:
            rgroupname = get_repgroup_name( rgroupid, con )
            if gid not in venn_data[rgroupname]:
                found_in_all = False
        if found_in_all:
            """Insert this union-gene pair into the DB, but only if it doesn't already exist."""
            cur.execute("SELECT COUNT(*) from UnionGenes where unionid=" + unionid.__str__() + " and geneid=" + gid.__str__() )
            if cur.fetchone()[0] == 0: # if this entry doesn't already exist...
                sql = "INSERT into UnionGenes (unionid, geneid) VALUES(" + unionid.__str__() + "," + gid.__str__() + ")"
                cur.execute(sql)
                con.commit()