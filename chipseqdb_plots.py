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

def get_xy_from_rscript(scriptpath):
    """Returns [x,y], where x and y are the data arrays plotted in the script."""
    x = None
    y = None
    fin = open(scriptpath, "r")
    for l in fin.xreadlines():
        if l.startswith("x<-c("):
            x = l
        elif l.startswith("y<-c("):
            y = l
    return [x,y]

def compute_summits_for_union(unionid, con):
    cur = con.cursor()
    sql = "DELETE from UnionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
    sql = "DELETE from UnionSummitStats where unionid=" + unionid.__str__()
    cur.execute(sql)
    con.commit()
    
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )

    print "\n. Comparing summits for union", unionname

    seen_genes = []

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
    count = 0
    for geneid in seen_genes:        
        nsummits = 0
        maxsummit = 0.0
        
        found_in_all = True # did we find this gene in all the repgroups?
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
                if x[2] > maxsummit:
                    maxsummit = x[2]
                nsummits += x[3]
            if x == None:
                found_in_all = False
                            
        if found_in_all:
            """If this gene is in all the repgroups for the union, then add values to the tables
            UnionGenes and UnionSummitStats."""
            sql = "INSERT into UnionGenes (unionid, geneid) VALUES(" + unionid.__str__() + "," + geneid.__str__() + ")"
            cur.execute(sql)
            #con.commit()
            
            sql = "INSERT into UnionSummitStats (unionid, geneid, maxsummit, nsummits)"
            sql += " VALUES(" + unionid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += maxsummit.__str__() + ","
            sql += (nsummits / float(rgroupids.__len__()) ).__str__() + ")"
            cur.execute(sql)
    con.commit()

def plot_summits_for_union(unionid, con):
    cur = con.cursor()
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )

    print "\n. Plotting summits for union", unionname

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
    xlpath = "summits." + unionname + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )        
        fout.write("max(" + rgroupname.__str__() + ")\t")
        fout.write("N(" + rgroupname.__str__() + ")\t")
    fout.write("\n")
    
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
            else:
                fout.write("0\t0\t")
        fout.write("\n")
    fout.close()
    
    add_unionfile(xlpath, unionid, "Excel table with summit stats for union " + unionname, con)
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with summit stats for union " + unionname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()

def compute_summits_for_speciesunion(uid, con):
    """A species-union is a meta union of unions.
    Each union contains multiple replicates from the same species.
    The Species-union compares these unions, across species boundaries."""
    
    cur = con.cursor()
    sql = "DELETE from SpeciesunionGenes where unionid=" + uid.__str__()
    cur.execute(sql)
    sql = "DELETE from SpeciesunionSummitStats where spunionid=" + uid.__str__()
    cur.execute(sql)
    con.commit()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    nuids = unionids.__len__()
    print "\n. Comparing summits for genes across species", spunionname

    seen_genes = []    
    geneid_summitstats = {} # key = translated gene ID, value = list of SQL results from UnionSummitStats

    count = 0
    unionid_genes = {}
    for unionid in unionids:
        unionname = get_unionname( unionid, con )
        
        sql = "SELECT * from UnionSummitStats where unionid=" + unionid.__str__()
        cur.execute(sql)
        results = cur.fetchall()
        for ii in results:
            geneid = ii[1]
            translated_id = get_geneid_from_aliasid( geneid, con )
            
            if geneid not in geneid_summistats:
                geneid_summitstats[geneid].append( ii )

    for gid in geneid_summitstats:    
        if geneid_summitstats[gid].__len__() == nuids:
            """Only deal with geneids that exist in all the unions of this speciesunion."""
            
            maxsummit = 0
            nsummits = 0
            for ii in geneid_summitstats[gid]:
                if ii[2] > maxsummit:
                    maxsummit = ii[2]
                nsummits += ii[3]
            
            sql = "INSERT INTO SpeciesunionGenes (unionid, geneid) VALUES(" + uid.__str__() + "," + gid.__str__() + ")"
            cur.execute(sql)
            #con.commit()

            sql = "INSERT INTO SpeciesunionSummitStats (spunionid, geneid, maxsummit, nsummits)"
            sql += " VALUES(" + uid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += maxsummit.__str__() + ","
            sql += nsummits.__str__() + ")"
            cur.execute(sql)
    con.commit()

def plot_summits_for_speciesunion(uid, con):
    """A species-union is a meta union of unions.
    Each union contains multiple replicates from the same species.
    The Species-union compares these unions, across species boundaries."""
    
    cur = con.cursor()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    print "\n. Plotting summits for genes across species", spunionname

    """seen_genes will be a list of un-aliased gene IDs for all genes that have
    a summit in their nearby regulatory regions, in all the replicates represented
    by at least one Union. Note that the gene IDs here are un-aliased; in other words,
    they represent the gene ID found from the pillars file. You may need to look up
    alias IDs for this gene's homolog to locate relevant summit scores in different species."""
    seen_genes = []

    """key = gene id, value = list of alias IDs."""
    gene_aliases = {}

    """This loop builds data for a Venn diagram of genes with/without summits in both replicates."""
    count = 0
    venn_data = {}
    for unionid in unionids:
        unionname = get_unionname( unionid, con )
        
        """genes is a list of all genes that have summits in nearby regulatory regions in all
        the replicates represented by the union."""
        genes = get_geneids_from_union(con, unionid)
        translated_genes = []
        for geneid in genes:
            
            """Print a period every 200 iterations, to indicate that this program is still alive."""
            count += 1
            if count%200 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            
            translated_id = get_geneid_from_aliasid( geneid, con )
            if translated_id != None:
                """Remember this alias, for later use in writing the Excel table."""
                if translated_id not in gene_aliases:
                    gene_aliases[translated_id] = []
                gene_aliases[translated_id].append( geneid )

                """translated_genes is the list of gene IDs that goes in the Venn diagram,
                so that homologous gene IDs can be compared between species."""
                translated_genes.append( translated_id ) 
                if translated_id not in seen_genes:
                    seen_genes.append( translated_id )
        venn_data[ unionname ] = translated_genes
    plot_venn_diagram( venn_data, "summits." + spunionname.__str__() )
    
    print "\n. The species-union", spunionname, "has", seen_genes.__len__(), "genes with summits in all replicates."
        
    #
    # Write an Excel table using data from UnionSummitStats
    #
    xlpath = "summits." + spunionname + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    for unionid in unionids:
        unionname = get_unionname( unionid, con )  
        fout.write("geneID(" + unionname.__str__() + ")\t")
        fout.write("name(" + unionname.__str__() + ")\t")      
        fout.write("max(" + unionname.__str__() + ")\t")
        fout.write("N(" + unionname.__str__() + ")\t")
    fout.write("\n")
    count = 0
    for gid in seen_genes: 
        fout.write(gid.__str__() + "\t" + get_genename(gid, con) + "\t" )
    
        for unionid in unionids:
            count += 1
            if count%100 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            
            foundit = False
            for geneid in gene_aliases[gid]:    

                if foundit == False:                    
                    sql = "SELECT * from UnionSummitStats where unionid=" + unionid.__str__()
                    sql += " and geneid=" + geneid.__str__()
                    cur.execute(sql)
                    x = cur.fetchone()
                    if x != None:
                        foundit = True
                        fout.write(geneid.__str__() + "\t")
                        fout.write(get_genename(geneid, con) + "\t")
                        fout.write(x[2].__str__() + "\t")
                        fout.write(x[3].__str__() + "\t")    
            if foundit == False:
                fout.write("---\t---\t0\t0\t")
        fout.write("\n")
    fout.close()
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with summit stats for speciesunion " + spunionname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + uid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()
    
    add_speciesunionfile( xlpath, uid, "Excel table with summit stats for speciesunion " + spunionname, con )
    
    #
    # to-do: scatter plots for SpeciesunionUnion data
    #

def compute_summits_for_reps_in_group(rgroupid, con):
    """This is a long method (sorry). 
    It correlates the summits from all the replicates in a group,
    using a variety of methods.
    It assumes that all the replicates come from the same species."""
     
    cur = con.cursor()
    sql = "DELETE from RepgroupGenes where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    sql = "DELETE from RepgroupSummitStats where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    con.commit()
         
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
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
    
    """For the following plots, we'll use these data structures."""   
    print "\n. Building a dictionary of summits & genes for", rep1name, "and", rep2name
    rep1_gene_summitscores = {} # key = geneid, value = summit scores from rep1
    rep2_gene_summitscores = {}
    count = 0
    
#     gene_result = {}
#     sql = "SELECT * from Summits where repid=" + repid1.__str__()
#     cur.execute(sql)
#     results = cur.fetchall()
#     for ii in results:
#         geneid = ii[0]
    
    genes = get_genes_for_species(con, species1)
    for g in genes:
        count += 1
        if count%50 == 0:
            sys.stdout.write(".")
            sys.stdout.flush() 
        
        gid = g[0]
        x = get_summit_scores_for_gene(gid, repid1, con)
        if x.__len__() > 0:
            rep1_gene_summitscores[gid] = x
        y = get_summit_scores_for_gene(gid, repid2, con)
        if y.__len__() > 0:
            rep2_gene_summitscores[gid] = y
    
    for g in genes:
        gid = g[0]
        """Ignore genes that aren't in at least on replicate."""
        if gid not in rep1_gene_summitscores and gid not in rep2_gene_summitscores:
            continue # skip to the next gene
        
        elif gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
            # This summit is in both replicates:
            sql = "INSERT into RepgroupGenes (repgroupid, geneid) VALUES(" + rgroupid.__str__() + "," + gid.__str__() + ")"
            cur.execute(sql)
            #con.commit()
        
            sql = "INSERT INTO RepgroupSummitStats(repgroupid, geneid, maxsummit, nsummits) "
            sql += " VALUES("
            sql += rgroupid.__str__() + ","
            sql += gid.__str__() + ","
            sql += max( rep1_gene_summitscores[gid] + rep2_gene_summitscores[gid] ).__str__() + ","
            sql += ( rep1_gene_summitscores[gid] + rep2_gene_summitscores[gid] ).__len__().__str__()
            sql += ")"
            cur.execute(sql)
    con.commit()                      

def plot_summits_for_reps_in_group(rgroupid, con):
    """This is a long method (sorry). 
    It correlates the summits from all the replicates in a group,
    using a variety of methods.
    It assumes that all the replicates come from the same species."""
    
    repgroupname = get_repgroup_name(rgroupid, con)
     
    cur = con.cursor()
    
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]
    
    print "\n. Plotting replicates", rep1name, "and", rep2name
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
    plot_venn_diagram( venn_data, "summits." + repgroupname)
    
 
    """For the following plots, we'll use these data structures."""   
    print "\n. Building a dictionary of summits & genes for", rep1name, "and", rep2name
    rep1_gene_summitscores = {} # key = geneid, value = summit scores from rep1
    rep2_gene_summitscores = {}
    count = 0
        
    genes = get_genes_for_species(con, species1)
    for g in genes:
        count += 1
        if count%50 == 0:
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
    xlpath = "summits." + repgroupname + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    fout.write("Nsites(" + rep1name.__str__() + ")\t")
    fout.write("mean(" + rep1name.__str__() + ")\t")
    fout.write("max(" + rep1name.__str__() + ")\t")
    fout.write("min(" + rep1name.__str__() + ")\t")
    fout.write("Nsites(" + rep2name.__str__() + ")\t")
    fout.write("mean(" + rep2name.__str__() + ")\t")
    fout.write("max(" + rep2name.__str__() + ")\t")
    fout.write("min(" + rep2name.__str__() + ")\n")

    for g in genes:
        gid = g[0]
        """Ignore genes that aren't in at least on replicate."""
        if gid not in rep1_gene_summitscores and gid not in rep2_gene_summitscores:
            continue # skip to the next gene
                
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
                fout.write("%d"%thisn + "\t")
                fout.write("%.3f"%thismean + "\t")
                fout.write("%.3f"%thismax + "\t")
                fout.write("%.3f"%thismin + "\t")
            else:
                fout.write("0\t0\t0\t0\t")               
    fout.close()
    
    add_repgroupfile(xlpath, rgroupid, "Excel table with summit stats for replicate group " + repgroupname, con)
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with summit stats for replicate group " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()

    """ Plot score correlations."""
    xvals = []
    yvals = []
    for g in genes:
        gid = g[0]
        if gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
            xvals.append( max(rep1_gene_summitscores[gid]) )
            yvals.append( max(rep2_gene_summitscores[gid]) )
    cranpath = scatter1(xvals,yvals,"summits.score." + repgroupname, xlab=rep1name+": max summit score for gene", ylab=rep2name+": max summit score for gene")
    
    
    add_repgroupfile(cranpath, rgroupid, "R script path to a scatterplot comparing gene-summit scores for replicate group " + repgroupname, con)
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing gene-summit scores for replicate group " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()


    """Plot rank correlations."""
    scoresx_genes = {} # bins of genes, organized by their score range
    scoresy_genes = {}
    chromids = get_chrom_ids(con, species1)
    for chrid in chromids:
        genes = get_genes_for_chrom(con, chrid)
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
    cranpath = scatter1(xvals,yvals,"summits.rank." + repgroupname, xlab=rep1name.__str__()+": rank order of genes by max summit score", ylab=rep2name.__str__()+": rank order of genes by max summit score")
    
    add_repgroupfile(cranpath, rgroupid, "R script path to a scatterplot comparing the rank of genes for replicate group " + repgroupname, con )
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the rank of genes for replicate group " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()

def compute_enrichments_for_union(unionid, con, keyword=None):
    """This method puts data into UnionEnrichmentStats"""
    
    cur = con.cursor()
    sql = "DELETE from UnionEnrichmentStats where unionid=" + unionid.__str__()
    cur.execute( sql )
    con.commit()
    
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    rgroupids.sort()
    nrgids = rgroupids.__len__()
    
    if nrgids < 1:
        print "\n. WARNING: The union", unionname, "contains only", nrgids,"replicate groups."
        print ". I'm skipping this union."
        return
    
    print "\n. Computing enrichments for union", unionid, unionname
    
    #
    # A better way?
    #
    sql = "SELECT * from GroupEnrichmentStats where"
    sqlbits = []
    for rgid in rgroupids:
        sqlbits.append( " rgroupid=" + rgid.__str__() )
    sql += " OR ".join( sqlbits )
    print "1084:", sql
    cur.execute(sql)
    results = cur.fetchall()
    
    gene_results = {}
    for ii in results:
        geneid = ii[1]
        if geneid not in results:
            gene_results[geneid] = []
        gene_results[geneid].append( ii )
    
    count = 0
    for geneid in gene_results.keys():        
        """Compute union-wide stats, but only if we have data from all the replicates in this union."""
        if gene_results[geneid].__len__() == nrgids:
            max_max = 0
            all_means = []
            all_sums = []
            
            for ii in gene_results[geneid]:
                rgid = ii[0]
                
                count += 1
                if count%20000 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
    
                if ii[2] > max_max:
                    max_max = ii[2]
                all_means.append( ii[3] )
                all_sums.append( ii[4] )
    
            
            sql = "INSERT INTO UnionEnrichmentStats (unionid, geneid, maxenrich, meanenrich, sumenrich)"
            sql += "VALUES (" + unionid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += max_max.__str__() + ","
            sql += mean(all_means).__str__() + ","
            sql += sum(all_sums).__str__() + ")"
            cur.execute(sql)
    con.commit() 
    

def plot_enrichment_4x4(unionid, con, keyword=None):
    cur = con.cursor()

    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    rgroupids.sort()
    repids = []
        
    for rgid in rgroupids:
        repname = get_repgroup_name(rgid, con)
        x = get_repids_in_group(rgid, con)
        for ii in x:
            repids.append( ii )
    repid_maxs = []
    repid_means = []
    repid_sums = []
    repid_name = []    
    speciesname = ""
    for repid in repids:
        sql = "SELECT name from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        repname = cur.fetchone()
        if repname == None:
            repname = repid.__str__()
        else:
            repname = repname[0].__str__()
        
        speciesid = get_speciesid_for_repid(repid, con)
        sql = "SELECT name from Species where id=" + speciesid.__str__()
        cur.execute(sql)
        speciesname = cur.fetchone()[0]
        
        if repname.__contains__(speciesname):
            repname = re.sub(speciesname, "", repname)
        
        (x_maxe, x_meane, x_sume, geneids) = get_plot_array_for_replicate(repid, speciesid, con)
        repid_means.append( x_meane )
        repid_maxs.append( x_maxe)
        repid_sums.append( x_sume)
        repid_name.append( repname )
    
    scatter_names = repid_name + repid_name + repid_name
    scatter12x4(repid_maxs + repid_means + repid_sums, scatter_names, "enrich.12x4." + unionname, title="Fold Enrichment " + speciesname, xlab="fold-enrichment", ylab="fold-enrichment")
    

def plot_enrichments_for_union(unionid, con, keyword=None):
    """This method makes the plots and excel tables relevant to the enrichment scores
    for all replicates in the union.
    This method assumes that data has been inserted into the table UnionEnrichmentStats"""

    cur = con.cursor()

    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    rgroupids.sort()
    
    print "\n. Plotting enrichments for union", unionname

    plot_enrichment_4x4(unionid, con)

    
    #
    # What gene IDs apply to this union?
    #
    sql = "SELECT * from GroupEnrichmentStats where"
    sqlbits = []
    for rgid in rgroupids:
        sqlbits.append( " rgroupid=" + rgid.__str__() )
    sql += " OR ".join( sqlbits )
    cur.execute(sql)
    results = cur.fetchall()
    geneid_results = {}
    for ii in results:
        geneid = ii[1]
        genename = get_genename(geneid, con)
        geneid, genename
        if genename == None:
            continue
        if geneid not in geneid_results:
            geneid_results[geneid] = []
        geneid_results[geneid].append( ii )
    
        
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
    fout.write("\n")    
    count = 0
    print "744:", geneid_results.keys().__len__()
    for geneid in geneid_results.keys():
        
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        rgid_foundvals = {}
        for rgid in rgroupids:
            
            if rgid not in rgid_foundvals:
                rgid_foundvals[rgid] = []
            
            count += 1
            if count%50 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            foundit = False
            
            for x in geneid_results[geneid]:
                if x[0] == rgid:
                    foundit = True
                    fout.write(x[2].__str__() + "\t")
                    fout.write(x[3].__str__() + "\t")
                    fout.write(x[4].__str__() + "\t")
    
                    rgid_foundvals[rgid].append(x[2])
    
                    rgid_maxvals[rgid].append(x[2])
                    rgid_meanvals[rgid].append(x[3])
                    rgid_sumvals[rgid].append(x[4])
 
               
            if foundit == False:
                fout.write("0\t0\t0\t")
                rgid_maxvals[rgid].append(0)
                rgid_meanvals[rgid].append(0)
                rgid_sumvals[rgid].append(0)

        #for rgid in rgroupids:
        #    x = rgid_foundvals[rgid].__len__()
        #    if x != 1:
        #        print "783:", x, rgid#, rgid_foundvals[rgid], rgid_maxvals[rgid]

        fout.write("\n")
    fout.close()
    
    add_unionfile(xlpath, unionid, "Excel table with enrichment stats for union " + unionname, con)
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with enrichment stats for union " + unionname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()
    
    #
    # Scatterplots
    #
#     for ii in range(0, rgroupids.__len__() ):
#         for jj in range(ii+1, rgroupids.__len__() ):
#             this_rgid = rgroupids[ii]
#             that_rgid = rgroupids[jj]
#             
#             if rgid_maxvals[this_rgid].__len__() == 0:
#                 continue
#             
#             this_rgname = get_repgroup_name(this_rgid, con)
#             that_rgname = get_repgroup_name(that_rgid, con)
#             
#             cranpath = scatter1(rgid_maxvals[this_rgid], rgid_maxvals[that_rgid], "enrich.max." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": max enrichment", ylab=that_rgname.__str__()+": max enrichment", force_square=True)
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the maximum fold-enrichment scores for replicates " + this_rgname + " to " + that_rgname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()
#         
#             cranpath = scatter1(rgid_meanvals[this_rgid], rgid_meanvals[that_rgid], "enrich.mean." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": mean enrichment", ylab=that_rgname.__str__()+": mean enrichment", force_square=True)
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the average fold-enrichment scores for replicates " + this_rgname + " to " + that_rgname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()
#             
#             cranpath = scatter1(rgid_sumvals[this_rgid], rgid_sumvals[that_rgid], "enrich.sum." + this_rgname.__str__() + "." + that_rgname.__str__(), xlab=this_rgname.__str__()+": sum enrichment", ylab=that_rgname.__str__()+": sum enrichment", force_square=True)
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the sum of fold-enrichment scores for replicates " + this_rgname + " to " + that_rgname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()
            
    values = {}
    values["max enrichment"] =  [ rgid_maxvals[rgroupids[0]], rgid_maxvals[rgroupids[1]]     ]
    values["mean enrichment"] = [ rgid_meanvals[rgroupids[0]], rgid_meanvals[rgroupids[1]]   ]
    values["sum enrichment"]=   [ rgid_sumvals[rgroupids[0]], rgid_sumvals[rgroupids[1]]     ]
    cranpath = scatter1xn(values, "enrich.1x3." + unionname, xlab=this_rgname, ylab=that_rgname, force_square=True )
    
    add_unionfile(cranpath, unionid,"R script path to a 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates " + this_rgname + " to " + that_rgname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()     

def compute_enrichments_for_reps_in_group(rgroupid, con):
    """Fills the table GroupEnrichmentStats"""
    
    cur = con.cursor()
    sql = "DELETE from GroupEnrichmentStats where rgroupid=" + rgroupid.__str__()
    cur.execute( sql )
    con.commit()
    
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
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
        if geneid in y.keys():
            geneids.append( geneid )
        
    for geneid in geneids:
        x_maxe.append( x[geneid][0] )
        y_maxe.append( y[geneid][0] )
        x_meane.append( x[geneid][1] )
        y_meane.append( y[geneid][1] )
        x_sume.append( x[geneid][2] )
        y_sume.append( y[geneid][2] )
    
    print "\n. Updating the table GroupEnrichmentStats"
    for ii in range(0, geneids.__len__() ):
        geneid = geneids[ii]
        sql = "INSERT into GroupEnrichmentStats (rgroupid, geneid, maxenrich, meanenrich, sumenrich)"
        sql += " VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ","
        sql += max( [x_maxe[ii], y_maxe[ii]] ).__str__() + ","
        sql += mean( [x_meane[ii], y_meane[ii]] ).__str__() + ","
        sql += sum( [x_sume[ii],y_sume[ii]] ).__str__() + ")"
        cur.execute(sql)
    con.commit()

def get_plot_array_for_replicate(repid, species, con):
    """Returns ( x_maxe[], x_meane[], x_sume[] ) for this replicate."""
    x_maxe = []
    x_meane = []
    x_sume = []
    x = get_enrichment_stats_for_replicate( repid, con )
    
    genes = get_genes_for_species(con, species)
    geneids = []
    for ii in genes:
        geneids.append( ii[0] )
    geneids.sort()
                
    for geneid in geneids:
        if geneid in x:
            x_maxe.append( x[geneid][0] )
            x_meane.append( x[geneid][1] )
            x_sume.append( x[geneid][2] )
        else:
            x_maxe.append(0)
            x_meane.append(0)
            x_sume.append(0)      
    return (x_maxe, x_meane, x_sume, geneids)

def plot_enrichments_for_reps_in_group(rgroupid, con, repgroupname=None, repids=None):
    cur = con.cursor()
        
    if repgroupname == None:
        repgroupname = get_repgroup_name(rgroupid, con)
    
    if repids == None:
        repids = get_repids_in_group(rgroupid, con)
        
    repid1 = repids[0]
    repid2 = repids[1]
    
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]
    
    speciesid1 = get_speciesid_for_repid(repid1, con)
    (x_maxe, x_meane, x_sume, geneids) = get_plot_array_for_replicate(repid1, speciesid1, con)
    speciesid2 = get_speciesid_for_repid(repid2, con)
    (y_maxe, y_meane, y_sume, geneids) = get_plot_array_for_replicate(repid2, speciesid2, con)

    
    #
    # Plot scatterplots
    #
#     cranpath = scatter1(x_maxe,y_maxe,"enrich.max." + repgroupname, xlab=rep1name.__str__()+": max enrichment", ylab=rep2name.__str__()+": max enrichment", force_square=True)
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'Scatterplot with maximum enrichments for replicate " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()
#     
#     
#     cranpath = scatter1(x_meane,y_meane,"enrich.mean." + repgroupname, xlab=rep1name.__str__()+": mean enrichment", ylab=rep2name.__str__()+": mean enrichment", force_square=True)
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'Scatterplot with mean enrichments for replicate " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()
#     
#     scatter1(x_sume,y_sume,"enrich.sum." + repgroupname, xlab=rep1name.__str__()+": sum enrichment", ylab=rep2name.__str__()+": sum enrichment", force_square=True)
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'Scatterplot with sum enrichments for replicate " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()

    values = {}
    values["max enrichment"] =  [x_maxe, y_maxe]
    values["mean enrichment"] = [x_meane, y_meane ]
    values["sum enrichment"]=   [x_sume, y_sume ]
    
    cranpath = scatter1xn(values, "enrich." + repgroupname, xlab=rep1name, ylab=rep2name, force_square=True )
    add_repgroupfile(cranpath, rgroupid, "R script for multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
    
#     sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script for multi-panel scatterplot with enrichment values for replicate group " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + cranpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit() 

    #
    # Write an Excel Table
    #
    xlpath = "enrich." + repgroupname + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for rgname in [rep1name,rep2name]:
        fout.write("max(" + rgname + ")\t")
        fout.write("mean(" + rgname + ")\t")
        fout.write("sum(" + rgname + ")\t")        
    fout.write("\n")
    for ii in range(0, geneids.__len__() ):
        genename = get_genename(geneids[ii], con)
        if genename == None:
            continue
        fout.write(geneids[ii].__str__() + "\t" + genename + "\t" )
        fout.write(x_maxe[ii].__str__() + "\t")
        fout.write(x_meane[ii].__str__() + "\t")
        fout.write(x_sume[ii].__str__() + "\t")
        fout.write(y_maxe[ii].__str__() + "\t")
        fout.write(y_meane[ii].__str__() + "\t")
        fout.write(y_sume[ii].__str__())
        fout.write("\n")
    fout.close()

    add_repgroupfile(xlpath,rgroupid,"Excel table with enrichment stats for replicate group " + repgroupname, con)

#     sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with enrichment stats for replicate group " + repgroupname + "')"
#     cur.execute(sql)
#     con.commit()
#     sql = "SELECT fileid from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     fileid = cur.fetchone()[0]
#     sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + rgroupid.__str__() + ","
#     sql += fileid.__str__() + ")"
#     cur.execute(sql)
#     con.commit()


def compute_enrichments_for_speciesunion(uid, con):
    cur = con.cursor()
    sql = "DELETE from SpeciesunionEnrichmentStats where unionid=" + uid.__str__()
    con.commit()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    nuids = unionids.__len__()
    
    print "\n. Comparing enrichments acorss species", spunionname, unionids
    sys.stdout.flush()

    """
    What gene IDs apply to this species union?
    In the following block, we lookup the genes that have enrichment
    scores in their nearby regulatory regions in all the replicates
    represented by at least one union. We then translate
    those gene IDs, using the table GeneHomology with data from
    the pillars file, in order to generate a list of genes
    which can be compared across species."""
    
    
    sql = "SELECT * from UnionEnrichmentStats where"
    sqlbits = []
    for unionid in unionids:
        sqlbits.append( " unionid=" + unionid.__str__() )
    sql += " OR ".join( sqlbits )
    cur.execute(sql)
    results = cur.fetchall()
    
    """genes will be a list of un-aliased genes. That is, gene IDs translated via pillars."""
    geneid_results = {}

    """A dictionary of aliased genes."""
    gene_aliases = {}
    for ii in results:       
        translated_id = get_geneid_from_aliasid( ii[1], con )
        if translated_id != None:
            if translated_id not in gene_aliases:
                gene_aliases[translated_id] = []
            gene_aliases[ translated_id ].append( ii[1]) 
            if translated_id not in genes:
                genes.append( translated_id )
        if translated_id not in geneid_results:
            geneid_results[ translated_id ] = []
        geneid_results[translated_id].append(ii)
        
    
    count = 0
    for gid in geneid_results.keys():
        """Only keep those genes that have enrichment scores in all the replicates in this union."""
        if geneid_results[gid].__len__() == nuids:    
            max_max = 0
            all_means = []
            all_sums = []
                   
            for ii in geneid_results[gid]:
                unionid = ii[0]
               
                count += 1
                if count%50 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
             
                if ii[2] > max_max:
                    max_max = ii[2]
                all_means.append( ii[3] )
                all_sums.append( ii[4] )

            sql = "INSERT INTO SpeciesunionEnrichmentStats (unionid, geneid, maxenrich, meanenrich, sumenrich)"
            sql += "VALUES (" + unionid.__str__() + ","
            sql += gid.__str__() + ","
            sql += max_max.__str__() + ","
            sql += mean(all_means).__str__() + ","
            sql += sum(all_sums).__str__() + ")"
            cur.execute(sql)
    con.commit()    

def plot_enrichments_for_speciesunion(uid, con):
    cur = con.cursor()
        
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    
    print "\n. Plotting enrichments acorss species", spunionname, unionids
    sys.stdout.flush()

    """
    What gene IDs apply to this species union?
    In the following block, we lookup the genes that have enrichment
    scores in their nearby regulatory regions in all the replicates
    represented by at least one union. We then translate
    those gene IDs, using the table GeneHomology with data from
    the pillars file, in order to generate a list of genes
    which can be compared across species."""
    
    sql = "SELECT geneid from UnionEnrichmentStats where"
    sqlbits = []
    for unionid in unionids:
        sqlbits.append( " unionid=" + unionid.__str__() )
    sql += " OR ".join( sqlbits )
    cur.execute(sql)
    results = cur.fetchall()
    
    geneid_results = {}

    """A dictionary of aliased genes."""
    gene_aliases = {}
    for ii in results:       
        translated_id = get_geneid_from_aliasid( ii[1], con )
        if translated_id != None:
            if translated_id not in gene_aliases:
                gene_aliases[translated_id] = []
            gene_aliases[ translated_id ].append( ii[1]) 
            if translated_id not in genes:
                genes.append( translated_id )
        if translated_id not in geneid_results:
            geneid_results[ translated_id ] = []
        geneid_results[translated_id].append(ii)
        
    #
    # Write an excel table,
    # while also filling data structures for the scatterplot
    #
    unionid_maxvals = {}
    unionid_meanvals = {}
    unionid_sumvals = {}
    
    xlpath = "enrich." + spunionname.__str__() + ".xls"
    print "\n. Writing an Excel table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for unionid in unionids:
        unionname = get_unionname(unionid, con)
        fout.write("GeneID(" + unionname + ")\t")
        fout.write("GeneName(" + unionname + ")\t")
        fout.write("max(" + unionname + ")\t")
        fout.write("mean(" + unionname + ")\t")
        fout.write("sum(" + unionname + ")\t")

        unionid_maxvals[unionid] = []
        unionid_meanvals[unionid] = []
        unionid_sumvals[unionid] = []
    
    count = 0
    fout.write("\n")
    for gid in geneid_results.keys():
        max_max = 0
        all_means = []
        all_sums = []
        fout.write(gid.__str__() + "\t" + get_genename(gid, con) + "\t" )        
        
        for unionid in unionids: # find the unionid that mathches result x
            foundit = False
            
            for x in geneid_results[gid]:    
                count += 1
                if count%50 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
        
                if x[0] == unionid:
                    foundit = True
                    fout.write(geneid.__str__() + "\t")
                    fout.write(get_genename(geneid, con) + "\t")
                    fout.write(x[2].__str__() + "\t")
                    fout.write(x[3].__str__() + "\t")
                    fout.write(x[4].__str__() + "\t")
    
                    unionid_maxvals[unionid].append(x[2])
                    unionid_meanvals[unionid].append(x[3])
                    unionid_sumvals[unionid].append(x[4])
                    
                    if x[2] > max_max:
                        max_max = x[2]
                    all_means.append( x[3] )
                    all_sums.append( x[4] )
                if foundit == True:
                    continue
                    
            #else:
            if foundit == False:
                fout.write("---\t---\t0\t0\t0\t")
                unionid_maxvals[unionid].append(0)
                unionid_meanvals[unionid].append(0)
                unionid_sumvals[unionid].append(0)
                #all_means.append( 0 )
        fout.write("\n")
    fout.close()

    add_speciesunionfile(xlpath, uid, "Excel table with enrichment stats for species union " + spunionname, con )

#     sql = "SELECT COUNT(*) from Files where path='" + xlpath + "'"
#     cur.execute(sql)
#     if cur.fetchone()[0] == 0:
#         sql = "INSERT INTO Files (path, note) VALUES('" + xlpath + "', 'Excel table with enrichment stats for species union " + spunionname + "')"
#         cur.execute(sql)
#         con.commit()
#         sql = "SELECT fileid from Files where path='" + xlpath + "'"
#         cur.execute(sql)
#         fileid = cur.fetchone()[0]
#         sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + uid.__str__() + ","
#         sql += fileid.__str__() + ")"
#         cur.execute(sql)
#         con.commit()
   
    #
    # Scatterplots
    #
    for ii in range(0, unionids.__len__() ):
        for jj in range(ii+1, unionids.__len__() ):
            #print "90:", ii, jj

            this_uid = unionids[ii]
            that_uid = unionids[jj]
            
            if unionid_maxvals[this_uid].__len__() == 0:
                continue
            
            this_uname = get_unionname(this_uid, con)
            that_uname = get_unionname(that_uid, con)
            
            cranpath = scatter1(unionid_maxvals[this_uid], unionid_maxvals[that_uid], "enrich.max." + this_uname.__str__() + "." + that_uname.__str__(), xlab=this_uname.__str__()+": max enrichment", ylab=that_uname.__str__()+": max enrichment", force_square=True)
            add_speciesunionfile(cranpath, uid, "R script path to a scatterplot comparing the maximum fold-enrichment scores for species union " + spunionname, con )
            
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the maximum fold-enrichment scores for species union " + spunionname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + uid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()
            
            cranpath = scatter1(unionid_meanvals[this_uid], unionid_meanvals[that_uid], "enrich.mean." + this_uname.__str__() + "." + that_uname.__str__(), xlab=this_uname.__str__()+": mean enrichment", ylab=that_uname.__str__()+": mean enrichment", force_square=True)
            add_speciesunionfile(cranpath, uid, "R script path to a scatterplot comparing the mean fold-enrichment scores for species union " + spunionname, con )
            
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the maximum fold-enrichment scores for species union " + spunionname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + uid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()
#             
            cranpath = scatter1(unionid_sumvals[this_uid], unionid_sumvals[that_uid], "enrich.sum." + this_uname.__str__() + "." + that_uname.__str__(), xlab=this_uname.__str__()+": sum enrichment", ylab=that_uname.__str__()+": sum enrichment", force_square=True)
            add_speciesunionfile(cranpath, uid, "R script path to a scatterplot comparing the sum fold-enrichment scores for species union " + spunionname, con )
            
#             sql = "INSERT INTO Files (path, note) VALUES('" + cranpath + "', 'R script path to a scatterplot comparing the maximum fold-enrichment scores for species union " + spunionname + "')"
#             cur.execute(sql)
#             con.commit()
#             sql = "SELECT fileid from Files where path='" + cranpath + "'"
#             cur.execute(sql)
#             fileid = cur.fetchone()[0]
#             sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + uid.__str__() + ","
#             sql += fileid.__str__() + ")"
#             cur.execute(sql)
#             con.commit()

    