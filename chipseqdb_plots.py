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

import scipy

def compute_summits_for_union(unionid, con):
    cur = con.cursor()
    sql = "DELETE from UnionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
    sql = "DELETE from UnionSummitStats where unionid=" + unionid.__str__()
    cur.execute(sql)
    con.commit()
    
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    nrepgroups = rgroupids.__len__()

    print "\n. Comparing summits for union", unionname

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
        print "The replicates in the union", unionname, "come from multiple species. Use the SPECIES_UNION command to compare across species boundaries."
        exit()    
        
    
    seen_genes = []    
    geneid_summitstats = {} # key = translated gene ID, value = list of SQL results from RepgroupSummitStats

    count = 1
    unionid_genes = {}
    for repgroupid in rgroupids:
        repgroupname = get_repgroup_name( repgroupid, con )
        
        sql = "SELECT * from RepgroupSummitStats where repgroupid=" + repgroupid.__str__()
        cur.execute(sql)
        results = cur.fetchall()
        for ii in results:
            geneid = ii[1]            
            if geneid not in geneid_summitstats:
                geneid_summitstats[geneid] = []
            geneid_summitstats[geneid].append( ii )

    """Calculate rank order of genes by max summit score."""
    maxscores_geneids = {}
    for gid in geneid_summitstats:
        maxsummit = 0
        for ii in geneid_summitstats[gid]:
            if ii[2] > maxsummit:
                maxsummit = ii[2]  
        if maxsummit not in maxscores_geneids:
            maxscores_geneids[ maxsummit ] = []
        maxscores_geneids[maxsummit].append(gid)
    scores = maxscores_geneids.keys()
    scores.sort(reverse=True)      

    gid_score = {}
    rank = 0
    for score in scores:
        newrank = rank + maxscores_geneids[score].__len__()/2.0 + 0.5
        for g in maxscores_geneids[score]:
            gid_score[g] = newrank
        rank += maxscores_geneids[score].__len__()

    for gid in geneid_summitstats:    
        
        """Only deal with geneids that exist in all the repgroups of this union."""
        if geneid_summitstats[gid].__len__() == nrepgroups:
            
            maxsummit = 0
            nsummits = 0
            for ii in geneid_summitstats[gid]:
                if ii[2] > maxsummit:
                    maxsummit = ii[2]
                nsummits += ii[3]
            
            sql = "INSERT INTO UnionGenes (unionid, geneid) VALUES(" + unionid.__str__() + "," + gid.__str__() + ")"
            cur.execute(sql)
            #con.commit()

            sql = "INSERT INTO UnionSummitStats (unionid, geneid, maxsummit, nsummits, ranksummit)"
            sql += " VALUES(" + unionid.__str__() + ","
            sql += gid.__str__() + ","
            sql += maxsummit.__str__() + ","
            sql += nsummits.__str__() + ","
            sql += gid_score[gid].__str__() + ")"
            cur.execute(sql)
    con.commit()
 
    sql = "SELECT COUNT(*) from UnionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
    ngenes = cur.fetchone()[0]
    print "\n    --> The union", unionname, "has", ngenes, "genes with summits in all members."

def plot_summits_for_union(unionid, con):
    cur = con.cursor()
    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )

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

    print "\n. Plotting summits for union", unionname

    seen_genes = [] # a list of genes that have data for this union
    
    """Scatterplot."""
    plot_summits_union_helper(unionid, con)
    
    """Venn diagram of genes with/without summits in both replicates."""
    venn_data = {}
    rgroupids
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )
        genes = get_geneids_from_repgroup(con, rgroupid)
        seen_genes += genes
        venn_data[ rgroupname ] = genes
    vennpath = plot_venn_diagram( venn_data, "summits." + unionname.__str__() )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for union " + unionname, con)

   
        
    """
    1. Update the table UnionGenes with genes that have a peak in all repgroups in this union.
    2. Write the Excel table listing the genes and their scores in the repgroups."""
    xlpath = "summits." + unionname + ".xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )        
        fout.write("max(" + rgroupname.__str__() + ")\t")
        fout.write("N(" + rgroupname.__str__() + ")\t")
    fout.write("\n")
    
    count = 1
    total_count = seen_genes.__len__() * rgroupids.__len__()
    
    for geneid in seen_genes:        
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        for rgroupid in rgroupids:
            
            count += 1
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
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

    count = 1
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
                geneid_summitstats[geneid] = []
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
    
    sql = "SELECT COUNT(*) from SpeciesunionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
    ngenes = cur.fetchone()[0]
    print "\n    --> The species-union", spunionname, "has", ngenes, "genes with summits in all members."

def plot_summits_for_speciesunion(uid, con):
    """A species-union is a meta union of unions.
    Each union contains multiple replicates from the same species.
    The Species-union compares these unions, across species boundaries."""
    
    cur = con.cursor()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    nunionids = unionids.__len__()
    print "\n. Plotting summits for genes across species", spunionname, unionids

    """seen_genes will be a list of un-aliased gene IDs for all genes that have
    a summit in their nearby regulatory regions, in all the repgroups represented by
    each union. Note that the gene IDs here are un-aliased; in other words,
    they represent the gene ID found from the pillars file. You may need to look up
    alias IDs for this gene's homolog to locate relevant summit scores in different species."""
    seen_genes = []

    """key = gene id, value = list of alias IDs."""
    gene_aliases = {}

    """This loop builds data for a Venn diagram of genes with/without summits in both replicates."""
    count = 1
    total_count = 0
    for unionid in unionids:
        total_count += get_geneids_from_union(con, unionid).__len__()
    
    unionid_genes = {}            
    venn_data = {}
    for unionid in unionids:
        unionname = get_unionname( unionid, con )
        
        """genes is a list of all genes that have summits in nearby regulatory regions in all
        the replicates represented by the union."""
        genes = get_geneids_from_union(con, unionid)
        unionid_genes[unionid] = genes
        translated_genes = []
        for geneid in genes:
            
            """Print a period every 200 iterations, to indicate that this program is still alive."""
            count += 1
            #if count%200 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
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
    vennpath = plot_venn_diagram( venn_data, "summits." + spunionname.__str__() )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for species-union " + spunionname, con)


    unionid_maxsummits = {}
    unionid_nsummits = {}
    unionid_ranksummits = {}
            
    #
    # Write an Excel table using data from UnionSummitStats
    #
    xlpath = "summits." + spunionname + ".xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    #fout.write("geneID\tname\t")
    for unionid in unionids:
        unionname = get_unionname( unionid, con )  
        fout.write("geneID(" + unionname.__str__() + ")\t")
        fout.write("name(" + unionname.__str__() + ")\t")      
        fout.write("max(" + unionname.__str__() + ")\t")
        fout.write("N(" + unionname.__str__() + ")\t")
        fout.write("rank(" + unionname.__str__() + ")\t")
    fout.write("\n")
    count = 1
    total_count = seen_genes.__len__() * unionids.__len__()
    for gid in seen_genes: 
        #fout.write(gid.__str__() + "\t" + get_genename(gid, con) + "\t" )
    
        for unionid in unionids:
            count += 1
            #if count%100 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            if unionid not in unionid_maxsummits:
                unionid_maxsummits[unionid] = []
                unionid_nsummits[unionid] = []
                unionid_ranksummits[unionid] = []
            
            """Find the appropriate gene alias for this union."""            
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
                        fout.write(x[4].__str__() + "\t")
                        
                        unionid_maxsummits[unionid].append( x[2] )
                        unionid_nsummits[unionid].append( x[3] )
                        unionid_ranksummits[unionid].append( x[4] )
                        
            if foundit == False:
                fout.write("---\t---\t0\t0\t--\t")
                unionid_maxsummits[unionid].append( 0 )
                unionid_nsummits[unionid].append( 0 )
                unionid_ranksummits[unionid].append( 0 )
        
        fout.write("\n")
    fout.close()
    
    add_speciesunionfile( xlpath, uid, "Excel table with summit stats for speciesunion " + spunionname, con )

    unionnames = []
    for unionid in unionids:
        unionnames.append(  get_unionname( unionid, con )  )
    

    maxes = []
    ns = []
    ranks = []
    for unionid in unionids:
        maxes.append( unionid_maxsummits[unionid] )
        ns.append( unionid_nsummits[unionid] )
        ranks.append( unionid_ranksummits[unionid] )
    
    scatterdata = maxes[0:nunionids] + ns[0:nunionids] + ranks[0:nunionids]
    scatter_names = unionnames + unionnames + unionnames
    #cranpath = scatter9x3(scatterdata, scatter_names, "summits.9x3." + spunionname, title="Summits " + spunionname)
    width = scatter_names.__len__()
    height = unionnames.__len__()
    filekeyword = "summits." + width.__str__() + "x" + height.__str__() + "."
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword + spunionname, title="Summits " + spunionname)
    add_speciesunionfile(cranpath, uid, "", con)
    add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid, "", con)


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
         
    repgroupname = get_repgroup_name(rgroupid, con)     
    
    repids = get_repids_in_group(rgroupid, con)
    repid1 = repids[0]
    repid2 = repids[1]
    
    cur.execute("SELECT name from Replicates where id=" + repid1.__str__())
    rep1name = cur.fetchone()[0]
    cur.execute("SELECT name from Replicates where id=" + repid2.__str__())
    rep2name = cur.fetchone()[0]
    
    #print "\n. Comparing replicates", rep1name, "and", rep2name
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

    
#     gene_result = {}
#     sql = "SELECT * from Summits where repid=" + repid1.__str__()
#     cur.execute(sql)
#     results = cur.fetchall()
#     for ii in results:
#         geneid = ii[0]

    count = 1    
    genes = get_genes_for_species(con, species1)
    total_count = genes.__len__()
    for g in genes:
        count += 1
        #if count%50 == 0:
            #sys.stdout.write(".")
            #sys.stdout.flush() 
        if count%10==0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
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
    
    sql = "SELECT COUNT(*) from RepgroupGenes where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    ngenes = cur.fetchone()[0]
    print "\n    --> The replicate group", repgroupname, "has", ngenes, "genes with summits in all replicates."                      

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
    vennpath = plot_venn_diagram( venn_data, "summits." + repgroupname)
    add_repgroupfile(vennpath, rgroupid, "Venn diagram comparing genes with summits for replicategroup " + repgroupname, con)
    
 
    """For the following plots, we'll use these data structures."""   
    print "\n. Building a dictionary of summits & genes for", rep1name, "and", rep2name
    rep1_gene_summitscores = {} # key = geneid, value = summit scores from rep1
    rep2_gene_summitscores = {}
    count = 1
        
    genes = get_genes_for_species(con, species1)
    ngenes = genes.__len__()
    for g in genes:
        count += 1
        if count%10 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*(count/float(ngenes) ) ) )
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
    print "\n. Writing a table to", xlpath
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
        fout.write("\n")               
    fout.close()
        
    add_repgroupfile(xlpath, rgroupid, "Excel table with summit stats for replicate group " + repgroupname, con)

    """ Plot score correlations."""
    q_xvals = []
    q_yvals = []
    for g in genes:
        gid = g[0]
        if gid in rep1_gene_summitscores and gid in rep2_gene_summitscores:
            q_xvals.append( max(rep1_gene_summitscores[gid]) )
            q_yvals.append( max(rep2_gene_summitscores[gid]) )
    
    #cranpath = scatter1(xvals,yvals,"summits.score." + repgroupname, xlab=rep1name+": max summit score for gene", ylab=rep2name+": max summit score for gene")
    #add_repgroupfile(cranpath, rgroupid, "R script path to a scatterplot comparing gene-summit scores for replicate group " + repgroupname, con)


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
    rank_xvals = []
    rank_yvals = []
    for gene in gene_xrank:
        rank_xvals.append( gene_xrank[gene] )
        rank_yvals.append( gene_yrank[gene] )   
        
    #cranpath = scatter1(xvals,yvals,"summits.rank." + repgroupname, xlab=rep1name.__str__()+": rank order of genes by max summit score", ylab=rep2name.__str__()+": rank order of genes by max summit score")
    #add_repgroupfile(cranpath, rgroupid, "R script path to a scatterplot comparing the rank of genes for replicate group " + repgroupname, con )

    values = {}
    values["Q score"] =  [q_xvals, q_yvals]
    values["summit rank"] = [rank_xvals, rank_yvals ]
    cranpath = scatter1xn(values, "summits.1x2." + repgroupname, title="Summit Scores", xlab=rep1name, ylab=rep2name, force_square=True )
    add_repgroupfile(cranpath, rgroupid, "R script for multi-panel scatterplot with summits scores and ranks for replicate group " + repgroupname, con)
    add_repgroupfile(re.sub("cran", "pdf", cranpath), rgroupid, "PDF multi-panel scatterplot with summits scores and ranks for replicate group " + repgroupname, con)

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
    #print "1084:", sql
    cur.execute(sql)
    results = cur.fetchall()
    
    gene_results = {}
    for ii in results:
        geneid = ii[1]
        if geneid not in gene_results:
            gene_results[geneid] = []
        gene_results[geneid].append( ii )
    
    meanfe_gids = {} #= bins of mean FE values, with lists of gene ids with that value.
    for geneid in gene_results.keys():
        if gene_results[geneid].__len__() == nrgids:
            max_max = 0
            all_means = []
            all_sums = []
            for ii in gene_results[geneid]:
                rgid = ii[0]    
                if ii[2] > max_max:
                    max_max = ii[2]    
            if max_max not in meanfe_gids:
                meanfe_gids[max_max] = []
            meanfe_gids[max_max].append( geneid )    
    
    meanfes = meanfe_gids.keys()
    meanfes.sort(reverse=True)
    gid_rankmeanfe = {}
    rank = 0
    for meanfe in meanfes:
        newrank = rank + meanfe_gids[meanfe].__len__()/2.0 + 0.5
        for g in meanfe_gids[meanfe]:
            gid_rankmeanfe[g] = newrank
        rank += meanfe_gids[meanfe].__len__()
    
#     count = 1
#     total_count = 0
#     for geneid in gene_results.keys():
#         for ii in gene_results[geneid]:
#             total_count += 1
    
    for geneid in gene_results.keys():        
        """Compute union-wide stats, but only if we have data from all the replicates in this union."""
        if gene_results[geneid].__len__() == nrgids:
            max_max = 0
            all_means = []
            all_sums = []
            
            for ii in gene_results[geneid]:
                rgid = ii[0]
                
                #count += 1
                #if count%20000 == 0:
                #if count%100 == 0:
                #    sys.stdout.write("\r    --> %.1f%%" % (100*(count/float(total_count) ) ) )
                #    sys.stdout.flush()
    
                if ii[2] > max_max:
                    max_max = ii[2]
                all_means.append( ii[3] )
                all_sums.append( ii[4] )
    
            
            sql = "INSERT INTO UnionEnrichmentStats (unionid, geneid, maxenrich, meanenrich, sumenrich, rankmeanenrich)"
            sql += "VALUES (" + unionid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += max_max.__str__() + ","
            sql += mean(all_means).__str__() + ","
            sql += sum(all_sums).__str__() + ","
            sql += gid_rankmeanfe[geneid].__str__() + ")"
            cur.execute(sql)
    con.commit() 
    


def plot_summits_union_helper(unionid, con, keyword=None):
    """This is a helper method for plot_summits_for_union."""
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
    repid_ns = []
    repid_rank = []
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
        
        (x_max, x_n, geneids, x_rank) = get_summit_plot_array_for_replicate_v2(repid, speciesid, con)
        repid_maxs.append( x_max )
        repid_ns.append( x_n)
        repid_rank.append( x_rank)
        repid_name.append( repname )

    #scatter_names = repid_name + repid_name
    scatter_names = []
    for r in repid_name:
        scatter_names.append( r + "-max")
    for r in repid_name:
        scatter_names.append( r + "-ns")
    scatter_data = repid_maxs + repid_ns
    width = scatter_names.__len__()
    height = repid_name.__len__()
    
    """IDR for summits is temporarily disabled."""
    #"""IDR scatters"""
    #filekeyword = "summits.idr." + width.__str__() + "x" + height.__str__() + "." + unionname
    #(cranpath, sinkpath, idr_stats) = scatter_idr_nxm(width, height, scatter_data, scatter_names, filekeyword, title="IDR for summit scores for " + speciesname, xlab="", ylab="")
    #add_unionfile( cranpath, unionid, "R script for a scatterplot with IDR statistics for the summit scores for union " + unionname, con)
    #add_unionfile( re.sub("cran", "pdf", cranpath), unionid, "PDF scatterplot with IDR statistics for the summit scores for union " + unionname, con)

    """Normal data: max, ns, and rank scatters."""
    scatter_names = repid_name + repid_name + repid_name
    scatter_data = repid_maxs + repid_ns + repid_rank
    width = scatter_names.__len__()
    height = repid_name.__len__()
    filekeyword = "summits." + width.__str__() + "x" + height.__str__() + "." + unionname
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="Summit scores for " + speciesname, xlab="", ylab="")
    add_unionfile( cranpath, unionid, "R script for a scatterplot with summit scores for union " + unionname, con)
    add_unionfile( re.sub("cran", "pdf", cranpath), unionid, "PDF scatterplot with summit scores for union " + unionname, con)


def plot_enrichments_for_union(unionid, con, keyword=None):
    """This method makes the plots and excel tables relevant to the enrichment scores
    for all replicates in the union.
    This method assumes that data has been inserted into the table UnionEnrichmentStats"""

    cur = con.cursor()

    unionname = get_unionname(unionid, con)
    rgroupids = get_repgroupids_in_union( unionid, con )
    rgroupids.sort()
    
    print "\n. Plotting enrichments for union", unionname
    plot_enrichment_union_helper(unionid, con)
    
    """Build geneid_results,
    where geneid_results[geneid] = SQL results which is a list of GroupEnrichmentStats objects"""
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
        
    """
    Write an excel table,
    while also filling data structures for the scatterplot.
    Ideally this little routine for writing an excel table would be
    its own function, but it presently relies on too many data
    structures that are intertwined with the scatterplot.
    """
    rgid_maxvals = {}
    rgid_meanvals = {}
    rgid_sumvals = {}
    
    meanvals_rgids = {}
    
    xlpath = "enrich." + unionname.__str__() + ".xls"
    print "\n. Writing a table to", xlpath
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
    count = 1
    total_count = geneid_results.__len__() * rgroupids.__len__()
    for geneid in geneid_results.keys():
        
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        rgid_foundvals = {}
        for rgid in rgroupids:
            
            if rgid not in rgid_foundvals:
                rgid_foundvals[rgid] = []
            
            count += 1
            #if count%50 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*(count/float(total_count) ) ) )
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
                    
                    if x[3] not in meanvals_rgids:
                        meanvals_rgids[ x[3] ] = []
                    meanvals_rgids[x[3]].append( rgid )
 
               
            if foundit == False:
                fout.write("0\t0\t0\t")
                rgid_maxvals[rgid].append(0)
                rgid_meanvals[rgid].append(0)
                rgid_sumvals[rgid].append(0)
                if 0 not in meanvals_rgids:
                    meanvals_rgids[ 0 ] = []
                meanvals_rgids[0].append( rgid )

        fout.write("\n")
    fout.close()
    add_unionfile(xlpath, unionid, "Excel table with enrichment stats for union " + unionname, con)

    
    """Plot a 1xn scatterplot using scores distilled from the replicates."""
    meanvals = meanvals_rgids.keys()
    meanvals.sort( reverse=True)
    rgid_rank = {}
    rank = 0
    for val in meanvals:
        newrank = rank + meanvals_rgids[val].__len__()/2.0 + 0.5
        for rgid in meanvals_rgids[val]:
            rgid_rank[ rgid ] = newrank
        rank += meanvals_rgids[val].__len__()
    values = {}
    values["max enrichment"] =  [ rgid_maxvals[rgroupids[0]], rgid_maxvals[rgroupids[1]]     ]
    values["mean enrichment"] = [ rgid_meanvals[rgroupids[0]], rgid_meanvals[rgroupids[1]]   ]
    values["sum enrichment"]=   [ rgid_sumvals[rgroupids[0]], rgid_sumvals[rgroupids[1]]     ]
    cranpath = scatter1xn(values, "enrich.1x3." + unionname, force_square=True )
    add_unionfile(cranpath, unionid,"R script path to a 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)
    add_unionfile(re.sub("cran", "pdf", cranpath), unionid,"PDF 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)   

    idr_scatter_values = [ rgid_maxvals[rgroupids[0]], rgid_maxvals[rgroupids[1]], rgid_meanvals[rgroupids[0]], rgid_meanvals[rgroupids[1]] ]
    idr_scatter_names = [ get_repgroup_name(rgroupids[0],con) + "-max", get_repgroup_name(rgroupids[1], con) + "-max", get_repgroup_name(rgroupids[0], con) + "-mean", get_repgroup_name(rgroupids[1], con) + "-mean" ]
    
    filekeyword = "enrich.idr.4x2." + unionname
    """idr_stats[gene number in the geneid list][rep ii][rep jj] = lidr"""
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(4, 2, idr_scatter_values, idr_scatter_names, filekeyword, title="", xlab="", ylab="", force_square=True)
    add_unionfile(cranpath, unionid,"" + unionname, con)
    add_unionfile(re.sub("cran", "pdf", cranpath), unionid,"" + unionname, con) 
    
    
    """Now insert the IDR stats into the SQL database."""
    
    """First, clear any previous entries for the pair of this replicate group."""
    for ii in range(0, rgroupids.__len__() ):
        for jj in range(0, rgroupids.__len__() ):
            sql = "DELETE from GeneUnionEnrichIdr where repgroupid1=" + rgroupids[ii].__str__() + " and repgroupid2=" + rgroupids[jj].__str__()
            cur.execute(sql)
            con.commit()
    
    """Next insert the data into SQL."""
    geneids = geneid_results.keys()
    for gg in range(0, geneid_results.__len__() ):
        geneid = geneids[gg]
        if gg in idr_stats:
            for ii in range(0, rgroupids.__len__() ):
                if ii in idr_stats[gg]:
                    for jj in range(0, rgroupids.__len__() ):
                        if jj in idr_stats[gg][ii]:
                            this_idr = idr_stats[gg][ii][jj]
                            ii_repgroupid = rgroupids[ii]
                            jj_repgroupid = rgroupids[jj]
                            sql = "INSERT into GeneUnionEnrichIdr(geneid,repgroupid1,repgroupid2,lidr)"
                            sql += " VALUES("
                            sql += geneid.__str__() + ","
                            sql += ii_repgroupid.__str__() + ","
                            sql += jj_repgroupid.__str__() + ","
                            sql += this_idr.__str__() + ")"
                            #print sql
                            cur.execute(sql)
                            
    con.commit()

def plot_enrichment_union_helper(unionid, con, keyword=None):
    """This is a helper method for plot_enrichments_for_union.
    All the repgroups in the union must be from the same species."""
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
    repid_ranks = []
    repid_names = []    
    speciesname = ""
    speciesid = None
    geneids = None
    for repid in repids:
        sql = "SELECT name from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        repname = cur.fetchone()
        if repname == None:
            repname = repid.__str__()
        else:
            repname = repname[0].__str__()
        
        this_speciesid = get_speciesid_for_repid(repid, con)
        if speciesid == None:
            speciesid = this_speciesid
        elif speciesid != this_speciesid:
            print "\n. Error in plot_enrichment_union_helper"
            print "All the replicate groups in union", unionname, "must be from the same species."
            print "Use a species_union to compare replicates across species."
            exit()
        
        sql = "SELECT name from Species where id=" + speciesid.__str__()
        cur.execute(sql)
        speciesname = cur.fetchone()[0]
        
        if repname.__contains__(speciesname):
            repname = re.sub(speciesname, "", repname)
        
        (x_maxe, x_meane, x_sume, x_ranke, geneids) = get_enrichment_plot_array_for_replicate(repid, speciesid, con)
        repid_means.append( x_meane )
        repid_maxs.append( x_maxe)
        repid_sums.append( x_sume)
        repid_ranks.append( x_ranke )
        repid_names.append( repname )
    
    scatter_names = []
    for repname in repid_names:
        scatter_names.append(repname + "-max")
    for repname in repid_names:
        scatter_names.append(repname + "-mean")

    scatter_data = repid_maxs + repid_means
    width = scatter_names.__len__()
    height = repid_names.__len__()    
    filekeyword = "enrich.idr." + width.__str__() + "x" + height.__str__() + "." + unionname
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_data, scatter_names, filekeyword, title="IDR for Fold Enrichment for " + unionname, xlab="", ylab="")
    add_unionfile(cranpath, unionid, "R script for a scatterplot comparing IDR statistics for fold enrichment for all replicates in the union " + unionname, con)
    add_unionfile(re.sub("cran", "pdf", cranpath), unionid, "PDF comparing IDR statistics for fold enrichment for all replicates in the union " + unionname, con)
    
    #scatter_names = repid_name + repid_name
    scatter_data = repid_maxs + repid_means
    width = scatter_names.__len__()
    height = repid_names.__len__()
    filekeyword = "enrich." + width.__str__() + "x" + height.__str__() + "." + unionname
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="Fold Enrichment for " + unionname, xlab="fold-enrichment", ylab="fold-enrichment")
    add_unionfile(cranpath, unionid, "R script for a scatterplot comparing fold enrichment for all replicates in the union " + unionname, con)
    add_unionfile(re.sub("cran", "pdf", cranpath), unionid, "PDF comparing fold enrichment for all replicates in the union " + unionname, con)

#     """Get a list of genes for this species."""
#     geneids.sort()
#     
#     idr_stats_keys = idr_stats.keys()
#     idr_stats_keys.sort()
#     
#     """Get a list of the comparison names."""
#     compnames = []
#     for ii in idr_stats_keys:
#         for cn in idr_stats[ii]:
#             if cn not in compnames:
#                 compnames.append( cn )
# 
#     return geneids, idr_stats

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
    
    #print "\n. Updating the table GroupEnrichmentStats"
    for ii in range(0, geneids.__len__() ):
        geneid = geneids[ii]
        sql = "INSERT into GroupEnrichmentStats (rgroupid, geneid, maxenrich, meanenrich, sumenrich)"
        sql += " VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ","
        sql += max( [x_maxe[ii], y_maxe[ii]] ).__str__() + ","
        sql += mean( [x_meane[ii], y_meane[ii]] ).__str__() + ","
        sql += sum( [x_sume[ii],y_sume[ii]] ).__str__() + ")"
        cur.execute(sql)
    con.commit()

#
# depricated
#
# def get_summit_plot_array_for_replicate(repid, species, con):
#     """Returns ( x_max[], x_n[] ) for this replicate."""
#     x_max = []
#     x_n = []
#     x_rank = []
#     
#     genes = get_genes_for_species(con, species)
#     geneids = []
#     for ii in genes:
#         geneids.append( ii[0] )
#     geneids.sort()
#     
#     count = 1
#     total_count = geneids.__len__()
#     for geneid in geneids:
#         count += 1
#         sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
#         sys.stdout.flush()
#         
#         x = get_summit_scores_for_gene(geneid, repid, con)
#     
#         if x.__len__() > 0:
#             x_max.append( x[0] )
#             x_n.append( x.__len__() )
#         else:
#             x_max.append(0)
#             x_n.append(0)      
#     return (x_max, x_n, geneids)

def get_summit_plot_array_for_replicate_v2(repid, species, con):
    cur = con.cursor()
    cur.execute("SELECT name from Replicates where id=" + repid.__str__())
    repname = cur.fetchone()[0]
    print "\n. Extracting summit scores for replicate", repname

    qvals = []
    geneids = []
    x_n = []
    x_rank = []

    maxsummit_genes = {} # bins of max summits values, each contains a list of geneIDs
            
    genes = get_genes_for_species(con, species)
    ngenes = genes.__len__()
    count = 0
    for g in genes:
        count += 1
        sys.stdout.write("\r    --> %.1f%%" % (100*(count/float(ngenes) ) ) )
        sys.stdout.flush()
        gid = g[0]
        geneids.append(gid)
        
        x = get_summit_scores_for_gene(gid, repid, con)
        if x.__len__() > 0:
            maxx = max(x)
            qvals.append( maxx )
            if maxx not in maxsummit_genes:
                maxsummit_genes[maxx] = []
            maxsummit_genes[maxx].append( g )
        else:
            qvals.append(0)
            if 0 not in maxsummit_genes:
                maxsummit_genes[0] = []
            maxsummit_genes[0].append( g )
        x_n.append( x.__len__() )
    
    """An extra loop here to rank the max summit scores and then build x_rank"""
    geneid_rank = {}
    maxs = maxsummit_genes.keys()
    maxs.sort(reverse=True)
    rank = 0
    for m in maxs:
        """newrank will be the rank of all the genes in this bin."""
        newrank = rank + maxsummit_genes[m].__len__()/2.0 + 0.5
        for gid in maxsummit_genes[m]:
            geneid_rank[gid] = newrank
            #x_rank.append(newrank)
        rank += maxsummit_genes[m].__len__()
    for g in genes:
        x_rank.append( geneid_rank[g] )
    return (qvals, x_n, geneids, x_rank)


def sort_gene_values():
    """Get the rank of genes according to an arbitrary score."""
    pass

def get_enrichment_plot_array_for_replicate(repid, species, con):
    """Returns ( x_maxe[], x_meane[], x_sume[] ) for this replicate."""
    x_maxe = []
    x_meane = []
    x_sume = []
    x_ranke = []
    x = get_enrichment_stats_for_replicate( repid, con )
    
    genes = get_genes_for_species(con, species)
    geneids = []
    for ii in genes:
        geneids.append( ii[0] )
    geneids.sort()
    
    meane_geneids = {}
    
    for geneid in geneids:
        if geneid in x:
            x_maxe.append( x[geneid][0] )
            x_meane.append( x[geneid][1] )
            x_sume.append( x[geneid][2] )
            
            if x[geneid][1] not in meane_geneids:
                meane_geneids[  x[geneid][1] ] = []
            meane_geneids[  x[geneid][1] ].append( geneid )
        else:
            x_maxe.append(0)
            x_meane.append(0)
            x_sume.append(0)
            if 0 not in meane_geneids:
                meane_geneids[ 0 ] = []
            meane_geneids[ 0 ].append( geneid )
    means = meane_geneids.keys()
    means.sort(reverse=True)
    rank = 0
    geneid_rank = {}
    for m in means:
        newrank = rank + meane_geneids[m].__len__()/2.0 + 0.5
        for geneid in meane_geneids[m]:
            geneid_rank[geneid] = newrank
        rank += meane_geneids[m].__len__()
    for geneid in geneids:
        x_ranke.append( geneid_rank[geneid] )
    return (x_maxe, x_meane, x_sume, x_ranke, geneids)

def get_enrichment_plot_array_for_union(unionid, gene_aliases, con):
    """Returns ( x_maxe[], x_meane[], x_sume[] ) averaged (or maxed) from all repgroups in this union."""
    x_maxe = []
    x_meane = []
    x_sume = []
    x_ranke = []
    x = get_enrichment_stats_for_union( unionid, con )
                    
    meane_geneids = {}
                    
    for geneid in x:
        translated_geneid = gene_aliases[geneid]
        
        if geneid in x:
            x_maxe.append( x[geneid][0] )
            x_meane.append( x[geneid][1] )
            x_sume.append( x[geneid][2] )
            
            if x[geneid][1] not in meane_geneids:
                meane_geneids[ x[geneid][1] ] = []
            meane_geneids[ x[geneid][1] ].append( geneid )
        else:
            x_maxe.append(0)
            x_meane.append(0)
            x_sume.append(0)
            if 0 not in meane_geneids:
                meane_geneids[ x[geneid][1] ] = []
            meane_geneids[ 0 ].append( geneid ) 
    means = meane_geneids.keys()
    means.sort(reverse=True)
    rank = 0
    geneid_rank = {}
    for m in means:
        newrank = rank + meane_geneids[m].__len__()/2.0 + 0.5
        for geneid in meane_geneids[m]:
            geneid_rank[geneid] = newrank
        rank += meane_geneids[m].__len__()
    for geneid in geneids:
        x_ranke.append(geneid_rank[geneid])
    return (x_maxe, x_meane, x_sume, x_ranke, genelist)

def plot_enrichments_for_reps_in_group(rgroupid, con, repgroupname=None, repids=None):
    """Plots a comparison of FE between replicates 1 and 2 of the replicates in the repgroup.
    If there are more than 2 replicates, then their data won't be plotted."""
    #
    # to-do: extend this method to work for N replicates
    
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
    
    """Build the scatterplot"""
    speciesid1 = get_speciesid_for_repid(repid1, con)
    (x_maxe, x_meane, x_sume, x_ranke, geneids) = get_enrichment_plot_array_for_replicate(repid1, speciesid1, con)
    speciesid2 = get_speciesid_for_repid(repid2, con)
    (y_maxe, y_meane, y_sume, y_ranke, geneids) = get_enrichment_plot_array_for_replicate(repid2, speciesid2, con)
    values = {}
    values["max enrichment"] =  [x_maxe, y_maxe]
    values["mean enrichment"] = [x_meane, y_meane ]
    values["sum enrichment"]=   [x_sume, y_sume ]
    values["rank enrichment"] = [x_ranke, y_ranke]
    cranpath = scatter1xn(values, "enrich.1x4." + repgroupname, xlab=rep1name, ylab=rep2name, force_square=True )
    add_repgroupfile(cranpath, rgroupid, "R script for multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
    add_repgroupfile(re.sub("cran", "pdf", cranpath), rgroupid, "PDF multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
    
    """Also plot the IDR stats"""
    idr_scatter_values = [ values["max enrichment"][0], values["max enrichment"][1], values["mean enrichment"][0], values["mean enrichment"][1] ]
    filekeyword = "enrich.idr.4x2." + repgroupname
    idr_scatter_names = ["rep1-max", "rep2-max","rep1-mean", "rep2-mean"]
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(4, 2, idr_scatter_values, idr_scatter_names, filekeyword, title="", xlab="", ylab="", force_square=True)
    add_repgroupfile(re.sub("cran", "pdf",cranpath),rgroupid,"PDF with IDR scatterplots for replicate group " + repgroupname, con)
    
    """Update the IDR stats into the database"""

    """First, clear any previous entries for the pair of this replicate group."""
    sql = "DELETE from GeneRepgroupEnrichIdr where repid1=" + repid1.__str__() + " and repid2=" + repid2.__str__()
    cur.execute(sql)
    con.commit()
    
    for gg in range(0, geneids.__len__() ):
        geneid = geneids[gg]
        if gg in idr_stats:
            for ii in range(0, repids.__len__() ):
                if ii in idr_stats[gg]:
                    for jj in range(0, repids.__len__() ):
                        if jj in idr_stats[gg][ii]:
                            this_idr = idr_stats[gg][ii][jj]
                            ii_repid = repids[ii]
                            jj_repid = repids[jj]
                            sql = "INSERT into GeneRepgroupEnrichIdr(geneid, repid1, repid2, lidr)"
                            sql += " VALUES(" + geneid.__str__() + ","
                            sql += ii_repid.__str__() + ","
                            sql += jj_repid.__str__() + ","
                            sql += this_idr.__str__() + ")"
                            #print sql
                            cur.execute(sql)
    con.commit()
            
    """Write an Excel Table"""
    xlpath = "enrich." + repgroupname + ".xls"
    print "\n. Writing a table with fold-enrichment and IDR to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for rgname in [rep1name,rep2name]:
        fout.write("max(" + rgname + ")\t")
        fout.write("mean(" + rgname + ")\t")
        fout.write("sum(" + rgname + ")\t")
        fout.write("rank(" + rgname + ")\t")
    for (ii,jj) in value_pairs:   
        fout.write("IDR(" + idr_scatter_names[ii].__str__() + "-" + idr_scatter_names[jj].__str__() + ")\t")     
    fout.write("\n")
    for gg in range(0, geneids.__len__() ):
        genename = get_genename(geneids[gg], con)
        if genename == None:
            continue
        fout.write(geneids[gg].__str__() + "\t" + genename + "\t" )
        fout.write(x_maxe[gg].__str__() + "\t")
        fout.write(x_meane[gg].__str__() + "\t")
        fout.write(x_sume[gg].__str__() + "\t")
        fout.write(x_ranke[gg].__str__() + "\t")
        fout.write(y_maxe[gg].__str__() + "\t")
        fout.write(y_meane[gg].__str__() + "\t")
        fout.write(y_sume[gg].__str__() + "\t")
        fout.write(y_ranke[gg].__str__() + "\t")
        """Write IDR stats if they exist"""
        if gg in idr_stats:
            for (ii,jj) in value_pairs:
                #print "1389:", ii, jj
                if ii in idr_stats[gg]:
                    if jj in idr_stats[gg][ii]:
                        fout.write(idr_stats[gg][ii][jj].__str__() + "\t")
                    else:
                        fout.write("\t")
                else:
                    fout.write("\t")
        fout.write("\n")
    fout.close()

    """Add an entry to the SQL database, listing the Excel file."""
    add_repgroupfile(xlpath,rgroupid,"Excel table with enrichment stats for replicate group " + repgroupname, con)


def compute_enrichments_for_speciesunion(uid, con):
    """Fills up the table SpeciesunionEnrichmentStats and also writes an excel tables."""
    
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
        
    
    count = 1
    total_count = 0
    for gid in geneid_results.keys():
        if geneid_results[gid].__len__() == nuids:
            total_count += nuids
    
    for gid in geneid_results.keys():
        """Only keep those genes that have enrichment scores in all the replicates in this union."""
        if geneid_results[gid].__len__() == nuids:    
            max_max = 0
            all_means = []
            all_sums = []
                   
            for ii in geneid_results[gid]:
                unionid = ii[0]
               
                count += 1
                #if count%50 == 0:
                sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
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
    
    print "\n. Plotting enrichments acorss species:", spunionname, unionids
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
    
    geneid_results = {}

    """A dictionary of aliased genes. This will be useful when we assign homology between
    enrichment values."""
    gene_aliases = {}
    for ii in results:       
        translated_id = get_geneid_from_aliasid( ii[1], con )
        if translated_id != None:
            """Have we seen this gene before?"""
            if translated_id not in gene_aliases:
                gene_aliases[translated_id] = []
            gene_aliases[ translated_id ].append( ii[1]) 

        """Remember this gene."""
        if translated_id not in geneid_results:
            geneid_results[ translated_id ] = []
            
        """Store the translation pairing"""
        geneid_results[translated_id].append(ii)
       
    #
    # Write an excel table,
    # while also filling data structures for the scatterplot
    #
    unionid_maxvals = {}
    unionid_meanvals = {}
    unionid_sumvals = {}
    unionid_rankmeanvals = {}
    
    xlpath = "enrich." + spunionname.__str__() + ".xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    #fout.write("GeneID\tGeneName\t")
    for unionid in unionids:
        unionname = get_unionname(unionid, con)
        fout.write("GeneID(" + unionname + ")\t")
        fout.write("GeneName(" + unionname + ")\t")
        fout.write("max(" + unionname + ")\t")
        fout.write("mean(" + unionname + ")\t")
        fout.write("sum(" + unionname + ")\t")
        fout.write("rankmean(" + unionname + ")\t")

        unionid_maxvals[unionid] = []
        unionid_meanvals[unionid] = []
        unionid_sumvals[unionid] = []
        unionid_rankmeanvals[unionid] = []
        
    count = 1
    total_count = 0
    for gid in geneid_results.keys():
        for unionid in unionids:
            for x in geneid_results[gid]:
                total_count += 1
       
    fout.write("\n")
    for gid in geneid_results.keys():
        max_max = 0
        all_means = []
        all_sums = []
        #fout.write(gid.__str__() + "\t" + get_genename(gid, con) + "\t" )        
        
        for unionid in unionids: # find the unionid that mathches result x
            
            """Did we find the SQL query that corresponds to this unionid?"""
            foundit = False
            
            for x in geneid_results[gid]:    
                count += 1
                sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
        
                if x[0] == unionid and foundit == False:
                    foundit = True
                    fout.write(x[1].__str__() + "\t")
                    fout.write(get_genename(x[1], con) + "\t")
                    fout.write(x[2].__str__() + "\t")
                    fout.write(x[3].__str__() + "\t")
                    fout.write(x[4].__str__() + "\t")
                    fout.write(x[5].__str__() + "\t")
    
                    unionid_maxvals[unionid].append(x[2])
                    unionid_meanvals[unionid].append(x[3])
                    unionid_sumvals[unionid].append(x[4])
                    unionid_rankmeanvals[unionid].append(x[5])
                    
                    if x[2] > max_max:
                        max_max = x[2]
                    #all_means.append( x[3] )
                    all_sums.append( x[4] )
                if foundit == True:
                    continue
                    
            #else:
            if foundit == False:
                fout.write("---\t---\t0\t0\t0\t--\t")
                unionid_maxvals[unionid].append(0)
                unionid_meanvals[unionid].append(0)
                unionid_sumvals[unionid].append(0)
                unionid_rankmeanvals[unionid].append(0)
                #all_means.append( 0 )
        fout.write("\n")
    fout.close()

    add_speciesunionfile(xlpath, uid, "Excel table with enrichment stats for species union " + spunionname, con )


    #plot_enrichment_12x4_speciesunion(uid, geneid_results, con )
    unionidnames = []
    for unionid in unionids:    
        unionidnames.append( get_unionname(unionid, con) )  
        
    scatterdata = []
    for unionid in unionids:
        scatterdata.append(  unionid_maxvals[unionid] )
    for unionid in unionids:
        scatterdata.append(  unionid_meanvals[unionid] )
    #for unionid in unionids:
    #    scatterdata.append(  unionid_rankmeanvals[unionid] )

    scatter_names = unionidnames + unionidnames
    #cranpath = scatter12x3(scatterdata, scatter_names, "enrich.12x3." + spunionname, title="Fold Enrichment " + spunionname, xlab="fold-enrichment", ylab="fold-enrichment") 
    width = scatter_names.__len__()
    height = unionidnames.__len__()
    filekeyword = "enrich." + width.__str__() + "x" + height.__str__() + "."
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword + spunionname, title="Fold Enrichment " + spunionname, xlab="fold-enrichment", ylab="fold-enrichment") 
    add_speciesunionfile(cranpath, uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )
    add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )    
     
    #print "\n1399: returning from plot_en.._spunion early."
    #return
   
    