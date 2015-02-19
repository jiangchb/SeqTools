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
    """This method extracts data from RepgroupSummits, compares those summits, and then
    writes data to UnionSummits."""
    cur = con.cursor()
    sql = "DELETE from UnionSummits where unionid=" + unionid.__str__()
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

    
    print "\n. Building a dictionary of summits & genes for union '" + unionname + "'"
    
    """What Genes Apply?
    
    
    
    
    """
    sql = "select distinct gene from GeneSummits where "
    sql_pieces = []
    for repgroupid in rgroupids:
        sql_pieces.append( " gene in (select distinct geneid from RepgroupSummits where repgroupid=" + repgroupid.__str__() + ") " )
    sql += " and ".join( sql_pieces )
    cur.execute(sql)
    x = cur.fetchall()
    seen_genes = []
    for ii in x:
        seen_genes.append( ii[0] )
    
    #print "\n. chipseqdb_plot 573 - repgroupid ", repgroupid, "has N genes:", seen_genes.__len__()
    
    for geneid in seen_genes:
        
        """Find the maximum scoring summit."""        
        sql = "select id, max(score) from Summits where "
        sql += " id in (select maxsummitid from RepgroupSummits where geneid=" + geneid.__str__()
        sql_pieces = []
        for repgroupid in rgroupids:
            sql_pieces.append(" repgroupid=" + repgroupid.__str__() )
        sql += " and (" + " or ".join( sql_pieces ) + ")"
        sql += ")"
        #print "76:", sql
        cur.execute(sql)
        max_summit_id = cur.fetchone()[0]
        
        """Find the nearest summit."""
        sql = "select summit, min( abs(distance) ) from GeneSummits where "
        sql += " summit in (select nearestsummitid from RepgroupSummits where geneid=" + geneid.__str__()
        sql_pieces = []
        for repgroupid in rgroupids:
            sql_pieces.append(" repgroupid=" + repgroupid.__str__() )
        sql += " and (" + " or ".join( sql_pieces ) + ")"
        sql += ")"
        #print "88:", sql
        cur.execute(sql)
        nearest_summit_id = cur.fetchone()[0] 
        
        """Find the mean of max peak scores, across replicates."""
        max_scores = []
        for repgroupid in rgroupids:        
            sql = "select max(score) from Summits where "
            sql += " id in (select maxsummitid from RepgroupSummits where geneid=" + geneid.__str__()
            sql += " and repgroupid=" + repgroupid.__str__()
            sql +=  ")"
            cur.execute(sql)
            max_score = cur.fetchone()[0]
            max_scores.append( max_score )

        sql = "INSERT INTO UnionSummits(unionid, geneid, maxsummitid, nearestsummitid, mean_maxsummitscore) "
        sql += " VALUES("
        sql += unionid.__str__() + ","
        sql += geneid.__str__() + ","
        sql += max_summit_id.__str__() + ","
        sql += nearest_summit_id.__str__()  + ","
        sql += mean( max_scores ).__str__()
        sql += ")"
        cur.execute(sql)
    con.commit()
 
    sql = "SELECT COUNT(*) from UnionSummits where unionid=" + unionid.__str__()
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

    """Scatterplot."""
    plot_summits_union_helper(unionid, con)
    
    """Venn diagram of genes with/without summits in both replicates."""
    venn_data = {}
    rgroupids
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )
        genes = get_geneids_from_repgroup(con, rgroupid)
        venn_data[ rgroupname ] = genes
    vennpath = plot_venn_diagram( venn_data, unionname.__str__() + ".summits" )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for union " + unionname, con)
        
    """Next we build some hashes that store data which will be used to write a big Excel table.
        These hashes prevent excessive SQL queries."""
    repgroup_summitid_score = {} # scores for all summits that may be relevant
    for rgroupid in rgroupids:
        repgroup_summitid_score[rgroupid] = {}
        sql = "select id, score from Summits where id in (select maxsummitid from RepgroupSummits where repgroupid=" + rgroupid.__str__() + ")"
        sql += " or id in (select nearestsummitid from RepgroupSummits where repgroupid=" + rgroupid.__str__() + ")"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            repgroup_summitid_score[rgroupid][ ii[0] ] = ii[1]
    
    repgroup_gene_maxsummitid = {}
    repgroup_gene_nearestsummitid = {}
    for rgroupid in rgroupids:
        repgroup_gene_maxsummitid[rgroupid] = {}
        repgroup_gene_nearestsummitid[rgroupid] = {}
        sql = "SELECT repgroupid, geneid, maxsummitid, nearestsummitid from RepgroupSummits where repgroupid=" + rgroupid.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            repgroup_gene_maxsummitid[rgroupid][ ii[1] ] = ii[2]
            repgroup_gene_nearestsummitid[rgroupid][ ii[1] ] = ii[3]
            
    summit_gene_mindistance = {}
    for rgroupid in rgroupids:
        sql = "select distance, gene, summit from GeneSummits where summit in "
        sql += "(SELECT maxsummitid from RepgroupSummits where repgroupid=" + rgroupid.__str__() + ") "
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            if ii[2] not in summit_gene_mindistance:
                summit_gene_mindistance[ ii[2] ] = {}
            summit_gene_mindistance[ ii[2] ][ ii[1] ] = ii[0]

        sql = "select distance, gene, summit from GeneSummits where summit in "
        sql += " (SELECT nearestsummitid from RepgroupSummits where repgroupid=" + rgroupid.__str__() + ")"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            if ii[2] not in summit_gene_mindistance:
                summit_gene_mindistance[ ii[2] ] = {}
            summit_gene_mindistance[ ii[2] ][ ii[1] ] = ii[0]
    
    """ seen_genes is a list of gene IDs that have summits in all the replicates in this union."""
    seen_genes = []
    sql = "select geneid from UnionSummits where unionid=" + unionid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        seen_genes.append( ii[0] )
    
    """ Write the Excel table listing the genes and their scores in the repgroups."""
    xlpath = unionname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("gene_ID\tgene_name\t")
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )        
        fout.write("max_summit_Qval("       + rgroupname.__str__() + ")\t")
        fout.write("max_summit_dist("       + rgroupname.__str__() + ")\t")
        fout.write("nearest_summit_Qval("   + rgroupname.__str__() + ")\t")
        fout.write("nearest_summit_dist("   + rgroupname.__str__() + ")\t")
    fout.write("\n")
    
    count = 1
    total_count = seen_genes.__len__() * rgroupids.__len__()

    """For each gene that exists in all replicate groups in the union."""    
    for geneid in seen_genes:        
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
    
        """For each replicate."""
        for rgroupid in rgroupids:
            
            count += 1
            if count%20 == 0:
                sys.stdout.write("\r    --> (Excel table) %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
            
            maxsummitid = repgroup_gene_maxsummitid[rgroupid][geneid]
            nearestsummtid = repgroup_gene_nearestsummitid[rgroupid][geneid]
            this_score = repgroup_summitid_score[rgroupid][ maxsummitid ]
            fout.write(this_score.__str__() + "\t")
            
            d = summit_gene_mindistance[maxsummitid][geneid]
            fout.write(d.__str__() + "\t")
            
            this_score = repgroup_summitid_score[rgroupid][ nearestsummtid ]
            fout.write(this_score.__str__() + "\t")
            
            d = summit_gene_mindistance[nearestsummtid][geneid]
            fout.write(d.__str__() + "\t")

        fout.write("\n")
    fout.close()
    
    add_unionfile(xlpath, unionid, "Excel table with summit stats for union " + unionname, con)


def compute_summits_for_speciesunion(uid, con):
    """A species-union is a meta union of unions.
        Each union contains multiple replicates from the same species.
        The Species-union compares these unions, across species boundaries."""
    
    cur = con.cursor()
    sql = "DELETE from SpeciesunionSummits where spunionid=" + uid.__str__()
    cur.execute(sql)
    con.commit()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    print "\n. Comparing summits for genes across species", spunionname

    """Get data from UnionSummits"""
    geneid_maxsummits = {} # key = translated gene ID, value = list of SQL results from RepgroupSummitStats
    geneid_nearestsummits = {} # key = translated gene ID, value = hash; key = summit id, value = distance
    translatedgene2ergene = {}
    gene_unions = {} # key = geneid, value = list of repgroup IDs with summit values for this gene.
    for unionid in unionids:
        unionname = get_unionname( unionid, con )
        
        sql = "SELECT * from UnionSummits where unionid=" + unionid.__str__()
        cur.execute(sql)
        results = cur.fetchall()

        for ii in results:
            erg_geneid = ii[1]
            """Translate the geneid to be species agnostic, using the pillars data"""
            geneid = get_geneid_from_aliasid( erg_geneid, con )
            max_summit_id = ii[2]        
            nearest_summit_id = ii[3]
            sql = "SELECT distance from GeneSummits where gene=" + erg_geneid.__str__() + " and summit=" + nearest_summit_id.__str__()
            cur.execute(sql)
            distance = cur.fetchone()[0]
            
            if geneid not in geneid_maxsummits:
                geneid_maxsummits[geneid] = []
            geneid_maxsummits[geneid].append( max_summit_id )
            
            if geneid not in geneid_nearestsummits:
                geneid_nearestsummits[geneid] = {}
            geneid_nearestsummits[geneid][nearest_summit_id] = distance

            if geneid not in gene_unions:
                gene_unions[geneid] = []
            if unionid not in gene_unions[geneid]:
                gene_unions[geneid].append( unionid )
                
    geneid_max = {}
    geneid_nearest = {} 
    for geneid in gene_unions:    
        if gene_unions[geneid].__len__() == unionids.__len__():
            """Only deal with geneids that exist in all the unions of this speciesunion."""
            """i.e., there is at least one summit at this gene in every replicate group."""
            max_max_summit_id = None
            max_max_summit_score = 0
            for summitid in geneid_maxsummits[geneid]:
                sql = "SELECT score from Summits where id=" + summitid.__str__()
                cur.execute(sql)
                this_score = cur.fetchone()[0]
                if this_score > max_max_summit_score:
                    max_max_summit_score = this_score
                    max_max_summit_id = summitid
            if max_max_summit_id == None:
                print "\n. Error. Checkpoint 86. there is no max summitid for gene", geneid, "in repgroup", repgroupname
                exit()
            geneid_max[geneid] = max_max_summit_id
        if geneid_nearestsummits[geneid].__len__() == unionids.__len__():
            nearest_nearest_summit_id = None
            nearest_nearest_summit_distance = None
            for summitid in geneid_nearestsummits[geneid]:
                this_distance = abs( geneid_nearestsummits[geneid][summitid] )
                #print "294:", this_distance
                if nearest_nearest_summit_distance == None:
                    nearest_nearest_summit_distance = this_distance
                    nearest_nearest_summit_id = summitid
                elif nearest_nearest_summit_distance > this_distance:
                    nearest_nearest_summit_distance = this_distance
                    nearest_nearest_summit_id = summitid
                    
            if nearest_nearest_summit_id == None:
                print "\n. Error, there is no nearest summitid for gene", geneid
                print "\n. chipseqdb_plot.py 103"
                exit()
            geneid_nearest[geneid] = nearest_nearest_summit_id
    
    """Insert the max and nearest summit stats into UnionSummits."""
    for geneid in geneid_max:    
        if geneid in geneid_nearest:
            
            sql = "INSERT INTO SpeciesunionSummits (spunionid, geneid, maxsummitid, nearestsummitid)"
            sql += " VALUES(" + uid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += geneid_max[geneid].__str__() + ","
            sql += geneid_nearest[geneid].__str__() + ")"
            #print sql
            cur.execute(sql)
    con.commit()
 
    sql = "SELECT COUNT(*) from SpeciesunionSummits where spunionid=" + uid.__str__()
    cur.execute(sql)
    ngenes = cur.fetchone()[0]
    print "\n    --> The species union", spunionname, "has", ngenes, "genes with summits in all species."

def plot_summits_for_speciesunion(uid, con):
    """A species-union is a meta union of unions.
    Each union contains multiple replicates from the same species.
    The Species-union compares these unions, across species boundaries."""
    
    cur = con.cursor()
    
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    nunionids = unionids.__len__()
    print "\n. Plotting summits for genes across species", spunionname, unionids

    """seen_genes will be a list of translated gene IDs for all genes that have
    a summit in their nearby regulatory regions, in all the repgroups represented by
    each union. Note that the gene IDs here are translated; in other words,
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
            if count%200 == 0:
                sys.stdout.write("\r    --> (venn diagram) %.1f%%" % (100*count/float(total_count)) )
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
    vennpath = plot_venn_diagram( venn_data, spunionname + ".summits" )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for species-union " + spunionname, con)
        
    """ Write an Excel table using data from UnionSummits
    """
    xlpath = spunionname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    #fout.write("geneID\tname\t")
    for unionid in unionids:
        unionname = get_unionname( unionid, con )  
        fout.write("gene_ID(" + unionname.__str__() + ")\t")
        fout.write("gene_name(" + unionname.__str__() + ")\t")      
        fout.write("max_summit_Qval(" + unionname + ")\t")
        fout.write("max_summit_dist(" + unionname + ")\t")
        fout.write("nearest_summit_Qval(" + unionname + ")\t")
        fout.write("nearest_summit_dist(" + unionname+ ")\t")

    fout.write("\n")
    count = 1
    total_count = seen_genes.__len__() * unionids.__len__()
    
    """These two hashes are used later for the scatterplots."""
    union_gene_qval = {}
    union_gene_nearestqval = {}
    for unionid in unionids:
        union_gene_qval[unionid] = {}
        union_gene_nearestqval[unionid] = {}

    """This hash is a shortcut to avoid excessive SQL queries."""
    unionid_summitid_score = {} # scores for all summits relevant to this union
    for unionid in unionids:
        unionid_summitid_score[unionid] = {}
        sql = "select id, score from Summits where id in (select maxsummitid from  UnionSummits where unionid=" + unionid.__str__() + ")"
        sql += " or id in (select nearestsummitid  from UnionSummits where unionid=" + unionid.__str__() + ")"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            unionid_summitid_score[unionid][ ii[0] ] = ii[1]
    
    """For each translated gene ID."""
    for gid in seen_genes: 
        #fout.write(gid.__str__() + "\t" + get_genename(gid, con) + "\t" )
        for unionid in unionids:
            
            
            
            count += 1
            if count%10 == 0:
                sys.stdout.write("\r    --> (Excel table) %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
                        
            #print "425:", gid, unionid, gene_aliases[gid].__len__()
            
            """Does this translated gene have a summit in the union?"""           
            foundgid = False
            for geneid in gene_aliases[gid]:    
                if foundgid == False:                                      
                    sql = "SELECT unionid, geneid, maxsummitid, nearestsummitid, mean_maxsummitscore from UnionSummits where unionid=" + unionid.__str__()
                    sql += " and geneid=" + geneid.__str__()
                    cur.execute(sql)
                    res = cur.fetchall()
                    if res == None or res.__len__() == 0:
                        continue
                    else:
                        foundgid = True
                    
                    x = res[0]
                    
                    foundit = True
                    fout.write(geneid.__str__() + "\t")
                    fout.write(get_genename(geneid, con) + "\t")
                    
                    #sql = "SELECT score from Summits where id=" + union_gene_summit[unionid][geneid].__str__()
                    #sql = "SELECT score from Summits where id=" + x[2].__str__()
                    #cur.execute(sql)
                    #this_score = cur.fetchone()[0]
                    this_score = unionid_summitid_score[unionid][ x[2] ]
                    fout.write(this_score.__str__() + "\t")
                    union_gene_qval[unionid][gid] = this_score
                    
                    sql = "SELECT min(distance) from GeneSummits where summit=" + x[2].__str__() + " and gene=" + geneid.__str__()
                    cur.execute(sql)
                    d = cur.fetchone()[0]
                    fout.write(d.__str__() + "\t")
                    
                    #sql = "SELECT score from Summits where id=" + x[3].__str__()
                    #cur.execute(sql)
                    #this_score = cur.fetchone()[0]
                    this_score = unionid_summitid_score[unionid][ x[3] ]
                    fout.write(this_score.__str__() + "\t")
                    union_gene_nearestqval[unionid][gid] = this_score

                    sql = "SELECT min(distance) from GeneSummits where summit=" + x[3].__str__() + " and gene=" + geneid.__str__()
                    cur.execute(sql)
                    d = cur.fetchone()[0]
                    fout.write(d.__str__() + "\t") 
            
            if foundgid == False:
                fout.write("---\t---\t---\t---\t---\t---\t")       
        fout.write("\n")
    fout.close()
    
    add_speciesunionfile( xlpath, uid, "Excel table with summit stats for speciesunion " + spunionname, con )

    """Find genes that have summits in all the unions."""
    sql = "select geneid from SpeciesunionSummits where spunionid=" + uid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    spgenes = []
    for ii in x:
        spgenes.append( ii[0] )
    
    """scatterdata is for the scatter_nxm function."""
    genes = []
    for unionid in unionids:
        for geneid in spgenes:
            if geneid not in union_gene_qval[unionid] or geneid not in union_gene_nearestqval[unionid]:
                pass
            else:
                genes.append( geneid )
    
    scatterdata = []
    for unionid in unionids:
        maxes = []
        for geneid in genes:
            maxes.append( union_gene_qval[unionid][geneid] )
        scatterdata.append( maxes )
    for unionid in unionids:
        nearest_maxes = []
        for geneid in genes:
            nearest_maxes.append( union_gene_nearestqval[unionid][geneid] )
        scatterdata.append( nearest_maxes )    
        
    scatter_names = []
    unit_labels = []
    for unionid in unionids:
        scatter_names.append( "max Qval (" + get_unionname( unionid, con ) + ")")
        unit_labels.append( "max. Qval")
    for unionid in unionids:
        scatter_names.append("nearest Qval (" +  get_unionname( unionid, con ) + ")")
        unit_labels.append( "nearest Qval")

    width = scatter_names.__len__()
    height = unionids.__len__()
    filekeyword = spunionname + ".summits"
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword, title="Summits " + spunionname, unit_labels=unit_labels)
    if cranpath != None:
        add_speciesunionfile(cranpath, uid, "", con)
        add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid, "", con)

def compute_summits_for_reps_in_group(rgroupid, con):
    """This is a long method (sorry). 
    It correlates the summits from all the replicates in a group,
    using a variety of methods.
    It assumes that all the replicates come from the same species.
    
    This method updates the table named RepgroupSummits with summits
    that are found in all the replicates within the replicate group with ID = rgroupid.
    """
     
    cur = con.cursor()
    sql = "DELETE from RepgroupSummits where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    con.commit()
         
    repgroupname = get_repgroup_name(rgroupid, con)     
    
    repids = get_repids_in_group(rgroupid, con)
    repid_repname = {}
    
    """Get display names for the replicates."""
    for repid in repids:
        sql = "SELECT name from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        repname = cur.fetchone()
        if repname == None:
            repname = repid.__str__()
        else:
            repname = repname[0].__str__()
        repid_repname[repid] = repname
    
    """Verify that all replicates belong to the same species."""
    repid_species = {}
    for repid in repids:
        cur.execute("SELECT species from Replicates where id=" + repid.__str__())
        data = cur.fetchone()
        if data != None:
            repid_species[repid] = data[0]
        else:
            repid_species[repid] = None
    for ii in range(0, repids.__len__()-1 ):
        if repid_species[ repids[ii] ] != repid_species[ repids[ii+1] ]:
            print "\n. Error while computing summits for replicates in the group", repgroupname
            print ". It appears that not all the replicates belong to the same species."
            exit()

    """Get the summit scores for every gene for every replicate."""
    print "\n. Building a dictionary of summits & genes for", repgroupname
         
    sql = "select distinct gene from GeneSummits where "
    sql_pieces = []
    for repid in repids:
        sql_pieces.append( " gene in (select distinct gene from GeneSummits where summit in (select id from Summits where replicate=" + repid.__str__() + ")) ")
    sql += " and ".join( sql_pieces )
    cur.execute(sql)
    x = cur.fetchall()
    seen_genes = []
    for ii in x:
        seen_genes.append( ii[0] )
    
    #print "\n. chipseqdb_plot 573 - repgroupid ", rgroupid, repgroupname, "has N genes:", seen_genes.__len__()
    
    for geneid in seen_genes:
        
        """Find the maximum scoring summit."""
        sql = "select id, max(score) from Summits where "
        sql += " id in (select summit from GeneSummits where gene=" + geneid.__str__() + ")"
        sql_pieces = []
        for repid in repids:
            sql_pieces.append(" replicate=" + repid.__str__() )
        sql += "and (" + " or ".join( sql_pieces ) + ")"
        cur.execute(sql)
        max_summit_id = cur.fetchone()[0]
        
        """Find the nearest summit."""
        sql = "select summit, min( abs(distance) ) from GeneSummits where gene=" + geneid.__str__() 
        sql += " and summit in (select summit from Summits where "
        sql_pieces = []
        for repid in repids:
            sql_pieces.append(" replicate=" + repid.__str__() )
        sql += " or ".join( sql_pieces )
        sql += " )"
        cur.execute(sql)
        nearest_summit_id = cur.fetchone()[0]
        
        """Find the mean of max peak scores, across replicates."""
        max_scores = []
        for repid in repids:        
            sql = "select max(score) from Summits where "
            sql += " id in (select summit from GeneSummits where gene=" + geneid.__str__() + ")"
            sql += " and replicate=" + repid.__str__()
            cur.execute(sql)
            max_score = cur.fetchone()[0]
            max_scores.append( max_score )

        sql = "INSERT INTO RepgroupSummits(repgroupid, geneid, maxsummitid, nearestsummitid, mean_maxsummitscore) "
        sql += " VALUES("
        sql += rgroupid.__str__() + ","
        sql += geneid.__str__() + ","
        sql += max_summit_id.__str__() + ","
        sql += nearest_summit_id.__str__() + ","
        sql += mean( max_scores ).__str__()
        sql += ")"
        cur.execute(sql)
    con.commit()
    
    sql = "SELECT COUNT(*) from RepgroupSummits where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    ngenes = cur.fetchone()[0]
    print "\n    --> The replicate group", repgroupname, "has", ngenes, "genes with summits in all replicates."                      

    """Write the Excel file."""
    
def plot_summits_for_reps_in_group(rgroupid, con):
    """This is a long method (sorry). 
    It correlates the summits from all the replicates in a group,
    using a variety of methods.
    It assumes that all the replicates come from the same species."""

    cur = con.cursor()
    
    repgroupname = get_repgroup_name(rgroupid, con) 
    repids = get_repids_in_group(rgroupid, con)
    repid_repname = {}
       
    """Get display names for the replicates."""
    for repid in repids:
        sql = "SELECT name from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        repname = cur.fetchone()
        if repname == None:
            repname = repid.__str__()
        else:
            repname = repname[0].__str__()
        repid_repname[repid] = repname

    """If the replicates are in different species, then this is an error.
    The user needs to fix their config file."""
    if False == are_replicates_same_species(repids, con):
        print "\n. Error, the replicates in", repgroupname, "come from different species."
        print ". Please fix your configuration file."
        exit()

    """Write an Excel table with genes and their scores in each replicate."""
    sql = "select geneid from RepgroupSummits where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    
    repid_qvals = {} # key = repid, value = list of Q values for summits
    for repid in repids:
        repid_qvals[repid] = []
    
    """Open the file and write the header row."""
    xlpath = repgroupname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("gene_id\tgene_name\t")
    for repid in repids:
        fout.write("NPeaks(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_Qval(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_dist(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_Qval(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_dist(" + repid_repname[repid].__str__() + ")\t")
    fout.write("\n")

    count = 1
    total_count = genes.__len__() * repids.__len__()

    # RepgroupSummits(repgroupid, geneid, maxsummitid, nearestsummitid, mean_maxsummitscore)
    for geneid in genes:
        genename = get_genename(geneid, con)      
        """Columns 1 and 2 are the geneID and the gene name."""
        fout.write(geneid.__str__() + "\t" + genename.__str__() + "\t") # Gene ID and name
    
        for repid in repids:
            
            count += 1
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            """N Peaks"""
            sql = "select count(id) from Summits where replicate=" + repid.__str__()
            sql += " and id in (select summit from GeneSummits where gene=" + geneid.__str__()
            sql += ")"
            cur.execute(sql)
            peakcount = cur.fetchone()[0]
            fout.write(peakcount.__str__() + "\t")
            
            """Max summit score"""
            sql = "select max(score), id from Summits where replicate=" + repid.__str__()
            sql += " and id in (select summit from GeneSummits where gene=" + geneid.__str__()
            sql += ")"
            cur.execute(sql)
            x = cur.fetchone()
            s = x[0]
            maxsummitid = x[1]
            fout.write(s.__str__() + "\t")
            
            repid_qvals[repid].append( s )
            
            """TSS distance for max summit"""
            sql = "select min(distance) from GeneSummits where summit=" + maxsummitid.__str__() + " and gene=" + geneid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            d = x[0]
            fout.write(d.__str__() + "\t")
        
            
            sql = "select min(distance), summit from GeneSummits where gene=" + geneid.__str__()
            sql += " and summit in (select id from Summits where replicate=" + repid.__str__() + ")"
            cur.execute(sql)
            x = cur.fetchone()
            mindistance = x[0]
            nearestsummtid = x[1]
            
            sql = "select score from Summits where id=" + nearestsummtid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            s = x[0]
            fout.write(s.__str__() + "\t")
            fout.write(mindistance.__str__() + "\t")
        
        """End of the row."""
        fout.write("\n")               
    fout.close()
    
    """Register the Excel file in the database."""
    add_repgroupfile(xlpath, rgroupid, "Excel table with summit stats for replicate group " + repgroupname, con)

    """If there's only one replicate in this group, then we're done."""
    if repids.__len__() < 2:
        return

    """Venn diagram of genes with/without summits in both replicates."""
    venn_data = {}
    for repid in repids:
        venn_data[ repid_repname[repid] ] = []
        x = get_geneids_with_summits(con, repid)
        for ii in x:
            venn_data[ repid_repname[repid] ].append(ii[0]) # NOTE: ii[0] is a gene ID
    vennpath = plot_venn_diagram( venn_data, repgroupname + ".summits")
    add_repgroupfile(vennpath, rgroupid, "Venn diagram comparing genes with summits for replicategroup " + repgroupname, con)
    

    """ Plot score and rank correlations."""
    scatter_values = []
    scatter_names = []
    plot_as_rank = []
    unit_labels = []
    count = 0
    for repid in repids:
        scatter_values.append( repid_qvals[repid] )
        scatter_names.append( "max Qval (" + repid_repname[repid] + ")" )
        unit_labels.append("max Qval")
        count += 1
    for repid in repids:
        plot_as_rank.append(count)
        scatter_values.append( repid_qvals[repid] )
        scatter_names.append( "rank Qval (" + repid_repname[repid] + ")" )
        unit_labels.append("rank Qval")
        count += 1
    width = scatter_values.__len__()
    height = repids.__len__()
    filekeyword = repgroupname + ".summits"
    cranpath = scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, plot_as_rank=plot_as_rank, unit_labels=unit_labels)
    if cranpath != None:
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
    sql = "SELECT rgroupid, geneid, maxenrich, meanenric, meanmaxenrich from GroupEnrichmentStats where"
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
    
            
            sql = "INSERT INTO UnionEnrichmentStats (unionid, geneid, maxenrich, meanenrich)"
            sql += "VALUES (" + unionid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += max_max.__str__() + ","
            sql += mean(all_means).__str__() +  ")"
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
        
        (x_max, x_n, geneids) = get_summit_plot_array_for_replicate_v2(repid, speciesid, con)
        repid_maxs.append( x_max )
        repid_ns.append( x_n)
        repid_name.append( repname )

    #scatter_names = repid_name + repid_name
    scatter_names = []
    unit_labels = []
    for r in repid_name:
        scatter_names.append("max Qval (" + r + ")")
        unit_labels.append("max Qval")
    for r in repid_name:
        scatter_names.append("N summits (" + r + ")" )
        unit_labels.append("count summits")
    scatter_data = repid_maxs + repid_ns
    width = scatter_names.__len__()
    height = repid_name.__len__()
    
    """Normal data: max, ns, and rank scatters."""
    scatter_names = repid_name + repid_name
    scatter_data = repid_maxs + repid_ns
    width = scatter_names.__len__()
    height = repid_name.__len__()
    filekeyword = unionname + ".summits"
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="Summit scores for " + speciesname, xlab="", ylab="", unit_labels=unit_labels)
    if cranpath != None:
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
    
    """This first plotting subfunction prints a maximalist plot comparing
    all the replicates in all replicate groups to each other.
    Later on, further in this function, a more refined plot function
    compares replicates groups rather than replicates."""
    #print "\n. Plotting enrichments for union", unionname
    plot_enrichment_union_helper(unionid, con)
    
    """Build geneid_results,
    where geneid_results[geneid] = SQL results which is a list of GroupEnrichmentStats objects"""
    sql = "SELECT rgroupid, geneid, maxenrich, meanenric, meanmaxenrich from GroupEnrichmentStats where"
    sqlbits = []
    for rgid in rgroupids:
        sqlbits.append( " rgroupid=" + rgid.__str__() )
    sql += " OR ".join( sqlbits )
    cur.execute(sql)
    results = cur.fetchall()
    if results.__len__() == 0:
        return
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
    rgid_maxvals = {} # rgid is a replicate group id
    rgid_meanvals = {}
    rgid_meanmaxvals = {}
    rgid_groupname = {}
    xlpath = unionname + ".enrich.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("gene_ID\tgene_name\t")
    for rgid in rgroupids:
        rgname = get_repgroup_name(rgid, con)
        rgid_groupname[rgid] = rgname
        fout.write("max_FE(" + rgname + ")\t")
        fout.write("mean_FE(" + rgname + ")\t")
        rgid_maxvals[rgid] = []
        rgid_meanvals[rgid] = []
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
            sys.stdout.write("\r    --> (Excel table) %.1f%%" % (100*(count/float(total_count) ) ) )
            sys.stdout.flush()
            foundit = False
            for x in geneid_results[geneid]:
                if x[0] == rgid:
                    foundit = True
                    fout.write(x[2].__str__() + "\t")
                    fout.write(x[3].__str__() + "\t")
    
                    rgid_foundvals[rgid].append(x[2])
                    rgid_maxvals[rgid].append(x[2])
                    rgid_meanvals[rgid].append(x[3])
                    rgid_meanmaxvals[rgid].append(x[4])
            if foundit == False:
                fout.write("0\t0\t")
                rgid_maxvals[rgid].append(0)
                rgid_meanvals[rgid].append(0)
                rgid_meanmaxvals[rgid].append(0)
        fout.write("\n")
    fout.close()
    add_unionfile(xlpath, unionid, "Excel table with enrichment stats for union " + unionname, con)

    
    """Plot a scatterplot using scores distilled from the replicates."""
    scatter_values = []
    scatter_names = []
    unit_labels = []
    for rgid in rgroupids:
        scatter_values.append( rgid_maxvals[rgid] )
        scatter_names.append( "max FE (" + rgid_groupname[rgid] + ")" )
        unit_labels.append("max FE")
    for rgid in rgroupids:
        scatter_values.append( rgid_meanvals[rgid] )
        scatter_names.append(  "mean FE (" + rgid_groupname[rgid] + ")" )
        unit_labels.append("mean FE")
    for rgid in rgroupids:
        scatter_values.append( rgid_meanmaxvals[rgid] )
        scatter_names.append( "mean of max FE (" + rgid_groupname[rgid] + ")" )
        unit_labels.append("mean-of-max FE")

    filekeyword = unionname + ".enrich.small"
    plot_as_rank = [] # a list of indices in scatter_data that should be plotted as rank, rather than raw.
    start = 3 * rgroupids.__len__()
    for ii in range(start, start + rgroupids.__len__() ):
        plot_as_rank.append( ii )
    width = scatter_values.__len__()
    height = rgroupids.__len__()
    cranpath = scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, plot_as_rank=plot_as_rank, unit_labels=unit_labels)
    if cranpath != None:
        add_unionfile(cranpath, unionid,"R script path to a 2-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)
        add_unionfile(re.sub("cran", "pdf", cranpath), unionid,"PDF 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)   

    return
#
# continue here
# to-do
#

#     """Also plot IDR."""
#     scatter_values = []
#     scatter_names = []
#     for rgid in rgroupids:
#         scatter_values.append( rgid_maxvals[rgid] )
#         scatter_names.append( rgid_groupname[rgid] + "-max")
#     for rgid in rgroupids:
#         scatter_values.append( rgid_meanvals[rgid] )
#         scatter_names.append( rgid_groupname[rgid] + "-mean")
#     width = scatter_values.__len__()
#     filekeyword = unionname + ".enrich.idr"
#     """idr_stats[gene number in the geneid list][rep ii][rep jj] = lidr"""
#     (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_values, scatter_names, filekeyword, title="", xlab="", ylab="", force_square=True)
#     if cranpath != None:
#         add_unionfile(cranpath, unionid,"" + unionname, con)
#         add_unionfile(re.sub("cran", "pdf", cranpath), unionid,"" + unionname, con) 
#     
#     
#     """Now insert the IDR stats into the SQL database."""
#     
#     """First, clear any previous entries for the pair of this replicate group."""
#     for ii in range(0, rgroupids.__len__() ):
#         for jj in range(0, rgroupids.__len__() ):
#             sql = "DELETE from GeneUnionEnrichIdr where repgroupid1=" + rgroupids[ii].__str__() + " and repgroupid2=" + rgroupids[jj].__str__()
#             cur.execute(sql)
#             con.commit()
#     
#     """Next insert the data into SQL."""
#     geneids = geneid_results.keys()
#     for gg in range(0, geneid_results.__len__() ):
#         geneid = geneids[gg]
#         if gg in idr_stats:
#             for ii in range(0, rgroupids.__len__() ):
#                 if ii in idr_stats[gg]:
#                     for jj in range(0, rgroupids.__len__() ):
#                         if jj in idr_stats[gg][ii]:
#                             this_idr = idr_stats[gg][ii][jj]
#                             ii_repgroupid = rgroupids[ii]
#                             jj_repgroupid = rgroupids[jj]
#                             sql = "INSERT into GeneUnionEnrichIdr(geneid,repgroupid1,repgroupid2,lidr)"
#                             sql += " VALUES("
#                             sql += geneid.__str__() + ","
#                             sql += ii_repgroupid.__str__() + ","
#                             sql += jj_repgroupid.__str__() + ","
#                             sql += this_idr.__str__() + ")"
#                             #print sql
#                             cur.execute(sql)
#                             
#     con.commit()

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
        
        x = get_enrichment_plot_array_for_replicate(repid, speciesid, con)
        if x == None:
            return
        (x_maxe, x_meane, geneids) = x
        repid_means.append( x_meane )
        repid_maxs.append( x_maxe)
        repid_names.append( repname )
            
    scatter_names = []
    for repname in repid_names:
        scatter_names.append("max FE (" + repname + ")")
    for repname in repid_names:
        scatter_names.append("mean FE (" + repname + ")")

    #scatter_names = repid_name + repid_name
    scatter_data = repid_maxs + repid_means
    width = scatter_names.__len__()
    height = repid_names.__len__()
    filekeyword = unionname + ".enrich.large"
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="Fold Enrichment for " + unionname, xlab="fold-enrichment", ylab="fold-enrichment")
    if cranpath != None:
        add_unionfile(cranpath, unionid, "R script for a scatterplot comparing fold enrichment for all replicates in the union " + unionname, con)
        add_unionfile(re.sub("cran", "pdf", cranpath), unionid, "PDF comparing fold enrichment for all replicates in the union " + unionname, con)

    """We don't process the IDR data because that already occurred in the plot_enrichment_for_reps_in_group."""
    

def compute_enrichments_for_reps_in_group(rgroupid, con):
    """Fills the table GroupEnrichmentStats with fold-enrichment stats
    for the genes that have FE data all the replicates within the replicate group
    whose ID is = rgroupid"""
    
    cur = con.cursor()
    sql = "DELETE from GroupEnrichmentStats where rgroupid=" + rgroupid.__str__()
    cur.execute( sql )
    con.commit()
    
    repids = get_repids_in_group(rgroupid, con)
    
    """Get replicate names."""
    repid_repname = {}
    for repid in repids:
        cur.execute("SELECT name from Replicates where id=" + repid.__str__())
        x = cur.fetchone()
        if x == None:
            repid_repname[repid] = repid.__str__()
        else:
            repid_repname[repid] = x[0]
    
    
    """Get FE data from the SQL database."""
    geneids = []
    repid_festats = {}
    for repid in repids:
        repid_festats[repid] = get_enrichment_stats_for_replicate( repid, con )
        """Record the seen gene IDs."""
        for geneid in repid_festats[repid]:
            if geneid not in geneids:
                geneids.append( geneid )
    
#     """Sanity Check: repid_festats[x] must be the same length for all x."""
#     for ii in range(0, repids.__len__()-2):
#         this_repid = repids[ii]
#         that_repid = repids[ii+1]
#         if repid_festats[this_repid].__len__() != repid_festats[that_repid].__len__():
#             print "\n. Error, chipseqdb_plot.py 1136, festats not the same length for", repid_repname[this_repid], "and", repid_repname[that_repid]
#             print "\.", repid_festats[this_repid].__len__().__str__() + " vs. " + repid_festats[that_repid].__len__().__str__()
#             exit()
    
    geneids_inall = []
    maxfe = []
    meanfe = []
    meanmaxfe = []
    sumfe = []

    for geneid in geneids:
        found_in_all = True
        for repid in repids:
            if geneid not in repid_festats[repid]:
                found_in_all = False
        if found_in_all:
            geneids_inall.append(geneid)
            this_max = 0
            all_means = [] 
            all_maxes = []
            for repid in repids:
                if geneid in repid_festats[repid]:
                    rmax = repid_festats[repid][geneid][0]
                    rmean = repid_festats[repid][geneid][1]
                    all_maxes.append( rmax )
                    all_means.append( rmean )
                if rmax > this_max:
                    this_max = rmax
            maxfe.append( this_max )
            meanfe.append( mean( all_means ) )
            meanmaxfe.append( mean(all_maxes) )
        else:
            pass
            
    for ii in range(0, geneids_inall.__len__() ):
        geneid = geneids_inall[ii]
        sql = "INSERT into GroupEnrichmentStats (rgroupid, geneid, maxenrich, meanenrich, meanmaxenrich)"
        sql += " VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ","
        sql += maxfe[ii].__str__() + ","
        sql += meanfe[ii].__str__() + ","
        sql += meanmaxfe[ii].__str__() + ")"
        cur.execute(sql)
    con.commit()


def get_summit_plot_array_for_replicate_v2(repid, species, con):
    """Returns a tuple with three lists: (qvals, x_n, geneids) """
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
        if count %100 == 0:
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
    
    return (qvals, x_n, geneids)


def sort_gene_values():
    """Get the rank of genes according to an arbitrary score."""
    pass

def get_enrichment_plot_array_for_replicate(repid, species, con):
    """Returns ( x_maxe[], x_meane[] ) for this replicate.
    There is an entry in these arrays for every gene in the species' genome.
    If something is wrong, then it returns None.
    """
    x_maxe = []
    x_meane = []
    x = get_enrichment_stats_for_replicate( repid, con )
    if x.__len__() == 0:
        return None
    
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
            
            if x[geneid][1] not in meane_geneids:
                meane_geneids[  x[geneid][1] ] = []
            meane_geneids[  x[geneid][1] ].append( geneid )
        else:
            x_maxe.append(0)
            x_meane.append(0)
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
        
    if x_maxe.__len__() == 0:
        print "\n. 1306"
        exit()
    if x_meane.__len__() == 0:
        print "\n. 1309"
        exit()
    return (x_maxe, x_meane, geneids)

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
    
    cur = con.cursor()
        
    if repgroupname == None:
        repgroupname = get_repgroup_name(rgroupid, con)
    
    if repids == None:
        repids = get_repids_in_group(rgroupid, con)
    
    """If there's only one replicate in this group, then we're done."""
    if repids.__len__() < 2:
        return
        
    """Get the names of the reps in the group."""
    repid_repname = {}
    for repid in repids:
        cur.execute("SELECT name from Replicates where id=" + repid.__str__())
        x = cur.fetchone()
        if x != None:
            repid_repname[repid] = x[0]
        else:
            repid_repname[repid] = repid.__str__()
    
    """Get the names of genes in the Repgroup"""
    sql = "select geneid from GroupEnrichmentStats where rgroupid=" + rgroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    geneids = []
    for ii in x:
        geneids.append( ii[0] )
        
    """This is the stuff we want to know for each gene. . . """
    repid_maxfe = {} # key = repid, value = list of max FEs, one for each gene
    repid_meanfe = {}
    repid_npeaks = {}
    repid_maxqval = {}
    repid_dist_maxsummit = {}
    for repid in repids:
        repid_maxfe[repid] = []
        repid_meanfe[repid] = []
        repid_maxqval[repid] = []
        repid_dist_maxsummit[repid] = []
        
    """Now iterate through all genes for all replicates (in repids), and get the data to
        fill the lists that were created above."""
    count = 0
    total_count = geneids.__len__() * repids.__len__()
    for geneid in geneids:
        for repid in repids:
            
            count += 1
            if count%10 == 0:
                sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
            
            sql = "select maxenrich, meanenrich from EnrichmentStats where geneid=" + geneid.__str__() + " and repid=" + repid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            repid_maxfe[repid].append( x[0] )
            repid_meanfe[repid].append( x[1] )
        
            count += 1
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            """N Peaks"""
            sql = "select count(id) from Summits where replicate=" + repid.__str__()
            sql += " and id in (select summit from GeneSummits where gene=" + geneid.__str__()
            sql += ")"
            cur.execute(sql)
            peakcount = cur.fetchone()[0]
            repid_npeaks[repid].append( peakcount )
            
            if peakcount == 0:
                """This gene in this replicate has no peaks."""
                repid_maxqval[repid].append(None)
                repid_dist_maxsummit[repid].append(None)
            elif peakcount > 0:
                """This gene in this replicate has at least one peak."""
                
                """Max summit score"""
                sql = "select max(score), id from Summits where replicate=" + repid.__str__()
                sql += " and id in (select summit from GeneSummits where gene=" + geneid.__str__()
                sql += ")"
                cur.execute(sql)
                x = cur.fetchone()
                s = x[0]
                maxsummitid = x[1]
                repid_maxqval[repid].append( s )
                
                """TSS distance for max summit"""            
                sql = "select min(distance), summit from GeneSummits where gene=" + geneid.__str__()
                sql += " and summit in (select id from Summits where replicate=" + repid.__str__() + ")"
                cur.execute(sql)
                x = cur.fetchone()
                mindistance = x[0]
                repid_dist_maxsummit[repid].append(  mindistance )
                #nearestsummtid = x[1]
            
    """Plot FE versus summits"""   
    for repid in repids:            
        repname = repid_repname[repid]
        filekeyword = rep
        scatter_data = [ repid_maxfe[repid], repid_maxqval[repid] ]
        scatter_names = [ "max FE (" + repname + ")", "max Qval (" + repname + ")" ]
        plot_as_rank = []
        width = 2
        height = 2
        filekeyword = repname + ".fe_x_qval"
        cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="", force_square=False)            

    """Plot FE stuff"""
    scatter_data = []
    scatter_names = []
    for repid in repids:
        scatter_data.append( repid_maxfe[repid] )
        scatter_names.append( "max FE (" + repid_repname[repid] + ")" )
    for repid in repids:
        scatter_data.append( repid_meanfe[repid] )
        scatter_names.append("mean FE (" + repid_repname[repid] + ")")
    width = scatter_data.__len__()
    height = repids.__len__()
    filekeyword = repgroupname + ".enrich"
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword)
    if cranpath != None:
        add_repgroupfile(cranpath, rgroupid, "R script for multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
        add_repgroupfile(re.sub("cran", "pdf", cranpath), rgroupid, "PDF multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
    
    """Also plot the IDR stats"""
    filekeyword = repgroupname + ".enrich.idr"
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_data, scatter_names, filekeyword, title="", xlab="", ylab="")
    if cranpath != None:
        add_repgroupfile(cranpath,rgroupid,"R script to make IDR scatterplots for replicate group " + repgroupname, con)
        add_repgroupfile(re.sub("cran", "pdf",cranpath),rgroupid,"PDF with IDR scatterplots for replicate group " + repgroupname, con)
    
    """Update the IDR stats into the database"""
    for (ii,jj) in value_pairs: 
        sql = "DELETE from GeneRepgroupEnrichIdr where repid1=" + repids[ii%repids.__len__()].__str__() + " and repid2=" + repids[jj%repids.__len__()].__str__()
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
                            cur.execute(sql)
    con.commit()

    """Write the Excel Table"""
    xlpath = repgroupname + ".enrich.xls"
    print "\n. Writing a table with fold-enrichment to", xlpath
    fout = open(xlpath, "w")
    
    """Header"""
    fout.write("gene_ID\tgene_name\t")
    for repid in repids:
        fout.write("max_FE(" + repid_repname[repid] + ")\t")
        fout.write("mean_FE(" + repid_repname[repid] + ")\t")

        """ Not all genes will have summit information... """
        fout.write("NPeaks(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_Qval(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_dist(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_Qval(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_dist(" + repid_repname[repid].__str__() + ")\t")
    fout.write("\n")
    
    """One row per gene"""
    for gg in range(0, geneids.__len__() ):
        geneid = geneids[gg]
        genename = get_genename(geneids[gg], con)
        if genename == None:
            continue
        fout.write(geneids[gg].__str__() + "\t" + genename + "\t" )
        for ii in range(0, repids.__len__() ):
            repid = repids[ii]
            fout.write(repid_maxfe[repid][gg].__str__() + "\t")
            fout.write(repid_meanfe[repid][gg].__str__() + "\t")
            fout.write(repid_npeaks[repid][gg].__str__() + "\t")
            fout.write(repid_maxqval[repid][gg].__str__() + "\t")
            fout.write(repid_dist_maxsummit[repid][gg].__str__() + "\t")
            
        fout.write("\n")
    fout.close()

    """Add an entry to the SQL database, listing the Excel file."""
    add_repgroupfile(xlpath,rgroupid,"Excel table with enrichment stats for replicate group " + repgroupname, con)

# depricated
# this is now folded into the plot enrichments method
# def plot_summits_vs_enrichments_for_replicates(repids, con):
#     cur = con.cursor()
#     
#     for repid in repids:
#         print "\n. Plotting summits vs. enrichments for replicate", repid
#         
#         cur.execute("SELECT name from Replicates where id=" + repid.__str__())
#         x = cur.fetchone()
#         repname = None
#         if x != None:
#             repname = x[0].__str__()
#         else:
#             print "\n. chipseqdb_plot.py 1477: Error, I can't find a name for the replicate", repid
#             exit()
#         
#         fearray = []
#         sarray = []
#         
#         """Get the genes with summits."""
#         #sql = "select gene from GeneSummits where summit in (select id from Summits where replicate=" + repid.__str__() + ")"
#         sql = "select id from Summits where replicate=" + repid.__str__()
#         cur.execute(sql)
#         x = cur.fetchall()
#         
#         sql = "select count(id) from Summits where replicate=" + repid.__str__()
#         cur.execute(sql)
#         count_total = cur.fetchone()[0]
#         
#         summits = []
#         count = 0
#         for ii in x:
#             """Progress indicator."""
#             count += 1
#             sys.stdout.write("\r    --> " + count.__str__() + " of " + count_total.__str__() + " summits." )
#             sys.stdout.flush()
#             
#             summitid = ii[0]
#             
# 
#             """Get the max FE at this gene."""
#             sql = "select max_enrichment from SummitsEnrichment where summit=" + summitid.__str__()
#             cur.execute(sql)
#             fe = cur.fetchone()[0]
#             
#             """Get the max summit score"""
#             sql = "select score from Summits where id=" + summitid.__str__()
#             cur.execute(sql)
#             score = cur.fetchone()[0]
#               
#             #"""Get the distance to the maximum summit."""
#             #sql = "select min(distance) from GeneSummits where summit=" + summitid.__str__()
#             #cur.execute(sql)
#             #d = cur.fetchone()[0]
#             
#             if fe != None and score != None:# and d != None:       
#                 summits.append(summitid)
#                 fearray.append( fe )
#                 sarray.append( score )
#                 #sdarray.append( d )
#             #else:
#             #    print "\n. This summit doesn't seem to map to a gene:", summitid, fe, score, d
#         
#         """Plot the data."""                
#         scatter_values = [fearray,sarray]
#         scatter_names = [ "max FE (" + repname + ")", "max Qval (" + repname + ")" ]
#         plot_as_rank = []
#         width = 2
#         height = 2
#         filekeyword = repname + ".enrich_x_summit"
#         cranpath = scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, title="", force_square=False)
#     
#         #
#         # to=do: add file to DB here
#         # continue here
#         #
#     
#         """Write an excel table."""
#         xlpath = repname + ".enrich_x_summit.xls"
#         fout = open(xlpath, "w")
#         fout.write("Summit_ID\tChrom.\tSite\tSummit_Q_Score\tFold-Enrich._at_Peak\n")
#         for ii in range( 0, summits.__len__() ): 
#             sql = "select chrom, site from Summits where id=" + summits[ii].__str__()
#             cur.execute(sql)
#             res = cur.fetchone()
#             chrom = res[0]
#             site = res[1]
#             fout.write(summits[ii].__str__() + "\t" + get_chrom_name(con, chrom).__str__() + "\t" + site.__str__() + "\t" + fearray[ii].__str__() + "\t" + sarray[ii].__str__() + "\n")
#         fout.close()

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
    
    
    sql = "SELECT unionid, geneid, maxenrich, meanenrich from UnionEnrichmentStats where"
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
            all_maxes = []
            all_means = []
                   
            for ii in geneid_results[gid]:
                unionid = ii[0]
               
                count += 1
                #if count%50 == 0:
                sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
             
                all_maxes.append( ii[2] )
                all_means.append( ii[3] )

            sql = "INSERT INTO SpeciesunionEnrichmentStats (unionid, geneid, maxenrich, meanenrich)"
            sql += "VALUES (" + uid.__str__() + ","
            sql += gid.__str__() + ","
            sql += max(all_maxes).__str__() + ","
            sql += mean(all_means).__str__() + ")"
            cur.execute(sql)
    con.commit()    

def plot_enrichments_for_speciesunion(uid, con):
    cur = con.cursor()
        
    spunionname = get_speciesunionname(uid, con)
    unionids = get_unionids_in_speciesunion( uid, con )
    
    print "\n. Plotting enrichments acorss species in the set '" + spunionname + "'"
    sys.stdout.flush()

    """
    What gene IDs apply to this species union?
    In the following block, we lookup the genes that have enrichment
    scores in their nearby regulatory regions in all the replicates
    represented by at least one union. We then translate
    those gene IDs, using the table GeneHomology with data from
    the pillars file, in order to generate a list of genes
    which can be compared across species."""
    
    sql = "SELECT unionid, geneid, maxenrich, meanenrich from SpeciesunionEnrichmentStats where unionid=" + uid.__str__()
    cur.execute(sql)
    results = cur.fetchall()
    
    #print "1766:", uid, results.__len__()
    
    """This will hold SQL data"""
    gene_aliases = {}
    geneid_results = {} # geneid_results[translated_geneid] = array of rows from UnionEnrichmentStats


    """Build a dictionary of aliased genes. This will be useful when we assign homology between
    enrichment values."""
    count = 1
    total_count = results.__len__()
    gene_aliases = {}
    for ii in results:       
        
        translated_id = ii[1]
        if translated_id not in gene_aliases:
            gene_aliases[translated_id] = []
        gene_aliases[translated_id] = get_aliasids_for_geneid(translated_id, con)
        #print "1820:", translated_id, gene_aliases[translated_id]
        
        """Remember this gene."""
        if translated_id not in geneid_results:
            geneid_results[ translated_id ] = []
            
        """Store the translation pairing"""
        geneid_results[translated_id].append(ii)
    
        count += 1
        if count %100:
            sys.stdout.write("\r    --> (translation) %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
        
    
    #print "\n. 1797:", geneid_results.__len__()
    
    """Write an excel table"""
    xlpath = spunionname + ".enrich.xls"
    print "\n. Writing a table to", xlpath
    
    fout = open(xlpath, "w")
    #fout.write("GeneID\tGeneName\t")
    for unionid in unionids:
        unionname = get_unionname(unionid, con)
        fout.write("gene_ID(" + unionname + ")\t")
        fout.write("gene_name(" + unionname + ")\t")
        fout.write("max_FE(" + unionname + ")\t")
        fout.write("max_FE(" + unionname + ")\t")
        
    count = 1
    total_count = 0
    for tgid in geneid_results.keys():
        for unionid in unionids:
            total_count += geneid_results[tgid].__len__()
    
    fout.write("\n")
    
    union_geneid_results = {}
    for unionid in unionids:
        union_geneid_results[unionid] = {}
        sql = "SELECT unionid, geneid, maxenrich, meanenrich  from UnionEnrichmentStats where unionid=" + unionid.__str__()
        x = cur.execute(sql)
        for ii in x:
            union_geneid_results[unionid][ ii[1] ] = ii
    
    #for unionid in unionids:
    #    print "1867", unionid, union_geneid_results[unionid].__len__()
    #exit()
    
    """One row for each translated gene_ID"""
    for tgid in geneid_results.keys():  
        for unionid in unionids: # find the unionid that mathches result x
            count += 1
            if count %10:
                sys.stdout.write("\r    --> (Excel table) %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
            
            foundit = False
            for alias in gene_aliases[tgid]:   
                if foundit == False:
                    if alias in union_geneid_results[unionid]:
                        x = union_geneid_results[unionid][alias]
                    
                    #sql = "SELECT unionid, geneid, maxenrich, meanenrich  from UnionEnrichmentStats where unionid=" + unionid.__str__()
                    #sql += " and geneid=" + alias.__str__()
                    #cur.execute(sql)
                    #x = cur.fetchone()
                    #if x != None:
                        """Did we find the SQL query that corresponds to this unionid?"""
                        foundit = True
                        fout.write(x[1].__str__() + "\t")
                        fout.write(get_genename(x[1], con) + "\t")
                        fout.write(x[2].__str__() + "\t")
                        fout.write(x[3].__str__() + "\t")
            if foundit == False:
                fout.write("---\t---\t0\t0\t")
        fout.write("\n")
    fout.close()

    add_speciesunionfile(xlpath, uid, "Excel table with enrichment stats for species union " + spunionname, con )


        
    scatterdata_part1 = []
    scatterdata_part2 = []

    for unionid in unionids:
        maxvals = []
        meanvals = []
        
        """For each gene that has FE data in all unions"""
        for translated_id in geneid_results:
            foundit = False
            for alias in gene_aliases[translated_id]:   
                if foundit == False:      
                    if alias in union_geneid_results[unionid]:
                        maxvals.append(union_geneid_results[unionid][alias][2] )
                        meanvals.append(union_geneid_results[unionid][alias][3] )
                        foundit = True
                        #print "1920:", translated_id, alias, unionid, union_geneid_results[unionid][alias][2], union_geneid_results[unionid][alias][3]
            if foundit == False:
                print "\n. Error 1922"
                exit()
        scatterdata_part1.append( maxvals )
        scatterdata_part2.append( meanvals )
    scatterdata = scatterdata_part1 + scatterdata_part2



    """Names for each union"""
    scatter_names = []
    for unionid in unionids:    
        scatter_names.append( "max FE (" + get_unionname(unionid, con) + ")" )
    for unionid in unionids:    
        scatter_names.append( "mean FE (" + get_unionname(unionid, con) + ")" )    
    width = scatter_names.__len__()
    height = unionids.__len__()
    filekeyword = spunionname + ".enrich"
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword, title="Fold Enrichment " + spunionname, xlab="fold-enrichment", ylab="fold-enrichment", skip_zeros=False) 
    if cranpath != None:
        add_speciesunionfile(cranpath, uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )
        add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )    
                            
    con.commit()

    