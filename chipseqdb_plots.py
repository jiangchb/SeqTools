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
    sql = "DELETE from UnionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
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
        
    """Get data from RepgroupSummits"""
    geneid_maxsummits = {} # key = translated gene ID, value = list of SQL results from RepgroupSummitStats
    geneid_nearestsummits = {}
    gene_repgroups = {} # key = geneid, value = list of repgroup IDs with summit values for this gene.
    for repgroupid in rgroupids:
        repgroupname = get_repgroup_name( repgroupid, con )
        
        sql = "SELECT * from RepgroupSummits where repgroupid=" + repgroupid.__str__()
        cur.execute(sql)
        results = cur.fetchall()

        for ii in results:
            geneid = ii[1]
            max_summit_id = ii[2]        
            nearest_summit_id = ii[3]
                        
            if geneid not in geneid_maxsummits:
                geneid_maxsummits[geneid] = []
            geneid_maxsummits[geneid].append( max_summit_id )
            
            if geneid not in geneid_nearestsummits:
                geneid_nearestsummits[geneid] = []
            geneid_nearestsummits[geneid].append( nearest_summit_id )

            if geneid not in gene_repgroups:
                gene_repgroups[geneid] = []
            if repgroupid not in gene_repgroups[geneid]:
                gene_repgroups[geneid].append( repgroupid )

    """Use the data from RepgroupSummits to find the max summit, and the nearest summit,
    for every gene."""
    geneid_max = {}
    geneid_nearest = {}
    for geneid in gene_repgroups:
        if gene_repgroups[geneid].__len__() == rgroupids.__len__():
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
                print "\n. Error, there is no max summitid for gene", geneid
                print "\n. chipseqdb_plot.py 86"
                exit()
            geneid_max[geneid] = max_max_summit_id
            
        if geneid_nearestsummits[geneid].__len__() == rgroupids.__len__():
            nearest_nearest_summit_id = None
            nearest_nearest_summit_distance = None
            for summitid in geneid_nearestsummits[geneid]:
                sql = "SELECT distance from GeneSummits where gene=" + geneid.__str__() + " and summit=" + summitid.__str__()
                cur.execute(sql)
                x = cur.fetchone()
                if x != None:  
                    this_distance = abs( x[0] )
                    if nearest_nearest_summit_distance == None:
                        nearest_nearest_summit_distance = this_distance
                        nearest_nearest_summit_id = summitid
                    elif nearest_nearest_summit_distance > this_distance:
                        nearest_nearest_summit_distance = this_distance
                        nearest_nearest_summit_id = summitid
                else:
                    print "Error: 106", summitid, geneid, repgroupid                    
                    exit()
                    
            if nearest_nearest_summit_id == None:
                print "\n. Error, there is no nearest summitid for gene", geneid
                print "\n. chipseqdb_plot.py 103"
                exit()
            geneid_nearest[geneid] = nearest_nearest_summit_id
    
    """Insert the max and nearest summit stats into UnionSummits."""
    for geneid in geneid_max:    
        if geneid in geneid_nearest:
            """UnionGenes is depricated, but we'll update it anyway."""            
            sql = "INSERT INTO UnionGenes (unionid, geneid) VALUES(" + unionid.__str__() + "," + geneid.__str__() + ")"
            cur.execute(sql)
            
            sql = "INSERT INTO UnionSummits (unionid, geneid, maxsummitid, nearestsummitid)"
            sql += " VALUES(" + unionid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += geneid_max[geneid].__str__() + ","
            sql += geneid_nearest[geneid].__str__() + ")"
            #print sql
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
    vennpath = plot_venn_diagram( venn_data, unionname.__str__() + ".summits" )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for union " + unionname, con)

        
    """
    1. Update the table UnionGenes with genes that have a peak in all repgroups in this union.
    2. Write the Excel table listing the genes and their scores in the repgroups."""
    xlpath = unionname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("geneID\tname\t")
    for rgroupid in rgroupids:
        rgroupname = get_repgroup_name( rgroupid, con )        
        fout.write("MaxSummit(" + rgroupname.__str__() + ")\t")
        fout.write("NearestSummit(" + rgroupname.__str__() + ")\t")
    fout.write("\n")
    
    count = 1
    total_count = seen_genes.__len__() * rgroupids.__len__()
    
    for geneid in seen_genes:        
        fout.write(geneid.__str__() + "\t" + get_genename(geneid, con) + "\t" )
        for rgroupid in rgroupids:
            
            count += 1
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            sql = "SELECT * from RepgroupSummits where repgroupid=" + rgroupid.__str__()
            sql += " and geneid=" + geneid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            if x != None:
                sql = "SELECT score from Summits where id=" + x[2].__str__()
                cur.execute(sql)
                this_score = cur.fetchone()[0]
                fout.write(this_score.__str__() + "\t")
                sql = "SELECT score from Summits where id=" + x[3].__str__()
                cur.execute(sql)
                this_score = cur.fetchone()[0]
                fout.write(this_score.__str__() + "\t")
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
                print "\n. Error, there is no max summitid for gene", geneid, "in repgroup", repgroupname
                print "\n. chipseqdb_plot.py 86"
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
            """SpeciesunionGenes is depricated, but we'll update it anyway."""            
            sql = "INSERT INTO SpeciesunionGenes (unionid, geneid) VALUES(" + uid.__str__() + "," + geneid.__str__() + ")"
            cur.execute(sql)
            
            sql = "INSERT INTO SpeciesunionSummits (spunionid, geneid, maxsummitid, nearestsummitid)"
            sql += " VALUES(" + uid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += geneid_max[geneid].__str__() + ","
            sql += geneid_nearest[geneid].__str__() + ")"
            #print sql
            cur.execute(sql)
    con.commit()
 
    sql = "SELECT COUNT(*) from SpeciesunionGenes where unionid=" + unionid.__str__()
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
    vennpath = plot_venn_diagram( venn_data, spunionname + ".summits" )
    add_unionfile(vennpath, unionid, "Venn diagram comparing genes with summits for species-union " + spunionname, con)


    unionid_maxsummits = {}
    unionid_nearestsummits = {}
            
    #
    # Write an Excel table using data from UnionSummitStats
    #
    xlpath = spunionname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    #fout.write("geneID\tname\t")
    for unionid in unionids:
        unionname = get_unionname( unionid, con )  
        fout.write("geneID(" + unionname.__str__() + ")\t")
        fout.write("name(" + unionname.__str__() + ")\t")      
        fout.write("MaxSummit(" + unionname.__str__() + ")\t")
        fout.write("NearestSummit(" + unionname.__str__() + ")\t")
        #fout.write("rank(" + unionname.__str__() + ")\t")
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
                unionid_nearestsummits[unionid] = []
            
            """Find the appropriate gene alias for this union."""            
            foundit = False
            for geneid in gene_aliases[gid]:    
                if foundit == False:                    
                    sql = "SELECT * from UnionSummits where unionid=" + unionid.__str__()
                    sql += " and geneid=" + geneid.__str__()
                    cur.execute(sql)
                    x = cur.fetchone()
                    if x != None:
                        foundit = True
                        fout.write(geneid.__str__() + "\t")
                        fout.write(get_genename(geneid, con) + "\t")
                        sql = "SELECT score from Summits where id=" + x[2].__str__()
                        cur.execute(sql)
                        this_score = cur.fetchone()[0]
                        fout.write(this_score.__str__() + "\t")
                        sql = "SELECT score from Summits where id=" + x[3].__str__()
                        cur.execute(sql)
                        this_score = cur.fetchone()[0]
                        fout.write(this_score.__str__() + "\t")
                        
                        unionid_maxsummits[unionid].append( x[2] )
                        unionid_nearestsummits[unionid].append( x[3] )
                        
            if foundit == False:
                fout.write("---\t---\t0\t0\t")
                unionid_maxsummits[unionid].append( 0 )
                unionid_nearestsummits[unionid].append( 0 )
        
        fout.write("\n")
    fout.close()
    
    add_speciesunionfile( xlpath, uid, "Excel table with summit stats for speciesunion " + spunionname, con )

    unionnames = []
    for unionid in unionids:
        unionnames.append(  get_unionname( unionid, con )  )
    

    maxes = []
    nearests = []
    ranks = []
    for unionid in unionids:
        maxes.append( unionid_maxsummits[unionid] )
        nearests.append( unionid_nearestsummits[unionid] )
    
    scatterdata = maxes[0:nunionids] + nearests[0:nunionids]
    scatter_names = unionnames + unionnames
    #cranpath = scatter9x3(scatterdata, scatter_names, "summits.9x3." + spunionname, title="Summits " + spunionname)
    width = scatter_names.__len__()
    height = unionnames.__len__()
    filekeyword = spunionname + ".summits"
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword, title="Summits " + spunionname)
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
    sql = "DELETE from RepgroupGenes where repgroupid=" + rgroupid.__str__()
    cur.execute(sql)
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
         
    """gene_repid_maxsummitid and gene_repid_nearestsummitid are intermediate hashtables 
    that we'll use to fill the table RepgroupSummits later in this method."""  
    gene_repid_maxsummitid = {} # key = geneid, value = hash; key = repid, value = (summit ID, summit score)
    gene_repid_nearestsummitid = {} # key = geneid, value = hash; key = repid, value = (summit ID, summit distance)
    seen_genes = [] # seen_genes is a list of genes that have summits in one or more replicates.
    gene_repcount = {} # key = geneid, value = count of replicates with a summit at that gene.
    for repid in repids:
        ssstats = get_summit_scores_for_replicate(repid, con)
        
        for geneid in ssstats:
            """Initialize our hashes outside this loop."""
            if geneid not in seen_genes:
                seen_genes.append( geneid )
            if geneid not in gene_repid_maxsummitid:
                gene_repid_maxsummitid[geneid] = {}
                gene_repid_nearestsummitid[geneid] = {}
            if repid not in gene_repid_maxsummitid[geneid]:
                gene_repid_maxsummitid[geneid][repid] = (None, 0, 0)
                gene_repid_nearestsummitid[geneid][repid] = (None, 0, 0)
                        
            for summitid in ssstats[geneid]:
                this_summit_score = ssstats[geneid][summitid][0]
                this_summit_distance = abs( ssstats[geneid][summitid][2] )
                
                """Does this summit have a higher score than other summits for geneid?"""
                if gene_repid_maxsummitid[geneid][repid][0] == None:
                    gene_repid_maxsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
                if gene_repid_maxsummitid[geneid][repid][1] < this_summit_score:
                    gene_repid_maxsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
            
                """Is this summit closer to the TSS of geneid?"""
                if gene_repid_nearestsummitid[geneid][repid][0] == None:
                    gene_repid_nearestsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
                if gene_repid_nearestsummitid[geneid][repid][2] > this_summit_distance:
                    gene_repid_nearestsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance) 
    
    for geneid in seen_genes:
        """Is a summit at this gene found in all replicates?
        If so, there should be at least one entry for every replicate."""
        if gene_repid_maxsummitid[geneid].__len__() == repids.__len__():

            """ This summit is in all replicates."""
            sql = "INSERT into RepgroupGenes (repgroupid, geneid) VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ")"
            cur.execute(sql)
        
            """Find the maximum summit for this gene."""
            max_summit_id = None # summit ID of the max summit for geneid
            max_summit_score = 0 # summit score for the maximum summit
            for repid in repids:
                this_summit_id = gene_repid_maxsummitid[geneid][repid][0]
                this_summit_score = gene_repid_maxsummitid[geneid][repid][1]
                if max_summit_score < this_summit_score:
                    max_summit_score = this_summit_score
                    max_summit_id = this_summit_id
                
            """Find the nearest summit for the gene"""
            nearest_summit_id = None # summit ID of the nearest summit to geneid
            nearest_summit_dist = None # distance from TSS for ths nearest summit
            for repid in repids:
                this_summit_id = gene_repid_nearestsummitid[geneid][repid][0]
                this_summit_distance =  gene_repid_nearestsummitid[geneid][repid][2]  
                if nearest_summit_dist == None:
                    neareest_summit_dist = this_summit_distance
                    nearest_summit_id = this_summit_id
                elif nearest_summit_dist > this_summit_distance:
                    neareest_summit_dist = this_summit_distance
                    nearest_summit_id = this_summit_id

            """How many summits are shared among replicates?"""
        
            sql = "INSERT INTO RepgroupSummits(repgroupid, geneid, maxsummitid, nearestsummitid) "
            sql += " VALUES("
            sql += rgroupid.__str__() + ","
            sql += geneid.__str__() + ","
            sql += max_summit_id.__str__() + ","
            sql += nearest_summit_id.__str__()
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

    seen_genes = []
    gene_repid_npeaks = {}
    gene_repid_maxsummitid = {}
    gene_repid_nearestsummitid = {}
    for repid in repids:
        """sssstats are the summit stats for all genes in replicate repid."""
        ssstats = get_summit_scores_for_replicate(repid, con)
        
        for geneid in ssstats:
            """Initialize our hashes outside this loop."""
            if geneid not in seen_genes:
                seen_genes.append( geneid )
            if geneid not in gene_repid_maxsummitid:    
                gene_repid_maxsummitid[geneid] = {}
                gene_repid_nearestsummitid[geneid] = {}
            if repid not in gene_repid_maxsummitid[geneid]:    
                gene_repid_maxsummitid[geneid][repid] = (None, 0, 0)
                gene_repid_nearestsummitid[geneid][repid] = (None, 0, 0)
             
            if geneid not in gene_repid_npeaks:
                gene_repid_npeaks[geneid] = {}
            if repid not in gene_repid_npeaks[geneid]:
                gene_repid_npeaks[geneid][repid] = 0
            gene_repid_npeaks[geneid][repid] = ssstats[geneid].__len__()
                        
            for summitid in ssstats[geneid]:
                this_summit_score = ssstats[geneid][summitid][0]
                this_summit_distance = ssstats[geneid][summitid][2]
                
                """Does this summit have a higher score than other summits for geneid?"""
                if gene_repid_maxsummitid[geneid][repid][0] == None:
                    gene_repid_maxsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
                if gene_repid_maxsummitid[geneid][repid][1] < this_summit_score:
                    gene_repid_maxsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
            
                """Is this summit closer to the TSS of geneid?"""
                if gene_repid_nearestsummitid[geneid][repid][0] == None:
                    gene_repid_nearestsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance)
                if gene_repid_nearestsummitid[geneid][repid][2] > this_summit_distance:
                    gene_repid_nearestsummitid[geneid][repid] = (summitid, this_summit_score, this_summit_distance) 
    
    if seen_genes.__len__() == 0: # is there at least one gene with a summit?
        """If we didn't find any summits, then return now."""
        return
    
    """Write an Excel table with genes and their scores in each replicate."""

    """Open the file and write the header row."""
    xlpath = repgroupname + ".summits.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for repid in repids:
        fout.write("NPeaks(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_score(" + repid_repname[repid].__str__() + ")\t")
        fout.write("max_summit_dist(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_score(" + repid_repname[repid].__str__() + ")\t")
        fout.write("nearest_summit_dist(" + repid_repname[repid].__str__() + ")\t")
    fout.write("\n")

    """Write one row for each gene that has a summit in at least one of the replicates."""
    for geneid in seen_genes:
        genename = get_genename(geneid, con)      
        """Columns 1 and 2 are the geneID and the gene name."""
        fout.write(geneid.__str__() + "\t" + genename.__str__() + "\t") # Gene ID and name
        
        for repid in repids:
            """Write column for N peaks."""
            if repid in gene_repid_npeaks[geneid]:
                fout.write( gene_repid_npeaks[geneid][repid].__str__() + "\t" )
            else:
                fout.write("0\t")
            
            """Write columns for max summit score."""
            foundit = False
            if geneid in gene_repid_maxsummitid:
                if repid in gene_repid_maxsummitid[geneid]:
                    foundit = True
                    this_summitid = gene_repid_maxsummitid[geneid][repid][0]
                    this_summitscore = gene_repid_maxsummitid[geneid][repid][1]
                    this_summit_distance = gene_repid_maxsummitid[geneid][repid][2] 
                    fout.write("%.3e"%this_summitscore + "\t")
                    fout.write(this_summit_distance.__str__() + "\t")
            if foundit == False:
                fout.write("0\t0\t")
            
            """Write columns for nearest summit score."""
            foundit = False
            if geneid in gene_repid_nearestsummitid:
                if repid in gene_repid_nearestsummitid[geneid]:
                    foundit = True
                    this_summitid = gene_repid_nearestsummitid[geneid][repid][0]
                    this_summitscore = gene_repid_nearestsummitid[geneid][repid][1]
                    this_summit_distance = gene_repid_nearestsummitid[geneid][repid][2] 
                    fout.write("%.3e"%this_summitscore + "\t")
                    fout.write(this_summit_distance.__str__() + "\t")
            if foundit == False:
                fout.write("0\t0\t")
        
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
    repid_qvals = {} # key = repid, value = list of Q values for summits
    for repid in repids:
        repid_qvals[repid] = []
    for geneid in seen_genes:
        for repid in repids:
            if repid in gene_repid_maxsummitid[geneid]:
                repid_qvals[repid].append( gene_repid_maxsummitid[geneid][repid][1] )
            else:
                repid_qvals[repid].append( 0.0 )

    scatter_values = []
    scatter_names = []
    plot_as_rank = []
    count = 0
    for repid in repids:
        scatter_values.append( repid_qvals[repid] )
        scatter_names.append( repid_repname[repid] + "-maxQ" )
        count += 1
    for repid in repids:
        plot_as_rank.append(count)
        scatter_values.append( repid_qvals[repid] )
        scatter_names.append( repid_repname[repid] + "-rank" )
        count += 1
    width = scatter_values.__len__()
    height = repids.__len__()
    filekeyword = repgroupname + ".summits"
    cranpath = scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, plot_as_rank=plot_as_rank)
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
    filekeyword = unionname + ".summits"
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, title="Summit scores for " + speciesname, xlab="", ylab="")
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
    rgid_sumvals = {}
    rgid_groupname = {}
    xlpath = unionname + ".enrich.xls"
    print "\n. Writing a table to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for rgid in rgroupids:
        rgname = get_repgroup_name(rgid, con)
        rgid_groupname[rgid] = rgname
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
            if foundit == False:
                fout.write("0\t0\t0\t")
                rgid_maxvals[rgid].append(0)
                rgid_meanvals[rgid].append(0)
                rgid_sumvals[rgid].append(0)
        fout.write("\n")
    fout.close()
    add_unionfile(xlpath, unionid, "Excel table with enrichment stats for union " + unionname, con)

    
    """Plot a scatterplot using scores distilled from the replicates."""
    scatter_values = []
    scatter_names = []
    for rgid in rgroupids:
        scatter_values.append( rgid_maxvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-max")
    for rgid in rgroupids:
        scatter_values.append( rgid_meanvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-mean")
    for rgid in rgroupids:
        scatter_values.append( rgid_sumvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-sum")
    for rgid in rgroupids:
        scatter_values.append( rgid_maxvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-rank")
    

    
    filekeyword = unionname + ".enrich.small"
    plot_as_rank = [] # a list of indices in scatter_data that should be plotted as rank, rather than raw.
    start = 3 * rgroupids.__len__()
    for ii in range(start, start + rgroupids.__len__() ):
        plot_as_rank.append( ii )
    width = scatter_values.__len__()
    height = rgroupids.__len__()
    cranpath = scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, plot_as_rank=plot_as_rank)
    if cranpath != None:
        add_unionfile(cranpath, unionid,"R script path to a 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)
        add_unionfile(re.sub("cran", "pdf", cranpath), unionid,"PDF 3-panel scatterplot showing max, mean, and sum fold-enrichment scores for replicates in the union " + unionname, con)   

    return
#
# continue here
# to-do
#


    """Also plot IDR."""
    scatter_values = []
    scatter_names = []
    for rgid in rgroupids:
        scatter_values.append( rgid_maxvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-max")
    for rgid in rgroupids:
        scatter_values.append( rgid_meanvals[rgid] )
        scatter_names.append( rgid_groupname[rgid] + "-mean")
    width = scatter_values.__len__()
    filekeyword = unionname + ".enrich.idr"
    """idr_stats[gene number in the geneid list][rep ii][rep jj] = lidr"""
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_values, scatter_names, filekeyword, title="", xlab="", ylab="", force_square=True)
    if cranpath != None:
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
        
        x = get_enrichment_plot_array_for_replicate(repid, speciesid, con)
        if x == None:
            return
        (x_maxe, x_meane, x_sume, x_ranke, geneids) = x
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
    filekeyword = unionname + ".enrich.idr.large"
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_data, scatter_names, filekeyword, title="IDR for Fold Enrichment for " + unionname, xlab="", ylab="")
    if cranpath != None:
        add_unionfile(cranpath, unionid, "R script for a scatterplot comparing IDR statistics for fold enrichment for all replicates in the union " + unionname, con)
        add_unionfile(re.sub("cran", "pdf", cranpath), unionid, "PDF comparing IDR statistics for fold enrichment for all replicates in the union " + unionname, con)
    
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
    sumfe = []

    for geneid in geneids:
        found_in_all = True
        for repid in repids:
            if geneid not in repid_festats[repid]:
                found_in_all = False
        if found_in_all:
            geneids_inall.append(geneid)
            this_max = 0
            this_mean = [] 
            this_sum = 0
            for repid in repids:
                if geneid in repid_festats[repid]:
                    rmax = repid_festats[repid][geneid][0]
                    rmean = repid_festats[repid][geneid][1]
                    rsum = repid_festats[repid][geneid][2]
                if rmax > this_max:
                    this_max = rmax
                this_mean.append(rmean)
                this_sum += rsum
            this_mean = mean( this_mean )
            maxfe.append( this_max )
            meanfe.append( this_mean )
            sumfe.append( this_sum )
        else:
            pass
            #print "1160: not found_in_all", geneid
            
    for ii in range(0, geneids_inall.__len__() ):
        geneid = geneids_inall[ii]
        sql = "INSERT into GroupEnrichmentStats (rgroupid, geneid, maxenrich, meanenrich, sumenrich)"
        sql += " VALUES(" + rgroupid.__str__() + "," + geneid.__str__() + ","
        sql += maxfe[ii].__str__() + ","
        sql += meanfe[ii].__str__() + ","
        sql += sumfe[ii].__str__() + ")"
        cur.execute(sql)
    con.commit()

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
    """Returns ( x_maxe[], x_meane[], x_sume[] ) for this replicate.
    There is an entry in these arrays for every gene in the species' genome.
    If something is wrong, then it returns None.
    """
    x_maxe = []
    x_meane = []
    x_sume = []
    x_ranke = []
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
        
    if x_maxe.__len__() == 0:
        print "\n. 1306"
        exit()
    if x_meane.__len__() == 0:
        print "\n. 1309"
        exit()
    if x_sume.__len__() == 0:
        print "\n. 1312"
        exit()
    if x_ranke.__len__() == 0:
        print "\n. 1315"
        exit()
    
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
    
    cur = con.cursor()
        
    if repgroupname == None:
        repgroupname = get_repgroup_name(rgroupid, con)
    
    if repids == None:
        repids = get_repids_in_group(rgroupid, con)
        
    repid_repname = {}
    for repid in repids:
        cur.execute("SELECT name from Replicates where id=" + repid.__str__())
        x = cur.fetchone()
        if x != None:
            repid_repname[repid] = x[0]
        else:
            repid_repname[repid] = repid.__str__()
    
    """Build the scatterplot"""
    scatter_max = []
    scatter_mean = []
    scatter_sum = []
    scatter_rank = []
    #scatter_geneids = []
    geneids = []
    for repid in repids:
        speciesid = get_speciesid_for_repid(repid, con)
        x = get_enrichment_plot_array_for_replicate(repid, speciesid, con)   
        if x == None:
            return 
        (x_maxe, x_meane, x_sume, x_ranke, geneids) = x
    
        scatter_max.append( x_maxe )
        scatter_mean.append( x_meane )
        scatter_sum.append( x_sume )
        scatter_rank.append( x_ranke )
        #scatter_geneid.append( geneids )

    """If there's only one replicate in this group, then we're done."""
    if repids.__len__() < 2:
        return
    
    scatter_data = scatter_max + scatter_mean + scatter_sum + scatter_max
    scatter_names = []
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-max")
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-mean")
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-sum")
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-rank(max)")
    plot_as_rank = [] # a list of indices in scatter_data that should be plotted as rank, rather than raw.
    start = 3 * repids.__len__()
    for ii in range(start, start + repids.__len__() ):
        plot_as_rank.append( ii )
    width = scatter_data.__len__()
    height = repids.__len__()
    filekeyword = repgroupname + ".enrich"
    cranpath = scatter_nxm(width, height, scatter_data, scatter_names, filekeyword, plot_as_rank=plot_as_rank)
    if cranpath != None:
        add_repgroupfile(cranpath, rgroupid, "R script for multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
        add_repgroupfile(re.sub("cran", "pdf", cranpath), rgroupid, "PDF multi-panel scatterplot with enrichment values for replicate group " + repgroupname, con)
    
        
    """Also plot the IDR stats"""
    scatter_data = scatter_max + scatter_mean
    scatter_names = []
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-max")
    for repid in repids:
        scatter_names.append(repid_repname[repid] + "-mean")
    width = scatter_data.__len__()
    height = repids.__len__()
    filekeyword = repgroupname + ".enrich.idr"
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, scatter_data, scatter_names, filekeyword, title="", xlab="", ylab="")
    if cranpath != None:
        add_repgroupfile(cranpath,rgroupid,"R script to make IDR scatterplots for replicate group " + repgroupname, con)
        add_repgroupfile(re.sub("cran", "pdf",cranpath),rgroupid,"PDF with IDR scatterplots for replicate group " + repgroupname, con)
    
    """Update the IDR stats into the database"""

    """First, clear any previous entries for the pair of this replicate group."""
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
                            #print sql
                            cur.execute(sql)
    con.commit()

    """Write an Excel Table"""
    xlpath = repgroupname + ".enrich.xls"
    print "\n. Writing a table with fold-enrichment and IDR to", xlpath
    fout = open(xlpath, "w")
    fout.write("GeneID\tGeneName\t")
    for repid in repids:
        fout.write("max(" + repid_repname[repid] + ")\t")
        fout.write("mean(" + repid_repname[repid] + ")\t")
        fout.write("sum(" + repid_repname[repid] + ")\t")
        #fout.write("rank(" + repid_repname[repid] + ")\t")
    for (ii,jj) in value_pairs:   
        fout.write("IDR(" + scatter_names[ii].__str__() + "-" + scatter_names[jj].__str__() + ")\t")     
    fout.write("\n")
    for gg in range(0, geneids.__len__() ):
        genename = get_genename(geneids[gg], con)
        if genename == None:
            continue
        fout.write(genename + "\t" )
        for ii in range(0, repids.__len__() ):
            fout.write(scatter_max[ii][gg].__str__() + "\t")
            fout.write(scatter_mean[ii][gg].__str__() + "\t")
            fout.write(scatter_sum[ii][gg].__str__() + "\t")
            
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

            
def plot_summits_vs_enrichments_for_replicates(repids, con):
    cur = con.cursor()
    
    for repid in repids:
        cur.execute("SELECT * from Replicates where id=" + repid.__str__())
        
        cur.execute("SELECT name from Replicates where id=" + repid.__str__())
        x = cur.fetchone()
        repname = None
        if x != None:
            repname = x[0].__str__()
        else:
            print "\n. chipseqdb_plot.py 1477: Error, I can't find a name for the replicate", repid
            exit()
        
        """Get the FE data"""
        festats = get_enrichment_stats_for_replicate(repid, con)
        if festats.keys().__len__() == 0:
            """Skip to the next replicate if there's no enrichment data for this rep."""
            continue
        
        """Get the summit data"""
        ssstats = get_summit_scores_for_replicate(repid, con)
        if ssstats.keys().__len__() == 0:
            """Skip to the next replicate if there's no summit data for this rep."""
            continue
        
        #print "\n. 1538 - ", festats
        print "\n. 159 - ", ssstats
    
        fearray = []
        fe_site_array = []
        sarray = []
        qarray = []
        sdarray = []
        geneids = []
        for geneid in festats:
            if geneid in ssstats:
                #print "1507:", geneid, festats[geneid][2], ssstats[geneid]
                sarray.append( ssstats[geneid][0] )
                qarray.append( ssstats[geneid][1] )
                sdarray.append( ssstats[geneid][2] )
            else:
                sarray.append( 0.0 )
                qarray.append( 0.0 )
                sdarray.append( 0.0 )
            fearray.append( festats[geneid][0] )
            fe_site_array.append(  festats[geneid][3] )
            geneids.append(geneid)
        
        scatter_values = [fearray,sarray,fe_site_array,sdarray]
        scatter_names = ["FE", "summit score", "d(max_FE)", "d(max_summit)"]
        plot_as_rank = []
        width = 4
        height = 2
        filekeyword = repname + ".enrich_x_summit"
        scatter_nxm(width, height, scatter_values, scatter_names, filekeyword, title="", force_square=False, plot_as_rank = plot_as_rank)
    
        #
        # to=do: add file to DB here
        # continue here
        #
    
        """Write an excel table."""
        xlpath = repname + ".enrich_x_summit.xls"
        fout = open(xlpath, "w")
        fout.write("GeneID\tGeneName\tmax_summit_score\td(max_summit)\tmax_fold_enrichment\td(max_fe)\n")
        for ii in range(0, geneids.__len__() ):
            geneid = geneids[ii]
            fout.write(geneid.__str__() + "\t")
            fout.write( get_genename(geneid, con) + "\t")
            score = sarray[ii]
            fout.write(score.__str__() + "\t")
            dist = sdarray[ii]
            fout.write(dist.__str__() + "\t")
            fe = fearray[ii]
            fout.write(fe.__str__() + "\t")
            dist = fe_site_array[ii]
            fout.write(dist.__str__() + "\n")
        fout.close()
    
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
    
    xlpath = spunionname + ".enrich.xls"
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

    scatter_names = unionidnames + unionidnames
    width = scatter_names.__len__()
    height = unionidnames.__len__()
    filekeyword = spunionname + ".enrich"
    cranpath = scatter_nxm(width, height, scatterdata, scatter_names, filekeyword, title="Fold Enrichment " + spunionname, xlab="fold-enrichment", ylab="fold-enrichment") 
    if cranpath != None:
        add_speciesunionfile(cranpath, uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )
        add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid, "Multipanel scatterplot showing fold enrichment data for the species-union " + spunionname, con )    
    
    idr_scatter_values = []  
    for unionid in unionids:
        idr_scatter_values.append( unionid_maxvals[unionid] )  
    width = unionidnames.__len__()
    height = unionidnames.__len__()
    filekeyword = spunionname + ".enrich.idr"
    """idr_stats[gene number in the geneid list][rep ii][rep jj] = lidr"""
    (cranpath, sinkpath, idr_stats, value_pairs) = scatter_idr_nxm(width, height, idr_scatter_values, unionidnames, filekeyword, title="Fold Enrichment IDR", xlab="", ylab="", force_square=True)
    if cranpath != None:
        add_speciesunionfile(cranpath, uid,"" + unionname, con)
        add_speciesunionfile(re.sub("cran", "pdf", cranpath), uid,"" + unionname, con)

    """Now insert the IDR stats into the SQL database."""
    
    """First, clear any previous entries for the pair of this replicate group."""
    for ii in range(0, unionids.__len__() ):
        for jj in range(0, unionids.__len__() ):
            sql = "DELETE from GeneSpeciesunionEnrichIdr where unionid1=" + unionids[ii].__str__() + " and unionid2=" + unionids[jj].__str__()
            cur.execute(sql)
            con.commit()
    
    """Next insert the data into SQL."""
    geneids = geneid_results.keys()
    for gg in range(0, geneid_results.__len__() ):
        geneid = geneids[gg]
        if gg in idr_stats:
            for ii in range(0, unionids.__len__() ):
                if ii in idr_stats[gg]:
                    for jj in range(0, unionids.__len__() ):
                        if jj in idr_stats[gg][ii]:
                            this_idr = idr_stats[gg][ii][jj]
                            ii_unionid = unionids[ii]
                            jj_unionid = unionids[jj]
                            sql = "INSERT into GeneSpeciesunionEnrichIdr(geneid,unionid1,unionid2,lidr)"
                            sql += " VALUES("
                            sql += geneid.__str__() + ","
                            sql += ii_unionid.__str__() + ","
                            sql += jj_unionid.__str__() + ","
                            sql += this_idr.__str__() + ")"
                            #print sql
                            cur.execute(sql)
                            
    con.commit()

    