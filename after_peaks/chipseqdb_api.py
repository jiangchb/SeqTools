import sqlite3 as lite
import os, re, sys

def get_species(con):
    cur = con.cursor()                
    cur.execute("SELECT * FROM Species")
    return cur.fetchall()

def get_species_ids(con):
    cur = con.cursor()                
    cur.execute("SELECT id FROM Species")
    return cur.fetchall()

def get_speciesid_for_rep(repid, con):
    cur = con.cursor()
    cur.execute("SELECT id FROM Species where id in (SELECT species from Replicates where id=" + repid.__str__() + ")")
    return cur.fetchone()[0]

def get_chrom_ids(con, speciesid):
    cur = con.cursor()
    cur.execute("SELECT id FROM Chromosomes where species=" + speciesid.__str__())
    x = cur.fetchall()
    if x == None:
        return None
    chromids = []
    for ii in x:
        chromids.append( ii[0] )
    return chromids

def get_chrom_id(con, name):
    cur = con.cursor()
    cur.execute("SELECT id from Chromosomes where name='" + name + "'")
    x = cur.fetchone()
    if x == None:
        return None
    else:
        return x[0]

def get_chrom_name(con, id):
    cur = con.cursor()
    cur.execute("SELECT name from Chromosomes where id=" + id.__str__())
    x = cur.fetchone()
    if x == None:
        return None
    else:
        return x[0]

def get_all_repids(con):
    cur = con.cursor()
    cur.execute("SELECT id from Replicates")
    x = cur.fetchall()
    if x != None:
        repids = []
        for ii in x:
            repids.append(ii[0])
        return repids
    else:
        return []

def get_repid(repname, speciesid, con):    
    cur = con.cursor()
    cur.execute("SELECT id from Replicates where name='" + repname.__str__() + "' and species=" + speciesid.__str__())
    return cur.fetchone()[0]

def get_repgroup_ids(con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups")
    return cur.fetchall()

def get_repgroup_id(groupname, con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups where name='" + groupname + "'")
    x = cur.fetchone()
    if x != None:
        return x[0]
    else:
        return None

def get_repgroup_name(rgroupid, con):
    cur = con.cursor()
    cur.execute("SELECT name from ReplicateGroups where id=" + rgroupid.__str__() )
    return cur.fetchone()[0]

def get_repids_in_group(rgroupid, con):
    cur = con.cursor()
    cur.execute("SELECT replicate from GroupReplicate where rgroup=" + rgroupid.__str__())
    x = cur.fetchall()
    if x == None:
        return None
    repids = []
    for ii in x:
        repids.append( ii[0] )
    return repids

def get_speciesid_for_repid(repid, con):
    cur = con.cursor()
    cur.execute("SELECT species from Replicates where id=" + repid.__str__() )
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]


def get_rgroupids_for_species(speciesid, con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups where id in (SELECT rgroup from GroupReplicate where replicate in (SELECT id from Replicates where species=" + speciesid.__str__() + "))")
    x = cur.fetchall()
    ids = []
    for ii in x:
        ids.append( ii[0] )
    return ids

def get_geneids_with_summits(con, repid):
    cur = con.cursor()
    sql = "SELECT gene FROM GeneSummits where summit in (SELECT id from Summits where replicate=" + repid.__str__() + ")"
    cur.execute(sql)
    return cur.fetchall()

def get_genes_for_species(con, speciesid):
    cur = con.cursor()
    sql = "SELECT * from Genes where chrom in (SELECT id from Chromosomes where species=" + speciesid.__str__() + ") order by id"
    cur.execute(sql)
    return cur.fetchall()

def get_genes_for_chrom(con, chromid):
    """Returns list of genes (represented as tuples), in order by their start sites"""
    cur = con.cursor()
    sql = "SELECT * FROM Genes where chrom=" + chromid.__str__() + " order by start ASC"
    #print "46:", sql
    cur.execute(sql)
    return cur.fetchall()

def get_genename(geneid, con):
    cur = con.cursor()
    sql = "SELECT name from Genes where id=" + geneid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    if x == None:
        return None
    else:
        return x[0]
    

def get_genename_for_aliasname(name, con):
    cur = con.cursor()
    sql = "SELECT realname from GeneAlias where alias='" + name + "' COLLATE NOCASE"
    cur.execute(sql)
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]

def get_aliasids_for_geneid(geneid, con):
     cur = con.cursor()
     sql = "SELECT aliasid from GeneHomology where geneid=" + geneid.__str__()
     cur.execute(sql)
     x = cur.fetchall()
     #print "x=", x
     ids = []
     for ii in x:
         ids.append(ii[0])
     return ids

def get_geneid_from_aliasid(aliasid, con):
    cur = con.cursor()
    sql = "SELECT geneid from GeneHomology where aliasid=" + aliasid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]
        
def get_geneids_from_repgroup(con, repgroupid):
    cur = con.cursor()
    sql = "SELECT geneid from RepgroupSummits where repgroupid=" + repgroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    return genes

def get_geneids_from_union(con, unionid):
    cur = con.cursor()
    sql = "SELECT geneid from UnionSummits where unionid=" + unionid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    return genes

def get_summits(con, repid, chromid):
     """Returns a list of Summits entires from the table Summits. Each element of the returned list
     is a tuple, containing all the parts of each Summit row."""
     cur = con.cursor()
     sql = "SELECT * FROM Summits where chrom=" + chromid.__str__() + " and replicate=" + repid.__str__()
     #print "53:", sql
     cur.execute(sql)
     x = cur.fetchall()
     if x == None:
         return []
     else:
         return x

def get_summit_scores_for_gene(geneid, repid, con):
    """Returns an array of tuples, one tuple for each summit at gene ID, with the tuples sorted in descending
    order based on their summit score."""
    cur = con.cursor()
    sql = "select score from Summits where replicate=" + repid.__str__() + " and id in (select summit from GeneSummits where gene=" + geneid.__str__() + ") order by score"
    cur.execute(sql)
    scores = []
    for s in cur.fetchall():
        scores.append( s[0] )
    return scores
                
def get_max_summit_score_for_gene(geneid, repid, con):
    scores = get_summit_scores_for_gene(geneid, repid, con)
    if scores.__len__() == 0:
        return None 
    return scores[0]

def get_enrichment_stats_for_gene(geneid, repid, con):
    cur = con.cursor()
    sql = "SELECT * from EnrichmentStats where repid=" + repid.__str__() + " and geneid=" + geneid.__str__()
    cur.execute(sql)
    return cur.fetchone()

def get_enrichment_stats_for_replicate(repid, con):
    """Returns gene_stats[geneid] = [stats]"""
    cur = con.cursor()
    sql = "SELECT * from EnrichmentStats where repid=" + repid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        return None
    gene_stats = {}
    for ii in x:
        geneid = ii[1]
        stats = ii[2:]
        gene_stats[ geneid ] = stats
    return gene_stats

def get_enrichment_stats_for_union( unionid, con ):
    """Returns gene_stats[geneid] = [stats]"""
    cur = con.cursor()
    sql = "SELECT * from UnionEnrichmentStats where unionid=" + unionid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        return None
    gene_stats = {}
    for ii in x:
        geneid = ii[1]
        stats = ii[2:]
        gene_stats[ geneid ] = stats
    return gene_stats

def does_species_exist(speciesname, con):
    cur = con.cursor()
    sql = "SELECT count(*) FROM Species where name='" + speciesname + "'"
    cur.execute(sql)
    if cur.fetchone()[0] > 0:
        return True
    return False

def add_species(speciesname, con):
    cur = con.cursor()
    sql = "INSERT INTO Species (name) VALUES('" + speciesname + "')"
    cur.execute( sql )
    con.commit()
    return con

def get_species_id(speciesname, con):
    cur = con.cursor()
    sql = "SELECT id FROM Species where name='" + speciesname + "'"
    cur.execute(sql)
    return cur.fetchone()[0]

def get_species_name(speciesid, con):
    cur = con.cursor()
    sql = "SELECT name FROM Species where id=" + speciesid.__str__()
    cur.execute(sql)
    return cur.fetchone()[0]

def does_repgroup_exist(rgroup, con):
    cur = con.cursor()
    sql = "SELECT count(*) from ReplicateGroups where name='" + rgroup.__str__() + "'"
    cur.execute(sql)
    if cur.fetchone()[0] > 0:
        return True
    return False
    
def add_repgroup(rgroup, con, note=None):
    cur = con.cursor()
    if note == None:
        note = "None"   
    sql = "INSERT INTO ReplicateGroups (name, note) VALUES('" + rgroup.__str__() + "','" + note + "')"
    cur.execute(sql)
    con.commit()
    return con        

def does_replicate_exist(repname, speciesid, con):
    cur = con.cursor()
    sql = "SELECT count(*) from Replicates where name='" + repname + "' and species=" + speciesid.__str__()
    cur.execute(sql)
    if cur.fetchone()[0] > 0:
        return True
    return False  
    

def add_replicate(repname, speciesid, con):
    cur = con.cursor()
    sql = "INSERT INTO Replicates (name,species) VALUES('" + repname.__str__() + "'," + speciesid.__str__() + ")"
    cur.execute( sql )
    con.commit()
    
    sql = "SELECT * from Replicates where name='" + repname.__str__() + "' and species=" + speciesid.__str__() + ""
    cur.execute(sql)
    print "\n. Adding the replicate", repname#, " [", cur.fetchall(), "]"
    
    return con

def add_rep2group(repid, rgroupid, con):
    cur = con.cursor()
    sql = "SELECT count(*) FROM GroupReplicate where rgroup=" + rgroupid.__str__() + " and replicate=" + repid.__str__()
    cur.execute(sql)
    if cur.fetchone()[0] > 0:
        return con
      
    sql = "INSERT INTO GroupReplicate (rgroup,replicate) VALUES(" + rgroupid.__str__() + "," + repid.__str__() + ")"
    cur.execute(sql)
    con.commit()
    
    return con

def clear_unions(con):
    cur = con.cursor()
    sql = "DROP TABLE IF EXISTS Unions"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS UnionRepgroups"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS UnionSummits"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS UnionSummitStats"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS UnionEnrichmentStats"
    cur.execute(sql)
    con.commit()
    
def clear_speciesunions(con):
    cur = con.cursor()
    sql = "DROP TABLE IF EXISTS Speciesunions"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS SpeciesunionUnions"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS SpeciesunionSummits"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS SpeciesunionSummitStats"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS SpeciesunionEnrichmentStats"
    cur.execute(sql)
    con.commit()
    

def add_union(unionname, repgroupnames, con):
    cur = con.cursor()
    sql = "SELECT count(*) from Unions where name='" + unionname + "'"
    cur.execute(sql)
    if cur.fetchone()[0] == 0:
        sql = "INSERT into Unions (name) VALUES('" + unionname + "')"
        cur.execute(sql)
        con.commit()
        
    # now get the union's id
    cur.execute("SELECT unionid from Unions where name='" + unionname + "'")
    unionid = cur.fetchone()[0]
    for repgroupname in repgroupnames:
        rgroupid = get_repgroup_id(repgroupname, con)
        if rgroupid == None:
            print "\n. The replicate group", repgroupname, "hasn't been imported yet."
            print " This union will be skipped."
            continue
        
        """Insert this union-repgroup pair, but only if we've never seen this pair before."""
        sql = "SELECT COUNT(*) from UnionRepgroups where unionid=" + unionid.__str__() + " and repgroupid=" + rgroupid.__str__()
        cur.execute(sql)
        if cur.fetchone()[0] == 0:
            sql = "INSERT INTO UnionRepgroups (unionid, repgroupid) VALUES(" + unionid.__str__() + "," + rgroupid.__str__() + ")"
            #print sql
            cur.execute(sql)
            con.commit()
    return con

def add_speciesunion(name, member_unions, con):
    cur = con.cursor()
    sql = "SELECT COUNT(*) from Speciesunions where name='" + name + "'"
    cur.execute(sql)
    if cur.fetchone()[0] == 0:
        sql = "INSERT into Speciesunions (name) VALUES('" + name + "')"
        cur.execute(sql)
        con.commit()
    
    # now get the union's id
    cur.execute("SELECT unionid from Speciesunions where name='" + name + "'")
    unionid = cur.fetchone()[0]
    
    for mu in member_unions:
        
        sql = "SELECT unionid from Unions where name='" + mu.__str__() + "'"
        cur.execute(sql)
        muid = cur.fetchone()
        if muid == None:            
            print "\n. The union group", mu, "hasn't been imported yet."
            print " This union will be skipped."
            continue
        muid = muid[0]
    
        sql = "SELECT COUNT(*) FROM SpeciesunionUnions where spunionid=" +unionid.__str__() + " and memunionid=" + muid.__str__()
        cur.execute(sql)
        if cur.fetchone()[0] == 0:
            """Insert this union-repgroup pair, but only if we've never seen this pair before."""
            sql = "INSERT INTO SpeciesunionUnions (spunionid, memunionid) VALUES(" + unionid.__str__() + "," + muid.__str__() + ")"
            #print sql
            cur.execute(sql)
            con.commit()
    return con

def get_unionids(con):
    cur = con.cursor()
    cur.execute("SELECT unionid from Unions")
    x = cur.fetchall()
    if x == None:
        return None
    unionids = []
    for ii in x:
        unionids.append(ii[0])
    return unionids

def get_unionname(unionid, con):
    cur = con.cursor()
    cur.execute("SELECT name from Unions where unionid=" + unionid.__str__())
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]

def get_speciesunionname(unionid, con):
    cur = con.cursor()
    cur.execute("SELECT name from Speciesunions where unionid=" + unionid.__str__())
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]

def get_repgroupids_in_union(unionid, con):
    """Returns the repgroupids sorted alphabetically by their group names."""
    cur = con.cursor()
    cur.execute("SELECT repgroupid from UnionRepgroups where unionid=" + unionid.__str__() )
    x = cur.fetchall()
    if x == None:
        return None
    repgroupname_id = {}
    for ii in x:
        repgroupname_id[  get_repgroup_name(ii[0],con)  ] = ii[0]
    names = repgroupname_id.keys()
    names.sort()
    sorted_repgroupids = []
    for name in names:
        sorted_repgroupids.append( repgroupname_id[name] )
    return sorted_repgroupids

def get_species_unionids(con):
    cur = con.cursor()
    sql = "SELECT unionid from Speciesunions"
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        return None
    unionids = []
    for ii in x:
        unionids.append( ii[0] )
    return unionids

def get_unionids_in_speciesunion(spunionid, con):
    """Returns the union IDs sorted alphabetically by their group names."""
    cur = con.cursor()
    cur.execute("SELECT memunionid from SpeciesunionUnions where spunionid=" + spunionid.__str__() )
    x = cur.fetchall()
    if x == None:
        return None
    ids = []
    for ii in x:
        ids.append( ii[0] )
    return ids

def remove_all_for_gene(gene, con):
    cur = con.cursor()
    geneid = gene[0]
    genename = gene[1]
        
    tables = ["GroupEnrichmentStats", "GroupEnrichmentStats", "RepgroupSummits", "GeneSummits", "GeneHomology", "GeneAlias", "Genes"]
    for t in tables:
        if t == "GeneSummits":
            sql = "DELETE FROM " + t + " where gene=" + geneid.__str__()
        elif t == "GeneAlias":
            sql = "DELETE FROM " + t + " where alias='" + genename.__str__() + "'"
        elif t == "Genes":
            sql = "DELETE FROM " + t + " where id=" + geneid.__str__()
        else:
            sql = "DELETE FROM " + t + " where geneid=" + geneid.__str__()
        #print sql
        cur.execute(sql)
    con.commit()
    
    # skip Summits for now

    
def reduce_db_for_test(con):
    """This method removes all entries in the DB related to Gene IDs that are not in the genelist."""
    
    print "\n. I'm building a test-sized database by eliminating entries."
    cur = con.cursor()
    
    sql = "SELECT * from GeneHomology"
    cur.execute(sql)
    results = cur.fetchall()
    import random
    random.shuffle(results)    
    if results.__len__() > 200:
        results = results[0:200]
    keep_genes = []
    for ii in results:
        geneid = ii[0]
        aliasid = ii[1]
        keep_genes.append( geneid )
        if geneid != aliasid:
            keep_genes.append( aliasid )

    count = 0

    sql = "SELECT * from Genes"
    cur.execute(sql)
    results = cur.fetchall()
    all_genes = []
    for ii in results:
        
        count += 1
        if count%50 == 0:
            print (count/float(results.__len__())).__str__() + "\%"
            #sys.stdout.write(".")
            #sys.stdout.flush()
        
        geneid = ii[0]
        if geneid not in keep_genes:
            remove_all_for_gene(ii, con)

def reset_files(con):
    print "\n. I'm removing several tables."
    cur = con.cursor()
    sqls = []
    sqls.append("DELETE from Files")
    sqls.append("DELETE from ReplicategroupFiles")
    sqls.append("DELETE from UnionFiles")
    sqls.append("DELETE from SpeciesunionFiles")
    for s in sqls:
        print s
        cur.execute(s)
    con.commit()

def add_file(filepath, note, con):
    """Returns the fileid of the newly inserted file."""
    cur = con.cursor()
    
    """Remove quotes from the string."""
    note = re.sub("'", "", note)
    note = re.sub("\"", "", note)
    
    sql = "SELECT count(*) from Files where path='" + filepath + "' and note='" + note + "'"
    cur.execute(sql)
    count = cur.fetchone()[0]
    if count == 0:  
        sql = "INSERT INTO Files (path, note) VALUES('" + filepath + "', '" + note + "')"
        cur.execute(sql)
        con.commit()
    sql = "SELECT fileid from Files where path='" + filepath + "'"
    cur.execute(sql)
    fileid = cur.fetchone()[0]
    return fileid

def add_repgroupfile(filepath, repgroupid, note, con):
    fileid = add_file(filepath, note, con)
    cur = con.cursor()
    sql = "INSERT INTO ReplicategroupFiles (repgroupid,fileid) VALUES(" + repgroupid.__str__() + ","
    sql += fileid.__str__() + ")"
    cur.execute(sql)
    con.commit()

def add_unionfile(filepath, unionid, note, con):
    fileid = add_file(filepath, note, con)
    cur = con.cursor()
    sql = "INSERT INTO UnionFiles (unionid,fileid) VALUES(" + unionid.__str__() + ","
    sql += fileid.__str__() + ")"
    cur.execute(sql)
    con.commit()

def add_speciesunionfile(filepath, spunionid, note, con):
    fileid = add_file(filepath, note, con)
    cur = con.cursor()
    sql = "INSERT INTO SpeciesunionFiles (spunionid,fileid) VALUES(" + spunionid.__str__() + ","
    sql += fileid.__str__() + ")"
    cur.execute(sql)
    con.commit()

def are_replicates_same_species(repids, con):
    """Verify that all replicates belong to the same species."""
    cur = con.cursor()
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
            return False
            #print "\n. Error while computing summits for replicates in the group", repgroupname
            #print ". It appears that not all the replicates belong to the same species."
    return True