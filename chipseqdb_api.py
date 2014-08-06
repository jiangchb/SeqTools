import sqlite3 as lite
import os, sys

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

def get_repid(repname, speciesid, con):    
    cur = con.cursor()
    cur.execute("SELECT id from Replicates where name='" + repname.__str__() + "' and species=" + speciesid.__str__())
    return cur.fetchone()[0]

def get_repgroup_ids(con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups")
    return cur.fetchall()

def get_repgroup_id(rgroup, con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups where name='" + rgroup + "'")
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
    sql = "SELECT * from Genes where chrom in (SELECT id from Chromosomes where species=" + speciesid.__str__() + ")"
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

# August 4 2014: This  method is not currently being used. It may be depricated.
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
    sql = "SELECT geneid from RepgroupGenes where repgroupid=" + repgroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    return genes

def get_geneids_from_union(con, unionid):
    cur = con.cursor()
    sql = "SELECT geneid from UnionGenes where unionid=" + unionid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    return genes

def get_summits(con, repid, chromid):
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
    cur = con.cursor()
    sql = "select score from Summits where replicate=" + repid.__str__() + " and id in (select summit from GeneSummits where gene=" + geneid.__str__() + ")"
    #print "60:", sql
    cur.execute(sql)
    scores = []
    for s in cur.fetchall():
        scores.append( s[0] )
    return scores

def get_max_summit_score_for_gene(geneid, repid, con):
    scores = get_summit_scores_for_gene(geneid, repid, con)
    if scores.__len__() == 0:
        return None 
    return max( scores )

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

def import_species(speciesname, con):
    cur = con.cursor()
    sql = "REPLACE INTO Species (name) VALUES('" + speciesname + "')"
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
    sql = "SELECT name FROM Species where id=" + speciesid
    cur.execute(sql)
    return cur.fetchone()[0]

def add_repgroup(rgroup, con, note=None):
    with con:
        cur = con.cursor()
        if note == None:
            note = "None"
        sql = "REPLACE INTO ReplicateGroups (name, note) VALUES('" + rgroup.__str__() + "','" + note + "')"
        cur.execute(sql)
        con.commit()
    return con        
    
def add_replicate(repname, speciesid, con):
    new_id = None
    with con:
        cur = con.cursor()
        sql = "REPLACE INTO Replicates (name,species) VALUES('" + repname.__str__() + "'," + speciesid.__str__() + ")"
        cur.execute( sql )
        con.commit()
    return con


def add_rep2group(repid, rgroupid, con):
    with con:
        cur = con.cursor()
        sql = "REPLACE INTO GroupReplicate (rgroup,replicate) VALUES(" + rgroupid.__str__() + "," + repid.__str__() + ")"
        cur.execute(sql)
        con.commit()
    return con

def clear_unions(con):
    cur = con.cursor()
    sql = "DROP TABLE IF EXISTS Unions"
    cur.execute(sql)
    con.commit()
    
    sql = "DROP TABLE IF EXISTS UnionRepgroups"
    cur.execute(sql)
    con.commit()
    
    sql = "DROP TABLE IF EXISTS UnionGenes"
    cur.execute(sql)
    con.commit()

    sql = "DROP TABLE IF EXISTS UnionSummitStats"
    cur.execute(sql)
    con.commit()

    sql = "DROP TABLE IF EXISTS UnionEnrichmentStats"
    cur.execute(sql)
    con.commit()
    
def clear_speciesunions(con):
    cur = con.cursor()
    sql = "DROP TABLE IF EXISTS Speciesunions"
    cur.execute(sql)
    con.commit()
    
    sql = "DROP TABLE IF EXISTS SpeciesunionUnions"
    cur.execute(sql)
    con.commit()
    
    sql = "DROP TABLE IF EXISTS SpeciesunionGenes"
    cur.execute(sql)
    con.commit()
    
    sql = "DROP TABLE IF EXISTS SpeciesunionSummitStats"
    cur.execute(sql)
    con.commit()

    sql = "DROP TABLE IF EXISTS SpeciesunionEnrichmentStats"
    cur.execute(sql)
    con.commit()


def build_unions(con):
    cur = con.cursor()
    # These tables describe which replicates are to be unioned.
    cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS UnionRepgroups(unionid INTEGER, repgroupid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS UnionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    cur.execute("CREATE TABLE IF NOT EXISTS UnionSummitStats(unionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
    con.commit()


def build_speciesunions(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Speciesunions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionUnions(spunionid INTEGER, memunionid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    
    # geneid in the Speciesunion tables point to translated gene IDs from the pillars.
    # This means that to find this gene ID in a particular Union, you may need to use
    # the table GeneHomology to find alias gene IDs for another species.
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionSummitStats(spunionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
    con.commit()   
    

def add_union(unionname, repgroupnames, con):
    cur = con.cursor()
    sql = "REPLACE into Unions (name) VALUES('" + unionname + "')"
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
        sql = "REPLACE INTO UnionRepgroups (unionid, repgroupid) VALUES(" + unionid.__str__() + "," + rgroupid.__str__() + ")"
        print sql
        cur.execute(sql)
        con.commit()
    return con

def add_speciesunion(name, member_unions, con):
    cur = con.cursor()
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
    
        """Insert this union-repgroup pair, but only if we've never seen this pair before."""
        sql = "REPLACE INTO SpeciesunionUnions (spunionid, memunionid) VALUES(" + unionid.__str__() + "," + muid.__str__() + ")"
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

def get_unionids_in_speciesunion(unionid, con):
    """Returns the union IDs sorted alphabetically by their group names."""
    cur = con.cursor()
    cur.execute("SELECT memunionid from SpeciesunionUnions where spunionid=" + unionid.__str__() )
    x = cur.fetchall()
    if x == None:
        return None
    ids = []
    for ii in x:
        ids.append( ii[0] )
    return ids



    