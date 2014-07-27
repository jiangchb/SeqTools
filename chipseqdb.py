############################################################
#
# A Python API / wrapper for the SQLLite ChIP-Seq database
#
############################################################

import sqlite3 as lite
import sys

def build_db(dbpath = None):
    """Initializes all the tables.
    Returns the DB connection object."""
    
    print "\n. Building the database. . ."
    
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
    con = lite.connect(dbpath)
    
    with con:
        cur = con.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT, start INT, stop INT, chrom INT, strand TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(gene INT, name TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species INT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key autoincrement, name TEXT)")

        cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, name TEXT unique, species INT)")        
        cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key autoincrement, name TEXT, note TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GroupReplicate(rgroup INTEGER, replicate INTEGER)")
        cur.execute("CREATE TABLE IF NOT EXISTS Peaks(id INTEGER primary key autoincrement, replicate INT, name TEXT, start INT, stop INT, chrom INT, enrich FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GenePeaks(gene INTEGER, peak INT, distance INT)") # Dist is distance to the trans. start site
        cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INT, distance INT)")

    return con

def get_species_ids(con):
    cur = con.cursor()                
    cur.execute("SELECT id FROM Species")
    return cur.fetchall()

def get_chrom_ids(con, speciesid):
    cur = con.cursor()
    cur.execute("SELECT id FROM Chromosomes where species=" + speciesid.__str__())
    return cur.fetchall()

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
    return cur.fetchone()[0]

def get_reps_in_group(rgroupid, con):
    cur = con.cursor()
    cur.execute("SELECT replicate from GroupReplicate where rgroup=" + rgroupid.__str__())
    return cur.fetchall()

def get_geneids(con, repid):
    cur = con.cursor()
    sql = "SELECT gene FROM GeneSummits where summit in (SELECT id from Summits where replicate=" + repid.__str__() + ")"
    cur.execute(sql)
    return cur.fetchall()

def get_genes(con, chromid):
     cur = con.cursor()
     sql = "SELECT * FROM Genes where chrom=" + chromid.__str__()
     #print "46:", sql
     cur.execute(sql)
     return cur.fetchall()

def get_summits(con, repid, chromid):
     cur = con.cursor()
     sql = "SELECT * FROM Summits where chrom=" + chromid.__str__() + " and replicate=" + repid.__str__()
     #print "53:", sql
     cur.execute(sql)
     return cur.fetchall()

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

def import_species(speciesname, con):
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


def import_gff(gffpath, speciesid, con, restrict_to_feature = "gene"):
    print "\n. Importing GFF from", gffpath
    
    """Returns chr_gene_sites"""
    cur = con.cursor()
    
    chr_gene_sites = {} # key = chromosome name, value = hash; key = gene name, value = (Start, stop)
    fin = open(gffpath, "r")
    for l in fin.xreadlines():
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            if tokens[2] != restrict_to_feature: # e.g., if restrict_to_feature == "gene", then we'll only import genes.
                continue
            chr = tokens[0]
            start = int( tokens[3] )
            stop = int( tokens[4] )
            strand = tokens[6]
            if strand == "-":
                x = start
                y = stop
                start = y
                stop = x
            gene = tokens[8].split(";")[0].split("=")[1] # orfName
            
            with con:
                cur.execute("SELECT COUNT(*) FROM Chromosomes WHERE name='" + chr + "'")
                data = cur.fetchone()
                if data[0] == 0: # that chromosome doesn't exist yet.
                    sql = "INSERT INTO Chromosomes (name,species) VALUES('" + chr + "'," + speciesid.__str__() + ")"
                    #print sql
                    cur.execute( sql )
                    con.commit()
            
                cur = con.cursor()                
                cur.execute("SELECT * FROM Chromosomes WHERE name='" + chr + "'")
                data = cur.fetchone()
                chrid = data[0]
                sql = "INSERT INTO Genes (name, start, stop, chrom,strand) VALUES('" + gene + "'," + start.__str__() + "," + stop.__str__() + "," + chrid.__str__() + ",'" + strand + "')"
                cur.execute(sql)
    fin.close() 
    
    cur.execute("SELECT COUNT(*) FROM Genes")
    count_genes = cur.fetchone()[0]
    #print "\n. The database now contains", count_genes, "total genes."
    
    return con

def add_repgroup(rgroup, con, note=None):
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM ReplicateGroups WHERE name='" + rgroup + "'")
        count = cur.fetchone()[0]
        if count == 0:
            if note == None:
                note = "None"
            sql = "INSERT INTO ReplicateGroups (name, note) VALUES('" + rgroup.__str__() + "','" + note + "')"
            cur.execute(sql)
            con.commit()
    return con        
    
def add_replicate(repname, speciesid, con):
    print "164:"
    new_id = None
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM Replicates WHERE name='" + repname.__str__() + "' and species=" + speciesid.__str__() )
        data = cur.fetchone()
        if data[0] == 0: # that chromosome doesn't exist yet.
            sql = "INSERT INTO Replicates (name,species) VALUES('" + repname.__str__() + "'," + speciesid.__str__() + ")"
            cur.execute( sql )
            con.commit()
    return con


def add_rep2group(repid, rgroupid, con):
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM GroupReplicate WHERE rgroup=" + rgroupid.__str__() + " and replicate=" + repid.__str__())
        count = cur.fetchone()[0]
        if count == 0:
            sql = "INSERT INTO GroupReplicate (rgroup,replicate) VALUES(" + rgroupid.__str__() + "," + repid.__str__() + ")"
            cur.execute(sql)
            con.commit()
    return con
        

def import_summits(summitpath, repid, con):
    """Reads a summit file (from MACS2, for example), and puts values into the Summits table."""
    
    print "\n. Importing summits from", summitpath,"for replicate", repid
    cur = con.cursor()
        
    #
    # Build a library of summits
    #
    chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
    fin = open(summitpath, "r")
    for l in fin.xreadlines():
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            chr = tokens[0]
            site = int( tokens[1] )
            name = tokens[3]
            score = float( tokens[4] )
            
            with con:
                cur = con.cursor()
                sql = "SELECT * FROM Chromosomes WHERE name='" + chr + "'"                
                cur.execute(sql)
                chrid = cur.fetchone()[0]
                sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chrid.__str__() + "," + score.__str__() + ")"
                cur.execute(sql)
    fin.close()
    
    cur.execute("SELECT COUNT(*) FROM Summits")
    count_genes = cur.fetchone()[0]
    #print "\n. I have", count_genes, "total summits."
    return con


def map_summits2genes(con, repid, speciesid=None, chromid=None):
    """This methods puts values into the table GeneSummits."""
    
    with con:
        cur = con.cursor()
        speciesids = []
        if speciesid:
            spids = [speciesid]
        else:
            spids = get_species_ids(con)
            
        for spid in spids:
            chroms = get_chrom_ids(con, spid)
            for chr in chroms:
                chrid = chr[0]
                genes = get_genes(con, chrid)
                summits = get_summits(con, repid, chrid)
                for s in summits:
                    sid = s[0]
                    sumsite = s[3]
                    score = s[5]
                    
                    min_up = None
                    closest_up = ""
                    min_down = None
                    closest_down = ""
                    for g in genes:
                        gid = g[0]
                        start = g[2]
                        stop = g[3]
                        d = start - sumsite
                        
                        """Sense direction, and upstream"""
                        if start < stop and start >= sumsite:
                            if min_up == None:
                                min_up = d
                                closest_up = gid
                            if min_up > d:
                                min_up = d
                                closest_up = gid
                                
                        elif start < stop and stop < sumsite:
                            """Sense direction and downstream"""
                            if min_down == None:
                                min_down = d
                                closest_down = None
                            if min_down < d:
                                min_down = d
                                closest_down = None
                            
                        elif start > stop and start <= sumsite:
                            """Antisense and downstream"""
                            if min_down == None:
                                min_down = d
                                closest_down = gid
                            if min_down < d: # be careful here, we're looking for the largest NEGATIVR number.
                                min_down = d
                                closest_down = gid
                                
                        elif start > stop and stop > sumsite:
                            """Antisense and upstream"""
                            if min_up == None:
                                min_up = d
                                closest_up = None
                            if min_up > d:
                                min_up = d
                                closest_up = None

                    if closest_up != None and min_up != None:
                        with con:
                            sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                            sql += " VALUES(" + closest_up.__str__() + "," 
                            sql += sid.__str__() + ","
                            sql += min_up.__str__() + ") "
                            #print sql         
                            cur.execute(sql) 
                    if closest_down != None and min_down != None:
                        with con:
                            sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                            sql += " VALUES(" + closest_down.__str__() + "," 
                            sql += sid.__str__() + ","
                            sql += min_down.__str__() + ") "  
                            #print sql           
                            cur.execute(sql) 
    return con
                
     