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
    
    if dbpath == None:
        dbpath = "test.db"
    con = lite.connect(dbpath)
    
    with con:
        cur = con.cursor()
        cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT, start INT, stop INT, chrom INT, strand TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(gene INT, name TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species ID)")
        cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key autoincrement, name TEXT)")
        
        cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, count INT, species INT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Peaks(id INTEGER primary key autoincrement, replicate INT, name TEXT, start INT, stop INT, chrom INT, enrich FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GenePeaks(gene INTEGER, peak INT, distance INT)") # Dist is distance to the trans. start site
        cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INT, distance INT)")

    return con

def get_species_ids(con):
    cur = con.cursor()                
    cur.execute("SELECT * FROM Species")
    return cur.fetchall()
     
def get_chrom_ids(con, speciesid):
    cur = con.cursor()
    cur.execute("SELECT * FROM Chromosomes where species=" + speciesid.__str__())
    return cur.fetchall()

def get_gene_ids(con, chromid):
     cur = con.cursor()
     cur.execute("SELECT * FROM Genes where chrom=" + chromid.__str__())
     return cur.fetchall()

def get_summit_ids(con, repid, chromid):
     cur = con.cursor()
     cur.execute("SELECT * FROM Summits where chrom=" + chromid.__str__() + " and replicate=" + repid.__str__())
     return cur.fetchall()

# Build a library of genes
#
def import_gff(gffpath, speciesid, con):
    print "\n. Importing GFF from", gffpath
    
    """Returns chr_gene_sites"""
    cur = con.cursor()
    
    chr_gene_sites = {} # key = chromosome name, value = hash; key = gene name, value = (Start, stop)
    fin = open(gffpath, "r")
    for l in fin.xreadlines():
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            if tokens[2] != "gene": # We're only interested in genes
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
    print "\n. I have ", count_genes, "total genes."
    
    return con

def import_summits(summitpath, repid, con):
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
                cur.execute("SELECT * FROM Chromosomes WHERE name='" + chr + "'")
                chrid = cur.fetchone()[0]
                sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chrid.__str__() + "," + score.__str__() + ")"
                cur.execute(sql)
    fin.close()
    
    cur.execute("SELECT COUNT(*) FROM Summits")
    count_genes = cur.fetchone()[0]
    print "\n. I have", count_genes, "total summits."
    return con


def map_summits2genes(con, repid, speciesid=None, chromid=None):
    
    with con:
        #
        # for each species:
        #
        cur = con.cursor()
        speciesids = []
        if speciesid:
            species = [speciesid]
        else:
            species = get_species_ids(con)
        for spid in species:
            chroms = get_chrom_ids(con, spid)
            for chr in chroms:
                chrid = chr[0]
                genes = get_gene_ids(con, chrid)
                summits = get_summit_ids(con, repid, chrid)
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
                            print sql         
                            cur.execute(sql) 
                    if closest_down != None and min_down != None:
                        with con:
                            sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                            sql += " VALUES(" + closest_down.__str__() + "," 
                            sql += sid.__str__() + ","
                            sql += min_down.__str__() + ") "  
                            print sql           
                            cur.execute(sql) 
    return con
                
     