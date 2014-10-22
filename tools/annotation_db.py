"""Read an annotation Excel file
See the NDT80 Analysis, October 2014, with Eugenio, for details."""

import os, sys, re
import sqlite3 as lite

def build_anno_db(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Settings(keyword, value)")
    cur.execute("CREATE TABLE IF NOT EXISTS Annotations(annoid INTEGER primary key autoincrement, sample TEXT, library_name TEXT, indexi INT, fastqpath TEXT, strain TEXT, species TEXT, tf TEXT, tag INT, media TEXT, condition TEXT, replicate INT, comment TEXT)")
    
    """All annotations will have an entry in BowtieOutput, but hybrid annotations will also
    have an entry in FilteredBowtieOutput."""
    cur.execute("CREATE TABLE IF NOT EXISTS BowtieOutput(annoid INTEGER primary key, sampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFF(speciesid TEXT, gffpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Hybrids(annoid INTEGER primary key, species1 TEXT, species2 TEXT)") # this means that the reads for annoid originally came from two species
    cur.execute("CREATE TABLE IF NOT EXISTS HybridPairs(annoid1 INTEGER, annoid2 INTEGER)") # these annos came from the same hybrid species
    cur.execute("CREATE TABLE IF NOT EXISTS Reads(readid INTEGER primary key autoincrement, readname TEXT, annoid INT, mismatch INT, order_seen INT)") # reads without mismatches
    cur.execute("CREATE TABLE IF NOT EXISTS UniqueReads(readid INTEGER primary key, annoid INT)") # reads that are unique to an annotation
    cur.execute("CREATE TABLE IF NOT EXISTS FilteredBowtieOutput(annoid INTEGER primary key, sampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReadStats(annoid INTEGER primary key, nperfect INT, ntotal INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UniqueReadStats(annoid INTEGER primary key, nunique INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsRun(exp_annoid INTEGER primary key, control_annoid INTEGER, name TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsPeakPaths(exp_annoid INTEGER primary key, treatment_pileup_path TEXT, control_lambda_path TEXT, peaks_path TEXT, summits_path TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsFE(exp_annoid INTEGER primary key, bdgpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS WigFiles(exp_annoid INTEGER primary key, org_bdgpath TEXT, wigpath TEXT)")
    
    """All annotations will be given an entry in SortedBamFiles."""
    cur.execute("CREATE TABLE IF NOT EXISTS SortedBamFiles(annoid INTEGER primary key, bampath TEXT)")
    con.commit()

def get_setting(keyword, con):
    cur = con.cursor()
    sql = "select value from Settings where keyword='" + keyword + "'"
    cur.execute(sql)
    result = cur.fetchone()
    if result.__len__() > 0:
        return result[0]
    else:
        return None

def get_setting_list(keyword, con):
    cur = con.cursor()
    sql = "select value from Settings where keyword='" + keyword + "'"
    cur.execute(sql)
    result = cur.fetchall()
    if result.__len__() > 0:
        ret_list = []
        for ii in result:
            ret_list.append( ii[0] )
        return ret_list
    else:
        return []

def get_name_for_macs(exp_annoid, control_annoid, con):
    cur = con.cursor()
    sql = "select count(*) from HybridPairs where annoid1=" + exp_annoid.__str__() + " or annoid2=" + exp_annoid.__str__()
    cur.execute(sql)
    count = cur.fetchone()[0]
    
    name = ""
    if count > 0:
        """This data is hybrid"""
        sql = "select library_name, species from Annotations where annoid=" + exp_annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        name = x[0] + "-" + x[1]
    else:
        """Not a hybrid read."""
        sql = "select library_name from Annotations where annoid=" + exp_annoid.__str__()
        cur.execute(sql)
        name = cur.fetchone()[0]
    return name

def get_macs_pairs(con):
    """Returns pairs of (experiment,control) annoids for MACS2 peak-calling."""
    pairs = [] #(treatment Id, control ID)
    
    cur = con.cursor()
    sql = "select distinct sample from Annotations"
    cur.execute(sql)
    for x in cur.fetchall():
        sample = x[0]
        sql = "select distinct replicate from Annotations where sample='" + sample + "'"
        cur.execute(sql)
        r = cur.execute(sql)
        for y in cur.fetchall():
            repid = y[0]     
            sql = "Select distinct species from Annotations where sample='" + sample + "' and replicate=" + repid.__str__()
            cur.execute(sql)
            s = cur.execute(sql)
            for z in cur.fetchall():
                species = z[0]
                #print sample, repid, species
                sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and tag=1"
                cur.execute(sql)
                treatments = cur.fetchall()
                sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and tag=0"
                cur.execute(sql)
                controls = cur.fetchall()                
                if controls.__len__() > 1:
                    print "\n. Error, I wasn't expecting to find multiple controls for", sample, species, repid
                    exit()
                control = controls[0][0]
                for t in treatments:
                    pairs.append( (t[0],control)  ) 
                    #print t[0], control, sample, repid, species
    return pairs

def get_db(dbpath):    
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
        print "\n. Creating a new database at", dbpath
    elif False == os.path.exists(dbpath):
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath
    con = lite.connect(dbpath)
    build_anno_db(con)
    return con


def import_annotations(apath, con):
    #
    # continue here -- it's not reading all the annotations correctly!
    #
    
    """Reads the annotation table and returns a Sqlite3 database object filled with data"""
    if False == os.path.exists(apath):
        print "\n. Error, I can't find your annotation file at", apath
        exit()
    
    cur = con.cursor()
    
    fin = open(apath, "r")
    lines = []
    for l in fin.readlines():
        ls = l.split("\r")
        for lt in ls:
            lines.append( lt )
    
    """We parse all the lines, EXCEPT the header line."""
    for l in lines[1:]:
        #print "142:", l
        tokens = l.split()
        if tokens.__len__() > 2:
            sample = tokens[0]
            """Restrict our import to a subset of rows, if the user specificied it."""
            if get_setting("restrict_to_sample",con) != None:
                if sample != get_setting("restrict_to_sample",con):
                    continue
            fastq = tokens[3]
            id = tokens[1]
            
            """Split the species token into (potentially) multiple tokens, seperated by underscores.
            This allows for hybrid species to be separated into two annocations, one for each parent species."""
            species = tokens[5]
            st = species.split("_")
            for species in st:
                if st.__len__() > 1:
                    id = tokens[1] + "-" + species
                    outsam = re.sub(".fastq", "-" + species + ".sam", fastq)
                else:
                    outsam = re.sub(".fastq", ".sam", fastq)
                tag = tokens[7]
                if tag == "YES":
                    tag = 1
                else:
                    tag = 0
                
                sql = "SELECT count(*) from Annotations where library_name='" + tokens[1] + "' and species='" + species + "'"
                cur.execute(sql)
                count = cur.fetchone()[0]
                if count == 0:
                    sql = "INSERT OR REPLACE INTO Annotations (sample, library_name, indexi, fastqpath, strain, species, tf, tag, media, condition, replicate, comment)"
                    sql += " VALUES('" + sample + "','" + tokens[1] + "'," + int(tokens[2]).__str__() + ",'" + fastq
                    sql += "','" + tokens[4] + "','" + species + "','" + tokens[6]
                    sql += "'," + tag.__str__() + ",'" + tokens[8] + "','" + tokens[9] + "'," + int(tokens[10]).__str__()
                    sql += ",'" + tokens[11] 
                    sql += "')"
                    #print sql
                    cur.execute(sql)
                    con.commit()
                
                """Get the annoid of the last record inserted."""
                sql = "SELECT annoid from Annotations where library_name='" + tokens[1] + "' and species='" + species + "'"
                cur.execute(sql)
                x = cur.fetchall()
                annoid = x[0][0]
                #print ". inserted annotation #", annoid, tokens[1], species
                
                """If it's a hybrid, then remember this for later."""
                if st.__len__() > 1:
                    sql = "insert or replace into Hybrids(annoid, species1, species2) "
                    sql += "VALUES('" + annoid.__str__() + "','" + st[0] + "','" + st[1] + "')"
                    cur.execute(sql)
                    #print sql
                    con.commit()           
    return con

def get_hybrid_pairs(con):
    cur = con.cursor()
    cur.execute("DELETE from HybridPairs")
    con.commit()
    
    hannoids = [] #list of hybrid annotations IDS
    sql = "select annoid from Hybrids"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        hannoids.append( ii[0] )
    seen_annoids = []
    for annoid in hannoids:
        
        if annoid in seen_annoids:
            continue # skip, we've already paired this one.
        
        sql = "select species1, species2 from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x.__len__() == 0:
            print "\n. An error occurred - 157 - ", annoid
            exit()
        species1 = x[0]
        species2 = x[1]
        
        sql = "select library_name from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x.__len__() == None:
            print "\n. An error occurred - 164 - ", annoid
            exit()
        library_name = x[0]
        
        sql = "select annoid from Annotations where library_name='" + library_name + "' and species='" + species2 + "'"
        cur.execute(sql)
        x = cur.fetchone()
        if x.__len__() == 0:
            print "\n. An error occurred - 171 -", annoid
            exit()
        pair_id = x[0]
        
        if pair_id in seen_annoids:
            continue
        
        seen_annoids.append( pair_id )
        seen_annoids.append( annoid )
        sql = "insert or replace into HybridPairs (annoid1, annoid2) VALUES("
        sql += annoid.__str__() + "," + pair_id.__str__() + ")"
        cur.execute(sql)
        con.commit()
                
            