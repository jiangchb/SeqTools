"""Read an annotation Excel file
See the NDT80 Analysis, October 2014, with Eugenio, for details."""

import os, sys, re
import sqlite3 as lite

def build_anno_db(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Settings(keyword, value)")
    
    """The Annotations table is structured based on the annotation Excel files generated by Eugenio Mancera.
        It would be good to normalize this table, however, and assign ID integers to species, tf, conditions, etc."""
    cur.execute("CREATE TABLE IF NOT EXISTS Annotations(annoid INTEGER primary key autoincrement, sample TEXT, library_name TEXT, indexi INT, fastqpath TEXT, strain TEXT, species TEXT, tf TEXT, tag INT, media TEXT, condition TEXT, replicate INT, comment TEXT)")
    
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesGenomepath(speciesname TEXT unique, genomepath TEXT)")
    
    """All annotations will have an entry in BowtieOutput, 
    but only hybrid annotations will also have an entry in FilteredBowtieOutput."""
    cur.execute("CREATE TABLE IF NOT EXISTS BowtieOutput(annoid INTEGER primary key, sampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS FilteredBowtieOutput(annoid INTEGER primary key, sampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFF(speciesid TEXT, gffpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Hybrids(annoid INTEGER primary key, species1 TEXT, species2 TEXT)") # this means that the reads for annoid originally came from two species
    cur.execute("CREATE TABLE IF NOT EXISTS HybridPairs(annoid1 INTEGER, annoid2 INTEGER)") # these annos came from the same hybrid species
    
    # A legacy hack to eliminate depricated tables (i.e. from the old schema).
    cur.execute("drop table if exists Reads")
    cur.execute("drop table if exists UniqueReads")
    
    cur.execute("CREATE TABLE IF NOT EXISTS ReadStats(annoid INTEGER primary key, nperfect INT, ntotal INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UniqueReadStats(annoid INTEGER primary key, nunique INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsRun(exp_annoid INTEGER primary key, control_annoid INTEGER, name TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsPeakPaths(exp_annoid INTEGER primary key, treatment_pileup_path TEXT, control_lambda_path TEXT, peaks_path TEXT, summits_path TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsFE(exp_annoid INTEGER primary key, bdgpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS FEWigFiles(exp_annoid INTEGER primary key, org_bdgpath TEXT, wigpath TEXT)")
    
    """All annotations will be given an entry in SortedBamFiles."""
    cur.execute("CREATE TABLE IF NOT EXISTS SortedBamFiles(annoid INTEGER primary key, bampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS BedgraphFiles(annoid INTEGER primary key, bedpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReadsWigFiles(annoid INTEGER primary key, wigpath TEXT)")
    
    cur.execute("create table if not exists Log(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    cur.execute("create table if not exists ErrorLog(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    con.commit()

def get_setting(keyword, con):
    cur = con.cursor()
    sql = "select value from Settings where keyword='" + keyword + "'"
    cur.execute(sql)
    result = cur.fetchone()
    if result == None:
        return None
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
    
    """Returns a list of tuples, each being (experiment,control) Annotation object ids for MACS2 peak-calling.
    Each pair shares the same unique value in the fields 'sample', 'species', and 'replicate' in the table Annotations.
    The control member of the pair is found by finding the Annotation with 'tag' = 0, 
    while the treatment member of the pair has 'tag'=1."""
    pairs = [] #(treatment Id, control ID)
    
    cur = con.cursor()
    sql = "select distinct tf from Annotations"
    cur.execute(sql)
    for dd in cur.fetchall():
        tf = dd[0]
        sql = "select distinct sample from Annotations where tf='" + tf.__str__() + "'"
        cur.execute(sql)
        for x in cur.fetchall():
            sample = x[0]
            sql = "select distinct replicate from Annotations where sample='" + sample + "' and tf='" + tf.__str__() + "'"
            cur.execute(sql)
            r = cur.execute(sql)
            for y in cur.fetchall():
                repid = y[0]     
                sql = "Select distinct species from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and  tf='" + tf.__str__() + "'"
                cur.execute(sql)
                s = cur.execute(sql)
                for z in cur.fetchall():
                    species = z[0]
                    sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and tag=1 and tf='" + tf.__str__() + "'"
                    cur.execute(sql)
                    treatments = cur.fetchall()
                    sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and tag=0 and tf='" + tf.__str__() + "'"
                    cur.execute(sql)
                    controls = cur.fetchall()               
                    if controls.__len__() > 1:
                        print "\n. Error, I wasn't expecting to find multiple controls for", sample, species, repid
                        exit()
                    control = controls[0][0]
                    for t in treatments:
                        pairs.append( (t[0],control)  ) 
                        print "Experiment:", t[0], "Control:", control, "Sample:", sample, "ReplicateID:",  repid, "Species:", species
    return pairs

def get_db(dbpath):    
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
        print "\n. Creating a new database at", dbpath
    elif False == os.path.exists(dbpath):
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath
    con = lite.connect(dbpath, timeout=1)
    build_anno_db(con)
    return con

def import_genome_list(gpath, con):
    """Reads the user-specificed genome list. Each line should contain two things:
    1. the species identified (a short string), and
    2. the path to the genome
    """
    
    cur = con.cursor()
    
    sql = "delete from SpeciesGenomepath"
    cur.execute(sql)
    con.commit()
    
    if False == os.path.exists(gpath.__str__()):
        print "\n. Error: I cannot find  your genome list file at " + gpath
        print os.getcwd()  
        print os.path.exists(gpath)      
        exit()
    
    fin = open(gpath, "r")
    lines = fin.readlines()
    fin.close()
    for l in lines:
        if l.startswith("#") or l.__len__() < 2:
            continue
        tokens = l.split()
        if tokens.__len__() < 2:
            continue
        speciesname = tokens[0]
        gpath = tokens[1]
        #if False == os.path.exists(gpath):
        #    print "\n. Error: your genome list includes a reference that I cannot find:"
        #    print l
        #    exit()
        sql = "insert into SpeciesGenomepath (speciesname, genomepath) VALUES("
        sql += "'" + speciesname + "','" + gpath + "')"
        cur.execute(sql)
        con.commit()
        
        print ". I found a genome for species", speciesname, "at", gpath
        
    return con
    
def import_annotations(apath, con):    
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
            if lt.startswith("#"):
                continue
            if lt.__len__() < 2:
                continue
            lines.append( lt )
    
    """We parse all the lines, skipping the first header line.
        Warning: if a header line is missing in the annotations file, then the first
        annotation will be skipped. To-do: add code to check if a header line exists."""
    for l in lines[1:]:
        tokens = l.split()
        if tokens.__len__() > 2:
            sample = tokens[0]
            """Restrict our import to a subset of rows, if the user specificied it."""
            if get_setting("restrict_to_sample",con) != None:
                if sample != get_setting("restrict_to_sample",con):
                    continue
            id = tokens[1] # id is also known as library_name
            indexi = int(tokens[2])
            fastq = tokens[3]
            strain = tokens[4]
            species = tokens[5]
            tf = tokens[6]
            tag = tokens[7]
            if tag == "YES":
                tag = 1
            else:
                tag = 0
            media = tokens[8]
            condition = tokens[9]
            replicate = int(tokens[10])
            comment = tokens[11]
            
            """Split the species token into (potentially) multiple tokens, each separated by underscores.
            This allows for hybrid species to be separated into two annocations, one for each parent species."""
            st = species.split("_") # st stands for species tokens
            for species in st:
                """If the species name is not a hybrid, then this loop will run just once."""
                if st.__len__() > 1:
                    id = tokens[1] + "-" + species
                    outsam = re.sub(".fastq", "-" + species + ".sam", fastq)
                else:
                    outsam = re.sub(".fastq", ".sam", fastq)

                
                """Does this Annotation's settings conflict with any previously-known Annotations?"""
                sql = "SELECT count(*) from Annotations where library_name='" + id + "' and species='" + species + "'"
                cur.execute(sql)
                count = cur.fetchone()[0]
                if count == 0:
                    """i.e. this Annotation is unique, so let's import it."""
                    sql = "INSERT OR REPLACE INTO Annotations (sample, library_name, indexi, fastqpath, strain, species, tf, tag, media, condition, replicate, comment)"
                    sql += " VALUES('" + sample + "',"
                    sql += "'" + id + "',"
                    sql += indexi.__str__() + ","
                    sql += "'" + fastq + "',"
                    sql += "'" + strain + "',"
                    sql += "'" + species + "',"
                    sql += "'" + tf + "',"
                    sql += tag.__str__() + ","
                    sql += "'" + media + "',"
                    sql += "'" + condition + "',"
                    sql += replicate.__str__() + ","
                    sql += "'" + comment + "'"
                    sql += ")"
                    #print sql
                    cur.execute(sql)
                    con.commit()
                
                """Get the annoid of the annotation we just imported."""
                sql = "SELECT annoid from Annotations where library_name='" + id + "' and species='" + species + "'"
                cur.execute(sql)
                x = cur.fetchall()
                annoid = x[0][0]
                
                """If it's a hybrid, then remember that we created TWO Annotation to represent it."""
                if st.__len__() > 1:
                    sql = "insert or replace into Hybrids(annoid, species1, species2) "
                    sql += "VALUES(" + annoid.__str__() + ",'" + st[0] + "','" + st[1] + "')"
                    cur.execute(sql)
                    #print sql
                    con.commit()           
    return con

def get_hybrid_pairs(con):
    cur = con.cursor()
    cur.execute("DELETE from HybridPairs")
    con.commit()
    
    hannoids = [] # hannoids is a list of hybrid annotations IDS
    sql = "select annoid from Hybrids"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        hannoids.append( ii[0] )
    seen_annoids = []
    
    print "285:", hannoids
    
    for annoid in hannoids:
        
        if annoid in seen_annoids:
            continue # skip, we've already paired this one.
        
        """Get the two species that form this hybrid"""
        sql = "select species1, species2 from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x == None:
            print "\n. An error occurred - 157 - ", annoid
            exit()
        species1 = x[0]
        species2 = x[1]
        
        print "302:", species1, species2
        
        """Get the library_name (i.e. ID)"""
        sql = "select library_name, strain from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x == None:
            print "\n. An error occurred - 164 - ", annoid
            exit()
        library_name = x[0]
        strain = x[1]
                
        """Get the annotation ID for the paired annotation."""
        sql = "select annoid from Annotations where species='" + species2.__str__() + "' and strain='" + strain.__str__() + "'"
        print sql
        cur.execute(sql)
        x = cur.fetchone()
        print "318:", x
        if x == None:
            print "\n. An error occurred - 171 -", annoid
            exit()
        the_other_annoiid = x[0]
        
        if the_other_annoiid in seen_annoids:
            continue
        
        seen_annoids.append( the_other_annoiid )
        seen_annoids.append( annoid )
        sql = "insert or replace into HybridPairs (annoid1, annoid2) VALUES("
        sql += annoid.__str__() + "," + the_other_annoiid.__str__() + ")"
        cur.execute(sql)
        con.commit()

def write_log(con, message, code=None):
    """
    Writes to the log file
    """
    cur = con.cursor()
    sql = "insert into Log (message"
    if code != None:
        sql += ",code"
    sql += ") values(\"" + message
    if code != None:
        sql += "\"," + code.__str__() + ")"
    else:
        sql += "\")"
    cur.execute(sql)
    con.commit()
    
    print "\n. " + message
    
def write_error(con, message, code=None):
    cur = con.cursor()
    sql = "insert into ErrorLog (message"
    if code != None:
        sql += ",code"
    sql += ") values(\"" + message
    if code != None:
        sql += "\"," + code.__str__() + ")"
    else:
        sql += "\")"
    cur.execute(sql)
    con.commit()
    print "\n. ERROR: " + message
                
            