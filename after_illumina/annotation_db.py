"""Read an annotation Excel file"""

import os, sys, re
import sqlite3 as lite

def build_anno_db(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Settings(keyword TEXT, value TEXT)")
    
    """The Annotations table is structured based on the annotation Excel files generated by Eugenio Mancera.
        It would be good to normalize this table, however, and assign ID integers to species, tf, conditions, etc."""
    cur.execute("CREATE TABLE IF NOT EXISTS Annotations(annoid INTEGER primary key, sample TEXT, library_name TEXT, indexi INT, fastqpath TEXT, strain TEXT, species TEXT, tf TEXT, tag INT, media TEXT, condition TEXT, replicate INT, comment TEXT)")
    
    """READS"""
    cur.execute("CREATE TABLE IF NOT EXISTS FastqFiles(id INTEGER primary key, filepath TEXT unique)")    
    cur.execute("CREATE TABLE IF NOT EXISTS Conditions(id INTEGER primary key, name TEXT unique)")    # growth condition
    cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key, name TEXT unique)")         # the tagged transcription factor
    cur.execute("CREATE TABLE IF NOT EXISTS Reads(id integer primary key, name TEXT unique, fastqid INT, speciesid INT, conditionid INT, geneid INT, tagged INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReadComments(readid integer primary key, comment TEXT)")
    
    """PAIRs"""
    cur.execute("CREATE TABLE IF NOT EXISTS Pairs(id integer primary key, name TEXT unique, controlid INT, taggedid INT)") # controlid and taggedid are IDs Read entries
    cur.execute("CREATE TABLE IF NOT EXISTS Shiftsizes(pairid integer, shiftsize float)") # controlid and taggedid are IDs Read entries
    
    cur.execute("CREATE TABLE IF NOT EXISTS Comparisons(id integer primary key, name TEXT unique)")
    cur.execute("CREATE TABLE IF NOT EXISTS PairsComparisons(compid INT, pairid INT)") # maps pairs into comparisons
    cur.execute("CREATE TABLE IF NOT EXISTS CompareComparisons(compid INT, targetid INT)")
    
    """Species and Genome"""
    cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key, name TEXT unique)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesGenomepath(speciesid INT unique, genomepath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesBowtieindex(speciesid TEXT unique, indexpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesGenomesize(speciesid INT unique, genomesize INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFF(speciesid INT, gffpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS IntergenicPaths(speciesid INT, intergenicpath TEXT)")
    
    """All annotations will have an entry in BowtieOutput, 
        but only hybrid annotations will also have an entry in FilteredBowtieOutput."""
    cur.execute("CREATE TABLE IF NOT EXISTS BowtieOutput(readid INTEGER primary key, sampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS FilteredBowtieOutput(readid INTEGER primary key, sampath TEXT)")

    cur.execute("CREATE TABLE IF NOT EXISTS HybridPairs(readid1 INTEGER, readid2 INTEGER)") # these readids came from the same hybrid species
    
    cur.execute("CREATE TABLE IF NOT EXISTS ReadStats(readid INTEGER primary key, nperfect INT, ntotal INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UniqueReadStats(readid INTEGER primary key, nunique INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsRun(id INTEGER primary key, pairid INT, name TEXT unique)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsPeakPaths(macsrunid INT primary key, treatment_pileup_path TEXT, control_lambda_path TEXT, peaks_path TEXT, summits_path TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS MacsFE(macsrunid INT primary key, bdgpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS FEWigFiles(macsrunid INT primary key, org_bdgpath TEXT, wigpath TEXT)")
    
    """All annotations will be given an entry in SortedBamFiles."""
    cur.execute("CREATE TABLE IF NOT EXISTS SortedBamFiles(readid INTEGER primary key, bampath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS BedgraphFiles(readid INTEGER primary key, bedpath TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReadsWigFiles(readid INTEGER primary key, wigpath TEXT)")
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

def import_generic(con, colval, table, colkeyword="name"):
    cur = con.cursor()
    sql = "select id, " + colkeyword + " from " + table + " where " + colkeyword + "='" + colval.__str__() + "'"
    cur.execute(sql)
    x = cur.fetchone()
    if x != None:
        return x[0]
    sql = "insert or ignore into " + table + " (" + colkeyword + ") values('" + colval.__str__() + "')"
    cur.execute(sql)
    con.commit()
    sql = "select id, " + colkeyword + " from " + table + " where " + colkeyword + "='" + colval.__str__() + "'"
    cur.execute(sql)
    x = cur.fetchone()
    if x != None:
        return x[0]
    else:
        msg = "An error occurred while importing into " + table + " " + colval.__str__()
        write_error(con, msg)
        print msg
        exit()

def import_fastq(con, fpath):
    return import_generic(con, fpath, "FastqFiles", colkeyword="filepath")

def import_species(con, name):
    return import_generic(con, name, "Species", colkeyword="name")

def import_condition(con, name):
    return import_generic(con, name, "Conditions", colkeyword="name")

def import_gene(con, name):
    return import_generic(con, name, "Genes", colkeyword="name")

def import_reads(con, name, fastqid, speciesid, conditionid, geneid, tagged):
    """Returns the Read ID upon success."""
    cur = con.cursor()
    sql = "insert or ignore into Reads (name, fastqid, speciesid, conditionid, geneid, tagged)"
    sql += " values('" + name + "',"
    sql += fastqid.__str__() + ","
    sql += speciesid.__str__() + ","
    sql += conditionid.__str__() + ","
    sql += geneid.__str__() + "," 
    sql += tagged.__str__() + ")"
    cur.execute(sql)
    con.commit()
    
    sql = "select id from Reads where name='" + name + "'"
    cur.execute(sql)
    x = cur.fetchone()
    if x != None:
        return x[0]
    else:
        msg = "An error occurred while importing the READ named " + name.__str__()
        write_error(con, msg)
        print msg
        exit()

#
# depricated
#
# def get_name_for_macs(exp_annoid, control_annoid, con):
#     cur = con.cursor()
#     sql = "select count(*) from HybridPairs where annoid1=" + exp_annoid.__str__() + " or annoid2=" + exp_annoid.__str__()
#     cur.execute(sql)
#     count = cur.fetchone()[0]
#     
#     name = ""
#     if count > 0:
#         """This data is hybrid"""
#         sql = "select library_name, species from Annotations where annoid=" + exp_annoid.__str__()
#         cur.execute(sql)
#         x = cur.fetchone()
#         name = x[0] + "-" + x[1]
#     else:
#         """Not a hybrid read."""
#         sql = "select library_name from Annotations where annoid=" + exp_annoid.__str__()
#         cur.execute(sql)
#         name = cur.fetchone()[0]
#     return name

def get_macs_pairs(con):
    """Returns a list of tuples, each being (experiment,control) Annotation object ids for MACS2 peak-calling.
    Each pair shares the same unique value in the fields 'tf', 'sample', 'species', and 'replicate' in the table Annotations.
    The control member of the pair is found by finding the Annotation with 'tag' = 0, 
    while the treatment member of the pair has 'tag'=1."""
    pairs = [] #(treatment Id, control ID)
    
    cur = con.cursor()
    sql = "select taggedid, controlid from Pairs"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        pairs.append( ii[0], ii[1] )
    return pairs
    
#     cur = con.cursor()
#     sql = "select distinct tf from Annotations"
#     cur.execute(sql)
#     for dd in cur.fetchall():
#         tf = dd[0]
#         sql = "select distinct sample from Annotations where tf='" + tf.__str__() + "'"
#         cur.execute(sql)
#         for x in cur.fetchall():
#             sample = x[0]
#             sql = "select distinct replicate from Annotations where sample='" + sample + "' and tf='" + tf.__str__() + "'"
#             cur.execute(sql)
#             for y in cur.fetchall():
#                 repid = y[0]     
#                 sql = "Select distinct species from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and  tf='" + tf.__str__() + "'"
#                 cur.execute(sql)
#                 s = cur.execute(sql)
#                 for z in cur.fetchall():
#                     species = z[0]
#                     sql = "select distinct media from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and  tf='" + tf.__str__() + "' and species='" + species.__str__() + "'"
#                     cur.execute(sql)
#                     mm = cur.fetchall()
#                     for media in mm:
#                         media = media[0]
#                         sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and media='" + media.__str__() + "' and tag=1 and tf='" + tf.__str__() + "'"
#                         cur.execute(sql)
#                         treatments = cur.fetchall()
#                         sql = "select annoid from Annotations where sample='" + sample + "' and replicate=" + repid.__str__() + " and species='" + species.__str__() + "' and media='" + media.__str__() + "' and tag=0 and tf='" + tf.__str__() + "'"
#                         cur.execute(sql)
#                         controls = cur.fetchall()               
#                         if controls.__len__() > 1:
#                             print "\n. Error, I wasn't expecting to find multiple controls for", sample, species, repid, tf
#                             print controls
#                             exit()
#                         control = controls[0][0]
#                         for t in treatments:
#                             pairs.append( (t[0],control)  ) 
#                             print "Experiment:", t[0], "Control:", control, "Sample:", sample, "ReplicateID:",  repid, "Species:", species
#     return pairs

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

#
# depricated
#
# def import_genome_list(gpath, con):
#     """Reads the user-specificed genome list. Each line should contain two things:
#     1. the species identified (a short string), and
#     2. the path to the genome
#     """
#     
#     cur = con.cursor()
#     
#     sql = "delete from SpeciesGenomepath"
#     cur.execute(sql)
#     con.commit()
#     
#     if False == os.path.exists(gpath.__str__()):
#         print "\n. Error: I cannot find  your genome list file at " + gpath
#         print os.getcwd()  
#         print os.path.exists(gpath)      
#         exit()
#     
#     fin = open(gpath, "r")
#     lines = fin.readlines()
#     fin.close()
#     for l in lines:
#         if l.startswith("#") or l.__len__() < 2:
#             continue
#         tokens = l.split()
#         if tokens.__len__() < 2:
#             continue
#         speciesname = tokens[0]
#         gpath = tokens[1]
#         #if False == os.path.exists(gpath):
#         #    print "\n. Error: your genome list includes a reference that I cannot find:"
#         #    print l
#         #    exit()
#         sql = "insert into SpeciesGenomepath (speciesid, genomepath) VALUES("
#         sql += "'" + speciesname + "','" + gpath + "')"
#         cur.execute(sql)
#         con.commit()
#         
#         print ". I found a genome for species", speciesname, "at", gpath
#         
#     return con
    
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
            """Restrict our import to a subset of rows, if the user specified it."""
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
            This allows for hybrid species to be separated into two annotations, one for each parent species."""
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
                    con.commit()           
    return con

# the configuration file is new for June 2015
def import_configuration(cpath, con):
    """Each line the configuration can start with one of the following keywords: 
        READS, PAIR, COMPARE, GENOMEPATH, BOWTIEINDEX, REDFLAGPATH
    
        The READS line has the following columns:
        * library name (unique for every FASTQ)
        * fastq path -- a filepath in the directory specified by the option "--datadir"
        * species name
        * tf name
        * tagged? YES/NO
        * condition/media
        * note (optional)
        
        The PAIR lines have these columns:
            <pair name> <experimental ID 1> <experimental ID 2>
            --> one of the experiments must be tagged YES, and the other must be tagged NO
            
        COMPARE <comparison name> = <ID 1> <ID 2> . . . <ID N>
            --> IDs can be pair names or comparison names 
    """
    if False == os.path.exists(cpath):
        print "\n. Error, I can't find your configuration file at", cpath
        exit()
    
    cur = con.cursor()
    
    print "\n. Reading the configuration file at ", cpath
    
    sql = "delete from SpeciesBowtieindex"
    cur.execute(sql)
    sql = "delete from SpeciesGenomepath"
    cur.execute(sql)
    sql = "delete from SpeciesGenomesize"
    cur.execute(sql)
    sql = "delete from GFF"
    cur.execute(sql)
    sql = "delete from IntergenicPaths"
    cur.execute(sql)
    sql = "delete from PairsComparisons"
    cur.execute(sql)
    sql = "delete from Comparisons"
    cur.execute(sql)
    sql = "delete from CompareComparisons"
    cur.execute(sql)
    sql = "delete from HybridPairs"
    cur.execute(sql)
    con.commit()
        
    fin = open(cpath, "r")
    """Get the lines with useful content."""
    lines = []
    for l in fin.readlines():
        ls = l.split("\r")
        """We split again on \r because files edited in OSX TextEdit use a different line
            break than normal Unix files."""
        for lt in ls:
            if lt.startswith("#"):
                continue
            if lt.__len__() < 2:
                continue
            tokens = lt.split()
            if tokens.__len__() < 2:
                continue
            lines.append( lt )
            
    """Parse lines for GENOMEPATH entries"""
    for ll in lines:
        if ll.startswith("GENOMEPATH"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            species = tokens[1]
            speciesid = import_species(con, species)
            gpath = tokens[2]
            sql = "insert or replace into SpeciesGenomepath (speciesid, genomepath) VALUES("
            sql += "'" + speciesid.__str__() + "','" + gpath + "')"
            cur.execute(sql)
            con.commit()
            #print ". I found a genome for species", speciesid, "at", gpath

    """Parse lines for BOWTIEINDEX entries"""
    for ll in lines:
        if ll.startswith("BOWTIEINDEX"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            species = tokens[1]
            speciesid = import_species(con, species)
            ipath = tokens[2]
            sql = "insert or replace into SpeciesBowtieindex (speciesid, indexpath) VALUES("
            sql += "'" + speciesid.__str__() + "','" + ipath + "')"
            cur.execute(sql)
            con.commit()
            #print ". I found a Bowtie index for species", speciesid, "at", ipath
    
    for ll in lines:
        if ll.startswith("REDFLAGPATH"):
            tokens = ll.split()
            if tokens.__len__() < 2:
                msg = "In your configuration file, this line doesn't make sense: " + ll
                write_error(con, msg)
                print msg
                exit()
            filepath = tokens[1]
            sql = "insert into Settings(keyword, value)"
            sql += " values('redflagpath','" + filepath + "')"
            cur.execute(sql)
            con.commit()
    
    for ll in lines:
        if ll.startswith("GENOMESIZE"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "In your configuration file, this line seems incomplete: " + ll
                write_error(con, msg)
                print msg
                exit()
            genomename = tokens[1]
            genomesize = int( tokens[2] )
            sql = "select id from Species where name='" + genomename + "'"
            cur.execute(sql)
            fetch = cur.fetchone()
            if fetch == None:
                msg = "Your configuration file includes the genome size for species " + genomename
                msg += ", but I cannot find that species in your database. Ensure that you included "
                msg += " a line for GENOMEPATH with the species " + genomename
                write_error(con, msg)
                print msg
                exit()
            speciesid = fetch[0]
            sql = "select count(*) from SpeciesGenomesize where speciesid = " + speciesid.__str__()
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count == 0:
                sql = "insert into SpeciesGenomesize(speciesid, genomesize) values("
                sql += speciesid.__str__() + "," + genomesize.__str__() + ")"
                cur.execute(sql)
                con.commit()
            
    """Parse lines for GFF entries"""
    for ll in lines:
        if ll.startswith("GFF"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()    
            species = tokens[1]
            speciesid = import_species(con, species)
            gffpath = tokens[2]
            sql = "insert or replace into GFF (speciesid, gffpath) VALUES("
            sql += "'" + speciesid.__str__() + "','" + gffpath + "')"
            cur.execute(sql)
            con.commit()
            #print ". I found an annotation for species", speciesid, "at", gffpath
    
    """Parse lines for GFF entries"""
    for ll in lines:
        if ll.startswith("INTERGENIC"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()    
            species = tokens[1]
            speciesid = import_species(con, species)
            intergenicpath = tokens[2]
            sql = "insert or replace into IntergenicPaths (speciesid, intergenicpath) VALUES("
            sql += "'" + speciesid.__str__() + "','" + intergenicpath + "')"
            cur.execute(sql)
            con.commit()

    """Parse lines for READS entries"""
    for ll in lines:
        if ll.startswith("READS"):
            tokens = ll.split()
            if tokens.__len__() < 7:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            libname = tokens[1]
            fastqpath = tokens[2]
            fastqid = import_fastq(con, fastqpath)
            species = tokens[3]
#             st = species.split("_") # split the species name, which will deal with hybrid species.
#             for species in st:
            speciesid = import_species(con, species)
            tf = tokens[4]
            tfid = import_gene(con, tf)
            tagged = tokens[5]
            if tagged not in ["YES", "NO", "yes", "no", "Y", "N", "y", "n"]:
                msg = "There is something wrong with the 'tagged?' column in the following configuration line: " + ll
                write_error(con, msg)
                print msg
                exit()
            elif tagged in ["YES", "yes", "Y", "y"]:
                tagged = 1
            elif tagged in ["NO", "no", "N", "n"]:
                tagged = 0
            condition = tokens[6]
            conditionid = import_condition(con, condition)
#             if st.__len__() > 1:
#                 """Special case for hybrid"""
#                 libname = libname + "-" + species
                    
                
            readid = import_reads(con, libname, fastqid, speciesid, conditionid, tfid, tagged)
            
            if tokens.__len__() > 7:
                note = tokens[7:]
                sql = "insert or replace into ReadComments (readid, comment) values(" + readid.__str__() + ",'" + (" ".join(note)) + "')"
                cur.execute(sql)
                con.commit()
                    
    """Parse lines for EXPERIMENT entries"""
    for ll in lines:
        if ll.startswith("EXPERIMENT"):
            tokens = ll.split()
            if tokens.__len__() < 5:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            pairname = tokens[1]
            liba = tokens[3] # the name of the first library
            libb = tokens[4] # the name of the second library
            shiftsize = 78 # the default shift size
            if tokens.__len__() > 5:
                shiftsize = tokens[5]
            ida = None # the Reads ID of the control library
            idb = None # the Reads ID of the experiment library
            for libname in [liba, libb]:
                sql = "select id, tagged from Reads where name='" + libname + "'"
                cur.execute(sql)
                x = cur.fetchone()
                if x == None:
                    msg = "I cannot find the library named " + libname + " in the line: " + ll
                    write_error(con, msg)
                    print msg
                    exit()
                if x[1] == 0:
                    ida = x[0]
                elif x[1] == 1:
                    idb = x[0]
            if ida == None or idb == None:
                msg = "Something is wrong with this EXPERIMENT line in your configuration: " + ll
                write_error(con, msg)
                print msg
                exit()
                
            """Validate that Reads ida and Reads idb come from the same species."""
            sql = "select speciesid from Reads where id=" + ida.__str__()
            cur.execute(sql)
            yy = cur.fetchone()
            sql = "select speciesid from Reads where id=" + idb.__str__()
            cur.execute(sql)
            zz = cur.fetchone()
            if yy[0] != zz[0]:
                msg = "Error: the controlled and tagged reads came from different species."
                print msg
                print ll
                exit()
                
            sql = "insert or ignore into Pairs (name, controlid, taggedid) VALUES"
            sql += " ('" + pairname.__str__() + "'," + ida.__str__() + "," + idb.__str__() + ")"
            cur.execute(sql)
            con.commit()
            
            sql = "SELECT last_insert_rowid()"
            cur.execute(sql)
            pairid = cur.fetchone()[0]
            
            sql = "insert or replace into Shiftsizes(pairid, shiftsize) values(" + pairid.__str__() + "," + shiftsize.__str__() + ")"
            cur.execute(sql)
            con.commit()
                
    """Parse lines for COMPARE entries"""
    for ll in lines:
        if ll.startswith("COMPARE"):
            tokens = ll.split()
            if tokens.__len__() < 5:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            name = tokens[1]
            sql = "insert or ignore into Comparisons (name) values('" + name.__str__() + "')"
            cur.execute(sql)
            con.commit()
            sql = "select id from Comparisons where name='" + name.__str__() + "'"
            cur.execute(sql)
            x = cur.fetchone()
            if x == None:
                msg = "An error occurred while parsing a COMPARE line in your configuration file: " + ll
                write_error(con, msg)
                print msg
                exit()
            compid = x[0]
                        
            for name in tokens[3:]:
                """Ensure the pair ID exists."""
                sql = "select id from Pairs where name='" + name.__str__() + "'"
                cur.execute(sql)
                targetid = None
                x = cur.fetchone()
                if x == None:
                    """Then maybe the comparison object is another comparison?"""
                    sql = "select id from Comparisons where name='" + name.__str__() + "'"
                    cur.execute(sql)
                    yy = cur.fetchall()
                    if yy.__len__() == 0:
                        msg = "There is something wrong with this COMPARE line: " + ll
                        write_error(con, msg)
                        print msg
                        exit()
                    if yy.__len__() > 3:
                        msg = "The ChIP-Seq Distillery does not yet support comparisons with more than three object."
                        msg += ": " + ll
                        write_error(con, msg)
                        print msg
                        exit()
                    targetid = int( yy[0][0] )
                    sql = "insert or replace into CompareComparisons (compid, targetid)"
                    sql += " values(" + compid.__str__() + "," + targetid.__str__() + ")"
                    cur.execute(sql)
                    con.commit()
                else:
                    targetid = int(x[0])
                    sql = "insert or replace into PairsComparisons (compid, pairid)"
                    sql += " values(" + compid.__str__() + "," + targetid.__str__() + ")"
                    cur.execute(sql)
                    con.commit()
                    
    """Parse lines for HYBRID entries"""
    for ll in lines:
        if ll.startswith("HYBRID"):
            tokens = ll.split()
            if tokens.__len__() < 3:
                msg = "This line in your configuration file doesn't have enough columns: " + ll
                write_error(con, msg)
                print msg
                exit()
            readname1 = tokens[1]
            readname2 = tokens[2]
            readids = []
            for rn in [readname1, readname2]:
                sql = "select id from Reads where name='" + rn.__str__() + "'"
                cur.execute(sql)
                x = cur.fetchall()
                if x == None:
                    msg = "The read name in this HYBRID line is causing problems:"
                    msg += ": " + ll
                    write_error(con, msg)
                    print msg
                    exit()
                if x.__len__() > 1:
                    msg = "The read name in this HYBRID line may not be unique, and is causing problems:"
                    msg += ": " + ll
                    write_error(con, msg)
                    print msg
                    exit()                    
                readids.append( x[0][0].__str__() )
            sql = "insert or replace into HybridPairs (readid1, readid2) values(" + ",".join(readids) + ")"
            cur.execute(sql)
            con.commit()
        
    return con 

def validate_configuration_import(con):
    print "\n. Summary of Analysis:"
    cur = con.cursor()
    
    sql = "select count(*) from Reads"
    cur.execute(sql)
    x = cur.fetchone()
    print "\t", x[0], "FASTQ files."
    
    sql = "select count(*) from HybridPairs"
    cur.execute(sql)
    x = cur.fetchone()
    print "\t", x[0], "pairs of reads from hybrid species."
    
    sql = "select id, name from Species"
    cur.execute(sql)
    count_good_species = 0
    x = cur.fetchall()
    for ii in x:
        sql = "select genomepath from SpeciesGenomepath where speciesid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y == None:
            msg = "Error, the species " + ii[1] + " doesn't have a GENOMEPATH entry."
            write_error(con, msg)
            print msg
            exit()
        if y.__len__() > 1:
            msg = "Error, the species " + ii[1] + " has multiple GENOMEPATH entries."
            write_error(con, msg)
            print msg
            exit()
            
        sql = "select indexpath from SpeciesBowtieindex where speciesid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y == None:
            msg = "Error, the species " + ii[1] + " doesn't have a BOWTIEINDEX entry."
            write_error(con, msg)
            print msg
            exit()
        if y.__len__() > 1:
            msg = "Error, the species " + ii[1] + " has multiple BOWTIEINDEX entries."
            write_error(con, msg)
            print msg
            exit()
        
        sql = "select genomesize from SpeciesGenomesize where speciesid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y == None:
            msg = "WARNING: the species " + ii[1] + " doesn't have a GENOMESIZE entry."
            print msg
        if y.__len__() > 1:
            msg = "Error, the species " + ii[1] + " has multiple GENOMESIZE entries."
            write_error(con, msg)
            print msg
            exit()
        
        sql = "select gffpath from GFF where speciesid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y == None:
            msg = "The species " + ii[1] + " doesn't have a GFF entry."
            write_error(con, msg)
            print msg
            exit()            
        if y.__len__() > 1:
            msg = "Error, the species " + ii[1] + " has multiple GFF entries."
            write_error(con, msg)
            print msg
            exit()
            
        """OK, at this point, the species has passed the tests."""
        count_good_species += 1
    print "\t", count_good_species, "annotated genomes"
    
    sql = "select count(*) from Pairs"
    cur.execute(sql)
    x = cur.fetchone()
    print "\t", x[0], "pairs of tagged/untagged experiments."
    
    sql = "select count(*) from Comparisons"
    cur.execute(sql)
    xx = cur.fetchone()
    print "\t", xx[0], "comparisons between experiments"


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
                
            