"""Read an annotation Excel file
See the NDT80 Analysis, October 2014, with Eugenio, for details."""

import os, sys, re
import sqlite3 as lite

def build_anno_db(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Annotation(annoid INTEGER primary key autoincrement, sample TEXT, library_name TEXT, index INT, fastqpath TEXT, strain TEXT, species TEXT, tf TEXT, tag INT, media TEXT, condition TEXT, replicate INT, comment TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFF(speciesid TEXT, gffpath TEXT)")
    con.commit()


def read_anno(apath):
    """Reads the annotation table and returns a Sqlite3 database object."""
    
    if False == os.path.exists(apath):
        print "\n. Error, I can't find your annotation file at", apath
        exit()
    
    dbpath = apath + ".db"
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath
    con = lite.connect(dbpath)
    build_anno_db(con)
    cur = con.cursor()
    
    fin = open(apath, "r")
    lines = []
    for l in fin.readlines():
        ls = l.split("\r")
        for lt in ls:
            lines.append( lt )
    
    for l in lines:
        tokens = l.split("\t")
        if tokens.__len__() > 2:
            if tokens[0] == sample:
                fastq = tokens[3]
                id = tokens[1]
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
                        tag = 0
                    else:
                        tag = 1
                     
                    sql = "INSERT INTO Annotation (sample, library_name, index, fastqpath, strain, species, tf, tag, media, condition, replicate, comment)"
                    sql += " VALUES(" + tokens[0] + "," + id + "," + int(tokens[2]).__str__() + "," + fastq
                    sql += "," + tokens[4] + "," + species + "," + tokens[6]
                    sql += "," + tag.__str__() + "," + tokens[8] + "," + tokens[9] + "," + int(tokens[10]).__str__()
                    sql += "," + tokens[11] 
                    sql += ")"
                    cur.execute(sql)
                    con.commit()
    
                    print id, fastq, species
    return con