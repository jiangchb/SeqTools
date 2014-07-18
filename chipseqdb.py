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
        cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INT, name TEXT, start INT, stop INT, chrom INT, strand TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Chromosome(id INT, name TEXT, species ID)")
        cur.execute("CREATE TABLE IF NOT EXISTS Species(id INT, name TEXT)")

        cur.execute("CREATE TABLE IF NOT EXISTS Peaks(id INT, name TEXT, start INT, stop INT, chrom INT, enrich FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Summit(id INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
        
        cur.execute("CREATE TABLE IF NOT EXISTS GenePeaks(gene INT, peak INT, distance INT)") # Dist is distance to the trans. start site
        cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INT, summit INT, distance INT)")
        
    return con
    
