############################################################
#
# A Python API / wrapper for the SQLLite ChIP-Seq database
#
############################################################

import sqlite3 as lite
import os, sys
from version import *
import traceback

from chipseqdb_api import *

#pillars

def print_db_stats(con):
    species = get_species(con)
    if species.__len__() > 0:
        print "\n. The database contains", species.__len__(), "species."
        for sp in species:
            print "    --> " + sp[1].__str__() + " (" + get_genes_for_species(con, sp[0]).__len__().__str__() + " genes)"
    return

def build_db(dbpath = None):
    """Initializes all the tables. Returns the DB connection object.
    If tables already exist, they will NOT be overwritten."""
        
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
        print "\n. Creating a new database at", dbpath
    else:
        print "\n. Restoring the existing database at", dbpath

    con = lite.connect(dbpath)

    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS About(version FLOAT)")
    cur.execute("DELETE from About")
    cur.execute("INSERT INTO About(version) VALUES(" + VERSION.__str__() + ")")
    
    # These data come from the GFF:
    cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key, name TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key, name TEXT COLLATE NOCASE, start INT, stop INT, chrom INT, strand TEXT)")
    cur.execute("create table if not exists Intergenics(id INTEGER primary key, chromid INT, start INT, stop INT, strand TEXT, geneid INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key, name TEXT, species INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFFs(id INTEGER primary key, species INT, filepath TEXT)")
    
    # This data comes from the pillars file:
    cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(realname TEXT COLLATE NOCASE, alias TEXT COLLATE NOCASE)")

    # This data from from the pillars file X GFF.
    cur.execute("CREATE TABLE IF NOT EXISTS GeneHomology(geneid INTEGER, aliasid INTEGER)")
    
    # genomic regions in which peaks should be ignored, for various reasons
    cur.execute("CREATE TABLE IF NOT EXISTS RedFlagRegions(chromid INT, start INT, stop INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key, name TEXT unique COLLATE NOCASE, species INT)")        
    
    # These data come from MACS2 output files
    cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SummitsEnrichment(summit INTEGER, max_enrichment FLOAT)")
    
    """GeneSummits maps summits to the their putative target gene."""
    cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INTEGER, distance INT)") # a mapping of Summits to nearby Genes    
    
    # The table Summits2Summits stores relationships about summits that are nearby, perhaps identical, between conditions and species.
    # Not all summits will have an entry in this table.
    # NOTE: The table Summits2Summits is not currently being used. Some work is required here.
    cur.execute("CREATE TABLE IF NOT EXISTS Summits2Summits(summitid1 INTEGER, summitid2 INTEGER, distance FLOAT)")
    
    """EnrichmentStats stores data about fold-enrichment (max, mean, etc.) for each gene."""
    cur.execute("CREATE TABLE IF NOT EXISTS EnrichmentStats(repid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT, maxenrichsite INT)")
    
    # The table Files contains paths to file generated during the analysis. These files will later be tar'd into
    # a package for results.
    cur.execute("CREATE TABLE IF NOT EXISTS Files(fileid INTEGER primary key, path TEXT, note TEXT)")
        
    cur.execute("create table if not exists Log(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    cur.execute("create table if not exists ErrorLog(id INTEGER primary key, time DATETIME DEFAULT CURRENT_TIMESTAMP,  message TEXT, code INT)")
    
    build_repgroups(con)
    build_unions(con)
    build_idr_tables(con)
    con.commit()
    
    if dbpath != None and dbpath != False:
        print_db_stats(con)
    
    return con

def clear_repgroups(con):
    cur = con.cursor()
    cur.execute("drop table RepgroupSummits")
    cur.execute("drop table GroupEnrichmentStats")
    con.commit()
    
def build_repgroups(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key, name TEXT COLLATE NOCASE, note TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GroupReplicate(rgroup INTEGER, replicate INTEGER, id INTEGER)")
    cur.execute("CREATE TABLE IF NOT EXISTS RepgroupSummits(repgroupid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT, mean_maxsummitscore FLOAT)") # maps summits that exist in all replicates in the group.
    cur.execute("CREATE TABLE IF NOT EXISTS GroupEnrichmentStats(rgroupid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, meanmaxenrich FLOAT, sumenrich FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicategroupFiles(repgroupid INTEGER, fileid INTEGER)")
    con.commit()

def clear_unions(con):
    cur = con.cursor()
    cur.execute("drop table UnionSummits")
    cur.execute("drop table UnionEnrichmentStats")
    con.commit()    

def build_unions(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key, name TEXT, type TEXT)") # defines a union set, type can be 'intersection' or 'union'
    cur.execute("CREATE TABLE IF NOT EXISTS UnionRepgroups(unionid INTEGER, repgroupid INTEGER, id INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS UnionSummits(unionid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT, mean_maxsummitscore FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, meanmaxenrich FLOAT, sumenrich FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionFiles(unionid INTEGER, fileid INTEGER)")
    con.commit()


def build_speciesunions(con):
    """Note: geneid in the Speciesunion tables point to translated gene IDs from the pillars.
        This means that to find this gene ID in a particular Union, you need to use
        the table GeneHomology to find alias gene IDs for another species."""
    
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Speciesunions(unionid INTEGER primary key, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionUnions(spunionid INTEGER, memunionid INTEGER, id INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionSummits(spunionid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT, mean_maxsummitscore FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, meanmaxenrich FLOAT, sumenrich FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionFiles(spunionid INTEGER, fileid INTEGER)")
    con.commit() 

def build_idr_tables(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS GeneRepgroupEnrichIdr(geneid INTEGER, repid1 INTEGER, repid2 INTEGER, lidr FLOAT, idr FLOAT)") # lidr is local IDR, idr is IDR for observations that are as irreproducible, or more, than this observation
    cur.execute("CREATE TABLE IF NOT EXISTS GeneUnionEnrichIdr(geneid INTEGER, repgroupid1 INTEGER, repgroupid2 INTEGER, lidr FLOAT, idr FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GeneSpeciesunionEnrichIdr(geneid INTEGER, unionid1 INTEGER, unionid2 INTEGER, lidr FLOAT, idr FLOAT)")
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


def import_gff(gffpath, speciesid, con, restrict_to_feature = "gene", filter_chrom = []):
    """filter_chrom is a list of chromosomes that should be ignored."""
    cur = con.cursor()
        
    print "\n. Importing genome features from", gffpath
    
    chromids = get_chrom_ids(con, speciesid)
    
    """Remove stale data."""
    sql = "DELETE FROM Chromosomes where species=" + speciesid.__str__()
    cur.execute(sql)
    try:
        for chromid in chromids:
            sql = "DELETE FROM Genes where chrom=" + chromid.__str__()
            cur.execute(sql)
    except:
        con.rollback()
    con.commit()
    
    chromname_id = {}
    
    count = 0
    total_count = estimate_line_count(gffpath)
    if total_count <= 1:          
        """If the line count is STILL 1 or 0, then something is wrong."""
        msg = "Something is wrong. There appears to be " + total_count.__str__() + " line(s) in the GFF "
        msg += gffpath.__str__()
        write_error(con, msg)
        print msg
        exit()
    
    fin = open(gffpath, "r")
    curr_chromname = None # the name of the last-seen chromosome.
    curr_chromid = None # the chromosome ID (from the table Chromosomes) of the last-seen chromosome.
    lines = fin.readlines()
    if lines.__len__() == 1:
        lines = lines[0].split("\r")
    try:
        for l in lines:
            count += 1
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
            
            if l.__len__() > 0 and False == l.startswith("#"):
                tokens = l.split()
                            
                if tokens.__len__() < 8:
                    print "\n. Error, something is wrong with the following line in the GFF:"
                    print l
                    exit()
                            
                if tokens[2] != restrict_to_feature: # e.g., if restrict_to_feature == "gene", then we'll only import genes.
                    continue
                chr = tokens[0]
                
                #print "179:", tokens[0], tokens[2]
                
                restrict = False
                for f in filter_chrom:
                    if False == chr.__contains__(f):
                        restrict = True
                if restrict == True:
                    #print "\n. Skipping", l
                    """Skip this entry"""
                    continue
    
                if chr != curr_chromname:
                    """Add (or potentially overwrite) this chromosome into the table Chromosomes."""
                    curr_chromname = chr
                    #sql = "SELECT COUNT(*) from Chromosomes where name='" + curr_chromname + "' and species=" + speciesid.__str__()
                    #cur.execute(sql)
                    #if cur.fetchone()[0] == 0:
                    if curr_chromname not in chromname_id:
                        sql = "SELECT COUNT(*) FROM Chromosomes where name='" + curr_chromname + "' and species=" + speciesid.__str__()
                        cur.execute(sql)
                        if cur.fetchone()[0] == 0:
                            """We've not seen this chrom yet... insert it!"""
                            sql = "INSERT INTO Chromosomes (name,species) VALUES('" + curr_chromname + "', " + speciesid.__str__() + ")"
                            cur.execute( sql )  
                            con.commit()
                            #print "209:", curr_chromname, speciesid 
                        """Get the row ID of the newly-added chromosome."""
                        cur.execute("SELECT id FROM Chromosomes WHERE name='" + curr_chromname + "' and species=" + speciesid.__str__())
                        curr_chromid = cur.fetchone()[0]
                        chromname_id[curr_chromname] = curr_chromid
                    """Remember the ID."""
                    curr_chromid = chromname_id[curr_chromname]
                
                #print "258", curr_chromname
                
                start = int( tokens[3] )
                stop = int( tokens[4] )
                strand = tokens[6]
                if strand == "-":
                    x = start
                    y = stop
                    start = y
                    stop = x
                #print "267a", tokens[8]
                #print "267b", tokens[8].split(";")[0]
                #print "267c", tokens[8].split(";")[0].split("=")
                gene = tokens[8].split(";")[0].split("=")[1] # orfName
                
                notetoks = tokens[8].split(";")
                if notetoks.__len__() > 2:
                    qrs = notetoks[2].split("=")
                    if qrs.__len__() > 1:
                        note = qrs[1]
                        if note.__contains__("orf"):
                            for t in note.split():
                                if t.startswith("orf"):
                                    gene = t # use this name instead of the orfName       
                if abs(stop - start) > 15000:
                    msg = "Warning, the gene named " + gene.__str__() + " is very long: " + abs(stop-start).__str__() +  " bp."
                    write_log(con, msg)
                    print "\n\t-->", msg
                                              
                sql = "INSERT INTO Genes (name, start, stop, chrom, strand) VALUES('" + gene.__str__() + "'," + start.__str__() + "," + stop.__str__() + "," + curr_chromid.__str__() + ",'" + strand + "')"
                cur.execute(sql) 
    except:
        msg = "An error occurred while parsing the GFF."
        write_error(con, msg)
        print msg
        traceback.print_exc()
        con.rollback()
        exit()
    fin.close()
    con.commit()
    
    sql = "delete from GFFs where filepath='" + gffpath + "'"
    cur.execute(sql)
    con.commit()
    
    sql = "insert or replace into GFFs (species, filepath) VALUES(" + speciesid.__str__() + ",'" + gffpath + "')"
    cur.execute(sql)
    con.commit() 
        
    print ""
    cur.execute("SELECT id,name from Chromosomes where species=" + speciesid.__str__())
    for ii in cur.fetchall():
        sql = "SELECT COUNT(*) from Genes where chrom=" + ii[0].__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        print "    --> Chromosome", ii[1], "contains", count, "genes."
    
    cur.execute("SELECT COUNT(*) FROM Genes")
    count_genes = cur.fetchone()[0]
    print "\n. The database contains", count_genes, "total genes."
    
    return con

def print_pillarsstats(con):
    cur = con.cursor()
    cur.execute("SELECT count(*) from GeneAlias")
    count = cur.fetchone()[0]
    print "    --> The pillars database contains", count, "entries."


def import_pillars(pillarspath, con):
    print "\n. Reading Pillars from", pillarspath
    
    """Remove any previous pillar definitions. We're going to rebuild the entire table."""
    cur = con.cursor()
    sql = "DELETE from GeneAlias"
    cur.execute(sql)
    con.commit()
    
    
    fin = open(pillarspath, "r")    
    count = 0
    try:
        for l in fin.xreadlines():
            if l.startswith("-") or l[0].isalnum():
                count += 1
                #if count%100 == 0:
                    #sys.stdout.write(".")
                    #sys.stdout.flush()            
                l = l.strip()
                tokens = l.split()
                orf_list = []
                for t in tokens:
                    if False == t.startswith("-"):
                        """Skip the --- entries."""
                        orf_list.append(t)
                if orf_list.__len__() > 1:
                    realname = orf_list[0]
                    for aliasname in orf_list:
                        sql = "INSERT INTO GeneAlias (realname, alias) VALUES('" + realname + "','" + aliasname + "')"
                        cur.execute(sql)
    except:
        print "\n. An error occurred while reading the pillars at", pillarspath
        print ". I'm rolling back all changes to the database."
        con.rollback()
    con.commit()
    fin.close()
    
    if count < 100:
        print "\n. Warning, your pillars file contained only", count, "rows, but I was expecting thousands of rows."
    else:
        print "\n. The pillars file contains", count, "gene sets."
    return con


def import_summits(summitpath, repid, con):
    """Reads a summit file (from MACS2, for example), and puts values into the Summits table."""
    
    print "\n. Importing summits from", summitpath,"for replicate", repid
    cur = con.cursor()
    
    if False == os.path.exists(summitpath):
        print "\n. Error: I cannot find the summit file", summitpath
        exit()
    
    """Clear existing DB entries that are stale."""
    sql = "DELETE FROM Summits where replicate=" + repid.__str__()
    cur.execute(sql)
    con.commit()
    
    """Get the species for this replicate."""
    sql = "select species from Replicates where id=" + repid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    speciesid = int( x[0] )
    #print "chipseqdb 398 - ", summitpath, repid, speciesid
    
    #
    # Build a library of summits
    #
    count = 0
    total_count = estimate_line_count(summitpath)
    chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
    fin = open(summitpath, "r")
    try:
        curr_chromname = None
        chromid = None
        
        for l in fin.xreadlines():
            count += 1
            if count%100==0:
                con.commit()
            
            if l.__len__() > 0 and False == l.startswith("#"):
                """If the summit line has content, and isn't a comment, then process it"""
                tokens = l.split()
                chrname = tokens[0]                    
                site = int( tokens[1] )
                name = tokens[3]
                score = float( tokens[4] )
                                
                if curr_chromname == None or curr_chromname != chrname:
                    curr_chromname = chrname
                    chromid = get_chrom_id(con, curr_chromname, speciesid, make_if_missing=True)
                
                """Verify if this fold-enrichment data is in a red-flag region."""
                sql = "select count(*) from RedFlagRegions where chromid=" + chromid.__str__()
                sql += " and (start < " + site.__str__() + " and stop > " + site.__str__() + ")"
                cur.execute(sql)
                countrf = cur.fetchone()[0]
                if countrf > 0:
                    print "Skipping the summit inside a red flag region:", chrname, site, name, score
                    continue
                
                
                if chromid != None:                    
                    sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chromid.__str__() + "," + score.__str__() + ")"
                    cur.execute(sql)

    except:
        traceback.print_exc()
        print "\n. Error (357) Something wrong occurred while importing summits from", summitpath
        print ". I'm rolling back the changes to the tables Summits"
        con.rollback()
        exit()
    fin.close()
    con.commit()
    
    cur.execute("SELECT COUNT(*) FROM Summits where replicate=" + repid.__str__())
    count_genes = cur.fetchone()[0]
    print "\t--> Found", count_genes, "summits in", summitpath
    return con


def blocks(files, size=65536):
    """THis is a helper method for estimate_line_count"""
    while True:
        b = files.read(size)
        if not b: break
        yield b

def estimate_line_count(filepath):
    """Returns an estimation of the number of lines in the file"""
    with open(filepath, "r") as f:
        return sum(bl.count("\n") for bl in blocks(f))

def convert_line_breaks(filepath):
    """Converts all \r linebreaks to \n linebreaks"""
    outlines = []
    fin = open(filepath, "r")
    for l in fin.xreadlines():
        ls = l.split("\r")
        for s in ls:
            outlines.append(s)
    fin.close()
    
    fout = open(filepath, "w")
    for l in outlines:
        fout.write(l + "\n")
    fout.close()

def import_foldenrichment(bdgpath, repid, con):
    """Imports a BDG file with enrichment scores tracked across genome sequence sites.
    The BDG must be sorted by site for each chromosome, ascending.
    For each FE value in the file, this method attempts to map that FE
    to any known summits, and to nearby genes.
    Upon success, this method inserts data into two tables: EnrichmentStats and SummitsEnrichment
    """ 
    print "\n. Importing fold-enrichment values from", bdgpath,"for replicate", repid
    cur = con.cursor()
    
    if False == os.path.exists(bdgpath):
        print "\n. Error: I cannot find the BedGraph file", bdgpath
        exit()
    
    """Remove any stale entries."""
    sql = "DELETE FROM EnrichmentStats where repid=" + repid.__str__()
    cur.execute(sql)
    con.commit()
    sql = "DELETE FROM SummitsEnrichment where summit in (select id from Summits where replicate=" + repid.__str__() + ")"
    cur.execute(sql)
    con.commit()
    
    """Get the species for this replicate."""
    sql = "select species from Replicates where id=" + repid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    speciesid = int( x[0] )
    
    total_count = estimate_line_count(bdgpath)
    
    """Open the BDG file"""
    fin = open(bdgpath, "r")
    curr_chromname = None # as we parse the BDG, this string will store the name of the chromosome we are currently evaluating
    curr_chromid = None   # the ID into the Chromosomes table for the chrom. we are currently evaluating
    last_start_site = 0 # an index used while parsing the BDG; the value of the start site of the last enrichment window
    
    """As we parse the BDG, we'll map enrichment values to genes,
        and store the following information. . . """
    geneid_sum = {} # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
    geneid_n = {} # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
    geneid_max = {}
    geneid_maxsite = {} # key = geneid, value = the distance from TSS of the maximum FE in the gene's regulatory region
    genes = None # an ordered list of gene objects, 
                    # it will be filled with data whenever we 
                    # encounter a new chromosome in the BDG file

    count = 0 # a count of how many lines we've analyzed -- used for displaying a progress indicator

    chromid_summitsites = {} # key = chromosome ID, value = list of summit sites on that chrom.
    count_found_summits = 0  # how many summits have we successfully matched with FE values?
    
    """
    For each line in the BDG file:
    """
    for l in fin.xreadlines():
        if l.__len__() <= 5:
            """Skip abnormally short lines."""
            continue
        
        """Print a progress indicator."""
        count += 1
        if count%10 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) + "     (" + count_found_summits.__str__() + " summits)" )
            sys.stdout.flush()
           
        """Split the line into pieces"""         
        tokens = l.split()
        
        """Parse the line"""
        chromname = tokens[0]
        festart = int(tokens[1]) # start site of this enrichment window
        festop = int(tokens[2])  # stop site of this enrichment window
        eval = float(tokens[3])  # enrichment value across this window
        
        """Is the chromosome in this line the same chromosome from the previous line?
            If not, then we need to retrieve information about this chromosome, including
            the genes on this chrom."""
        if curr_chromname == chromname:
            chromid = curr_chromid
        else:
            """We need to look up the ID for this chromosome,
                and get this chrom's list of gene pairs and summits."""
            curr_chromname = chromname
            chromid = get_chrom_id(con, curr_chromname, speciesid, make_if_missing = True)
            
            """Reset the red flag index pointer"""
            last_rfptr = 0
                        
            if chromid == None:
                """We don't know anything about this chromosome; skip to the next FE window."""
                msg = "The chromosome named " + curr_chromname.__str__() + " exists in your FE file, but not in your GFF."
                write_error(con, msg)
                continue
            
            curr_chromid = chromid
            
            if curr_chromid not in chromid_summitsites:
                """Get the list of summits for this chromosome."""
                chromid_summitsites[curr_chromid] = {}
                summits = get_summits(con, repid, chromid)
                for ss in summits:
                    chromid_summitsites[curr_chromid][ ss[3] ] = ss
            
            last_start_site = 0
            print "\n\t", curr_chromname

        """Check if this fold-enrichment data is in a red-flag region."""
        sql = "select count(*) from RedFlagRegions where chromid=" + chromid.__str__()
        sql += " and ( (start < " + festart.__str__() + " and stop > " + festart.__str__() + ")"
        sql += "      or (start < " + festop.__str__() + " and stop > " + festop.__str__() + ") )"
        cur.execute(sql)
        countrf = cur.fetchone()[0]
        if countrf > 0:
            print "Skipping FE data inside a red flag region:", chromname, festart, festop
            continue
  
        """Check for discontinuous data in the BDG file"""
        if last_start_site < festart and last_start_site != 0:
            msg = "Warning: the BDG file may skip some sites, at site: " + festart.__str__() + " for chrom " + curr_chromname.__str__() + " for BDG " + bdgpath.__str__()
            write_log(msg)
            print msg
        last_start_site = festop
    
        """Can we map this enrichment site to a summit?"""        
        for fesite in xrange(festart, festop):
            
            """For each site in the enrichment window, is there a known summit at this site?"""
            if fesite in chromid_summitsites[curr_chromid]:              
                summitid = chromid_summitsites[curr_chromid][fesite][0] 
                count_found_summits += 1
                sql = "insert into SummitsEnrichment (summit, max_enrichment) "
                sql += " VALUES(" + summitid.__str__() + ","
                sql += eval.__str__() + ")"
                cur.execute(sql)
                con.commit()
        
        """ Can we amp this enrichment values to a gene?"""
        for fesite in xrange(festart, festop):            
            sql = "select id, geneid from Intergenics where chromid=" + chromid.__str__()
            sql += " and start > " + fesite.__str__() + " and stop < " + fesite.__str__()
            cur.exectue(sql)
            x = cur.fetchall()
            for ii in x:
                this_intergenic_id = ii[0]
                this_gene_id = ii[1]
                if this_gene_id not in geneid_sum:
                    geneid_sum[ this_gene_id ] = 0.0
                geneid_sum[ this_gene_id ] += eval
                
                if this_gene_id not in geneid_max:
                    geneid_max[ this_gene_id ] = eval
                    geneid_maxsite[ this_gene_id ] = fesite
                elif geneid_max[ this_gene_id ] < eval:
                    geneid_max[ this_gene_id ] = eval
                    geneid_maxsite[ this_gene_id ] = fesite   
                
                if this_gene_id not in geneid_n:
                    geneid_n[this_gene_id] = 0
                geneid_n[this_gene_id] += 1
    fin.close()
        
    """Finally, write all our findings into the table EnrichmentStats."""
    count = 0
    total_count = geneid_sum.__len__()
    try:
        for geneid in geneid_sum:
            if geneid_n[geneid] > 0:
                count += 1
                if count%5000==0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                    con.commit()
                sql = "INSERT INTO EnrichmentStats (repid, geneid, maxenrich, meanenrich, sumenrich, maxenrichsite)  "
                sql += "VALUES(" + repid.__str__() + "," + geneid.__str__()
                sql += "," + geneid_max[geneid].__str__()
                sql += "," + (geneid_sum[geneid]/float(geneid_n[geneid])).__str__()
                sql += "," + geneid_sum[geneid].__str__()
                sql += "," + geneid_maxsite[geneid].__str__()
                sql += ")"
                cur.execute(sql)
    except:
        print "\n. An error occurred updating the table EnrichmentStats."
        print ". I'm rolling back the changes."
        con.rollback()
    con.commit()        

    return con

def validate_summits_fe(repid, con):
    """This method validates the data in the tables SummitsEnrichment and EnrichmentStats."""    
       
    print "\n. Validating Summits and FE values for replicate", repid
    
    fout = open("summit_validation.log", "w")
    
    cur = con.cursor()
    speciesid = get_speciesid_for_repid(repid, con)
    if speciesid == None:
        print "\n. An error occurred, I cannot find the species for replicate", repid
        exit()
    
    """Look for any summits that don't have a corresponding FE value."""
    sql = "select count(*) from Summits where replicate=" + repid.__str__()
    sql += " and id not in (select summit from SummitsEnrichment)"
    cur.execute(sql)
    count_missing = int( cur.fetchone()[0] )
    if count_missing > 0:
        print "\n. Error: I could not find FE values for", count_missing, "of the summits:"

        sql = "select id, site, chrom from Summits where replicate=" + repid.__str__()
        sql += " and id not in (select summit from SummitsEnrichment)"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            print "Summit", ii[0], "at site:", ii[1], "on chrom:", get_chrom_name(con, ii[2])

        sql = "select name from Chromosomes where id in ("
        sql += "select distinct chrom from Summits where replicate=" + repid.__str__()
        sql += " and id not in (select summit from SummitsEnrichment)"        
        sql += " ) "
        cur.execute(sql)
        x = cur.fetchall()
        print ". Overall, there are missing summits on the following chromosomes:", x
        exit()
    
    """Look for any genes that lack a max FE value."""
    sql = "select * from Genes where "
    sql += "id in (select id from Genes where chrom in (select id from Chromosomes where species=" + speciesid.__str__() + "))"
    sql += "and id not in (select geneid from EnrichmentStats where repid=" + repid.__str__() + ") "
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() > 0:
        print "\n. Warning: I found no FE values for", x.__len__().__str__(), "genes."
        for ii in x:
            msg = "I found no FE values for gene " + ii[0].__str__() + ", replicate " + repid.__str__()
            write_error(con, msg)
    
    fout.close()
    

def resolve_aliasids(con):
    """This method inserts data into the table GeneHomology."""
    print "\n. Indexing homologous gene IDs between species."
    
    """Remove any prior GeneHomology entires. We're going to rebuild the entire table."""
    cur = con.cursor()
    sql = "DELETE from GeneHomology"
    cur.execute(sql)
    con.commit()
    
    """Build a dictionary of gene names to gene IDs."""
    genename2id = {}
    sql = "SELECT * FROM Genes"
    cur.execute(sql)
    for ii in cur.fetchall():
        name = ii[1]
        id = ii[0]
        genename2id[ name ] = id
    
    sql = "SELECT count(*) from Genes"
    cur.execute(sql)
    total_count = cur.fetchone()[0]
    
    count = 0
    sql = "SELECT * from Genes"
    cur.execute(sql)
    try:
        for g in cur.fetchall():
            count += 1
            if count%5 == 0:
                sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                sys.stdout.flush()
            if count%2000 == 0:
                """Note the commit happens here"""
                con.commit()
            
            this_id = g[0]
            this_name = g[1]
            
            """If there is an alias reference for this gene, then get the reference name for this gene."""
            realname = get_genename_for_aliasname( this_name, con )
            if realname == None:
                realname = this_name
                realid = this_id
            
            if realname in genename2id:
                realid = genename2id[realname]
                    
            sql = "INSERT INTO GeneHomology (geneid, aliasid) VALUES(" + realid.__str__() + "," + this_id.__str__() + ")"
            cur.execute(sql)
    except:
        con.rollback()
    
    """Save the commit for the end."""
    con.commit()

def get_geneorder(con, chromid):
    """Returns a list of tuples, each tuple containing two adjacent genes.
        In each tuple, the 0th element is a 5' gene and the 1st element is a 3' gene.
        This list of tuples can be iterated over in order to examine all intergenic regions
        on a chromosome.
        The genes on the ends of the chromosome are paired with None, i.e. (None, gene) for the 5'
        and (gene, None) for the 3'.
    """
    genes = get_genes_for_chrom(con, chromid)
    if genes.__len__() == 0:
        return []
                    
    """Build a list of gene pairs, corresponding to upstream-downstream neighbors
    in the genome. We'll iterate through these pairs in order to map FE values
    to intergenic regions."""
    genepairs = []
    for ii in xrange(0, genes.__len__()):                    
        if ii == 0:
            pair = (None,ii)
            genepairs.append( pair )
        elif ii == genes.__len__()-1:
            pair = (ii,None)
            genepairs.append( pair )
            pair = (ii-1,ii)
            genepairs.append( pair )
        else:
            pair = (ii-1,ii)
            genepairs.append( pair )                            
    return genepairs


def get_genes4site(con, repid, site, chromid, speciesid=None):
    """Given a site on a chromosome, this method finds the nearest upstream and downstream gene.
        Returns a tuple (upstream gene ID, distance to upstream TSS, downstream gene ID, distance to downstream TSS), 
        If a nearby gene cannot be found, then the gene IDs and the distance values will be None."""
    cur = con.cursor()
    
    if speciesid == None:
        sql = "select species from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x == None:
            print "\n. Error 913: there is no species defined for replicate ID", repid
            exit()
        speciesid = int( x[0] )
    
    try:
        """Iterate through all the genes on this chromosome and find the closest upstream
            and downstream gene from the given 'site'."""
        genes = get_genes_for_chrom(con, chromid)        
        min_up = None       # the minimum, so far, distance to an upstream TSS from the given site
        closest_up = None
        min_down = None     # the minimum, so far, distance to the a downstream TSS from the given site
        closest_down = None
        for g in genes:            
            gid = g[0]          # gene ID
            start = g[2]        # gene start site
            stop = g[3]         # gene stop site
            d = start - site    # distance to the TSS from the given site
            
            """Sense direction, and upstream"""
            if start < stop and start >= site:
                if min_up == None:
                    min_up = d
                    closest_up = gid
                if min_up > d:
                    min_up = d
                    closest_up = gid
                    
            elif start < stop and stop < site:
                """Sense direction and downstream"""
                if min_down == None:
                    min_down = d
                    closest_down = None
                if min_down < d:
                    min_down = d
                    closest_down = None
                
            elif start > stop and start <= site:
                """Antisense and downstream"""
                if min_down == None:
                    min_down = d
                    closest_down = gid
                if min_down < d: # be careful here, we're looking for the largest NEGATIVR number.
                    min_down = d
                    closest_down = gid
                    
            elif start > stop and stop > site:
                """Antisense and upstream"""
                if min_up == None:
                    min_up = d
                    closest_up = None
                if min_up > d:
                    min_up = d
                    closest_up = None
        return (closest_up, min_up, closest_down, min_down)
    except:
        print "\n. An error occurred (968) while looking for a nearby gene."
        exit()
    return (None, None, None, None)
    
    
def map_intergenic_regions(con, speciesid, intergenic_path):
    """This methods fills the DB table named Intergenics
        by reading the intergenic definitions in the file at
        at the path intergenic_path""" 
    
    print "\n. Mapping intergenic regions for species", get_species_name(speciesid, con), " using intergenic file", intergenic_path
      
    cur = con.cursor()
   
    fin = open(intergenic_path, "r")
    for l in fin.xreadlines():
        tokens = l.split()
        if tokens.__len__() < 5:
            continue
        chrom_name = tokens[0]
        chromid = get_chrom_id(con, chrom_name, speciesid, make_if_missing = False)
        
        start = int(tokens[1])
        stop = int(tokens[2])
        gene_name = tokens[3]
        strand = tokens[4]
    
        geneid = None
        sql = "select id from Genes where name='" + gene_name + "' and chrom=" + chromid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x == None or x.__len__() == 0:
            print "\n. I cannot find a Gene entry for the gene named", gene_name, "in your intergenic file."
            continue
        geneid = x[0]
        
        sql = "insert or replace into Intergenics (chromid, start, stop, strand, geneid)"
        sql += " values(" + chromid.__str__() + "," + start.__str__() + ","
        sql += stop.__str__() + ",'" + strand.__str__() + "'," + geneid.__str__()
        sql += ")"
        cur.execute(sql)
        con.commit()
           
    chromids = get_chrom_ids(con, speciesid)
    for chromid in chromids:       
        sql = "select count(*) from Intergenics where chromid=" + chromid.__str__()
        cur.execute(sql)
        count_inserted = cur.fetchone()[0]
        print ". Found", count_inserted, "intergenic regions for chromosome", get_chrom_name(con, chromid)
            
            
def map_summits2genes(con, repid, speciesid=None, chroms=None):
    """This methods puts values into the table GeneSummits.
        speciesid can be the integer ID of a species whose genome will be used to map summits.
        chroms can be a list of chromosome IDs to restrict the mapping (i.e. map summits only to
            genes on these chromosomes)."""

    cur = con.cursor()    
    print "\n. Mapping summits to nearby genes. . ."
    if speciesid == None:
        sql = "select species from Replicates where id=" + repid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        if x == None:
            print "\n. Error 857: there is no species defined for replicate ID", repid
            exit()
        speciesid = int( x[0] )
    if chroms == None:
        chroms = get_chrom_ids(con, speciesid)
    
    for chrid in chroms:
        summits = get_summits(con, repid, chrid)
            
        for s in summits:
            sid = s[0] # summit ID
            sumsite = s[3] # summit site in the genome
            score = s[5] # summit score
            
            summit_name = None
            sql = "select name from Summits where id=" + sid.__str__()
            cur.execute(sql)
            summit_name = cur.fetchone()[0]
            
            sql = "Select id, geneid, strand, start, stop from Intergenics where chromid=" + chrid.__str__()
            sql += " and start >" + sumsite.__str__() + " and stop < " + sumsite.__str__()
            cur.execute(sql)
            x = cur.fetchall()
            for ii in x:
                this_intergenic_id = ii[0]
                this_geneid = ii[1]
                this_strand = ii[2]
                this_start = ii[3]
                this_stop = ii[4]
                
                this_genename = get_genename(this_geneid, con)
                
                if this_strand == "+":
                    distance = this_start - sumsite
                elif this_start == "-":
                    distance = sumsite - this_stop
            
                sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                sql += " VALUES(" + this_geneid.__str__() + "," 
                sql += sid.__str__() + ","
                sql += distance.__str__() + ") "         
                cur.execute(sql) 
                
                print "\n. Mapping summit", summit_name, " to gene", this_genename

    con.commit()
    return con
                
     