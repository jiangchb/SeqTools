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
    cur.execute("create table if not exists Intergenics(id INTEGER primary key, downstreamgeneid INT, upstreamgeneid INT, chromid INT, start INT, stop INT)")
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
    print "chipseqdb 398 - ", summitpath, repid, speciesid
    
    #
    # Build a library of summits
    #
    count = 0
    total_count = estimate_line_count(summitpath)
    chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
    fin = open(summitpath, "r")
    try:
        for l in fin.xreadlines():
            count += 1
            #sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            #sys.stdout.flush()
            if count%100==0:
                con.commit()
            
            if l.__len__() > 0 and False == l.startswith("#"):
                """If the summit line has content, and isn't a comment, then process it"""
                tokens = l.split()
                chr = tokens[0]
                site = int( tokens[1] )
                name = tokens[3]
                score = float( tokens[4] )
                chrid = get_chrom_id(con, chr, speciesid)
                
                if chrid != None:                    
                    sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chrid.__str__() + "," + score.__str__() + ")"
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
    print "chipseqdb.py 496 - ", repid, speciesid
    
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
        
    chromid_genepairs = {} # the pair of genes before and after this enrichment window
    pairi = 0 # working index into chromid_genepairs[curr_chromid]
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
            print "554 - ", curr_chromname, speciesid
            chromid = get_chrom_id(con, curr_chromname, speciesid)
            if chromid == None:
                """We don't know anything about this chromosome; skip to the next FE window."""
                msg = "The chromosome named " + curr_chromname.__str__() + " exists in your FE file, but not in your GFF."
                write_error(con, msg)
                print msg
                print "559:", curr_chromname, speciesid
                print l
                continue
            curr_chromid = chromid
            if curr_chromid not in chromid_genepairs:
                """Get the list of ordered gene pair tuples for this chrom."""
                chromid_genepairs[curr_chromid] = get_geneorder(con, repid, curr_chromid)
                pairi = 0 # reset the pair index
                genes = get_genes_for_chrom(con, chromid)
                if genes.__len__() == 0:
                    msg = "There are no genes on chromosome " + chromid.__str__()
                    print msg
                    write_error(msg)
                    return None
            if curr_chromid not in chromid_summitsites:
                """Get the list of summits for this chromosome."""
                chromid_summitsites[curr_chromid] = {}
                summits = get_summits(con, repid, chromid)
                for ss in summits:
                    chromid_summitsites[curr_chromid][ ss[3] ] = ss
            last_start_site = 0
            print "\n\t", curr_chromname
                                         
        """Check for discontinuous data in the BDG file"""
        if last_start_site < festart and last_start_site != 0:
            msg = "Warning: the BDG file may skip some sites, at site: " + festart.__str__() + " for chrom " + curr_chromname.__str__() + " for BDG " + bdgpath.__str__()
            write_log(msg)
            print msg
        last_start_site = festop
    
        """Can we map this enrichment site to a summit?"""        
        for fesite in range(festart, festop):

#             if fesite == 1720362:
#                 print "566:", l
#                 if fesite in chromid_summitsites[curr_chromid]:
#                     print "Summit on chrom", curr_chromid
#                     print chromid_summitsites[curr_chromid][fesite]
#                 else:
#                     print "No summit on chrom", curr_chromid

#             # A test of some genes we know about:
#             test_flag1 = False
#             if curr_chromname.__contains__( "Chr_1" ):
#                 if ii < 221489 and ii >= 221489 - 30:
#                     #81|PICST_37571|221489|223006|1|+
#                     print l
#                     test_flag1 = True
#                 elif ii < 1101228 and ii >= 1101228 - 30:
#                     #409|PICST_28392|1101228|1103089|1|+
#                     print l
#                     test_flag1 = True

            """For each site in the enrichment window, is there a known summit at this site?"""
            if fesite in chromid_summitsites[curr_chromid]:              
                summitid = chromid_summitsites[curr_chromid][fesite][0] 
                count_found_summits += 1
                sql = "insert into SummitsEnrichment (summit, max_enrichment) "
                sql += " VALUES(" + summitid.__str__() + ","
                sql += eval.__str__() + ")"
#                 if fesite == 1720362:
#                     print "593:", sql
                cur.execute(sql)
                con.commit()
           
           #361
                  
            """If the current enrichment window ('festart') is outside the intergenic region defined by the
                current gene pair, then we need to advance to the next gene pair.
            """
            this_gene_pair = chromid_genepairs[ curr_chromid ][pairi]
            while this_gene_pair[1] != None and (genes[ this_gene_pair[1] ][2] < festart and genes[ this_gene_pair[1] ][3] < festart):
                pairi += 1
                this_gene_pair = chromid_genepairs[ curr_chromid ][pairi]
                
            """Can we map enrichment to both upstream and downstream genes?"""
            down_ok = False # is there a downstream gene?
            up_ok = False   # is there an upstream gene?
            ups_ii = chromid_genepairs[curr_chromid][pairi][0] # the ID of the downstream gene
            down_ii = chromid_genepairs[curr_chromid][pairi][1]   # the ID of the upstream gene
            
            if ups_ii != None:
                if genes[ups_ii][2] < festart and genes[ups_ii][3] < festart:
                    if genes[ups_ii][5] == "-":
                        """Yes, map scores to the downstream gene."""
                        up_ok = True
            if down_ii != None:
                if genes[down_ii][2] > festart and genes[down_ii][3] > festart:
                    if genes[down_ii][5] == "+":
                        """Yes, map scores to the upstream gene."""
                        down_ok = True
    
            #if test_flag1 == True:
            #    print ". 643 - ", up_ok, ups_ii, genes[ups_ii][0], down_ok, down_ii, genes[down_ii][0]
    
            if up_ok:     
                geneid = genes[ups_ii][0]
                """Initialize some data structures about this gene."""
                if geneid not in geneid_sum:
                    geneid_sum[geneid] = 0
                if geneid not in geneid_n:
                    geneid_n[geneid] = 0
                if geneid not in geneid_max:
                    geneid_max[geneid] = -1
                if geneid not in geneid_maxsite:
                    geneid_maxsite[geneid] = 0
                
                geneid_sum[geneid] += eval
                geneid_n[geneid] += 1
                
                """Is this fe value larger than we've seen before?"""
                if eval > geneid_max[geneid]:
                    geneid_max[geneid] = eval
                    geneid_maxsite[geneid] = fesite - genes[ups_ii][2]
            
            if down_ok:     
                geneid = genes[down_ii][0]
                """Initialize some data structures about this gene."""
                if geneid not in geneid_sum:
                    geneid_sum[geneid] = 0
                if geneid not in geneid_n:
                    geneid_n[geneid] = 0
                if geneid not in geneid_max:
                    geneid_max[geneid] = -1
                if geneid not in geneid_maxsite:
                    geneid_maxsite[geneid] = 0

                geneid_sum[geneid] += eval
                geneid_n[geneid] += 1
                
                """Is this fe value larger than we've seen before?"""
                if eval > geneid_max[geneid]:
                    geneid_max[geneid] = eval
                    geneid_maxsite[geneid] = genes[down_ii][2] - fesite
                        
#             if fesite == 1720362:
#                 print "666:", up_ok, ups_ii, down_ok, down_ii
            # The following print statement is too noisy:
            #if up_ok == False and down_ok == False:
            #    print "\n. FE data at site", ii, "doesn't map to any regulatory regions. ( Chrom:", curr_chromname, ")"
            
    fin.close()
        
    """Finally, write all our findings into the table EnrichmentStats."""
    count = 0
    total_count = geneid_sum.__len__()
    try:
        for geneid in geneid_sum:
            if geneid_n[geneid] > 0:
                count += 1
                if count%5000==0:
                    #sys.stdout.write(".")
                    #sys.stdout.flush()
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

def get_geneorder(con, repid, chromid):
    """Returns a list of tuples, each tuple containing two adjacent genes.
        In each tuple, the 0th element is a 5' gene and the 1st element is a 3' gene.
        This list of tuples can be iterated over in order to examine all intergenic regions
        on a chromosome.
        The genes on the ends of the chromosome are paired with None, i.e. (None, gene) for the 5'
        and (gene, None) for the 3'.
    """
    genes = get_genes_for_chrom(con, chromid)
    if genes.__len__() == 0:
        msg = "There are no genes on chromosome " + chromid.__str__()
        print msg
        write_error(msg)
        return None
                    
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
    
def map_intergenic_regions(con, repid, speciesid=None, chroms=None):
    """This methods fills the DB table IntergenicRegions""" #(id INTEGER primary key, downstreamgeneid INT, upstreamgeneid INT, chromid INT, start INT, stop INT)
    cur = con.cursor()
        
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
    count = 0    
    #try:
    for chrid in chroms:
        genes = get_genes_for_chrom(con, chrid) #genes is a list = id, name, start, stop, chrom, strand
        genepairs = get_geneorder(con, repid, chrid)
        for count, pair in enumerate(genepairs):
            leftgeneid = pair[0]
            rightgeneid = pair[1]
            downstreamgeneid = "'NULL'"
            upstreamgeneid = "'NULL'"
            start = 1
            stop = "'NULL'"
            if leftgeneid != None:
                leftgene = genes[leftgeneid]
                start = max( leftgene[2], leftgene[3])
                if leftgene[5] == "-" and leftgene[2] < leftgene[3]:
                    downstreamgeneid = leftgeneid.__str__()
                
            if rightgeneid != None:
                rightgene = genes[rightgeneid]
                stop = max(rightgene[2], rightgene[3])
                if rightgene[5] == "+" and rightgene[3] > rightgene[2]:
                    upstreamgeneid = rightgeneid.__str__()
                
            sql = "insert or replace into Intergenics (downstreamgeneid, upstreamgeneid, chromid, start, stop)"
            sql += " values(" + downstreamgeneid + "," + upstreamgeneid + ","
            sql += chrid.__str__() + "," + start.__str__() + "," + stop.__str__()
            sql += ")"
            cur.execute(sql)
        con.commit()
            
            
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
    
    count = 0    
    #try:
    for chrid in chroms:
        genes = get_genes_for_chrom(con, chrid)
        summits = get_summits(con, repid, chrid)
        
        genepairs = get_geneorder(con, repid, chrid)
        pairi = 0
        
        for s in summits:
            sid = s[0] # summit ID
            sumsite = s[3] # summit site in the genome
            score = s[5] # summit score
            
            this_gene_pair = genepairs[pairi]
            while this_gene_pair[1] != None and (genes[ this_gene_pair[1] ][2] < sumsite and genes[ this_gene_pair[1] ][3] < sumsite):
                pairi += 1
                this_gene_pair = genepairs[pairi]
        
            """Can we map enrichment to both upstream and downstream genes?"""
            down_ok = False # is there a downstream gene?
            up_ok = False   # is there an upstream gene?
            up_ii = this_gene_pair[0] # the ID of the downstream gene
            down_ii = this_gene_pair[1]   # the ID of the upstream gene
            
            if up_ii != None:
                if genes[up_ii][2] < sumsite and genes[up_ii][3] < sumsite:
                    if genes[up_ii][5] == "-":
                        """Yes, map scores to the downstream gene."""
                        up_ok = True
            if down_ii != None:
                if genes[down_ii][2] > sumsite and genes[down_ii][3] > sumsite:
                    if genes[down_ii][5] == "+":
                        """Yes, map scores to the upstream gene."""
                        down_ok = True

            #(closest_up, min_up, closest_down, min_down) = get_genes4site(con, repid, sumsite, chrid, speciesid=None)

            if up_ok and up_ii != None:
                distance = abs(sumsite - genes[up_ii][2])
                sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                sql += " VALUES(" + genes[up_ii][0].__str__() + "," 
                sql += sid.__str__() + ","
                sql += distance.__str__() + ") "         
                cur.execute(sql) 
            if down_ok and down_ii != None:
                distance = abs(genes[down_ii][2] - sumsite)
                sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                sql += " VALUES(" + genes[down_ii][0].__str__() + "," 
                sql += sid.__str__() + ","
                sql += distance.__str__() + ") "           
                cur.execute(sql) 
            
            if down_ok == None and up_ok == None:
                print "\n. I cannot find a nearby gene for summit ID", sid, " in replicate", repid
#     except:
#         print "\n. An error occurred (980). I'm rolling-back changes to the table GeneSummits."
#         print sys.exc_info()[0]
#         con.rollback()
    con.commit()
    return con
                
     