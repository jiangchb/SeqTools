############################################################
#
# A Python API / wrapper for the SQLLite ChIP-Seq database
#
############################################################

import sqlite3 as lite
import os, sys
from version import *

from chipseqdb_api import *

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
    cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key autoincrement, name TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, start INT, stop INT, chrom INT, strand TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species INT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GFFs(id INTEGER primary key autoincrement, species INT, filepath TEXT)")
    
    # This data comes from the pillars file:
    cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(realname TEXT COLLATE NOCASE, alias TEXT COLLATE NOCASE)")

    # This data from from the pillars file X GFF.
    cur.execute("CREATE TABLE IF NOT EXISTS GeneHomology(geneid INTEGER, aliasid INTEGER)")
    
    # genomic regions in which peaks should be ignored, for various reasons
    cur.execute("CREATE TABLE IF NOT EXISTS RedFlagRegions(chromid INT, start INT, stop INT)")
    
    cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, name TEXT unique COLLATE NOCASE, species INT)")        
    
    # These data come from MACS2 output files
    cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
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
    cur.execute("CREATE TABLE IF NOT EXISTS Files(fileid INTEGER primary key autoincrement, path TEXT, note TEXT)")
    
#     cur.execute("CREATE TABLE IF NOT EXISTS Compgroups(groupid INTEGER primary key autoincrement, name TEXT, note TEXT)")
#     cur.execute("CREATE TABLE IF NOT EXISTS CompgroupMembership(groupid INTEGER, repgroupid INTEGER, series INTEGER)") #series is 0 or 1. 
#     cur.execute("CREATE TABLE IF NOT EXISTS CompgroupSummits(groupid INTEGER)")
#     cur.execute("CREATE TABLE IF NOT EXISTS CompgroupFiles(groupid INTEGER, fileid INTEGER)")
    
    
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
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, note TEXT)")
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
    cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key autoincrement, name TEXT, type TEXT)") # defines a union set, type can be 'intersection' or 'union'
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
    cur.execute("CREATE TABLE IF NOT EXISTS Speciesunions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
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
    fin = open(gffpath, "r")
    curr_chromname = None # the name of the last-seen chromosome.
    curr_chromid = None # the chromosome ID (from the table Chromosomes) of the last-seen chromosome.
    try:
        for l in fin.xreadlines():
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
                
                start = int( tokens[3] )
                stop = int( tokens[4] )
                strand = tokens[6]
                if strand == "-":
                    x = start
                    y = stop
                    start = y
                    stop = x
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
                                                      
                sql = "INSERT INTO Genes (name, start, stop, chrom, strand) VALUES('" + gene + "'," + start.__str__() + "," + stop.__str__() + "," + curr_chromid.__str__() + ",'" + strand + "')"
                cur.execute(sql) 
                #print "240:", gene, start, stop, curr_chromid, strand
    except:
        con.rollback()
    fin.close()
    con.commit()
    
    sql = "insert or replace into GFFs (species, filepath) VALUES(" + speciesid.__str__() + ",'" + gffpath + "')"
    cur.execute(sql)
    con.commit() 
    
    #GFFs(id INTEGER primary key autoincrement, species INT, filepath TEXT)
    
    print ""
    cur.execute("SELECT id,name from Chromosomes where species=" + speciesid.__str__())
    for ii in cur.fetchall():
        sql = "SELECT COUNT(*) from Genes where chrom=" + ii[0].__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        print "    --> Chromosome", ii[1], "contains", count, "genes."
    
    #cur.execute("SELECT COUNT(*) FROM Genes")
    #count_genes = cur.fetchone()[0]
    #print "\n. The database now contains", count_genes, "total genes."
    
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
            if l.startswith("orf"):
                count += 1
                #if count%100 == 0:
                    #sys.stdout.write(".")
                    #sys.stdout.flush()            
                l = l.strip()
                tokens = l.split()
                #print tokens
                orf_list = []
                for t in tokens:
                    if False == t.startswith("-"):
                        orf_list.append(t)
                orf_list = orf_list
                #print orf_list
    
                realname = orf_list[0]
                for aliasname in orf_list:
                    sql = "INSERT INTO GeneAlias (realname, alias) VALUES('" + realname + "','" + aliasname + "')"
                    cur.execute(sql)
    except:
        con.rollback()
    con.commit()
    fin.close()
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
        print "\n. Error (357) Something wrong occurred while importing summits from", summitpath
        print ". I'm rolling back the changes to the tables Summits"
        con.rollback()
        exit()
    fin.close()
    con.commit()
    
    cur.execute("SELECT COUNT(*) FROM Summits where replicate=" + repid.__str__())
    count_genes = cur.fetchone()[0]
    print "\n. Found", count_genes, "summits in", summitpath
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

def import_foldenrichment(bdgpath, repid, con):
    """Imports a BDG file with enrichment scores tracked across genome sequence sites.
    The BDG must be sorted by site for each chromosome, ascending.
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
    
    """Open the BDG file"""
    fin = open(bdgpath, "r")
    curr_chromname = None
    curr_chromid = None
    last_start_site = 0
    
    """As we parse the BDG, we'll map enrichment values to genes,
        and store the following information. . . """
    geneid_sum = {} # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
    geneid_n = {} # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
    geneid_max = {}
    geneid_maxsite = {} # key = geneid, value = the distance from TSS of the maximum FE in the gene's regulatory region
    genes = None # an ordered list of gene objects, 
                    # it will be filled with data whenever we 
                    # encounter a new chromosome in the BDG file
    
#     allgenes = get_genes_for_species(con, speciesid)
#     genes = {}
#     for ii in allgenes:
#         genes[ ii[0] ] = ii
    
    chromid_genepairs = {} # the pair of genes before and after this enrichment window
    pairi = 0 # index into chromid_genepairs[curr_chromid]
    count = 0
    total_count = estimate_line_count(bdgpath)
    
    chromid_summitsites = {} # key = chromosome ID, value = list of summit sites on that chrom.
    count_found_summits = 0
    
    """For each line in the BDG file"""
    for l in fin.xreadlines():
        if l.__len__() <= 5:
            """Skip to the next line."""
            continue
        
        """Print a progress indicator."""
        count += 1
        if count%10 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) + "     (" + count_found_summits.__str__() + " summits)" )
            sys.stdout.flush()
                    
        tokens = l.split()
        
        """Is the chromosome in this line the same chromosome from the previous line?
            If not, then we need to retrieve information about this chromosome, including
            the genes on this chrom."""
        chromname = tokens[0]
        if curr_chromname == chromname:
            """If this is the same chrom. we saw on the previous line,
                then we already know the chrom. ID"""
            chromid = curr_chromid
        else:
            """We need to look up the chrom. ID"""
            chromid = get_chrom_id(con, chromname, speciesid)
            if chromid == None:
                """We don't know anything about this chromosome; skip to the next FE window."""
                continue

            if curr_chromid != chromid:
                """genes will be sorted by the start sites of the genes."""
                genes = get_genes_for_chrom(con, chromid)
                if genes.__len__() == 0:
                    print "\n. An error occurred at chipseqdb.py 322"
                    exit()

            curr_chromname = chromname
            curr_chromid = chromid
            last_start_site = 0
                                         
            if curr_chromid not in chromid_summitsites:
                chromid_summitsites[curr_chromid] = []
                sql = "select site from Summits where replicate=" + repid.__str__() + " and chrom=" + curr_chromid.__str__()
                sql += " order by site ASC"
                cur.execute(sql)
                x = cur.fetchall()
                for ii in x:
                    chromid_summitsites[curr_chromid].append( ii[0] )
                chromid_summitsites[curr_chromid].sort()
             
            if curr_chromid not in chromid_genepairs:
                chromid_genepairs[curr_chromid] = []
                for ii in xrange(0, genes.__len__()):
                    if ii == 0:
                        pair = (None,ii)
                        chromid_genepairs[curr_chromid].append( pair )
                    if ii == genes.__len__()-1:
                        pair = (ii,None)
                        chromid_genepairs[curr_chromid].append( pair )
                        pair = (ii-1,ii)
                        chromid_genepairs[curr_chromid].append( pair )
                    else:
                        pair = (ii-1,ii)
                        chromid_genepairs[curr_chromid].append( pair )                        
                pairi = 0
                                        
        """Get the fold-enrichment values for this line"""
        start = int(tokens[1]) # start of this enrichment window
        stop = int(tokens[2])
        eval = float(tokens[3]) # enrichment value across this window
        
        if last_start_site < start and last_start_site != 0:
            print ". Warning: the BDG file may skip some sites, at site:", start, "for chrom", curr_chromname, "for BDG", bdgpath
        last_start_site = stop
    
        """Can we map this enrichment site to a summit?"""        
        summit_here = False
        for ii in range(start, stop):
            """For each site in the enrichment window, is there a known summit at this site?"""
            if ii in chromid_summitsites[curr_chromid]:
                summit_here = ii
        if summit_here != False:
            """Yes, there is a summit in this enrichment window."""
            sql = "select id, score from Summits where replicate=" + repid.__str__() + " and chrom=" + curr_chromid.__str__()
            sql += " and site=" + summit_here.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            if x == None:
                print "\n. Error, I can't find the summit at site", summit_here.__str__()," for replicate", repid.__str__(), "on chrom", curr_chromid.__str__()
            
            if x != None:
                count_found_summits += 1
                sql = "insert into SummitsEnrichment (summit, max_enrichment) "
                sql += " VALUES(" + x[0].__str__() + ","
                sql += eval.__str__() + ")"
                cur.execute(sql)
                con.commit()
        """NOTE: if summit_here == False, then this enrichment window has no summit"""
              
        """Does this enrichment site map to a known gene?"""
#         (closest_up, min_up, closest_down, min_down) = get_genes4site(con, repid, start, curr_chromid, speciesid=speciesid)
#  
#         up_ok = False
#         down_ok = False
#         if closest_up != None:
#             down_ii = closest_up
#             up_ok = True
#         if closest_down != None:
#             ups_ii = closest_down
#             down_ok = True
        
        """Ensure that pairi points to correct intergenic region."""
        this_gene_pair = chromid_genepairs[ curr_chromid ][pairi]
        
        """If the current enrichment window ('start') is beyond the intergenic region defined by the
            current gene pair, then we need to advance to the next gene pair.
        """
        while this_gene_pair[1] != None and (genes[ this_gene_pair[1] ][2] < start and genes[ this_gene_pair[1] ][3] < start):
            pairi += 1
            
        """Can we map enrichment to both upstream and downstream genes?"""
        down_ok = False # is there a downstream gene?
        up_ok = False   # is there an upstream gene?
        ups_ii = chromid_genepairs[curr_chromid][pairi][0] # the ID of the downstream gene
        down_ii = chromid_genepairs[curr_chromid][pairi][1]   # the ID of the upstream gene
        
        if ups_ii != None:
            if genes[ups_ii][2] < start and genes[ups_ii][3] < start:
                if genes[ups_ii][5] == "-":
                    """Yes, map scores to the downstream gene."""
                    up_ok = True
        if down_ii != None:
            if genes[down_ii][2] > start and genes[down_ii][3] > start:
                if genes[down_ii][5] == "+":
                    """Yes, map scores to the upstream gene."""
                    down_ok = True

        if up_ok:     
            geneid = genes[ups_ii][0]
            """Initialize some data structures about this gene."""
            if geneid not in geneid_sum:
                geneid_sum[geneid] = 0
            if geneid not in geneid_n:
                geneid_n[geneid] = 0
            if geneid not in geneid_max:
                geneid_max[geneid] = 0
            if geneid not in geneid_maxsite:
                geneid_maxsite[geneid] = 0
            
            for ii in range(start, stop):
                geneid_sum[geneid] += eval
                geneid_n[geneid] += 1
                
                """Is this fe value larger than we've seen before?"""
                if eval > geneid_max[geneid]:
                    geneid_max[geneid] = eval
                    geneid_maxsite[geneid] = ii - genes[ups_ii][2]
        
        if down_ok:     
            geneid = genes[down_ii][0]
            """Initialize some data structures about this gene."""
            if geneid not in geneid_sum:
                geneid_sum[geneid] = 0
            if geneid not in geneid_n:
                geneid_n[geneid] = 0
            if geneid not in geneid_max:
                geneid_max[geneid] = 0
            if geneid not in geneid_maxsite:
                geneid_maxsite[geneid] = 0
            
            for ii in range(start, stop):
                geneid_sum[geneid] += eval
                geneid_n[geneid] += 1
                
                """Is this fe value larger than we've seen before?"""
                if eval > geneid_max[geneid]:
                    geneid_max[geneid] = eval
                    geneid_maxsite[geneid] = genes[down_ii][2] - ii
            
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
    print sql
    cur.execute(sql)
    x = cur.fetchall()
    if x.__len__() > 0:
        print "\n. Error: I found no FE values for", x.__len__().__str__(), "genes."
        print ". Replicate: ", repid
        for ii in x:
            print "Gene:", ii[0].__str__()
        exit()
    

    #     cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, start INT, stop INT, chrom INT, strand TEXT)")
    #cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species INT)
    #EnrichmentStats(repid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT, maxenrichsite INT)

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
                #sys.stdout.write(".")
                #sys.stdout.flush()
                """Note the commit happens here"""
                con.commit()
            
            this_id = g[0]
            this_name = g[1]
            
            """If there is an alias reference for this gene, then get the reference name for this gene."""
            realname = get_genename_for_aliasname( this_name, con )
            if realname == None:
                #print "No pillar entry for", this_name
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
    
    
def map_summits2genes(con, repid, speciesid=None, chroms=None):
    """This methods puts values into the table GeneSummits.
        speciesid can be the integer ID of a species whose genome will be used to map summits.
        chroms can be a list of chromosome IDs to restrict the mapping (i.e. map summits only to
            genes on these chromosomes)."""

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
    try:
        for chrid in chroms:
            genes = get_genes_for_chrom(con, chrid)
            summits = get_summits(con, repid, chrid)
                    
            for s in summits:
                sid = s[0] # summit ID
                sumsite = s[3] # summit site in the genome
                score = s[5] # summit score
                
                (closest_up, min_up, closest_down, min_down) = get_genes4site(con, repid, sumsite, chrid, speciesid=None)
                
#                 min_up = None
#                 closest_up = ""
#                 min_down = None
#                 closest_down = ""
#                 for g in genes:
#                     if count%100==0:
#                         #sys.stdout.write(".")
#                         #sys.stdout.flush()
#                         con.commit()
                                   
#                     gid = g[0]
#                     start = g[2]
#                     stop = g[3]
#                     d = start - sumsite
#                     
#                     """Sense direction, and upstream"""
#                     if start < stop and start >= sumsite:
#                         if min_up == None:
#                             min_up = d
#                             closest_up = gid
#                         if min_up > d:
#                             min_up = d
#                             closest_up = gid
#                             
#                     elif start < stop and stop < sumsite:
#                         """Sense direction and downstream"""
#                         if min_down == None:
#                             min_down = d
#                             closest_down = None
#                         if min_down < d:
#                             min_down = d
#                             closest_down = None
#                         
#                     elif start > stop and start <= sumsite:
#                         """Antisense and downstream"""
#                         if min_down == None:
#                             min_down = d
#                             closest_down = gid
#                         if min_down < d: # be careful here, we're looking for the largest NEGATIVR number.
#                             min_down = d
#                             closest_down = gid
#                             
#                     elif start > stop and stop > sumsite:
#                         """Antisense and upstream"""
#                         if min_up == None:
#                             min_up = d
#                             closest_up = None
#                         if min_up > d:
#                             min_up = d
#                             closest_up = None
    
                if closest_up != None and min_up != None:
                    sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                    sql += " VALUES(" + closest_up.__str__() + "," 
                    sql += sid.__str__() + ","
                    sql += min_up.__str__() + ") "
                    #print sql         
                    cur.execute(sql) 
                if closest_down != None and min_down != None:
                    sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                    sql += " VALUES(" + closest_down.__str__() + "," 
                    sql += sid.__str__() + ","
                    sql += min_down.__str__() + ") "  
                    #print sql           
                    cur.execute(sql) 
                
                if min_up == None and min_down == None:
                    print "\n. I cannot find a nearby gene for summit", sid, "(species", speciesid, ")"
    except:
        print "\n. An error occurred (980). I'm rolling-back changes to the table GeneSummits."
        con.rollback()
    con.commit()
    return con
                
     