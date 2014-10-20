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
    print "\n. The database contains", species.__len__(), "species:"
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
    
    # This data comes from the pillars file:
    cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(realname TEXT COLLATE NOCASE, alias TEXT COLLATE NOCASE)")
    
    # This data from from the pillars file X GFF.
    cur.execute("CREATE TABLE IF NOT EXISTS GeneHomology(geneid INTEGER, aliasid INTEGER)")
    
    # genomic regions in which peaks should be ignored, for various reasons
    cur.execute("CREATE TABLE IF NOT EXISTS RedFlagRegions(chromid INT, start INT, stop INT)")
    
    cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, name TEXT unique COLLATE NOCASE, species INT)")        
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, note TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GroupReplicate(rgroup INTEGER, replicate INTEGER)")

    # to-do:
    # the table RepgroupGenes can be removed and its uses
    # can be replaced by using the genes found in the table RepgroupSummits
    cur.execute("CREATE TABLE IF NOT EXISTS RepgroupGenes(repgroupid INTEGER, geneid INTEGER)") # genes that have summits in all replicates of a repgroup 

    # These data come from MACS2 output files
    cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INTEGER, distance INT)") # a mapping of Summits to nearby Genes
    cur.execute("CREATE TABLE IF NOT EXISTS RepgroupSummits(repgroupid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT)") # maps summits that exist in all replicates in the group.
    
    cur.execute("CREATE TABLE IF NOT EXISTS EnrichmentStats(repid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT, maxenrichsite INT)") #geneid is the canonical geneID from pillars
    cur.execute("CREATE TABLE IF NOT EXISTS GroupEnrichmentStats(rgroupid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
    
    cur.execute("CREATE TABLE IF NOT EXISTS Files(fileid INTEGER primary key autoincrement, path TEXT, note TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicategroupFiles(repgroupid INTEGER, fileid INTEGER)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionFiles(unionid INTEGER, fileid INTEGER)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionFiles(spunionid INTEGER, fileid INTEGER)")
        
    build_unions(con)
    build_idr_tables(con)
    con.commit()
    
    if dbpath != None and dbpath != False:
        print_db_stats(con)
    
    return con

def build_unions(con):
    cur = con.cursor()
    # These tables describe which replicates are to be unioned.
    cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS UnionRepgroups(unionid INTEGER, repgroupid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS UnionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    cur.execute("CREATE TABLE IF NOT EXISTS UnionSummits(unionid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT)")
    
    # UnionSummitStats is depricated:
    cur.execute("CREATE TABLE IF NOT EXISTS UnionSummitStats(unionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT, ranksummit FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT, rankmeanenrich FLOAT)")
    con.commit()


def build_speciesunions(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Speciesunions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionUnions(spunionid INTEGER, memunionid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionSummits(spunionid INTEGER, geneid INTEGER, maxsummitid INT, nearestsummitid INT)")
    
    # geneid in the Speciesunion tables point to translated gene IDs from the pillars.
    # This means that to find this gene ID in a particular Union, you need to use
    # the table GeneHomology to find alias gene IDs for another species.

    # SpeciesunionSummitStats is depricated
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionSummitStats(spunionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
    con.commit() 

def build_idr_tables(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS GeneRepgroupEnrichIdr(geneid INTEGER, repid1 INTEGER, repid2 INTEGER, lidr FLOAT, idr FLOAT)") # lidr is local IDR, idr is IDR for observations that are as irreproducible, or more, than this observation
    cur.execute("CREATE TABLE IF NOT EXISTS GeneUnionEnrichIdr(geneid INTEGER, repgroupid1 INTEGER, repgroupid2 INTEGER, lidr FLOAT, idr FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS GeneSpeciesunionEnrichIdr(geneid INTEGER, unionid1 INTEGER, unionid2 INTEGER, lidr FLOAT, idr FLOAT)")
    con.commit()

def import_gff(gffpath, speciesid, con, restrict_to_feature = "gene"):
    cur = con.cursor()
    
    print "\n. Importing genome features from", gffpath
    
    chromids = get_chrom_ids(con, speciesid)
    
    """Remove stale data."""
    sql = "DELETE FROM Chromosomes where species=" + speciesid.__str__()
    cur.execute(sql)
    for chromid in chromids:
        sql = "DELETE FROM Genes where chrom=" + chromid.__str__()
        cur.execute(sql)
    con.commit()
    
    chromname_id = {}
    
    count = 0
    total_count = estimate_line_count(gffpath)
    fin = open(gffpath, "r")
    curr_chromname = None # the name of the last-seen chromosome.
    curr_chromid = None # the chromosome ID (from the table Chromosomes) of the last-seen chromosome.
    for l in fin.xreadlines():
        count += 1
        sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
        sys.stdout.flush()
        
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            if tokens[2] != restrict_to_feature: # e.g., if restrict_to_feature == "gene", then we'll only import genes.
                continue
            chr = tokens[0]
            
            if chr != curr_chromname:
                """Add (or potentially overwrite) this chromosome into the table Chromosomes."""
                curr_chromname = chr
                #sql = "SELECT COUNT(*) from Chromosomes where name='" + curr_chromname + "' and species=" + speciesid.__str__()
                #cur.execute(sql)
                #if cur.fetchone()[0] == 0:
                if curr_chromname not in chromname_id:
                    sql = "SELECT COUNT(*) FROM Chromosomes where name='" + curr_chromname + "'"
                    cur.execute(sql)
                    if cur.fetchone()[0] == 0:
                        """We've not seen this chrom yet... insert it!"""
                        sql = "INSERT INTO Chromosomes (name,species) VALUES('" + curr_chromname + "', " + speciesid.__str__() + ")"
                        cur.execute( sql )  
                        con.commit()
                    """Get the row ID of the newly-added chromosome."""
                    cur.execute("SELECT id FROM Chromosomes WHERE name='" + curr_chromname + "'")
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
            
            if tokens[8].split(";").__len__() > 2:
                note = tokens[8].split(";")[2].split("=")[1]
                if note.__contains__("orf"):
                    for t in note.split():
                        if t.startswith("orf"):
                            gene = t # use this name instead of the orfName
                                                  
            sql = "INSERT INTO Genes (name, start, stop, chrom, strand) VALUES('" + gene + "'," + start.__str__() + "," + stop.__str__() + "," + curr_chromid.__str__() + ",'" + strand + "')"
            cur.execute(sql) 
    fin.close()
    con.commit()
    
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
    con.commit()
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
    
    #
    # Build a library of summits
    #
    count = 0
    total_count = estimate_line_count(summitpath)
    chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
    fin = open(summitpath, "r")
    for l in fin.xreadlines():
        count += 1
        sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
        sys.stdout.flush()
        
        if count%100==0:
            con.commit()
        
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            chr = tokens[0]
            site = int( tokens[1] )
            name = tokens[3]
            score = float( tokens[4] )
            
            cur = con.cursor()
            sql = "SELECT * FROM Chromosomes WHERE name='" + chr + "'"                
            cur.execute(sql)
            x = cur.fetchone()
            if x == None:
                print "\n. Warning: Your summit file includes summits on chromosome", chr, "but your GFF file lacks this chromosome."
                print ". This is most likely because there are no annoted genes located on", chr
                print ". You should verify if this is true."
                print ". In the meantime, I am not importing the summits on", chr
            else:
                chrid = x[0]
                sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chrid.__str__() + "," + score.__str__() + ")"
                cur.execute(sql)
    fin.close()
    con.commit()
    
    cur.execute("SELECT COUNT(*) FROM Summits")
    count_genes = cur.fetchone()[0]
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

def import_bdg(bdgpath, repid, con):
    """Imports a BDG file with enrichment scores tracked across genome sequence sites.
    The BDG must be sorted by site for each chromosome, ascending.""" 
    print "\n. Importing fold-enrichment values from", bdgpath,"for replicate", repid

    cur = con.cursor()
    
    if False == os.path.exists(bdgpath):
        print "\n. Error: I cannot find the BedGraph file", bdgpath
        exit()
    
    """Remove any stale entries."""
    sql = "DELETE FROM EnrichmentStats where repid=" + repid.__str__()
    cur.execute(sql)
    con.commit()
    
    fin = open(bdgpath, "r")
    curr_chromname = None
    curr_chromid = None
    last_start_site = -1
    geneid_sum = {} # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
    geneid_n = {} # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
    geneid_max = {}
    geneid_maxsite = {} # key = geneid, value = the site number of the maximum FE in the gene's regulatory region
    genes = None # an ordered list of gene objects, 
                    # it will be filled with data whenever we 
                    # encounter a new chromosome in the BDG file
    gene_pairs = []
    pairi = 0 # index into gene_pairs
    count = 0
    total_count = estimate_line_count(bdgpath)
    for l in fin.xreadlines():
        if l.__len__() <= 5:
            """Skip to the next line."""
            continue
        
        """Progress indicator."""
        count += 1
        if count%10 == 0:
            sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
            sys.stdout.flush()
                    
        tokens = l.split()
        
        """Get a list of genes for this chromosome."""
        chromname = tokens[0]
        if curr_chromname == chromname:
            chromid = curr_chromid
        else:
            """Found a new chromosome."""
            chromid = get_chrom_id( con, chromname )
            if chromid == None:
                """Skip to the next FE window."""
                continue
            
            curr_chromname = chromname
            curr_chromid = chromid
            last_start_site = -1
            
            """genes will be sorted by the start sites of the genes."""
            genes = get_genes_for_chrom(con, chromid)
            if genes.__len__() == 0:
                print "\n. An error occurred at chipseqdb.py 322"
                exit()
            gene_pairs = []
            for ii in xrange(0, genes.__len__()):
                if ii == 0:
                    pair = (None,ii)
                if ii == genes.__len__()-1:
                    pair = (ii,None)
                else:
                    pair = (ii-1,ii)
                gene_pairs.append( pair )
            pairi = 0
                
        """Get the fold-enrichment values in this window."""
        start = int(tokens[1]) # start of this enrichment window
        stop = int(tokens[2])
        eval = float(tokens[3]) # enrichment value across this window
    
        """Ensure that pairi points to correct intergenic region."""
        #print "404:", gene_pairs[pairi]
        while gene_pairs[pairi][1] != None and (genes[ gene_pairs[pairi][1] ][2] < start and genes[ gene_pairs[pairi][1] ][3] < start):
            #print "405:", gene_pairs[pairi]
            pairi += 1
            
        """Can we map enrichment to both upstream and downstream genes?"""
        down_ok = False
        up_ok = False
        down_ii = gene_pairs[pairi][0]
        up_ii = gene_pairs[pairi][1]
        if down_ii != None:
            if genes[down_ii][2] < start and genes[down_ii][3] < start:
                if genes[down_ii][5] == "-":
                    """Yes, map scores to the downstream gene."""
                    down_ok = True
        if up_ii != None:
            if genes[up_ii][2] > start and genes[up_ii][3] > start:
                if genes[up_ii][5] == "+":
                    """Yes, map scores to the upstream gene."""
                    up_ok = True

        if down_ok:     
            geneid = genes[down_ii][0]
            #print "426:", start, pairi, geneid
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
                    geneid_maxsite[geneid] = ii - genes[down_ii][2]
        
        if up_ok:     
            geneid = genes[up_ii][0]
            #print "448:", start, pairi, geneid
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
                    geneid_maxsite[geneid] = genes[up_ii][2] - ii
    fin.close()
    
    
    """Finally, write all our findings into the table EnrichmentStats."""
    count = 0
    total_count = geneid_sum.__len__()
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
    con.commit()        

    return con

# def import_bdg_original(bdgpath, repid, con):
#     """Imports a BDG file with enrichment scores tracked across genome sequence sites.
#     The BDG must be sorted by site for each chromosome, ascending.""" 
#     print "\n. Importing fold-enrichment values from", bdgpath,"for replicate", repid
# 
#     cur = con.cursor()
#     
#     if False == os.path.exists(bdgpath):
#         print "\n. Error: I cannot find the BedGraph file", bdgpath
#         exit()
#     
#     """Remove any stale entries."""
#     sql = "DELETE FROM EnrichmentStats where repid=" + repid.__str__()
#     cur.execute(sql)
#     con.commit()
#     
#     fin = open(bdgpath, "r")
#     curr_chromname = None
#     curr_chromid = None
#     last_start_site = -1
#     geneid_sum = {} # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
#     geneid_n = {} # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
#     geneid_max = {}
#     geneid_maxsite = {} # key = geneid, value = the site number of the maximum FE in the gene's regulatory region
#     genes = None # an ordered list of gene objects, 
#                     # it will be filled with data whenever we 
#                     # encounter a new chromosome in the BDG file
#     nearest_up_gene_i = 0; # an index into genes, the nearest gene going up in order
#     nearest_down_gene_i = 0 # an index into genes, the nearest gene going down in order
#     count = 0
#     total_count = estimate_line_count(bdgpath)
#     for l in fin.xreadlines():
#         if l.__len__() <= 5:
#             """Skip to the next line."""
#             continue
#         
#         """Progress indicator."""
#         count += 1
#         if count%10 == 0:
#             sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
#             sys.stdout.flush()
#                     
#         tokens = l.split()
#         
#         """Get a list of genes for this chromosome."""
#         chromname = tokens[0]
#         if curr_chromname == chromname:
#             chromid = curr_chromid
#         else:
#             chromid = get_chrom_id( con, chromname )
#             if chromid == None:
#                 """Skip to the next FE window."""
#                 continue
#             
#             curr_chromname = chromname
#             curr_chromid = chromid
#             last_start_site = -1
#             
#             """genes will be sorted by the start sites of the genes."""
#             genes = get_genes_for_chrom(con, chromid)
#             if genes.__len__() == 0:
#                 print "\n. An error occurred at chipseqdb.py 322"
#                 exit()
#                 
#                 
#             """Reset our gene counters"""
#             nearest_up_gene_i = 0
#             nearest_down_gene_i = 0
#         
#         """Get the fold-enrichment values in this window."""
#         start = int(tokens[1]) # start of this enrichment window
#         stop = int(tokens[2])
#         eval = float(tokens[3]) # enrichment value across this window
#     
#         """Sanity check:"""
#         if last_start_site >= start:
#             """Then the BDG probably isn't sorted, because we're seeing sites out of order."""
#             print "\n. Error importing BDG. I don't think your BDG file is sorted:", bdgpath
#             exit()
#         last_start_site = start
#     
#         """Is nearest_up_gene_i actually the nearest upstream gene?"""
#         while nearest_up_gene_i < genes.__len__() and (start > genes[nearest_up_gene_i][2] and start > genes[nearest_up_gene_i][3] ):
#             """ we've moved beyond this genes."""
#             nearest_up_gene_i += 1
#     
#         skip_up = False
#         """If this FE window is inside a coding region, then skip it."""
#         if (start >= genes[nearest_up_gene_i][2] and start <= genes[nearest_up_gene_i][3] ):
#             """ we're in the middle of a + gene"""
#             skip_up = True
#         elif (start <= genes[nearest_up_gene_i][2] and start >= genes[nearest_up_gene_i][3] ):
#             """we're in the middle of a - gene"""
#             skip_up = True
#         
#         """A normal regulatory region for + genes"""
#         if skip_up == False and (start < genes[nearest_up_gene_i][2] and start < genes[nearest_up_gene_i][3] ):
#             """Get a tuple representing the gene's meta information"""
#             this_gene = genes[nearest_up_gene_i]
#             
#             if this_gene[5] != "+":
#                 """Then this gene is negative strand. Continue with the next loop iteration."""
#                 continue
#             
#             geneid = this_gene[0]
#             """Initialize some data structures about this gene."""
#             if geneid not in geneid_sum:
#                 geneid_sum[geneid] = 0
#             if geneid not in geneid_n:
#                 geneid_n[geneid] = 0
#             if geneid not in geneid_max:
#                 geneid_max[geneid] = 0
#             if geneid not in geneid_maxsite:
#                 geneid_maxsite[geneid] = 0
#             
#             for ii in range(start, stop):
#                 geneid_sum[geneid] += eval
#                 geneid_n[geneid] += 1
#                 
#                 """Is this fe value larger than we've seen before?"""
#                 if eval > geneid_max[geneid]:
#                     geneid_max[geneid] = eval
#                     geneid_maxsite[geneid] = genes[nearest_up_gene_i][2] - ii
# 
#         """Repeat for downstream_gene..."""
#         while nearest_down_gene_i != 0 and nearest_down_gene_i < genes.__len__() and (start < genes[nearest_down_gene_i][2] and start < genes[nearest_down_gene_i][3] ):
#             """ we've moved beyond this gene."""
#             nearest_down_gene_i += 1
#     
#         """If this FE window is inside a coding region, then skip it."""
#         if (start >= genes[nearest_down_gene_i][2] and start <= genes[nearest_down_gene_i][3] ):
#             """ we're in the middle of a + gene"""
#             continue
#         elif (start <= genes[nearest_down_gene_i][2] and start >= genes[nearest_down_gene_i][3] ):
#             """we're in the middle of a - gene"""
#             continue
#         
#         """A normal regulatory region for + genes"""
#         if (start < genes[nearest_up_gene_i][2] and start < genes[nearest_up_gene_i][3] ):
#             """Get a tuple representing the gene's meta information"""
#             this_gene = genes[nearest_up_gene_i]
#             
#             if this_gene[5] != "+":
#                 """Then this gene is negative strand. Continue with the next loop iteration."""
#                 continue
#             
#             geneid = this_gene[0]
#             """Initialize some data structures about this gene."""
#             if geneid not in geneid_sum:
#                 geneid_sum[geneid] = 0
#             if geneid not in geneid_n:
#                 geneid_n[geneid] = 0
#             if geneid not in geneid_max:
#                 geneid_max[geneid] = 0
#             if geneid not in geneid_maxsite:
#                 geneid_maxsite[geneid] = 0
#             
#             for ii in range(start, stop):
#                 geneid_sum[geneid] += eval
#                 geneid_n[geneid] += 1
#                 
#                 """Is this fe value larger than we've seen before?"""
#                 if eval > geneid_max[geneid]:
#                     geneid_max[geneid] = eval
#                     geneid_maxsite[geneid] = genes[nearest_up_gene_i][2] - ii    
# 
#                                 
#             
#     count = 0
#     total_count = geneid_sum.__len__()
#     for geneid in geneid_sum:
#         if geneid_n[geneid] > 0:
#             count += 1
#             if count%5000==0:
#                 #sys.stdout.write(".")
#                 #sys.stdout.flush()
#                 con.commit()
#             sql = "INSERT INTO EnrichmentStats (repid, geneid, maxenrich, meanenrich, sumenrich, maxenrichsite)  "
#             sql += "VALUES(" + repid.__str__() + "," + geneid.__str__()
#             sql += "," + geneid_max[geneid].__str__()
#             sql += "," + (geneid_sum[geneid]/float(geneid_n[geneid])).__str__()
#             sql += "," + geneid_sum[geneid].__str__()
#             sql += "," + geneid_maxsite[geneid].__str__()
#             sql += ")"
#             cur.execute(sql)
#     con.commit()        
#     fin.close()
#     
#     #print "384:", geneid_sum # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
#     #print "385:", geneid_n # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
#     #print "386:", geneid_max     
# 
#     #con.commit()
#     return con


def resolve_aliasids(con):
    """This method inserts data into the table GeneHomology."""
    print "\n. Resolving homologous gene IDs between species."
    
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
    
    """Save the commit for the end."""
    con.commit()

def map_summits2genes(con, repid, speciesid=None, chromid=None):
    """This methods puts values into the table GeneSummits."""

    cur = con.cursor()
    speciesids = []
    if speciesid:
        spids = [speciesid]
    else:
        spids = get_species_ids(con)
    
    count = 0
    for spid in spids:
        chroms = get_chrom_ids(con, spid)
        for chrid in chroms:
            genes = get_genes_for_chrom(con, chrid)
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
                    if count%100==0:
                        #sys.stdout.write(".")
                        #sys.stdout.flush()
                        con.commit()
                    
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
    con.commit()
    return con
                
     