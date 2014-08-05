############################################################
#
# A Python API / wrapper for the SQLLite ChIP-Seq database
#
############################################################

import sqlite3 as lite
import os, sys

from chipseqdb_api import *

def build_db(dbpath = None):
    """Initializes all the tables. Returns the DB connection object.
    If tables already exist, they will NOT be overwritten."""
    
    print "\n. Building the database. . ."
    
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
    con = lite.connect(dbpath)
    
    with con:
        cur = con.cursor()
        # These data come from the GFF:
        cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key autoincrement, name TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, start INT, stop INT, chrom INT, strand TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species INT)")
        
        # This data comes from the pillars file:
        cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(realname TEXT COLLATE NOCASE, alias TEXT COLLATE NOCASE)")
        
        # This data from from the pillars file X GFF.
        cur.execute("CREATE TABLE IF NOT EXISTS GeneHomology(geneid INTEGER, aliasid INTEGER)")
        
        cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, name TEXT unique COLLATE NOCASE, species INT)")        
        cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key autoincrement, name TEXT COLLATE NOCASE, note TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GroupReplicate(rgroup INTEGER, replicate INTEGER)")

        cur.execute("CREATE TABLE IF NOT EXISTS RepgroupGenes(repgroupid INTEGER, geneid INTEGER)") # genes that have summits in all replicates of a repgroup 

        # These data come from MACS2 output files
        cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INT, distance INT)") # a mapping of Summits to nearby Genes
        cur.execute("CREATE TABLE IF NOT EXISTS RepgroupSummitStats(repgroupid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
        
        cur.execute("CREATE TABLE IF NOT EXISTS EnrichmentStats(repid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)") #geneid is the canonical geneID from pillars
        cur.execute("CREATE TABLE IF NOT EXISTS GroupEnrichmentStats(rgroupid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
        
        build_unions(con)
    return con

def import_gff(gffpath, speciesid, con, restrict_to_feature = "gene"):
    print "\n. Importing genome features from", gffpath
    
    """Returns chr_gene_sites"""
    cur = con.cursor()
    count = 0
    chr_gene_sites = {} # key = chromosome name, value = hash; key = gene name, value = (Start, stop)
    fin = open(gffpath, "r")
    for l in fin.xreadlines():
        if l.__len__() > 0 and False == l.startswith("#"):
            tokens = l.split()
            if tokens[2] != restrict_to_feature: # e.g., if restrict_to_feature == "gene", then we'll only import genes.
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
            
            if tokens[8].split(";").__len__() > 2:
                note = tokens[8].split(";")[2].split("=")[1]
                if note.__contains__("orf"):
                    for t in note.split():
                        if t.startswith("orf"):
                            print "342:", gene, t
                            gene = t # use this name instead of the orfName
            
            count += 1
            if count%100 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                          
            with con:
                cur.execute("SELECT COUNT(*) FROM Chromosomes WHERE name='" + chr + "'")
                data = cur.fetchone()
                if data[0] == 0: # that chromosome doesn't exist yet.
                    sql = "INSERT INTO Chromosomes (name,species) VALUES('" + chr + "'," + speciesid.__str__() + ")"
                    #print sql
                    cur.execute( sql )
                    con.commit()
            
                cur = con.cursor()                
                cur.execute("SELECT * FROM Chromosomes WHERE name='" + chr + "'")
                data = cur.fetchone()
                chrid = data[0]
                cur.execute("SELECT COUNT(*) from Genes where name='" + gene + "' and chrom=" + chrid.__str__())
                if cur.fetchone()[0] == 0:                
                    sql = "INSERT INTO Genes (name, start, stop, chrom, strand) VALUES('" + gene + "'," + start.__str__() + "," + stop.__str__() + "," + chrid.__str__() + ",'" + strand + "')"
                    cur.execute(sql)    
    fin.close() 
    
    cur.execute("SELECT COUNT(*) FROM Genes")
    count_genes = cur.fetchone()[0]
    #print "\n. The database now contains", count_genes, "total genes."
    
    return con



def import_pillars(pillarspath, con):
    print "\n. Reading Pillars from", pillarspath
    cur = con.cursor()
    fin = open(pillarspath, "r")
    count = 0
    for l in fin.xreadlines():
        if l.startswith("orf"):
            count += 1
            if count%100 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()            
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
    cur.execute("SELECT count(DISTINCT realname) from GeneAlias")
    x = cur.fetchone()
    print "\n. The pillars file contains", x[0], "orthologs."
    return con

def import_summits(summitpath, repid, con):
    """Reads a summit file (from MACS2, for example), and puts values into the Summits table."""
    
    print "\n. Importing summits from", summitpath,"for replicate", repid
    cur = con.cursor()
    
    if False == os.path.exists(summitpath):
        print "\n. Error: I cannot find the summit file", summitpath
        exit()
    
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
                sql = "SELECT * FROM Chromosomes WHERE name='" + chr + "'"                
                cur.execute(sql)
                x = cur.fetchone()
                if x == None:
                    print "Warning: Your summit file includes the chromosome", chr, "but your GFF file lacks this chromosome. I'm skipping it."
                else:
                    chrid = x[0]
                    sql = "INSERT INTO Summits (replicate,name,site,chrom,score) VALUES(" + repid.__str__() + ",'" + name + "'," + site.__str__() + "," + chrid.__str__() + "," + score.__str__() + ")"
                    cur.execute(sql)
    fin.close()
    
    cur.execute("SELECT COUNT(*) FROM Summits")
    count_genes = cur.fetchone()[0]
    return con

def import_bdg(bdgpath, repid, con):
    """Imports a BedGraph file with enrichment scores tracked across genome sequence sites."""
    print "\n. Importing enrichment values from", bdgpath,"for replicate", repid

    cur = con.cursor()
    
    if False == os.path.exists(bdgpath):
        print "\n. Error: I cannot find the BedGraph file", bdgpath
        exit()
    
    fin = open(bdgpath, "r")
    curr_chromname = None
    curr_chromid = None
    geneid_sum = {} # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
    geneid_n = {} # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
    geneid_max = {}
    genes = None # an ordered list of gene objects, 
                    # it will be filled with data whenever we 
                    # encounter a new chromosome in the BDG file
    nearest_up_gene_i = 0; # an index into genes, the nearest gene going up in order
    nearest_down_gene_i = 0 # an index into genes, the nearest gene going down in order
    nearest_up_gene = None
    nearest_down_gene = None
    for l in fin.xreadlines():
        if l.__len__() > 5:
            tokens = l.split()
            chromname = tokens[0]
            if curr_chromname == chromname:
                chromid = curr_chromid
            else:
                chromid = get_chrom_id( con, chromname )
                if chromid == None:
                    #print "\n. I can't find chromosome", chromname, "in the descriptions provided by your GFF."
                    #print "\n. Your bedgraph file", bdgpath,"references chromosome", chromname, ". I'm skipping this BedGraph entry."
                    continue
                curr_chromname = chromname
                curr_chromid = chromid
                
                genes = get_genes_for_chrom(con, chromid) #an ordered list of gene objects
                if genes.__len__() == 0:
                    print "\n. An error occurred at chipseqdb.py 322"
                    exit()
                nearest_up_gene_i = 0
                nearest_down_gene_i = 0
                nearest_up_gene = genes[ nearest_up_gene_i ][0]
                nearest_down_gene = genes[ nearest_down_gene_i ][0]
                geneid_sum[nearest_up_gene] = 0.0
                geneid_sum[nearest_down_gene] = 0.0
                geneid_n[nearest_up_gene] = 0
                geneid_n[nearest_down_gene] = 0
                geneid_max[nearest_up_gene] = 0
                geneid_max[nearest_down_gene] = 0
                
            start = int(tokens[1]) # start of this enrichment window
            stop = int(tokens[2])
            eval = float(tokens[3]) # enrichment value across this window
            
            """Print a period every 30K sites to indicate to the user that this program is still alive."""
            for ii in range(start, stop):
                if ii%30000 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
            
            """Add this score to the tally for the upstream gene."""
            if nearest_up_gene_i != None:
                if start < genes[nearest_up_gene_i][2] and start < genes[nearest_up_gene_i][3]:
                    if genes[nearest_up_gene_i][5] == "+":
                        # everything is good
                        for ii in range(start, stop):
                            geneid_sum[nearest_up_gene] += eval
                            geneid_n[nearest_up_gene] += 1
                            if eval > geneid_max[nearest_up_gene]:
                                geneid_max[nearest_up_gene] = eval
                else:
                    """The nearest_up_gene_i counter is stale. Update it."""
                    while nearest_up_gene_i != None and nearest_down_gene_i != None and False == (start < genes[nearest_up_gene_i][2] and start < genes[nearest_up_gene_i][3]):
                        nearest_up_gene_i += 1
                        if nearest_up_gene_i - nearest_down_gene_i > 1:
                            nearest_down_gene_i += 1
                            nearest_down_gene = genes[nearest_down_gene_i][0]
                        
                        if nearest_up_gene_i >= genes.__len__():
                            nearest_up_gene_i = None
                        else:
                            nearest_up_gene = genes[nearest_up_gene_i][0]
                            geneid_sum[nearest_up_gene] = 0.0
                            geneid_n[nearest_up_gene] = 0
                            geneid_max[nearest_up_gene] = 0
                        
                        if nearest_down_gene_i >= genes.__len__():
                            nearest_down_gene_i = None
                        else:
                            geneid_sum[nearest_down_gene] = 0.0
                            geneid_n[nearest_down_gene] = 0
                            geneid_max[nearest_down_gene] = 0
                        #print nearest_up_gene_i, nearest_down_gene_i
                        
            """Add to the tally for the downstream gene."""
            #print "Gene:", nearest_down_gene_i
            #print genes[nearest_down_gene_i]
            
            if nearest_down_gene_i != None:
                if start > genes[nearest_down_gene_i][2] and start > genes[nearest_down_gene_i][3]:
                    if genes[nearest_down_gene_i][5] == "-":
                        # everything is good
                        for ii in range(start, stop):
                            geneid_sum[nearest_down_gene] += eval
                            geneid_n[nearest_down_gene] += 1
                            if eval > geneid_max[nearest_down_gene]:
                                geneid_max[nearest_down_gene] = eval
            
    for geneid in geneid_sum:
        if geneid_n[geneid] > 0:
            sql = "INSERT INTO EnrichmentStats (repid, geneid, maxenrich, meanenrich, sumenrich)  "
            sql += "VALUES(" + repid.__str__() + "," + geneid.__str__()
            sql += "," + geneid_max[geneid].__str__()
            sql += "," + (geneid_sum[geneid]/float(geneid_n[geneid])).__str__()
            sql += "," + geneid_sum[geneid].__str__()
            sql += ")"
            cur.execute(sql)
            con.commit()
            
    fin.close()
    
    #print "384:", geneid_sum # key = geneid, value = sum of enrichment scores in its nearby regulatory regions
    #print "385:", geneid_n # key = geneid, value = number of sites with enrichment scores in its regulatory regions 
    #print "386:", geneid_max     

    #con.commit()
    return con

def resolve_aliasids(con):
    """This method inserts data into the table GeneHomology."""
    print "\n. Resolving homologous gene IDs."
    
    count = 0
    cur = con.cursor()
    sql = "SELECT * from Genes"
    cur.execute(sql)
    for g in cur.fetchall():
        count += 1
        if count%50 == 0:
            sys.stdout.write(".")
            sys.stdout.flush()
        
        this_id = g[0]
        this_name = g[1]
        realname = get_genename_for_aliasname( this_name, con )
        if realname == None:
            #print "No pillar entry for", this_name
            realname = this_name
        
        sql = "SELECT id from Genes where name='" + realname + "'"
        cur.execute(sql)
        x = cur.fetchone()
        if x != None:
            realid = x[0]
            sql = "INSERT INTO GeneHomology (geneid, aliasid) VALUES(" + realid.__str__() + "," + this_id.__str__() + ")"
            cur.execute(sql)
            con.commit()

def map_summits2genes(con, repid, speciesid=None, chromid=None):
    """This methods puts values into the table GeneSummits."""
    
    with con:
        cur = con.cursor()
        speciesids = []
        if speciesid:
            spids = [speciesid]
        else:
            spids = get_species_ids(con)
            
        for spid in spids:
            chroms = get_chrom_ids(con, spid)
            for chr in chroms:
                chrid = chr[0]
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
                            #print sql         
                            cur.execute(sql) 
                    if closest_down != None and min_down != None:
                        with con:
                            sql = "INSERT INTO GeneSummits (gene,summit,distance)" 
                            sql += " VALUES(" + closest_down.__str__() + "," 
                            sql += sid.__str__() + ","
                            sql += min_down.__str__() + ") "  
                            #print sql           
                            cur.execute(sql) 
    return con
                
     