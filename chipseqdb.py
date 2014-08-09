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
    
    cur.execute("CREATE TABLE IF NOT EXISTS Files(fileid INTEGER primary key autoincrement, path TEXT, note TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS ReplicategroupFiles(repgroupid INTEGER, fileid INTEGER)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionFiles(unionid INTEGER, fileid INTEGER)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionFiles(spunionid INTEGER, fileid INTEGER)")
        
    build_unions(con)
    con.commit()
    return con

def build_unions(con):
    cur = con.cursor()
    # These tables describe which replicates are to be unioned.
    cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS UnionRepgroups(unionid INTEGER, repgroupid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS UnionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    cur.execute("CREATE TABLE IF NOT EXISTS UnionSummitStats(unionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS UnionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
    con.commit()


def build_speciesunions(con):
    cur = con.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS Speciesunions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionUnions(spunionid INTEGER, memunionid INTEGER)") # puts repgroups into union sets
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
    
    # geneid in the Speciesunion tables point to translated gene IDs from the pillars.
    # This means that to find this gene ID in a particular Union, you may need to use
    # the table GeneHomology to find alias gene IDs for another species.
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionSummitStats(spunionid INTEGER, geneid INTEGER, maxsummit FLOAT, nsummits FLOAT)")
    cur.execute("CREATE TABLE IF NOT EXISTS SpeciesunionEnrichmentStats(unionid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")
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
    fin = open(gffpath, "r")
    curr_chromname = None # the name of the last-seen chromosome.
    curr_chromid = None # the chromosome ID (from the table Chromosomes) of the last-seen chromosome.
    for l in fin.xreadlines():
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
            
            count += 1
            if count%50 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                                       
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
        print ". Chromosome", ii[1], "contains", count, "genes."
    
    #cur.execute("SELECT COUNT(*) FROM Genes")
    #count_genes = cur.fetchone()[0]
    #print "\n. The database now contains", count_genes, "total genes."
    
    return con



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
    chr_site_score = {} # key = chromosome name, value = hash; key = site of summit, value = score for summit
    fin = open(summitpath, "r")
    for l in fin.xreadlines():
        count += 1
        if count%100==0:
            sys.stdout.write(".")
            sys.stdout.flush()
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

def import_bdg(bdgpath, repid, con):
    """Imports a BedGraph file with enrichment scores tracked across genome sequence sites."""
    print "\n. Importing enrichment values from", bdgpath,"for replicate", repid

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
    count = 0
    for l in fin.xreadlines():
        if l.__len__() > 5:
            count += 1
            if count%20000==0:
                sys.stdout.write(".")
                sys.stdout.flush()
            
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
            
#            """Print a period every 30K sites to indicate to the user that this program is still alive."""
#             for ii in range(start, stop):
#                 if ii%10000 == 0:
#                     sys.stdout.write(".")
#                     sys.stdout.flush()
            
            """Add this score to the tally for the upstream gene."""
            if nearest_up_gene_i != None:
                
                """If this enrichment site is lower than the nearest upstream gene:"""
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
                """If the enrichment site is above the nearest downstream gene."""
                if start > genes[nearest_down_gene_i][2] and start > genes[nearest_down_gene_i][3]:
                    if genes[nearest_down_gene_i][5] == "-":
                        # everything is good
                        for ii in range(start, stop):
                            geneid_sum[nearest_down_gene] += eval
                            geneid_n[nearest_down_gene] += 1
                            if eval > geneid_max[nearest_down_gene]:
                                geneid_max[nearest_down_gene] = eval
            
    count = 0
    for geneid in geneid_sum:
        if geneid_n[geneid] > 0:
            count += 1
            if count%1000==0:
                sys.stdout.write(".")
                sys.stdout.flush()
                con.commit()
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
    
    count = 0
    sql = "SELECT * from Genes"
    cur.execute(sql)
    for g in cur.fetchall():
        
        count += 1
        if count%50 == 0:
            sys.stdout.write(".")
            sys.stdout.flush()
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
                        sys.stdout.write(".")
                        sys.stdout.flush()
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
                
     