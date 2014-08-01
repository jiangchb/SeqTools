############################################################
#
# A Python API / wrapper for the SQLLite ChIP-Seq database
#
############################################################

import sqlite3 as lite
import os, sys

def build_db(dbpath = None):
    """Initializes all the tables.
    Returns the DB connection object."""
    
    print "\n. Building the database. . ."
    
    if dbpath == None or dbpath == False:
        dbpath = "test.db"
    con = lite.connect(dbpath)
    
    with con:
        cur = con.cursor()
        # These data come from the GFF:
        cur.execute("CREATE TABLE IF NOT EXISTS Species(id INTEGER primary key autoincrement, name TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Genes(id INTEGER primary key autoincrement, name TEXT, start INT, stop INT, chrom INT, strand TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS Chromosomes(id INTEGER primary key autoincrement, name TEXT, species INT)")
        
        # This data comes from the pillars file:
        cur.execute("CREATE TABLE IF NOT EXISTS GeneAlias(realname TEXT, alias TEXT)")
        
        cur.execute("CREATE TABLE IF NOT EXISTS Replicates(id INTEGER primary key autoincrement, name TEXT unique, species INT)")        
        cur.execute("CREATE TABLE IF NOT EXISTS ReplicateGroups(id INTEGER primary key autoincrement, name TEXT, note TEXT)")
        cur.execute("CREATE TABLE IF NOT EXISTS GroupReplicate(rgroup INTEGER, replicate INTEGER)")

        # These data come from MACS2 output files
        cur.execute("CREATE TABLE IF NOT EXISTS Summits(id INTEGER primary key autoincrement, replicate INT, name TEXT, site INT, chrom INT, score FLOAT, pvalue FLOAT, qvalue FLOAT)")
        cur.execute("CREATE TABLE IF NOT EXISTS EnrichmentStats(repid INTEGER, geneid INTEGER, maxenrich FLOAT, meanenrich FLOAT, sumenrich FLOAT)")

        # These tables describe which replicates are to be unioned.
        cur.execute("CREATE TABLE IF NOT EXISTS Unions(unionid INTEGER primary key autoincrement, name TEXT)") # defines a union set
        cur.execute("CREATE TABLE IF NOT EXISTS UnionRepgroups(unionid INTEGER, repgroupid INTEGER)") # puts repgroups into union sets
        cur.execute("CREATE TABLE IF NOT EXISTS UnionGenes(unionid INTEGER, geneid INTEGER)") # genes that have summits in all the repgroups in this union
        
        # These data come from mapping summits to genes
        cur.execute("CREATE TABLE IF NOT EXISTS GeneSummits(gene INTEGER, summit INT, distance INT)")
        cur.execute("CREATE TABLE IF NOT EXISTS RepgroupGenes(repgroupid INTEGER, geneid INTEGER)") # genes that have summits in all replicates of a repgroup
        cur.execute("CREATE TABLE IF NOT EXISTS SpeciesGenes(speciesid INTEGER, geneid INTEGER)") # genes that have summits in all repgroups of a species
        
        
    return con

def get_species_ids(con):
    cur = con.cursor()                
    cur.execute("SELECT id FROM Species")
    return cur.fetchall()

def get_speciesid_for_rep(repid, con):
    cur = con.cursor()
    cur.execute("SELECT id FROM Species where id in (SELECT species from Replicates where id=" + repid.__str__() + ")")
    return cur.fetchone()[0]

def get_chrom_ids(con, speciesid):
    cur = con.cursor()
    cur.execute("SELECT id FROM Chromosomes where species=" + speciesid.__str__())
    return cur.fetchall()

def get_chrom_id(con, name):
    cur = con.cursor()
    cur.execute("SELECT id from Chromosomes where name='" + name + "'")
    x = cur.fetchone()
    if x == None:
        return None
    else:
        return x[0]

def get_repid(repname, speciesid, con):    
    cur = con.cursor()
    cur.execute("SELECT id from Replicates where name='" + repname.__str__() + "' and species=" + speciesid.__str__())
    return cur.fetchone()[0]

def get_repgroup_ids(con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups")
    return cur.fetchall()

def get_repgroup_id(rgroup, con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups where name='" + rgroup + "'")
    x = cur.fetchone()
    if x != None:
        return x[0]
    else:
        return None

def get_repgroup_name(rgroupid, con):
    cur = con.cursor()
    cur.execute("SELECT name from ReplicateGroups where id=" + rgroupid.__str__() )
    return cur.fetchone()[0]

def get_reps_in_group(rgroupid, con):
    cur = con.cursor()
    cur.execute("SELECT replicate from GroupReplicate where rgroup=" + rgroupid.__str__())
    return cur.fetchall()

def get_rgroupids_for_species(speciesid, con):
    cur = con.cursor()
    cur.execute("SELECT id from ReplicateGroups where id in (SELECT rgroup from GroupReplicate where replicate in (SELECT id from Replicates where species=" + speciesid.__str__() + "))")
    x = cur.fetchall()
    ids = []
    for ii in x:
        ids.append( ii[0] )
    return ids

def get_geneids(con, repid):
    cur = con.cursor()
    sql = "SELECT gene FROM GeneSummits where summit in (SELECT id from Summits where replicate=" + repid.__str__() + ")"
    cur.execute(sql)
    return cur.fetchall()

def get_genes_for_species(con, speciesid):
    cur = con.cursor()
    sql = "SELECT * from Genes where chrom in (SELECT id from Chromosomes where species=" + speciesid.__str__() + ")"
    cur.execute(sql)
    return cur.fetchall()

def get_genes_for_chrom(con, chromid):
    """Returns list of genes (represented as tuples), in order by their start sites"""
    cur = con.cursor()
    sql = "SELECT * FROM Genes where chrom=" + chromid.__str__() + " order by start ASC"
    #print "46:", sql
    cur.execute(sql)
    return cur.fetchall()
 
def get_geneid_for_aliasname(aliasname, con):
    """Returns the geneID that is canonical for the given gene name alias.
    The mapping of alias names to gene ID comes from the pillars file."""
    realname = None
    cur = con.cursor()
    sql = "SELECT realname from GeneAlias where alias='" + aliasname + "'"
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        print "\n. Warning: the gene name", aliasname, "has no alias mappings."
        realname = aliasname
    else:
        realname = x[0]
    sql = "SELECT id from Genes where name='" + realname + "'"
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        print "\n. Error: The gene named", realname, "cannot be found."
        return None
    else:
        return x[0]
        
def get_geneid_for_aliasid(aliasid, con):
    cur = con.cursor()
    sql = "SELECT name from Genes where id=" + aliasid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    if x == None:
        print "\n. Warning: the gene with ID", aliasid, "cannot be found."
        exit()
    aliasname = x[0]
    return get_geneid_for_aliasname(aliasname, con)
 
def get_geneids_from_repgroup(con, repgroupid):
    cur = con.cursor()
    sql = "SELECT geneid from RepgroupGenes where repgroupid=" + repgroupid.__str__()
    cur.execute(sql)
    x = cur.fetchall()
    genes = []
    for ii in x:
        genes.append( ii[0] )
    return genes

def get_summits(con, repid, chromid):
     cur = con.cursor()
     sql = "SELECT * FROM Summits where chrom=" + chromid.__str__() + " and replicate=" + repid.__str__()
     #print "53:", sql
     cur.execute(sql)
     x = cur.fetchall()
     if x == None:
         return []
     else:
         return x

def get_summit_scores_for_gene(geneid, repid, con):
    cur = con.cursor()
    sql = "select score from Summits where replicate=" + repid.__str__() + " and id in (select summit from GeneSummits where gene=" + geneid.__str__() + ")"
    #print "60:", sql
    cur.execute(sql)
    scores = []
    for s in cur.fetchall():
        scores.append( s[0] )
    return scores

def get_max_summit_score_for_gene(geneid, repid, con):
    scores = get_summit_scores_for_gene(geneid, repid, con)
    if scores.__len__() == 0:
        return None 
    return max( scores )

def get_enrichment_stats_for_gene(geneid, repid, con):
    cur = con.cursor()
    sql = "SELECT * from EnrichmentStats where repid=" + repid.__str__() + " and geneid=" + geneid.__str__()
    cur.execute(sql)
    return cur.fetchone()
    

def import_species(speciesname, con):
    cur = con.cursor()
    sql = "INSERT INTO Species (name) VALUES('" + speciesname + "')"
    cur.execute( sql )
    con.commit()
    return con

def get_species_id(speciesname, con):
    cur = con.cursor()
    sql = "SELECT id FROM Species where name='" + speciesname + "'"
    cur.execute(sql)
    return cur.fetchone()[0]

def get_species_name(speciesid, con):
    cur = con.cursor()
    sql = "SELECT name FROM Species where id=" + speciesid
    cur.execute(sql)
    return cur.fetchone()[0]

def import_gff(gffpath, speciesid, con, restrict_to_feature = "gene"):
    print "\n. Importing genome features from", gffpath
    
    """Returns chr_gene_sites"""
    cur = con.cursor()
    
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
                sql = "INSERT INTO Genes (name, start, stop, chrom, strand) VALUES('" + gene + "'," + start.__str__() + "," + stop.__str__() + "," + chrid.__str__() + ",'" + strand + "')"
                cur.execute(sql)
    fin.close() 
    
    cur.execute("SELECT COUNT(*) FROM Genes")
    count_genes = cur.fetchone()[0]
    #print "\n. The database now contains", count_genes, "total genes."
    
    return con

def add_repgroup(rgroup, con, note=None):
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM ReplicateGroups WHERE name='" + rgroup + "'")
        count = cur.fetchone()[0]
        if count == 0:
            if note == None:
                note = "None"
            sql = "INSERT INTO ReplicateGroups (name, note) VALUES('" + rgroup.__str__() + "','" + note + "')"
            cur.execute(sql)
            con.commit()
    return con        
    
def add_replicate(repname, speciesid, con):
    new_id = None
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM Replicates WHERE name='" + repname.__str__() + "' and species=" + speciesid.__str__() )
        data = cur.fetchone()
        if data[0] == 0: # that chromosome doesn't exist yet.
            sql = "INSERT INTO Replicates (name,species) VALUES('" + repname.__str__() + "'," + speciesid.__str__() + ")"
            cur.execute( sql )
            con.commit()
    return con


def add_rep2group(repid, rgroupid, con):
    with con:
        cur = con.cursor()
        cur.execute("SELECT COUNT(*) FROM GroupReplicate WHERE rgroup=" + rgroupid.__str__() + " and replicate=" + repid.__str__())
        count = cur.fetchone()[0]
        if count == 0:
            sql = "INSERT INTO GroupReplicate (rgroup,replicate) VALUES(" + rgroupid.__str__() + "," + repid.__str__() + ")"
            cur.execute(sql)
            con.commit()
    return con

def add_union(unionname, repgroupnames, con):
    cur = con.cursor()
    cur.execute("SELECT COUNT(*) FROM Unions where name='" + unionname + "'")
    count = cur.fetchone()[0]
    if count == 0:
        sql = "INSERT into Unions (name) VALUES('" + unionname + "')"
        cur.execute(sql)
        con.commit()
        
    # now get the union's id
    cur.execute("SELECT unionid from Unions where name='" + unionname + "'")
    unionid = cur.fetchone()[0]
    for repgroupname in repgroupnames:
        rgroupid = get_repgroup_id(repgroupname, con)
        print "288:", rgroupid
        if rgroupid != None:
            sql = "INSERT INTO UnionRepgroups (unionid, repgroupid) VALUES(" + unionid.__str__() + "," + rgroupid.__str__() + ")"
            print sql
            cur.execute(sql)
            con.commit()
    return con

def get_unionids(con):
    cur = con.cursor()
    cur.execute("SELECT unionid from Unions")
    x = cur.fetchall()
    if x == None:
        return None
    unionids = []
    for ii in x:
        unionids.append(ii[0])
    return unionids

def get_unionname(unionid, con):
    cur = con.cursor()
    cur.execute("SELECT name from Unions where unionid=" + unionid.__str__())
    x = cur.fetchone()
    if x == None:
        return None
    return x[0]

def get_repgroupds_in_union(unionid, con):
    cur = con.cursor()
    cur.execute("SELECT repgroupid from UnionRepgroups where unionid=" + unionid.__str__() )
    x = cur.fetchall()
    if x == None:
        return None
    repgroupids = []
    for ii in x:
        repgroupids.append( ii[0] )
    return repgroupids

def import_pillars(pillarspath, con):
    print "\n. Reading Pillars from", pillarspath
    cur = con.cursor()
    fin = open(pillarspath, "r")
    count = 0
    for l in fin.xreadlines():
        if l.startswith("orf"):
            count += 1
            if count%20 == 0:
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
                    print "\n. I can't find chromosome", chromname, "in the descriptions provided by your GFF."
                    print "\n. Your bedgraph file", bdgpath,"references chromosome", chromname, ". I'm skipping this BedGraph entry."
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
                
     