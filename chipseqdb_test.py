#
# A test harness for chipseqdb.py
#

from read_config import *
from chipseqdb import *
from chipseqdb_plots import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

def splash():
    print "==================================="
    print ". Project Guacamole"
    print ". . . by Victor Hanson-Smith"
    print ". . . victorhansonsmith@gmail.com"
    print "==================================="
#
# Build, or add to, the database
#
def import_gffs(con):
    for sp in ap.params["species"]:
        speciesname = ap.params["species"][sp]["name"]
        
        if False == does_species_exist(speciesname, con):
            con = import_species(speciesname, con)
        speciesid = get_species_id(speciesname, con)
                
        gffpath = ap.params["species"][sp]["gff"]
        con = import_gff(gffpath, speciesid, con)

def import_data(con):
    """If this analysis uses pillars to translate gene names, then those values
    must be imported prior to this method."""
    
    """Remove  stale data."""
    cur = con.cursor()
    sql = "DELETE FROM GeneSummits"
    cur.execute(sql)
    con.commit()
    
    for sp in ap.params["species"]:
        print "\n. Importing data for species", sp
        speciesname = ap.params["species"][sp]["name"]
        
        if False == does_species_exist(speciesname, con):
            con = import_species(speciesname, con)
        speciesid = get_species_id(speciesname, con)
                        
        ii = -1
        for groupname in ap.params["species"][sp]["rgroups"]:
            if False == does_repgroup_exist(groupname, con):
                con = add_repgroup(groupname, con)
            rgroupid = get_repgroup_id(groupname, con)
                        
            for jj in ap.params["species"][sp]["rgroups"][groupname]["reps"]:
                #print "Replicate:", jj, ", Species:", sp, speciesname, speciesid, ", RepGroup:", groupname, rgroupid
                ii += 1

                repname = groupname + "-" + jj.__str__()
                if False == does_replicate_exist(repname, speciesid, con):
                    con = add_replicate(repname, speciesid, con)
                repid = get_repid(repname, speciesid, con)
                con = add_rep2group(repid, rgroupid, con)
                
                if "summitpath" in ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]:
                    summitpath = ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]["summitpath"]
                    con = import_summits(summitpath, repid, con)
                if "bdgpath" in ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]:
                    bdgpath = ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]["bdgpath"]
                    con = import_bdg(bdgpath, repid, con)
                
                con = map_summits2genes(con, repid, speciesid=speciesid)
    return con

def setup_unions(con):    
    clear_unions(con)
    build_unions(con)
    #print "67:", ap.params["species"]
    for sp in ap.params["species"]:
        #print "68:", ap.params["species"][sp]
        for unionname in ap.params["species"][sp]["unions"]:
            repgroupnames = ap.params["species"][sp]["unions"][unionname]
            print "\n. Adding a repgroup union:", unionname, repgroupnames
            con = add_union(unionname, repgroupnames, con)

def setup_speciesunions(con):
    clear_speciesunions(con)
    build_speciesunions(con)  
    if "unionssp" in ap.params:
        for u in ap.params["unionssp"]:
            member_unions = ap.params["unionssp"][u]
            print "\n. Adding a species union:", u, member_unions
            con = add_speciesunion(u, member_unions, con)
    
def correlate_replicates(con):
    """Correlates all the replicates within each REPGROUP."""
    rgroups = get_repgroup_ids(con)
    for rg in rgroups:
        rgroupid = rg[0]
        compute_summits_for_reps_in_group(rgroupid, con)
        compute_enrichments_for_reps_in_group(rgroupid, con)
    return con

def plot_replicates(con):
    rgroups = get_repgroup_ids(con)
    for rg in rgroups:
        rgroupid = rg[0]
        plot_summits_for_reps_in_group(rgroupid, con)
        plot_enrichments_for_reps_in_group(rgroupid, con)
    return con

def correlate_unions(con):
    """Correlates a set of replicate groups, as defined by UNION commands in the configuration file."""
    unions = get_unionids(con)
    for uid in unions:
        compute_summits_for_union( uid, con )
        compute_enrichments_for_union(uid, con)
        print "107:", uid
        
def plot_unions(con):
    unions = get_unionids(con)
    for uid in unions:
        plot_summits_for_union( uid, con )
        print "113", uid
        plot_enrichments_for_union(uid, con)
    print "114"

def correlate_species(con):
    """Correlates the unions for each species-union, that is, sets of unions that cross species boundaries."""
    species_unions = get_species_unionids(con)
    for spuid in species_unions:
        compute_summits_for_speciesunion( spuid, con)
        compute_enrichments_for_speciesunion( spuid, con)

def plot_species(con):
    """Correlates the unions for each species-union, that is, sets of unions that cross species boundaries."""
    species_unions = get_species_unionids(con)
    for spuid in species_unions:
        plot_summits_for_speciesunion( spuid, con)
        plot_enrichments_for_speciesunion( spuid, con)

######################################################################
#
# main
#
splash()

"""If the user gives a configpath, then it's referenced
data will be imported into the database."""
configpath = ap.getOptionalArg("--configpath")

if configpath != False:
    ap.params = read_config( configpath )
    
    if False == validate_config(ap.params):
        print "\n. Error: something is wrong with your configuration file."
        print ". Look at previous errors for more detail."
        exit(0)
    
"""If the user didn't give a configpath, then we at
least need a dbpath (to load a prior database)."""
dbpath = ap.getOptionalArg("--dbpath")

if configpath == False and dbpath == False:
    print "\n. Error, you need to specify either:"
    print "1. a configuration path, using --configpath, that describes data that will be imported into the DB."
    print "2. an existing database, using --dbpath, whose contents will be read and imported."
    print ""
    exit()
    
con = build_db(dbpath=dbpath)


#
# 
#
if True == ap.getOptionalToggle("--make_testdb"):
    reduce_db_for_test(con)
    exit()

#
# IMPORT
#
pillarspath = ap.getArg("--pillarspath")
if False == ap.getOptionalToggle("--skip_pillars_import"):
    con = import_pillars(pillarspath, con)

if configpath != False:
    if False == ap.getOptionalToggle("--skip_gff"):
        import_gffs(con)
        resolve_aliasids(con)
    
    if False == ap.getOptionalToggle("--skip_import"):
        con = import_data(con)

#
# ANALYSIS
#
if False == ap.getOptionalToggle("--skip_reps"):
    correlate_replicates(con)

if False == ap.getOptionalToggle("--skip_unions"):
    setup_unions(con)
    
    # for debugging:
    cur = con.cursor()
    sql = "SELECT * from UnionRepgroups"
    cur.execute(sql)
    print "194:", cur.fetchall()

    correlate_unions(con)
    
    # for debugging:
    cur = con.cursor()
    sql = "SELECT * from UnionRepgroups"
    cur.execute(sql)
    print "204:", cur.fetchall()
    
if False == ap.getOptionalToggle("--skip_species"):
    setup_speciesunions(con)
    correlate_species(con)    

#
# PLOTS
#
if False == ap.getOptionalToggle("--skip_repsplot"):
    print "204"
    plot_replicates(con)
if False == ap.getOptionalToggle("--skip_unionsplot"):
    print "207"
    plot_unions(con)
if False == ap.getOptionalToggle("--skip_speciesplot"):
    print "210"
    plot_species(con)