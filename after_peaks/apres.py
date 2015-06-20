#
# APRES ChIP-Seq
# (Analaysis Pipeline for Replicate ExperimentS)
#

from read_config import *
from chipseqdb import *
from chipseqdb_plots import *
from overlap_stats import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

def splash():
    print "============================================================================"
    print "."
    print ". APRES: Analysis Pipeline for Replicate Experiments across Species."
    print "."
    print ". [Part of the ChIP-Seq Distillery Suite]"
    print "."
    print ". Last Updated January 2014"
    print "."
    print ". Written by Victor Hanson-Smith -- victorhansonsmith@gmail.com"
    print "============================================================================="
#
# Build, or add to, the database
#
def import_gffs(con):
    for sp in ap.params["species"]:
        speciesname = ap.params["species"][sp]["name"]
        
        if False == does_species_exist(speciesname, con):
            con = add_species(speciesname, con)
        speciesid = get_species_id(speciesname, con)
                
        gffpath = ap.params["species"][sp]["gff"]
        if False == os.path.exists(gffpath):
            print "\n. Error: I can't find your genome feature file at", gffpath
            exit()
        
        chrom_restrict_list = ap.getOptionalList("--skip_chroms")
        con = import_gff(gffpath, speciesid, con, filter_chrom=chrom_restrict_list)

def import_data(con):    
    """If this analysis uses pillars to translate gene names, then those values
    must be imported prior to this method."""
    cur = con.cursor()
    
    """Remove any output files, because they are now stale."""
    sql = "DELETE FROM Files"
    cur.execute(sql)
    con.commit()
    
    """Remove  stale data."""
    sql = "DELETE FROM GeneSummits"
    cur.execute(sql)
    con.commit()
        
    for sp in ap.params["species"]:
        #print "\n. Importing data for species", sp
        speciesname = ap.params["species"][sp]["name"]
        
        if False == does_species_exist(speciesname, con):
            con = add_species(speciesname, con)
        speciesid = get_species_id(speciesname, con)
                        
        ii = -1
        for groupname in ap.params["species"][sp]["rgroups"]:
            if False == does_repgroup_exist(groupname, con):
                con = add_repgroup(groupname, con)
            rgroupid = get_repgroup_id(groupname, con)
                        
            for jj in ap.params["species"][sp]["rgroups"][groupname]["reps"]:
                print "Replicate:", jj, ", Species:", sp, speciesname, speciesid, ", RepGroup:", groupname, rgroupid
                ii += 1

                repname = groupname + "-" + jj.__str__()
                con = add_replicate(repname, speciesid, con)
                repid = get_repid(repname, speciesid, con)
                con = add_rep2group(repid, rgroupid, con)
                
                if "summitpath" in ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]:
                    summitpath = ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]["summitpath"]
                    con = import_summits(summitpath, repid, con)
                    con = map_summits2genes(con, repid, speciesid=speciesid)
                
                if "bdgpath" in ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]:
                    bdgpath = ap.params["species"][sp]["rgroups"][groupname]["reps"][jj]["bdgpath"]
                    #validate_enrichment(bdgpath, repid, con)
                    #print "91:", speciesname, speciesid, rgroupid, repname, repid
                    con = import_foldenrichment(bdgpath, repid, con)                
                
                """Check that all summits have a corresponding FE value,
                    and that all genes have an FE value."""
                validate_summits_fe(repid, con)     
    return con

def setup_unions(con):    
    clear_unions(con)
    build_unions(con)
    for sp in ap.params["species"]:
        if "unions" in ap.params["species"][sp]:
            for unionname in ap.params["species"][sp]["unions"]:
                repgroupnames = ap.params["species"][sp]["unions"][unionname]
                union_type = ap.params["union_type"][unionname]
                con = add_union(unionname, union_type, repgroupnames, con)

def setup_speciesunions(con):
    clear_speciesunions(con)
    build_speciesunions(con)  
    if "unionssp" in ap.params:
        for u in ap.params["unionssp"]:
            member_unions = ap.params["unionssp"][u]
            #print "\n. Adding a species union:", u, member_unions
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
    cur = con.cursor()
        
    rgroups = get_repgroup_ids(con)
    for rg in rgroups:
        rgroupid = rg[0]
        sql = "DELETE FROM Files where fileid in (SELECT fileid FROM ReplicategroupFiles WHERE repgroupid=" + rgroupid.__str__() + ")"
        cur.execute(sql)
        sql = "DELETE FROM ReplicategroupFiles where repgroupid=" + rgroupid.__str__()
        cur.execute(sql)
        con.commit()
        #plot_summits_for_reps_in_group(rgroupid, con)
        plot_fexfe_replicates(rgroupid, con)
        plot_enrichments_for_reps_in_group(rgroupid, con)
    return con

def correlate_unions(con):
    """Correlates a set of replicate groups, as defined by UNION commands in the configuration file."""
    unions = get_unionids(con)
    for uid in unions:
        compute_summits_for_union( uid, con )
        compute_enrichments_for_union(uid, con)
        
def plot_unions(con):
    cur = con.cursor()
    unions = get_unionids(con)
    for uid in unions:
        sql = "DELETE FROM Files where fileid in (SELECT fileid FROM UnionFiles WHERE unionid=" + uid.__str__() + ")"
        cur.execute(sql)
        sql = "DELETE FROM UnionFiles where unionid=" + uid.__str__()
        cur.execute(sql)
        con.commit()
        plot_summits_for_union( uid, con )
        plot_enrichments_for_union(uid, con)

def correlate_species(con):
    """Correlates the unions for each species-union, that is, sets of unions that cross species boundaries."""
    species_unions = get_species_unionids(con)
    for spuid in species_unions:
        compute_summits_for_speciesunion( spuid, con)
        compute_enrichments_for_speciesunion( spuid, con)

def plot_species(con):
    """Correlates the unions for each species-union, that is, sets of unions that cross species boundaries."""
    cur = con.cursor()
    species_unions = get_species_unionids(con)
    for spuid in species_unions:
        sql = "DELETE FROM Files where fileid in (SELECT fileid FROM SpeciesunionFiles WHERE spunionid=" + spuid.__str__() + ")"
        cur.execute(sql)
        sql = "DELETE FROM SpeciesunionFiles where spunionid=" + spuid.__str__()
        cur.execute(sql)
        con.commit()
        plot_summits_for_speciesunion( spuid, con)
        plot_enrichments_for_speciesunion( spuid, con)
        

def package_results(con):
    path_note = {}
    cur = con.cursor()
    sql = "select * from Files"
    cur.execute(sql)
    x = cur.fetchall()
    if x != None:
        for ii in x:
            path = ii[1]
            note = ii[2]
            if path not in path_note:
                path_note[path] = note

    files = path_note.keys()
    files.sort()
    fout = open("INVENTORY.xls", "w")
    for f in files:
        if os.path.exists(f):
            fout.write(f + "\t" + path_note[f]  + "\n")
    fout.close()
    files.append("INVENTORY.xls")

    print "\n. Packaging the results into a ZIP file."
    print "    ", files.__len__(), "files."
        
    command = "tar czf results.tar.gz "
    for f in files:
        if os.path.exists(f):
            command += f + " "
    os.system(command)
    print "\n. Results are packaged into results.tar.gz"

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

"""Regardless of which options are being executed, the DB gets
built, or rebuilt, depending on its status."""
con = build_db(dbpath=dbpath)

"""--make_testdb is a dev-only option"""
if True == ap.getOptionalToggle("--make_testdb"):
    reduce_db_for_test(con)
    exit()

"""--reset_files is a dev-only option"""
if True == ap.getOptionalToggle("--reset_files"):
    reset_files(con)
    exit()

# 
# IMPORT
#
pillarspath = ap.getArg("--pillarspath")
if False == ap.getOptionalToggle("--skip_pillars_import") and False == ap.getOptionalToggle("--skip_import"):
    con = import_pillars(pillarspath, con)
else:
    print "\n. I'm skipping pillars import. I will use the existing pillars in the database."
print_pillarsstats(con)

if configpath != False:
    if False == ap.getOptionalToggle("--skip_gff") and False == ap.getOptionalToggle("--skip_import"):
        import_gffs(con)
        resolve_aliasids(con)
    
    if False == ap.getOptionalToggle("--skip_import"):
        con = import_data(con)

#
# ANALYSIS
#
if False == ap.getOptionalToggle("--skip_analysis"):
    if False == ap.getOptionalToggle("--skip_reps"):
        clear_repgroups(con)
        build_repgroups(con)
        correlate_replicates(con)
    
    if False == ap.getOptionalToggle("--skip_unions"):
        clear_unions(con)
        build_unions(con)
        
        setup_unions(con)
        correlate_unions(con)
         
    if False == ap.getOptionalToggle("--skip_species"):
        setup_speciesunions(con)
        correlate_species(con)    

#IDR

#
# continue here:
# plot FE versus summit score.
#

#
# PLOTS
#
if False == ap.getOptionalToggle("--skip_plots"):
    if False == ap.getOptionalToggle("--skip_repsplot"):
        plot_replicates(con)        
    if False == ap.getOptionalToggle("--skip_unionsplot"):
        plot_unions(con)
    if False == ap.getOptionalToggle("--skip_speciesplot"):
        plot_species(con)

#
# More Stats
#
#test_species_overlap(con)

#
# PACKAGE the results    
#
if False == ap.getOptionalToggle("--skip_zip"):
    add_file(dbpath, "The Sqlite3 database file for this project.", con)
    if configpath != False:
        add_file(configpath, "The configuration file for this project.", con)
    package_results(con)
    