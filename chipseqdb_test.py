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
def import_data(con):
    for sp in ap.params["species"]:
        print "\n. Importing data for species", sp
        speciesname = ap.params["species"][sp]["name"]
        con = import_species(speciesname, con)
        speciesid = get_species_id(speciesname, con)
                
        gffpath = ap.params["species"][sp]["gff"]
        con = import_gff(gffpath, speciesid, con)
        
        ii = -1
        for groupname in ap.params["species"][sp]["rgroups"]:
            con = add_repgroup(groupname, con)
            rgroupid = get_repgroup_id(groupname, con)
                        
            for jj in ap.params["species"][sp]["rgroups"][groupname]["reps"]:
                #print "Replicate:", jj, ", Species:", sp, speciesname, speciesid, ", RepGroup:", groupname, rgroupid
                ii += 1

                repname = groupname + "-" + jj.__str__()
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
    for sp in ap.params["species"]:
        for unionname in ap.params["species"][sp]["unions"]:
            repgroupnames = ap.params["species"][sp]["unions"][unionname]
            print "\n. Adding a repgroup union:", unionname, repgroupnames
            con = add_union(unionname, repgroupnames, con)

def correlate_replicates(con):
    """Correlates all the replicates within each REPGROUP."""
    rgroups = get_repgroup_ids(con)
    for rg in rgroups:
        rgroupid = rg[0]
        con = correlate_summits_for_reps_in_group(rgroupid, con)

def correlate_unions(con):
    """Correlates a set of replicate groups, as defined by UNION commands in the configuration file."""
    unions = get_unionids(con)
    for uid in unions:
        con = correlate_summits_for_union( uid, con )

def correlate_enrichments_species(speciesid1, speciesid2, con):
    series_gene_eval = {}
    chromids = get_chrom_ids(con, speciesid1)
    for chromid in chromids:
        geneids = get_geneids(con, chromid)
        #
        # continue here
        #
        for geneid in geneids:
            x = get_enrichment_stats_for_gene(geneid, repid, con)

    correlate_enrichments(series_gene_eval, con, keyword=None)

    #correlate_enrichments(series_gene_eval, con, keyword=None)

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
# IMPORT
#
pillarspath = ap.getOptionalArg("--pillarspath")
if pillarspath != False:
    con = import_pillars(pillarspath, con)

if configpath != False:
    skip_import = ap.getOptionalArg("--skip_import")
    if skip_import == False:
        con = import_data(con)
    # Setup union relationships, regardless
    con = setup_unions(con)

#
# ANALYSIS
#
if False == ap.getOptionalToggle("--skip_repcorr"):
    correlate_replicates(con)

if False == ap.getOptionalArg("--skip_groupcorr"):
    correlate_unions(con)
    
if False == ap.getOptionalArg("--skip_speciescorr"):
    correlate_enrichments_species(con)    
