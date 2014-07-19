#
# A test harness for chipseqdb.py
#

from chipseqdb import *
from chipseqdb_plots import *

from argParser import ArgParser
ap = ArgParser(sys.argv)


con = build_db("test.db")
gffpath = ap.getArg("--gff")
summits = ap.getList("--summits")
species = ["C. albicans"]

for sp in species:
    speciesid = 1
    repid = 0
    con = import_gff(gffpath, speciesid, con)
    for s in summits:
        con = add_replicate(repid, speciesid, con)
        con = import_summits(s, repid, con)
        con = map_summits2genes(con, repid, speciesid=1)
        repid += 1

print get_max_summit_score_for_gene(644, 1, con)
correlate_two_reps(0, 1, con)