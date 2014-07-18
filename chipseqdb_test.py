#
# A test harness for chipseqdb.py
#

from chipseqdb import *

from argParser import ArgParser
ap = ArgParser(sys.argv)


con = build_db("test.db")
gffpath = ap.getArg("--gff")
summits = ap.getList("--summits")
species = ["C. albicans"]

for sp in species:
    speciesid = 1
    con = import_gff(gffpath, speciesid, con)
    repid = 0
    for s in summits:
        con = import_summits(s, repid, con)
        con = map_summits2genes(con, repid, speciesid=1)
        repid += 1
