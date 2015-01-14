#
# bedgraph2wig.py
#
# Written by Victor Hanson-Smith
# victorhansonsmith@gmail.com
#
# This script will convert BEDgraph files into WIG files.
#

import os, re, sys

def usage():
    print "\n"
    print "====================================================="
    print "bedgraph2wig.py"
    print ""
    print "written by Victor Hanson-Smith"
    print "victorhansonsmith@gmail.com"
    print ""
    print "USAGE:"
    print "python bedgraph2wig.py FILEIN FILEOUT"
    print ""
    print "... where FILEIN is the filepath to a bedgraph"
    print "... where FILEOUT is the filepath to a WIG file,"
    print "      that will be created (or overwritten) by this script."
    print ""
    print "This script will parse the bedgraph and create a new"
    print "WIG file for every chromosome found within the bedgraph."
    print ""
    print "====================================================="

if len( sys.argv ) < 3:
    usage()
    exit()

def splash():
    print "\n"
    print "======================================"
    print "bedgraph2wig.py"
    print ""
    print "written by Victor Hanson-Smith"
    print "victorhansonsmith@gmail.com"
    print "======================================"

bedgraphpath = sys.argv[1]
if False == os.path.exists(bedgraphpath):
    print "I can't find you bedgraph file at " + bedgraphpath
    exit()

wigpath = sys.argv[2]


printspan = 100000 # print an update every N sites
count = 0
last_seen_chrom = None
fin = open(bedgraphpath, "r")
print "\n. Converting BDG->WIG:", bedgraphpath, "-->", wigpath
fout = open(wigpath, "w")
fout.write("track type=WIG\n")
for l in fin.xreadlines():
    if l.__len__() > 2:
        tokens = l.split()
        chromname = tokens[0]

        if last_seen_chrom != None and chromname != last_seen_chrom:
            last_seen_chrom = None


        if last_seen_chrom == None:            
            fout.write("variableStep chrom=" + chromname + "\n")
            last_seen_chrom = chromname
        start = int( tokens[1] )
        stop = int( tokens[2] )
        value = tokens[3]
        for ii in range(start, stop):
            if count%printspan == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            count += 1
            fout.write(ii.__str__() + "\t" + value + "\n")
fin.close()
fout.close()
