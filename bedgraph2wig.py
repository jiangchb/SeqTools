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
    print "python bedgraph2wig.py FILE"
    print ""
    print "... where FILE is the filepath to a bedgraph"
    print ""
    print "This script will parse the bedgraph and create a new"
    print "WIG file for every chromosome found within the bedgraph."
    print ""
    print "The new WIG files will be named FILE.CHROM.wig"
    print "...where FILE is the path to the original bedgraph,"
    print "and CHROM is the name of the chromosome."
    print "====================================================="

if len( sys.argv ) < 2:
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

# The output data will be built into the following hashtable:
#site_value = {}

# We'll keep a second list of the seen sites. I realize this list could
# also be acquired from site_value.leys(), but by keeping this second
# list we avoid the need to enumerate through a *large* number of hash
# keys.
#sites = []

printspan = 100000 # print an update every N sites
count = 0
last_seen_chrom = None
fin = open(bedgraphpath, "r")
wigpath = re.sub(".bdg", ".wig", bedgraphpath)
print "\n. Converting BDG->WIG:", bedgraphpath, "-->", wigpath
fout = open(wigpath, "w")
fout.write("track type=WIG\n")
for l in fin.xreadlines():
    if l.__len__() > 2:
        tokens = l.split()
        #print tokens
        chromname = tokens[0]

        if last_seen_chrom != None and chromname != last_seen_chrom:
            #fout.close()
            last_seen_chrom = None


        if last_seen_chrom == None:
            #outpath = bedgraphpath + "." + chromname + ".wig"
            #fout = open( outpath, "w")
            
            fout.write("variableStep chrom=" + chromname + "\n")
            last_seen_chrom = chromname
            #print "\n\n.Writing a new WIG file for chromosome " + chromname + "\n -> " + outpath

        start = int( tokens[1] )
        stop = int( tokens[2] )
        value = tokens[3]
        for ii in range(start, stop):
            #sites.append( ii )
            #site_value[ii] = value
            if count%printspan == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
            count += 1
            fout.write(ii.__str__() + "\t" + value + "\n")
fin.close()
fout.close()
