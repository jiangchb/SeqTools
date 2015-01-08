import os, sys
from argParser import ArgParser
ap = ArgParser(sys.argv)

gffpath = ap.getArg("--gff")
if False == os.path.exists(gffpath):
    print "\n. Error. I can't find your GFF at", gffpath
    exit()

bdgpath = ap.getArg("--bdg")
if bdgpath != False:
    if False == os.path.exists(bdgpath):
        print "\n. Error. I can't find your BDG at", bdgpath
        exit()

genename = ap.getArg("--genename")

regsize = int( ap.getArg("--regsize") )

chrom = None
start = None
stop = None
strand = None
fin = open(gffpath, "r")
for l in fin.xreadlines():
    if l.__len__() > 0 and False == l.startswith("#"):
        tokens = l.split()
        this_gene = tokens[8].split(";")[0].split("=")[1] # orfName 
        if this_gene == genename:
            chrom = tokens[0]
            start = int( tokens[3] )
            stop = int( tokens[4] )
            strand = tokens[6]
            if strand == "-":
                x = start
                y = stop
                start = y
                stop = x


print "\n. I found Gene", genename, "chrom:", chrom
print ". start:", start, "stop:", stop, "strand:", strand

site_fe = {}
max_fe = None
max_fe_site = None
min_fe = None
fin = open(bdgpath, "r")
for l in fin.xreadlines():
    if l.__len__() < 2:
        continue
    elif False == l.startswith(chrom):
        continue
    tokens = l.split()
    this_site_start = int(tokens[1])
    this_site_stop = int(tokens[2])
    if (strand == "+" and start - this_site_start <= regsize and start > this_site_start) or (strand == "-" and this_site_start - start <= regsize and start < this_site_start):
        this_score = float(tokens[3])
        for ii in range(this_site_start, this_site_stop):
            site_fe[ ii ] = this_score
            if max_fe == None:
                max_fe = this_score
                max_fe_site = ii
                min_fe = this_score
            if max_fe < this_score:
                max_fe = this_score
                max_fe_site = ii
            if min_fe > this_score:
                min_fe = this_score

"""Display the data on the command-line"""
fe_width = 50.0
fe_step = (max_fe - min_fe) / fe_width

count = 0
sites = site_fe.keys()
sites.sort()
for site in sites:
    if count%40 == 0:
        header = "Site # \tD(tss)\t"
        for ii in range(0, int(fe_width) ):
            header += "-"
        print header
    count += 1
    
    l = site.__str__() + "\t" + (site-start).__str__() + "\t"
    for ii in range(0, int( (site_fe[site] - min_fe)/fe_step) ):
        l += " "
    l += "+ (" + site_fe[site].__str__() + ")"
    print l
print "============================="
print "max FE = ", max_fe, "at site", max_fe_site, (max_fe_site-start)



