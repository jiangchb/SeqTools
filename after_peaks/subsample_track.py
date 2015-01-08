#
# This script is for generating toy-sized test cases
#
# Subsample a track file, such as BDGs and BEDs
#

import sys, os
from argParser import ArgParser
ap = ArgParser(sys.argv)

bdgpath = ap.getArg("--in")
outpath = ap.getArg("--out") # output path of subsampled BDG

# Lines in the BDG that contain the keyword will be sampled.
# all other lines will be discarded
keyword = ap.getArg("--keyword")

fout = open(outpath, "w")

fin = open(bdgpath, "r")
for l in fin.xreadlines():
    if l.__contains__(keyword):
        fout.write( l )
fin.close()
fout.close()