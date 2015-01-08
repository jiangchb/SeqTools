#
# Draw a subsample of reads from a FASTQ file.
#
import sys, os
from argParser import ArgParser
ap = ArgParser(sys.argv)

inpath = ap.getArg("--in")
outpath = ap.getArg("--out")
n = int( ap.getArg("--n") ) 

fin = open(inpath, "r")
fout = open(outpath, "w")

count_reads = 0
lc = 0
for l in fin.xreadlines():
    if lc == 0:
        count_reads += 1

    if count_reads <= n:
        fout.write( l )
    else:
        break
    
    lc += 1
    if lc == 4:
        lc = 0
    
    if count_reads%10 == 0:
        sys.stdout.write("\r    --> %.1f%%" % (100*count_reads/float(n)) )
        sys.stdout.flush()
    

fin.close()
fout.close()
