#
# This script processes ChIP-Seq results from FASTQ to (mapped) SAM to (sorted) BAM
# to fold-enrichment track and summit list.
#

import re, os, sys
from read_annotation import *

"""In what directory do the FASTQ files exist?"""
DATADIR = "/Network/Servers/udp015817uds.ucsf.edu/Users/eugeniom/Documents/ChIP-seq_analysis/Ndt80/"

"""In what directory should output be written?"""
OUTDIR = "./"

"""Path to annocation file, describining FASTQ files."""
APATH = "mancera_chipseq_sample_annotation.txt"

"""Analyze rows with this sample ID from the annotation file."""
sample = "EM14918"


# id is library_name
id_fastq = {}
id_sampleid = {}
id_repid = {}
id_tag = {} # True/False
id_species = {}
id_sam = {}

con = read_anno(DATADIR + "/" + APATH)



fin = open(DATADIR + "/" + APATH, "r")

lines = []
for l in fin.readlines():
    ls = l.split("\r")
    for lt in ls:
        lines.append( lt )

for l in lines:
    tokens = l.split("\t")
    if tokens.__len__() > 2:
        if tokens[0] == sample:
            fastq = tokens[3]
            id = tokens[1]
            species = tokens[5]
            st = species.split("_")

            for species in st:
                if st.__len__() > 1:
                    id = tokens[1] + "-" + species
                    outsam = re.sub(".fastq", "-" + species + ".sam", fastq)
                else:
                    outsam = re.sub(".fastq", ".sam", fastq)
                id_sam[id] = outsam
                id_fastq[id] = fastq
                id_species[id] = species
                id_sampleid[id] = tokens[0]
                id_repid[id] = tokens[10]
                id_tag[id] = tokens[7]
                if id_tag[id] == "YES":
                    id_tag[id] = True
                else:
                    id_tag[id] = False
                 

                print id, fastq, species

#print id_sam
#print id_fastq
#print id_species
#print id_sampleid
#print id_repid
#print id_tag
exit()

"""Make a list of commands to launch the bowtie script (a.k.a. hybrid_runme.py)"""
commands = []
for sample in id_sampleid.values():
    for repid in id_repid



for (f,s) in fastq_sam:

    if False == os.path.exists(DATADIR + f):
        print "\n. I can't find " + DATADIR + f
        continue
    #if False == os.path.exists(DATADIR + s):
    #    print "\n. I can't find " + DATADIR + s
    #print s

    c = "python ~/Applications/SeqTools/hybrid_runme.py"
    if s.__contains__("Cdub"):
        c += " --iprefix /Users/Shared/sequencing_analysis/indexes/06-Nov-2013C_dubliniensis_CD36"
    if s.__contains__("Calb"):
        c += " --iprefix /Users/Shared/sequencing_analysis/indexes/06-Apr-2014C_albicans_SC5314"
    if s.__contains__("Ctro"):
        c += " --iprefix /Users/Shared/sequencing_analysis/indexes/11-Dec-2013C_tropicalis_MYA-3404_corrected"
    c += " --fastqpath " + DATADIR + f 
    c += " --samoutpath " + s
    #c += " --jump 1"
    #print c
    #os.system(c)
    commands.append(c)

fout = open("bowtie_commands.sh", "w")
for c in commands:
    fout.write(c + "\n")
fout.close()

os.system("source bowtie_commands.sh")

#os.system("mpirun -np 8 -machinefile hosts.txt /common/REPOSITORY/mpi_dispatch/mpi_dispatch bowtie_commands.sh")