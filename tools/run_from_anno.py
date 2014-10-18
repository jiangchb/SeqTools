#
# This script processes ChIP-Seq results from FASTQ to (mapped) SAM to (sorted) BAM
# to fold-enrichment track and summit list.
#

import re, os, sys, time
from annotation_db import *
from hybrid_tools import *
from tools import *
from argParser import ArgParser
ap = ArgParser(sys.argv)

"""In what directory do the FASTQ files exist?"""
DATADIR = "/Network/Servers/udp015817uds.ucsf.edu/Users/eugeniom/Documents/ChIP-seq_analysis/Ndt80/"
x = ap.getOptionalArg("--datadir")
if x != False:
    DATADIR = x

"""In what directory should output be written?"""
OUTDIR = "./"
x = ap.getOptionalArg("--outdir")
if x != False:
    OUTDIR = x

"""Path to annotation file, describining FASTQ files."""
APATH = "mancera_chipseq_sample_annotation.txt"
x = ap.getOptionalArg("--annopath")
if x != False:
    APATH = x
DBPATH = ap.getOptionalArg("--dbpath")
if DBPATH == False:
    DBPATH = APATH + ".db"
    
PROJECT_NAME = ap.getArg("--project_name")

"""Software Paths:"""
BOWTIE2 = "bowtie2"
SAMTOOLS = "samtools"
MPIRUN = "mpirun -np 6 --machinefile hosts.txt /common/bin/mpi_dispatch"
MACS2 = "/common/REPOSITORY/MACS2-2.0.10.07132012/bin/macs2"
SEQTOOLSDIR = "~/Applications/SeqTools"

def splash():
    print "============================================="
    print " run_from_anno.py"
    print " October 2014"
    print " A script to analyze ChIP-Seq data."
    print "=============================================" 

splash()


"""restricts analysis to only those annotation lines whose sample column is this value."""
restrict_to_sample = ap.getOptionalArg("--restrict_to_sample")

PRACTICE_MODE = ap.getOptionalToggle("--practice_mode")

"""Jump allows for some steps to be skipped."""
jump = ap.getOptionalArg("--jump")
if jump == False:
    jump = 0
else:
    jump = float(jump)
    
"""Build or restore the DB"""
con = get_db(DBPATH)
cur = con.cursor()

sql = "delete from Settings"
cur.execute(sql)
con.commit()
"""Write the settings into the DB"""
sql = "insert or replace into Settings (keyword, value) VALUES('bowtie2','" + BOWTIE2 + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('datadir','" + DATADIR + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('outdir','" + OUTDIR + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('samtools','" + SAMTOOLS + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('macs2','" + MACS2 + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('seqtoolsdir','" + SEQTOOLSDIR + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('mpirun','" + MPIRUN + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('restrict_to_sample','" + restrict_to_sample + "')"
cur.execute(sql)
sql = "insert or replace into Settings (keyword, value) VALUES('project_name','" + PROJECT_NAME + "')"
cur.execute(sql)
x = 0
if PRACTICE_MODE:
    print "\n. Practice Mode Enabled."
    x = 1
sql = "insert or replace into Settings (keyword, value) VALUES('practice_mode','" + x.__str__() + "')"
cur.execute(sql)
con.commit()

con = import_annotations(DATADIR + "/" + APATH, con, sample_restrict=get_setting("restrict_to_sample",con) )

sql = "select count(*) from Annotations"
cur.execute(sql)
count = cur.fetchone()[0]
print "\n. I found", count, "annotations in the database."

"""Run Bowtie2 for each FASTQ path"""
if jump <= 1:
    """Clear previous Bowtie output paths."""
    sql = "delete from BowtieOutput"
    cur.execute(sql)
    con.commit()
    
    bowtie_commands = []
    sql = "select annoid, fastqpath, species from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        c = get_setting("bowtie2",con)
        annoid = ii[0]
        fastq = ii[1]
        species = ii[2]
        
        """Sanity check:"""
        if False == os.path.exists( get_setting("datadir",con) + "/" + fastq ):
            print "\n. Error, I can't find your FASTQ file at", fastq
            print ". (run_from_anno.py 43)"
            exit()
        c += " -U " + fastq
        
        """Is it a hybrid?"""
        sql = "select count(*) from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        samoutpath = re.sub(".fastq", ".sam", fastq)
        if count > 0:
            """Hybrids get a special SAM path"""
            samoutpath = re.sub(".fastq", "-" + species + ".sam", fastq)
        c += " -S " + get_setting("outdir",con) + "/" + samoutpath
        c += " --no-unal "
        
        """Path to directory with genome sequences"""
        if species == "Cdub":
            c += " -x /Users/Shared/sequencing_analysis/indexes/06-Nov-2013C_dubliniensis_CD36"
        if species == "Calb":
            c += " -x /Users/Shared/sequencing_analysis/indexes/06-Apr-2014C_albicans_SC5314"
        if species == "Ctro":
            c += " -x /Users/Shared/sequencing_analysis/indexes/11-Dec-2013C_tropicalis_MYA-3404_corrected"
        bowtie_commands.append(c)
        #print annoid, species, fastq, samoutpath
        
        sql = "insert or replace into BowtieOutput (annoid, sampath) VALUES(" + annoid.__str__() + ",'" + samoutpath + "')"
        cur.execute(sql)
        con.commit()
        
    fout = open("bowtie_commands.sh", "w")
    for c in bowtiw_commands:
        fout.write(c + "\n")
    fout.close()
    
    print "\n. Launching Bowtie2 with the following command:"
    for c in bowtie_commands:
        os.system(c)

"""
    Step 2: Hybrid-related stuff
"""

if jump <= 2:
    sql = "select * from Hybrids"
    cur.execute(sql)
    x = cur.fetchall()
    annoids = []
    for ii in x:
        annoids.append( ii[0] )
    
    sql = "drop table Reads"
    cur.execute(sql)
    build_anno_db(con)
    for annoid in annoids:
        chrom_filter = ap.getOptionalArg("--chrom_filter")
        if chrom_filter == False:
            chrom_filter = None
        extract_perfect_reads(annoid, con, chrom_filter=chrom_filter)

"""Build a map of which annotations map to which species."""
if jump <= 2.1:
    get_hybrid_pairs(con)

"""Remove hybrid reads that map to both parental genomes."""    
if jump <= 2.3:
    find_hybrid_unique_reads(con)

if jump <= 2.4:
    """Write SAM files for hybrid reads, containing only those reads that
    are 100% match and unique to each parent genome."""
    write_filtered_sam(con)

if jump <= 3:
    """Convert SAM files to sorted BAM files."""
    write_sorted_bam(con)

if jump <= 3.1:
    """Verify the SAM -> sorted BAM occurred correctly."""
    check_bams(con)

if jump <= 4:
    """Launch MACS2"""
    run_peak_calling(con)

#
# check macs output
if jump <= 4.1:
    check_peaks(con)

#
# make BDG for F-E
if jump <= 5:
    calculate_fe(con)

if jump <= 5.1:
    check_fe(con)

if jump <= 6:
    bed2wig(con)

if jump <= 6.1:
    check_wig(con)

#
# write config
if jump <= 7:
    write_viz_config(con)

#
# launch APRES

if jump <= 100:
    print_read_stats(con)

exit()

