#
# This script processes ChIP-Seq results from FASTQ to (mapped) SAM to (sorted) BAM
# to fold-enrichment track and summit list.
#

import re, os, sys, time
from annotation_db import *
from config import *
from hybrid_tools import *
from tools import *
from argParser import ArgParser
ap = ArgParser(sys.argv)

def splash():
    print "============================================="
    print " run_from_anno.py"
    print " October 2014"
    print " A script to analyze ChIP-Seq data."
    print "=============================================" 

splash()

"""read_cli reads the command-line, the annotation file, and creates a SQLite3 database."""
con = read_cli(ap)

cur = con.cursor()
sql = "select count(*) from Annotations"
cur.execute(sql)
count = cur.fetchone()[0]
print "\n. I found", count, "annotations in the database."

"""Jump allows for some steps to be skipped."""
jump = ap.getOptionalArg("--jump")
if jump == False:
    jump = 0
else:
    jump = float(jump)

"""Run Bowtie2 for each FASTQ path"""
if jump <= 1:
    run_bowtie(con)

if jump <= 1.1:
    check_bowtie_output(con)

"""Hybrid only: find non-mismatch reads."""
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
        extract_perfect_reads(annoid, con)

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

if jump <= 4.1:
    check_peaks(con)

if jump <= 5:
    """Make BDG for FE"""
    calculate_fe(con)

if jump <= 5.1:
    check_fe(con)

if jump <= 6:
    """Make WIG for FE BDG"""
    bed2wig(con)

if jump <= 6.1:
    check_wig(con)

# write config
if jump <= 7:
    write_viz_config(con)
#
# launch APRES

if jump <= 100:
    print_read_stats(con)

exit()

