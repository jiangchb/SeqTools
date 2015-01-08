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
    print "======================================================================"
    print "."
    print ". The ChIP-Seq Distillery"
    print "."
    print ". A pipeline to process your FASTQ reads."
    print "."
    print ". Written by Victor Hanson-Smith -- victorhansonsmith@gmail.com"
    print "."
    print ". Last Updated: January 2014"
    print "."
    print ". This program does the following steps:"
    print "\t1. Maps reads to reference genome (bowtie)"
    print "\t2. Separates reads from hybrid genomes (optional)"
    print "\t3. Sorts all reads and creates BAM files (samtools)"
    print "\t4. Finds peaks (macs2)"
    print "\t5. Calculates fold-enrichment (macs2)"
    print "\t6. Compares peaks and fold-enrichments across "
    print "\t\treplicates conditions, and species."
    print ""
    print "======================================================================" 

splash()

"""read_cli reads the command-line, the annotation file, and creates a SQLite3 database."""
con = read_cli(ap)

print_settings(con)

cur = con.cursor()
sql = "select count(*) from Annotations"
cur.execute(sql)
count = cur.fetchone()[0]
print "\n\t--> I found", count, "annotations in the database."

# sql = "select * from Annotations"
# cur.execute(sql)
# x = cur.fetchall()
# for ii in x:
#     print ii
# exit()

"""Jump allows for some steps to be skipped."""
jump = ap.getOptionalArg("--jump")
stop = ap.getOptionalArg("--stop")

if jump == False:
    jump = 0
else:
    jump = float(jump)
    
if stop == False:
    stop = 10000000
else:
    stop = float(stop)

# macs_pairs = get_macs_pairs(con)
# for p in macs_pairs:
#     print p
# exit()

"""Run Bowtie2 for each FASTQ path"""
if jump <= 1 and stop > 1:
    run_bowtie(con)

if jump <= 1.1 and stop > 1.1:
    check_bowtie_output(con)

"""Hybrid only: find non-mismatch reads."""
if jump <= 2 and stop > 2:
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
if jump <= 2.1 and stop > 2.1:
    get_hybrid_pairs(con)

"""Remove hybrid reads that map to both parental genomes."""    
if jump <= 2.3 and stop > 2.3:
    find_hybrid_unique_reads(con)

if jump <= 2.4 and stop > 2.4:
    """Write SAM files for hybrid reads, containing only those reads that
    are 100% match and unique to each parent genome."""
    write_filtered_sam(con)

if jump <= 3 and stop > 3:
    """Convert SAM files to sorted BAM files."""
    write_sorted_bam(con)

if jump <= 3.1 and stop > 3.1:
    """Verify the SAM -> sorted BAM occurred correctly."""
    check_bams(con)

if jump <= 3.2 and stop > 3.2:
    """"Make WIGs from sorted BAMs"""
    bam2bedgraph(con)

if jump <= 3.3 and stop > 3.3:
    """Check the bedgraphs containing raw read densities (not peaks)"""
    check_bedgraphs(con)

if jump <= 4 and stop > 4:
    """Launch MACS2"""
    run_peak_calling(con)

if jump <= 4.1 and stop > 4.1:
    check_peaks(con)

if jump <= 5 and stop > 5:
    """Make BDG for FE"""
    calculate_fe(con)

if jump <= 5.1 and stop > 5.1:
    check_fe(con)

if jump <= 6 and stop > 6:
    """Make WIGs from bedgraph files containing raw reads (step 3.2)
    and bedgraph files containing fold-enrichment values (step 5)"""
    bed2wig(con)

if jump <= 6.2 and stop > 6.2:
    """Checks all the wigs written in the previous step."""
    check_wig(con)

# write config
if jump <= 7 and stop > 7:
    """Write a *.config file for the apres.py script"""
    write_viz_config(con)
#
# launch APRES
if jump <= 8 and stop > 8:
    """Launch the apres.py script."""
    launch_viz(con)

if jump <= 100 and stop > 100:
    print_read_stats(con)

print "\n. ChIP-Seq distillation is complete.  Goodbye."
exit()

