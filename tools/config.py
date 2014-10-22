from annotation_db import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

def read_cli(ap):
    """Returns the database."""
    
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
    x = ap.getOptionalArg("--bowtie2")
    if x != False:
        BOWTIE2 = x
    
    SAMTOOLS = "samtools"
    x = ap.getOptionalArg("--samtools")
    if x != False:
        SAMTOOLS = x
    
    MPIRUN = "mpirun -np 11 --machinefile hosts.txt /common/bin/mpi_dispatch"
    
    MACS2 = "/common/REPOSITORY/MACS2-2.0.10.07132012/bin/macs2"
    x = ap.getOptionalArg("--macs2")
    if x != False:
        MACS2 = x
    
    SEQTOOLSDIR = "~/Applications/SeqTools"
    x = ap.getOptionalArg("--seqtoolsdir")
    if x != False:
        SEQTOOLSDIR = x
    
    """The pillars file provides alias names for orthologous genes in different species.
    For example, a pillars file for several model yeast species can be downloaded from the Candida Genome Database.
    This information is used later in the script apres.py"""
    PILLARSPATH = "/Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/Pillars.tab.txt"
    x = ap.getOptionalArg("--pillarspath")
    if x != False:
        PILLARSPATH = x
        
    """restricts analysis to only those annotation lines whose sample column is this value."""
    restrict_to_sample = ap.getOptionalArg("--restrict_to_sample")
    
    PRACTICE_MODE = ap.getOptionalToggle("--practice_mode")
        
    """Build or restore the DB"""
    con = get_db(DBPATH)
    cur = con.cursor()
    
    sql = "delete from Settings"
    cur.execute(sql)
    con.commit()
    
    """Write all the settings into the DB.
    Later in the program, lookups of settings will use the values in the DB, 
    rather than the values in the ArgParser object."""
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
    sql = "insert or replace into Settings (keyword, value) VALUES('pillars_path','" + PILLARSPATH + "')"
    cur.execute(sql)
    
    """chrom_filter will restrict the analysis to only those reads whose chromosome name contains the value of --chrom_filter"""
    chrom_filters = ap.getOptionalList("--chrom_filter")
    for cf in chrom_filters:
        sql = "insert or replace into Settings (keyword, value) VALUES('chrom_filter','" + cf + "')"
        cur.execute(sql)

    """If practice mode is enabled, then major programs will not be launched, such as Bowtie, MACS2, etc.
    However, regardless of practice mode, shell scripts for these programs will be written.
    Essentially, practice mode provides a means to test if the pipeline will succeed."""
    x = 0
    if PRACTICE_MODE:
        print "\n. Practice Mode Enabled."
        x = 1
    sql = "insert or replace into Settings (keyword, value) VALUES('practice_mode','" + x.__str__() + "')"
    cur.execute(sql)
    con.commit()
    
    con = import_annotations(APATH, con)
    return con