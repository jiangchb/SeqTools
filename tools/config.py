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
    SAMTOOLS = "samtools"
    MPIRUN = "mpirun -np 11 --machinefile hosts.txt /common/bin/mpi_dispatch"
    MACS2 = "/common/REPOSITORY/MACS2-2.0.10.07132012/bin/macs2"
    SEQTOOLSDIR = "~/Applications/SeqTools"
        
    """restricts analysis to only those annotation lines whose sample column is this value."""
    restrict_to_sample = ap.getOptionalArg("--restrict_to_sample")
    
    PRACTICE_MODE = ap.getOptionalToggle("--practice_mode")
        
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
    
    chrom_filters = ap.getOptionalList("--chrom_filter")
    for cf in chrom_filters:
        sql = "insert or replace into Settings (keyword, value) VALUES('chrom_filter','" + cf + "')"
        cur.execute(sql)
    x = 0
    if PRACTICE_MODE:
        print "\n. Practice Mode Enabled."
        x = 1
    sql = "insert or replace into Settings (keyword, value) VALUES('practice_mode','" + x.__str__() + "')"
    cur.execute(sql)
    con.commit()
    
    con = import_annotations(APATH, con)
    return con