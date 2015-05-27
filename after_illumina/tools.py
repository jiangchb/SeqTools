import re, os, sys, time
from annotation_db import *

def run_bowtie(con):
    """Clear previous Bowtie output paths."""
    cur = con.cursor()
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
        full_fastq_path = get_setting("datadir",con) + "/" + fastq
        species = ii[2]
        
        """Sanity check:"""
        if False == os.path.exists( full_fastq_path ):
            print "\n. Error, I can't find your FASTQ file at", full_fastq_path
            print ". (tools.py 43)"
            exit()
        c += " -U " + full_fastq_path
        
        """Is it a hybrid?"""
        sql = "select count(*) from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        samoutpath = re.sub(".fastq", ".sam", fastq)
        if count > 0:
            """Hybrids get a special SAM path"""
            samoutpath = re.sub(".fastq", "-" + species + ".sam", fastq)
            samoutpath = get_setting("outdir",con) + "/" + samoutpath
        c += " -S " + samoutpath
        c += " --no-unal "
        
        """Path to directory with genome sequences"""
        sql = "select genomepath from SpeciesGenomepath where speciesname='" + species + "'"
        cur.execute(sql)
        x = cur.fetchall()
        if x.__len__() == 0:
            print "\n. Error: I cannot find a genome file reference for the species " + species
            print ". Check your genome list, specified by the --genome_list argument."
            exit()
        else:
            gpath = x[0][0]
        c += " -x " + gpath
        
        bowtie_commands.append(c)
        
        sql = "insert or replace into BowtieOutput (annoid, sampath) VALUES(" + annoid.__str__() + ",'" + samoutpath + "')"
        cur.execute(sql)
        con.commit()
        
    print "\n. Launching the following Bowtie2 commands:"
    fout = open("bowtie_commands.sh", "w")
    for c in bowtie_commands:
        print c
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode", con) == "0":
        """MPI for Bowtie is temporarily disabled because Bowtie is using the -p option."""
        #if get_setting("use_mpi", con) == "1":
         #   os.system( get_setting("mpirun",con) + " bowtie_commands.sh")
        #else:
        time_start_bowtie = time.time()
        os.system("source bowtie_commands.sh")
        bowtie_runtime = time.time() - time_start_bowtie
        write_log(con, "The script bowtie_commands.sh completed in " + bowtie_runtime.__str__() + " seconds.")

def check_bowtie_output(con):
    cur = con.cursor()
    sql = "select annoid, sampath from BowtieOutput"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        sampath = ii[1]
        if False == os.path.exists( sampath ):
            print "\n. Error, I can't find your BowtieOutput file at", sampath
            exit()
    print "\n. Bowtie Output is OK."
    return


def extract_matched_reads(annoid, con, chrom_filter = None):
    """This method parses the SAM output from Bowtie and extracts those reads that pass a gauntlet of test:
        1. Is the read's mismatch level <= the user-specified mismatch level?
        2. Is the read aligned to multiple locations in the genome? (Default: if yes, then discard the read)
        3. Did the user specify the parameter --chrom_filter to process reads only from a subset of chromosomes?
        Output: Writes a new (shorted) SAM file, and also returns a hash of read IDs.
    
        If the user specified the parameter --eliminate_multialign, then
        this method will also cull reads mapped to multiple locations within the same genome.
    
    """
    cur = con.cursor()
    
    sql = "create table if not exists Reads" + annoid.__str__() + " (readid INTEGER primary key autoincrement, readname TEXT, mismatch INT, order_seen INT)" # reads without mismatches
    #sql = "delete from Reads where annoid=" + annoid.__str__()
    cur.execute(sql)
    con.commit()
        
    sql = "delete from Reads" + annoid.__str__()
    cur.execute(sql)
    con.commit()
    
    
    sampath = None
    sql = "select sampath from BowtieOutput where annoid=" + annoid.__str__()
    cur.execute(sql)
    x = cur.fetchone()
    if x.__len__() > 0:
        sampath = x[0]
    else:
        print "\n. An error occurred; I can't find a sampath for annotation", annoid
        exit()
    outsampath = re.sub(".sam", ".perfect.sam", sampath)
    
    """Get the mismatch threshold:"""
    MTHRESH = int( get_setting("mismatch_thresh", con) )
    if MTHRESH < 0:
        MTHRESH = 100000000 # effectively, gathers all reads regardless of their mismatch level.
        print "\n. Extracting all reads from", sampath
    else:
        print "\n. Extracting reads with mismatch <= " + MTHRESH.__str__() + ", from", sampath
    
    """Open the SAM output of Bowtie""" 
    fin = open(sampath, "r")
    total = 0 # the total count of all reads seen
    inserted = 0
    for line in fin.xreadlines():
        """Parse each line of the SAM"""
        try:
            """The var outline will be a printed statement when we're done processing this line."""
            outline = ""
            if line.startswith("@"):
                continue
            if line.__len__() < 5:
                continue
            
            """Display a progress indicator."""
            total += 1
            if total%10000 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                con.commit()
            
            """Parse this line."""
            tokens = line.split()
            readname = tokens[0]
            mismatchlevel = None
            multilocscore = None
            
            for t in tokens[10:]:
                """The tokens from index 10+ contains useful information about this read."""
                if False == t.startswith("X"):
                    """Informative tokens begin with an X, e.g. XM:..., XS:...."""
                    continue
                elif t.startswith("XM:i:"): #0":
    
                    """Skip reads with a mismatch score higher than the
                        mismatch threshold."""
                    tparts = t.split(":")
                    mismatchlevel = int(tparts[2])
                elif t.startswith("XS:i:"):
                    tparts = t.split(":")
                    multilocscore = int(tparts[2])
            
            if mismatchlevel > MTHRESH:
                """Too many mismatches."""
                continue
            if get_setting("eliminate_multialign", con) != None and multilocscore != None:
                """Reads with multiple aligned locations should be skipped."""
                continue
                 
            """The chrom_filter allows for reads to be considered
            only if they match the filter, such "Chr1".
            This is useful for generating toy-sized test databases
            on a restricted number of reads."""
            cflist = get_setting_list("chrom_filter",con)
            if cflist.__len__() > 0:
                chrom = tokens[2]
                if chrom not in cflist:
                    continue
            
            """The read passed all the validation checks. Let's import it."""
            inserted += 1
            sql = "insert into Reads" + annoid.__str__() + "(readname,mismatch,order_seen) VALUES("
            sql += "'" + readname + "'," + mismatchlevel.__str__() + "," + inserted.__str__() + ")"
            cur.execute(sql)
        except:
            con.rollback() # Roll back the stuff we didn't commit.
        
    con.commit()
    fin.close()
    
    ratio = 100.0*float(inserted)/total
    ratio = "%.1f"%ratio
    print "\n\t--> I found " + inserted.__str__() + " reads out of " + total.__str__() + " reads (" + ratio.__str__() + "%) that satisfied all validation checks."  

    """How many reads were perfect?"""
    sql = "select count(*) from Reads" + annoid.__str__() + " where mismatch=0"
    cur.execute(sql)
    count_perfect = cur.fetchone()[0]
    
    sql = "insert or replace into ReadStats(annoid, nperfect, ntotal) VALUES("
    sql += annoid.__str__() + "," + count_perfect.__str__() + "," + total.__str__() + ")"
    cur.execute(sql)
    con.commit()

def write_sorted_bam(con):
    print "\n. Writing sorted BAM files."
    
    """Writes sorted BAM files for all SAM files.
    if delete_sam == True, then the original SAM file will be deleted."""
    pairs = [] # pairs is a list of (annoid, sampath) pairs.
    
    cur = con.cursor()
    sql = "delete from SortedBamFiles"
    cur.execute(sql)
    con.commit()
    
    """First, get the special filtered SAM files for hybrids."""
    sql = "select annoid, sampath from FilteredBowtieOutput"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        sampath = ii[1]
        pairs.append( (annoid,sampath) )
    
    """Next, get the normal SAM files for all other annotations."""
    sql = "select annoid, sampath from BowtieOutput where annoid not in (select annoid from FilteredBowtieOutput)"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        sampath = ii[1]
        pairs.append( (annoid,sampath) )
    
    samtools_commands = []
    for pair in pairs:
        annoid = pair[0]
        sampath = pair[1]
        sorted_bampath = re.sub(".sam", ".sort", sampath)
        c = get_setting("samtools",con) + " view -bS " + sampath
        c += " | " + get_setting("samtools",con) + " sort - " + sorted_bampath
        samtools_commands.append(c)
        
        """Update the DB"""
        sql = "insert or replace into SortedBamFiles(annoid, bampath) VALUES("
        sql += annoid.__str__() + ",'" + sorted_bampath + ".bam" + "')"
        cur.execute(sql)
        con.commit()

    fout = open("run_sam2bam.sh", "w")
    for c in samtools_commands:
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode",con) == "0":
        if get_setting("use_mpi", con) == "1":
            os.system(get_setting("mpirun", con) + " run_sam2bam.sh")
        else:
            os.system("source run_sam2bam.sh")
    
def check_bams(con, delete_sam = True):
    """Checks for the presence of sorted BAM files. If they all exist,
        then the SAM files (from which the BAMs were created) will be deleted.
        If the files aren't found, then the program will terminate."""
    cur = con.cursor()        
    sql = "select annoid from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        sql = "select bampath from SortedBamFiles where annoid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y.__len__() < 1:
            emsg = "Error, the table SortedBamFiles has no records for annotation #" + ii[0].__str__()
            print "\n.", emsg
            write_error(con, emsg)
            exit()
        bampath = y[0][0]
        if False == os.path.exists( bampath ):
            emsg = "Error, I can't find your BAM file at " + bampath
            print "\n.", emsg
            write_error(con, emsg)
            exit()
    
    """At this point, all the BAMs look good."""
    if delete_sam == True and get_setting("practice_mode",con) == "0":
        """Delete the original SAM files."""
        sql = "select sampath from BowtieOutput"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            sampath = ii[0]
            os.system("rm " + sampath)
#         sql = "select sampath from FilteredBowtieOutput"
#         cur.execute(sql)
#         x = cur.fetchall()
#         for ii in x:
#             sampath = ii[0]
#             os.system("rm " + sampath)
    write_log(con, "The sorted BAMs appear OK")
    print "\n. The sorted BAMs appear OK."
    return


def run_peak_calling(con):
    """Runs MACS2 to call peaks."""
    cur = con.cursor()
    
    sql = "delete from MacsRun"
    cur.execute(sql)
    con.commit()
    
    """Verify that all the BAM files exist."""
    sql = "select annoid, bampath from SortedBamFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        bampath = ii[1]
        if False == os.path.exists(bampath):
            emsg = "Error, I can't find your BAM at " + bampath.__str__()
            print "\n.", emsg
            write_error(con, emsg)
            exit()
    
    macs_commands = []
    
    macs_pairs = get_macs_pairs(con)
    
    print "\n. I found the following MACS pairs:"
    print macs_pairs
    
    for pair in macs_pairs:       
        exp_annoid = pair[0]
        control_annoid = pair[1]
        
        sql = "select bampath from SortedBamFiles where annoid=" + exp_annoid.__str__()
        cur.execute(sql)
        exp_bampath = cur.fetchone()[0]
    
        sql = "select bampath from SortedBamFiles where annoid=" + control_annoid.__str__()
        cur.execute(sql)
        control_bampath = cur.fetchone()[0]
        
        macs_cmd = get_setting("macs2", con) + " callpeak "
        macs_cmd += " -t " + exp_bampath
        macs_cmd += " -c " + control_bampath
        macs_cmd += " --gsize 1.43e+07 --nomodel --shiftsize 78 -B --SPMR "
        macs_cmd += " --name " + get_name_for_macs(exp_annoid, control_annoid, con)
        macs_commands.append( macs_cmd )
    
        sql = "insert or replace into MacsRun(exp_annoid, control_annoid, name) VALUES("
        sql += exp_annoid.__str__() + "," + control_annoid.__str__()
        sql += ",'" + get_name_for_macs(exp_annoid, control_annoid, con) + "')"
        cur.execute(sql)
        con.commit()
    
    fout = open("macs_commands.sh", "w")
    for c in macs_commands:
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode", con) == "0":
        if get_setting("use_mpi", con) == "1":
            os.system( get_setting("mpirun",con) + " macs_commands.sh" )
        else:
            os.system("source macs_commands.sh")


def check_peaks(con):
    cur = con.cursor()
    
    sql = "delete from MacsPeakPaths"
    cur.execute(sql)
    con.commit()
    
    sql = "select exp_annoid, control_annoid, name from MacsRun"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        exp_annoid = ii[0]
        tbdg = ii[2] + "_treat_pileup.bdg"
        cbdg = ii[2] + "_control_lambda.bdg"
        peaks = ii[2] + "_peaks.bed"
        summits = ii[2] + "_summits.bed"
    
        outpaths = [tbdg, cbdg, peaks, summits]
        for f in outpaths:
            if False == os.path.exists(f):
                print "\n. Error, I can't find the MACS2 output file", tbdg
                exit()
                
        """Update the SQL db."""
        sql = "insert or replace into MacsPeakPaths(exp_annoid, treatment_pileup_path, control_lambda_path, peaks_path, summits_path) VALUES("
        sql += exp_annoid.__str__() + ",'" + tbdg + "','" + cbdg + "','" + peaks + "','" + summits + "')"
        cur.execute(sql)
        con.commit()
    return


def calculate_fe(con):
    cur = con.cursor()
    sql = "delete from MacsFE"
    cur.execute(sql)
    con.commit()
    
    sql = "select exp_annoid, control_annoid, name from MacsRun"
    cur.execute(sql)
    x = cur.fetchall()
    fe_commands = []
    for ii in x:
        exp_annoid = ii[0]
        tbdg = ii[2] + "_treat_pileup.bdg"
        cbdg = ii[2] + "_control_lambda.bdg"
        outbdg = ii[2] + "_output_FE.bdg"
        c = get_setting("macs2", con) + " bdgcmp "
        c += " -t " + tbdg
        c += " -c " + cbdg
        c += " -o " + outbdg
        c += " -m FE"
        fe_commands.append( c )
        
        sql = "insert or replace into MacsFE(exp_annoid, bdgpath) VALUES("
        sql += exp_annoid.__str__()
        sql += ",'" + outbdg + "')"
        cur.execute(sql)
        con.commit()
        
    fout = open("fe_commands.sh", "w")
    for c in fe_commands:
        print c
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode", con) == "0":
        if get_setting("use_mpi", con) == "1":
            os.system( get_setting("mpirun",con) + " fe_commands.sh" )
        else:
            os.system("source fe_commands.sh")


def check_fe(con):
    cur = con.cursor()
    sql = "select exp_annoid, bdgpath from MacsFE"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        bdgpath = ii[1]
        if False == os.path.exists(bdgpath):
            print "\n. Error, I can't find the BDG file at", bdgpath
            exit()
    print "\n. The fold-enrichment bedgraph files are OK."
    return

def bam2bedgraph(con):
    """Converts sorted BAM files to Bedgraph files."""
    print "\n. Converting sorted BAMs to bedgraph files."
    
    cur = con.cursor()        
    
    """Clear any prior entries in the table Bedgraphfiles"""
    sql = "delete from BedgraphFiles"
    cur.execute(sql)
    con.commit()
    
    sql = "select annoid from Annotations"
    cur.execute(sql)
    x = cur.fetchall()

    annoid_bampath = {}
    for ii in x:
        sql = "select annoid, bampath from SortedBamFiles where annoid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y.__len__() < 1:
            print "\n. Error, the table SortedBamFiles has no records for annotation #", ii[0].__str__()
            exit()
        annoid = y[0][0]
        bampath = y[0][1]
        if False == os.path.exists( bampath ):
            print "\n. Error, I can't find your BAM file at", bampath
            exit()
        annoid_bampath[annoid] =  bampath
    
    commands = []
    for annoid in annoid_bampath:
        bedpath = re.sub( ".sort.bam", ".bed", annoid_bampath[annoid] )
        bedpath = re.sub( ".unique", "", bedpath)
        c = get_setting("gcb", con) + " -ibam " + annoid_bampath[annoid]
        c += " -bga "
        c += " > " + bedpath
        commands.append(c)
    
        sql = "insert or replace into BedgraphFiles (annoid, bedpath) VALUES("
        sql += annoid.__str__() + ",'" + bedpath + "')"
        cur.execute(sql)
        con.commit()
    
    """Write a shell script with all the genomCoveredBed commands,
    and then launch the script with mpi_dispatch."""
    scriptpath = get_setting("outdir", con) + "/run_gcb.sh"
    fout = open(scriptpath, "w")
    for c in commands:
        print c
        fout.write(c + "\n")
    fout.close()
    if get_setting("practice_mode", con) == "0":
        if get_setting("use_mpi", con) == "1":
            os.system(get_setting("mpirun",con) + " " + scriptpath)
        else:
            os.system("source " + scriptpath)

def check_bedgraphs(con):
    cur = con.cursor()
    sql = "select annoid, bedpath from BedgraphFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        bedpath = ii[1]
        if False == os.path.exists(bedpath):
            emsg = " Error, I can't find your bedgraph file at " + bedpath.__str__()
            print "\n.", emsg
            write_error(con, emsg)
            exit()
        fin = open(bedpath, "r")
        firstline = fin.readline()
        tokens = firstline.split()
        fin.close()
        if tokens.__len__() != 4:
            emsg = "Error, the first line of your bedgraph doesn't contain 4 columns. It appears to not be a bedgraph file."
            print "\n.", emsg
            write_error(con, emsg)
            exit()
    return


def bed2wig_helper(bedgraphpath, wigpath):
    printspan = 100000 # print an update every N sites
    count = 0
    last_seen_chrom = None
    fin = open(bedgraphpath, "r")
    msg = "Converting BDG->WIG:" + bedgraphpath.__str__() + " --> " + wigpath.__str__()
    print "\n.", msg
    write_log(con, msg)
    fout = open(wigpath, "w")
    fout.write("track type=WIG\n")
    for l in fin.xreadlines():
        if l.__len__() > 2:
            tokens = l.split()
            chromname = tokens[0]
    
            if last_seen_chrom != None and chromname != last_seen_chrom:
                last_seen_chrom = None
    
            if last_seen_chrom == None:            
                fout.write("variableStep chrom=" + chromname + "\n")
                last_seen_chrom = chromname
            start = int( tokens[1] )
            stop = int( tokens[2] )
            value = tokens[3]
            for ii in range(start, stop):
                if count%printspan == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                count += 1
                fout.write(ii.__str__() + "\t" + value + "\n")
    fin.close()
    fout.close()


def bed2wig(con):
    """Converts Bedgraph files to WIG files, for viewing in MochiView.
    Bedgraph filepaths are drawn from the tables ReadsWigFiles and FEWigFiles."""
    
    cur = con.cursor()
    sql = "delete from ReadsWigFiles"
    cur.execute(sql)
    sql = "delete from FEWigFiles"
    cur.execute(sql)
    con.commit()

    #commands = []

    """First convert the reads wig files."""
    sql = "select annoid, bedpath from BedgraphFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        bedpath = ii[1]
        wigpath = re.sub(".bed", ".wig", bedpath)
        
        """THis is a sanity check to prevent the BED being overwritten."""
        if wigpath == bedpath: # wtf.
            wigpath = bedpath + ".wig"
        
        #
        # continue here
        #
        #c = "python " + get_setting("seqtoolsdir", con) + "/bedgraph2wig.py " + bedpath + "  " + wigpath
        #commands.append(c)
        bed2wig_helper(bedpath, wigpath)
        
        sql = "insert or replace into ReadsWigFiles(annoid, wigpath) VALUES("
        sql += annoid.__str__()
        sql += ",'" + wigpath + "')"
        cur.execute(sql)
        con.commit()
    
    """Then convert the FE wig files."""
    sql = "select exp_annoid, bdgpath from MacsFE"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        bdgpath = ii[1]
        wigpath = re.sub(".bdg", ".wig", bdgpath)
        
        """This is a sanity check to prevent the BED being overwritten."""
        if wigpath == bdgpath: # wtf.
            wigpath = bdgpath + ".wig"
        
        #c = "python " + get_setting("seqtoolsdir", con) + "/bedgraph2wig.py " + bdgpath + "  " + wigpath
        #commands.append(c)
        bed2wig_helper(bdgpath, wigpath)
        
        sql = "insert or replace into FEWigFiles(exp_annoid, org_bdgpath, wigpath) VALUES("
        sql += annoid.__str__()
        sql += ",'" + bdgpath
        sql += "','" + wigpath + "')"
        cur.execute(sql)
        con.commit()

    #scriptpath = "run_bed2wig.sh"
    #fout = open(scriptpath, "w")
    #for c in commands:
    #    print c
    #    fout.write(c + "\n")
    #fout.close()
    
    #if get_setting("practice_mode", con) == "0":
    #    if get_setting("use_mpi", con) == "1":
    #        os.system( get_setting("mpirun",con) + " " + scriptpath ) 
    #    else:
    #        os.system("source " + scriptpath)  
    
def check_wig(con):
    cur = con.cursor()

    sql = "select annoid, wigpath from ReadsWigFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        wigpath = ii[1]
        if False == os.path.exists( wigpath ):
            print "\n. Error, I can't find your WIG file ", wigpath
            exit()
        print wigpath

    sql = "select exp_annoid, wigpath from FEWigFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        wigpath = ii[1]
        if False == os.path.exists( wigpath ):
            print "\n. Error, I can't find your WIG file ", wigpath
            exit()
        fin = open(wigpath, "r")
        header = fin.readline()
        fin.close()
        if False == header.__contains__("track type=WIG"):
            print "\n. Error, I cannot find a valid header in your WIG file", wigpath
            exit()
    
    print "\n. WIG files OK."
    return



def write_viz_config(con):
    """This method writes the configuration file for the 'after_peaks' pipeline, which is primarily
        concerned with comparing replicates, comparing groups of replicates, and creating visulizations."""
    
    cur = con.cursor()
    
    configpath = get_setting("outdir", con) + "/" + get_setting("project_name", con) + ".config"
    sql = "insert or replace into Settings (keyword, value) VALUES('viz_configpath','" + configpath + "')"
    cur.execute(sql)
    con.commit()
        
    species = []
    sql = "select distinct species from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        species.append( ii[0] )
        
    fout = open(configpath, "w")
    
    for s in species:
        fout.write("SPECIES " + s + "\n")
        fout.write("NAME = " + s + "\n")
        
        """Here are some genome feature files that are known."""
        if s == "Cdub":
            fout.write("GFF = /Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/gff/C_dubliniensis_CD36_version_s01-m02-r08_features.gff\n")
        if s == "Calb":
            fout.write("GFF = /Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/gff/C_albicans_SC5314_A21_current_features.gff\n")
        if s == "Ctro":
            fout.write("GFF = /Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/gff/C_tropicalis_MYA-3404_features.gff\n")
        if s == "Scer":
            fout.write("GFF = /Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/gff/saccharomyces_cerevisiae.gff\n")
        if s == "Pstip":
            fout.write("GFF = /Network/Servers/udp015817uds.ucsf.edu/Users/Shared/sequencing_analysis/gff/pstip.gff\n")
        
        repgroups = []
        sql = "select distinct strain from Annotations where species='" + s + "'"
        sql += " and Annotations.annoid in (select exp_annoid from MacsRun)"
        
        #sql = "select distinct strain from Annotations where species='" + s + "'"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            strain = ii[0]
            repgroups.append( strain )
        for strain in repgroups:
            fout.write("\tREPGROUP " + strain + "-" + s + "\n")
            
            replicates = []
            sql = "select distinct replicate from Annotations where species='" + s + "' and strain='" + strain + "'"
            cur.execute(sql)
            x = cur.fetchall()
            for ii in x:
                replicates.append( ii[0] )
            
            for ii in range(0, replicates.__len__() ):
                fout.write("\t\tREPLICATE " + (ii+1).__str__() + "\n")
                annoids = []
                #sql += "where "
                sql = "select exp_annoid from MacsRun where exp_annoid in (select annoid from Annotations where species='" + s + "' and strain='" + strain + "' and replicate=" + replicates[ii].__str__() + ")"
                cur.execute(sql)
                x = cur.fetchall()
                for ii in x:
                    annoids.append( ii[0] )
                
                for id in annoids:
                    sql = "select summits_path from MacsPeakPaths where exp_annoid=" + id.__str__()
                    cur.execute(sql)
                    x = cur.fetchall()
                    if x.__len__() < 1:
                        print "\n. An error occurred. Checkpoint 413"
                        exit()
                    summitpath = x[0][0]
                    fout.write("\t\tSUMMITS = " + summitpath + "\n")

                    sql = "select bdgpath from MacsFE where exp_annoid=" + id.__str__()
                    cur.execute(sql)
                    x = cur.fetchall()
                    if x.__len__() < 1:
                        print "\n. An error occurred. Checkpoint 413"
                        exit()
                    bdgpath = x[0][0]
                    fout.write("\t\tENRICHMENTS = " + bdgpath + "\n")
                
                pass
    fout.close()

def launch_viz(con):
    cur = con.cursor()
    pname = get_setting("project_name", con)
    vcpath = get_setting("viz_configpath", con)
    
    vizdbpath = get_setting("outdir", con) + "/" + pname + ".viz.db"
    if os.path.exists(vizdbpath):
        os.system("rm -rf " + vizdbpath)
    
    scriptpath = get_setting("outdir", con) + "/run_viz.sh"
    fout = open(scriptpath, "w")
    
    c = "python /common/REPOSITORY/SeqTools/after_peaks/apres.py "
    c += "--dbpath " + vizdbpath
    c += " --pillarspath " + get_setting("pillars_path", con)
    c += " --configpath " + vcpath
    fout.write(c + "\n")
    fout.close()
    if get_setting("practice_mode", con) == "0":
        os.system("source " + scriptpath)
    
    