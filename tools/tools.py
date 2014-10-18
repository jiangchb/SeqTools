import re, os, sys
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
        species = ii[2]
        
        """Sanity check:"""
        if False == os.path.exists( get_setting("datadir",con) + "/" + fastq ):
            print "\n. Error, I can't find your FASTQ file at", fastq
            print ". (run_from_anno.py 43)"
            exit()
        c += " -U " + get_setting("datadir",con) + "/" + fastq
        
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
        
        sql = "insert or replace into BowtieOutput (annoid, sampath) VALUES(" + annoid.__str__() + ",'" + samoutpath + "')"
        cur.execute(sql)
        con.commit()
        
    fout = open("bowtie_commands.sh", "w")
    for c in bowtie_commands:
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode", con) == "0":
        print "\n. Launching Bowtie2 with the following command:"
        for c in bowtie_commands:
            os.system(c)

def write_sorted_bam(con, delete_sam = True):
    """Writes sorted BAM files for all SAM files.
    if delete_sam == True, then the original SAM file will be deleted."""
    pairs = [] # list of (annoid, sampath) pairs.
    
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

    if get_setting("practice_mode",con) == "0":
        fout = open("run_sam2bam.sh", "w")
        for c in samtools_commands:
            fout.write(c + "\n")
        fout.close()

        mpirun = get_setting("mpirun", con)
        os.system(mpirun + " run_sam2bam.sh")
    
    if delete_sam == True and get_setting("practice_mode",con) == "0":
        """Delete the original SAM files."""
        sql = "select sampath from BowtieOutput"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            sampath = ii[0]
            os.system("rm " + sampath)
        sql = "select sampath from FilteredBowtieOutput"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            sampath = ii[0]
            os.system("rm " + sampath)

def check_bams(con):
    """Will exit on error."""
    cur = con.cursor()        
    sql = "select annoid from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        sql = "select bampath from SortedBamFiles where annoid=" + ii[0].__str__()
        cur.execute(sql)
        y = cur.fetchall()
        if y.__len__() < 1:
            print "\n. Error, the table SortedBamFiles has no records for annotation #", ii[0].__str__()
            exit()
        bampath = y[0][0]
        if False == os.path.exists( bampath ):
            print "\n. Error, I can't find your BAM file at", bampath
            exit()
    return

def run_peak_calling(con):
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
            print "\n. Error, I can't find your BAM at", bampath
            exit()
    
    macs_commands = []
    
    macs_pairs = get_macs_pairs(con)
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
        os.system( get_setting("mpirun",con) + " macs_commands.sh" )

def check_peaks(con):
    cur = con.cursor()
    sql = "select exp_annoid, control_annoid, name from MacsRun"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        tbdg = ii[2] + "_treat_pileup.bdg"
        cbdg = ii[2] + "_control_lambda.bdg"
        peaks = ii[2] + "_peaks.bed"
        summits = ii[2] + "_summits.bed"
        
        outpaths = [tbdg, cdbg, peaks, summits]
        for f in outpaths:    
            #outbdg = ii[2] + "_output_FE.bdg"
            if False == os.path.exists(f):
                print "\n. I can't find the MACS2 output file", tbdg
                exit()
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
        os.system( get_setting("mpirun",con) + " fe_commands.sh" )
    exit()

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
    print "\n. FE bedgraphs OK."
    return

def bed2wig(con):
    cur = con.cursor()
    sql = "delete from WigFiles"
    cur.execute(sql)
    con.commit()
    
    sql = "select exp_annoid, bdgpath from MacsFE"
    cur.execute(sql)
    x = cur.fetchall()
    
    commands = []
    for ii in x:
        annoid = ii[0]
        bdgpath = ii[1]
        wigpath = re.sub(".bdg", ".wig", bdgpath)
        
        c = "python " + get_setting("seqtoolsdir", con) + "/bedgraph2wig.py " + bdgpath
        commands.append(c)
        
        sql = "insert or replace into WigFiles(exp_annoid, org_bdgpath, wigpath) VALUES("
        sql += annoid.__str__()
        sql += ",'" + bdgpath
        sql += "','" + wigpath + "')"
        cur.execute(sql)
        con.commit()

    fout = open("run_bed2wig.sh", "w")
    for c in commands:
        print c
        fout.write(c + "\n")
    fout.close()
    
    if get_setting("practice_mode", con) == "0":
        for c in commands:
            os.system(c)
    
def check_wig(con):
    cur = con.cursor()
    sql = "select exp_annoid, wigpath from WigFiles"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        annoid = ii[0]
        wigpath = ii[1]
        if False == os.path.exists( wigpath ):
            print "\n. Error, I can't find your WIG file ", wigpath
            exit()
        print wigpath
    print "\n. WIG files OK."
    return

def write_viz_config(con):
    cur = con.cursor()
    
    configpath = get_setting("project_name", con) + ".config"
    
    species = []
    sql = "select distinct species from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        species.append( ii[0] )
        
    for s in species:
        repgroups = []
        
        # A hack for Eugenio's data, in which the library_names are maybe wrong?
#         sql = "select library_name from Annotations where species='" + s + "'"
#         cur.execute(sql)
#         x = cur.fetchall()
#         for ii in x:
#             strain_name = ii[0].split("_")[2]
#             repgroups.append(strain_name)
        
        # But normally do this instead...
        sql = "select distinct strain from Annotations where species='" + s + "'"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            strain = ii[0]
            repgroups.append( strain )
        
        for strain in repgroups:
            replicates = []
            sql = "select distinct replicate from Annotations where species='" + s + "' and strain='" + strain + "'"
            cur.execute(sql)
            x = cur.fetchall()
            for ii in x:
                replicates.append( ii[0] )
            
            annoids = []
            for repid in replicates:
                sql = "select annoid from Annotations where species='" + s + "' and strain='" + strain + "' and replicate=" + repid.__str__()
                cur.execute(sql)
                x = cur.fetchall()
                for ii in x:
                    annoids.append( ii[0] )
                
            for annoid in annoids:
                print s, strain, repid, annoids
                #
                # continue here
                #
                
    