from annotation_db import *
import re, os, sys
from sets import Set

def find_hybrid_unique_reads(con):
    """For hybrid READS X and Y, find the reads that uniquely
        map to X or Y, but not both. Before this method was called,
        reads were mapped to X and Y using Bowtie and the method
        'extract_matched_reads'.
        
        Output: fills the table UniqueReads, UniqueReadsStats,
        and writes a new SAM file containing only those reads that
        were uniquely mapped."""
    cur = con.cursor()
    
    """pairs is a list of hybrid pairs."""
    pairs = []
    sql = "select * from HybridPairs"
    cur.execute(sql)
    x = cur.fetchall()
    for ii in x:
        pairs.append( (ii[0],ii[1]) )
    
    for pair in pairs:
        readid1 = pair[0]
        readid2 = pair[1]

        sql = "CREATE TABLE IF NOT EXISTS UniqueReads" + readid1.__str__() + "(readid INTEGER primary key)" # reads that are unique to an annotation
        cur.execute(sql)
        con.commit()
        sql = "CREATE TABLE IF NOT EXISTS UniqueReads" + readid2.__str__() + "(readid INTEGER primary key)" # reads that are unique to an annotation
        cur.execute(sql)
        con.commit()

        sql = "delete from UniqueReads" + readid1.__str__()
        cur.execute(sql)
        sql = "delete from UniqueReads" + readid2.__str__()
        cur.execute(sql)
        con.commit()

        print "\n. Finding unique reads between", readid1, readid2
           
        sql = "select readname, readid from Reads" + readid1.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        set1 = {}
        for ii in x:
            set1[ ii[0] ] = ii[1]
        
        sql = "select readname, readid from Reads" + readid2.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        set2 = {}
        for ii in x:
            set2[ ii[0] ] = ii[1]
        
        """Convert to Set, so that we can call the difference method"""
        names1 = set( set1.keys() ) # read names in annotation 1
        names2 = set( set2.keys() ) # read names in annotation 2
        
        """Process unique reads for the first species."""
        unique_read_names = names1.difference(names2)
        print "\n. N unique reads in", readid1,"=", unique_read_names.__len__(), "of", names1.__len__(), "total reads."

#         """Sanity Check"""
        if len(unique_read_names) == 0:
            print "\n. Warning: There are no reads that map uniquely to READ ", readid1
            write_error(con, "There are no reads that map uniquely to READ " + readid1)       
        else:
            """Write the Unique reads to the UniqueReads table."""
            print "\n. Updating the unique reads table"
            count = 0
            total_count = names1.__len__()
            for name in unique_read_names:
                count += 1
                if count%10000 == 0:
                    sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                    sys.stdout.flush()
                    con.commit()
                
                readid = set1[name]
                sql = "insert into UniqueReads" + readid1.__str__() + "(readid) VALUES("
                sql += readid.__str__() + ")"
                cur.execute(sql)  
        
        sys.stdout.write("\r    --> %100.0\n")
        sys.stdout.flush()      
        
        """Repeat for the second species."""
        unique_read_names = names2.difference(names1)
        print "\n. N unique reads in", readid2,"=", unique_read_names.__len__(), "of", names2.__len__(), "total reads."
        
        if len(unique_read_names) == 0:
            write_error(con, "There are no reads that map uniquely to the genome for READ " + readid2)
            print "\nWarning: There are no reads that map uniquely to READ ", readid2
        else:
            unique_readids = []
            count = 0
            total_count = names2.__len__()
            try:
                for name in unique_read_names:
                    count += 1
                    if count%10000 == 0:
                        sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                        sys.stdout.flush()
                        con.commit()
                        
                    readid = set2[name]
                    sql = "insert into UniqueReads" + readid2.__str__() + "(readid) VALUES("
                    sql += readid.__str__() + ")"
                cur.execute(sql)            
            except:
                con.rollback()
            con.commit()
        
        sys.stdout.write("\r    --> %100.0f\n")
        sys.stdout.flush()   
        
        """Now check what we inserted."""
        sql = "select count(*) from UniqueReads" + readid1.__str__()
        cur.execute(sql)
        count1 = cur.fetchone()[0]
        
        sql = "select count(*) from UniqueReads" + readid2.__str__()
        cur.execute(sql)
        count2 = cur.fetchone()[0]
        
        sql = "insert or replace into UniqueReadStats(readid, nunique) VALUES("
        sql += readid1.__str__() + "," + count1.__str__() + ")"
        cur.execute(sql)
        sql = "insert or replace into UniqueReadStats(readid, nunique) VALUES("
        sql += readid2.__str__() + "," + count2.__str__() + ")"
        cur.execute(sql) 
        con.commit()

def print_read_stats(con):
    """This method prints basic stats about the extent to which reads mapped to genomes."""    
    cur = con.cursor()
    sql = "select id, name from Reads"
    cur.execute(sql)
    x = cur.fetchall()
    readids = {}
    for ii in x:
        readids[ ii[0] ] = ii[1]

    fout = open("read_stats.xls", "w")        
    print "\nlibrary_name\tN_total\tN_perfect\tN_unique"
    fout.write("\nlibrary_name\tN_total\tN_perfect\tN_unique\n")
     
    for readid in readids:
        sql = "select nperfect, ntotal from ReadStats where readid=" + readid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        nperfect = x[0]
        ntotal = x[1] # total reads mapped to this genome
        
        """Is this annotation part of a hybrid genome?"""
        sql = "select count(*) from HybridPairs where readid1=" + readid.__str__() + " or readid2=" + readid.__str__()
        cur.execute(sql)
        x = cur.fetchone()[0]
        if x > 0:
            """Yes, it's part of a hybrid."""
            sql = "select nunique from UniqueReadStats where readid=" + readid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            nunique = x[0]
        else:
            """Nope, not a hybrid"""
            nunique = ntotal
        
        rowname = readids[ readid ]
        l = rowname + "\t" + ntotal.__str__() + "\t" + nperfect.__str__() + "\t" + nunique.__str__() + "\t" + "%.3f"%(nperfect/float(ntotal)) + "\t" + "%.3f"%(nunique/float(ntotal))
        print l
        fout.write(l + "\n")
    fout.close()

def print_read_histograms(con):
    print "\n. Plotting read matchcount histograms...."
    
    cur = con.cursor()

    #sql = "select * from Hybrids"
    sql = "select id from Reads"
    cur.execute(sql)
    x = cur.fetchall()
    readids = []
    for ii in x:
        readids.append( ii[0] )

    max_mismatch = 18

    readid_totalreads = {}
    readid_bars = {} # key = annoid, value = list of bargraph Y values
    max_count = 0 # the maximum Y value
    for readid in readids:
        print "\n. Counting mismatches for READ", readid
        
        readid_totalreads[readid] = 0
        readid_bars[readid] = []

        for ii in range(0, max_mismatch+1):
            sql = "select count(*) from Reads" + readid.__str__() + " where mismatch=" + ii.__str__()
            cur.execute(sql)
            this_count = cur.fetchone()[0]
            readid_totalreads[readid] += this_count
            readid_bars[readid].append( this_count )
            
            if max_count < this_count:
                max_count = this_count
    
        print "\n. READ", readid, "has", readid_totalreads[readid], "total reads."
    

    """This next loop is for hybrids only."""
    readid_countunique = {} # key = readid, value count of unique reads
    readid_unique_bars = {} # data for uniquely-mapped reads
    for readid in readids:
        print "\n. Counting unique reads for hybrid read", readid
        sql = "select count(*) from HybridPairs where readid1=" + readid.__str__() + " or readid2=" + readid.__str__()
        cur.execute(sql)
        if cur.fetchone()[0] > 0: # i.e., is this annotation a hybrid?
            readid_countunique[readid] = 0
            readid_unique_bars[readid] = []
            for ii in range(0, max_mismatch+1):
                sql = "select count(*) from Reads" + readid.__str__() + " where mismatch=" + ii.__str__() + " and readid in (select readid from UniqueReads" + readid.__str__() + ")"
                cur.execute(sql)
                this_count = cur.fetchone()[0]
                readid_countunique[readid] += this_count
                readid_unique_bars[readid].append( this_count )        

    bins = range(0, max_mismatch+1)

    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    pdfpath = get_setting("project_name", con) + ".read_histograms.pdf"
    pdf_pages = PdfPages(pdfpath)
    
    for ii in range(0, readids.__len__()):
        readid = readids[ii]
        sql = "select name, speciesid from Reads where id=" + readid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        this_library_name = x[0]
        speciesid = x[1]
        sql = "Select name from Species where id=" + speciesid.__str__()
        cur.execute(sql)
        yy = cur.fetchone()
        this_species = yy[0]
        
        fig = plt.figure(figsize=(8,4))
        p1 = plt.bar(bins, readid_bars[readid], 0.75, log=1, align="center", color="#99CCFF")
        if readid in readid_countunique: # i.e., is this annotation a hybrid?
            p2 = plt.bar(bins, readid_unique_bars[readid], 0.75, align="center", log=1, color="#0066CC")
        plt.xlim( -1,max_mismatch+1 )
        plt.ylim( 1,max_count    )
        plt.xlabel('Mismatched Sites')
        plt.ylabel('Reads')
        plt.title('Histogram of Read Mismatches -- ' + this_library_name)
        if readid in readid_countunique: # i.e., is this annotation a hybrid?
            plt.legend( (p1[0], p2[0]), ('Reads Mapped to ' + this_species + ' (' + readid_totalreads[readid].__str__() + ')', 'Reads Unique to ' + this_species + ' (' + readid_countunique[readid].__str__() + ')') )
        plt.tight_layout()
        pdf_pages.savefig(fig)
    pdf_pages.close()


def write_filtered_sam(con):
    """Writes a new SAM file containing only those reads whose mismatch level is
    equal or below the mismatch threshold,
    and (for hybrids) which are uniquely mapped to only one of the hybrid parent species."""
    cur = con.cursor()
    sql = "select readid1, readid2 from HybridPairs"
    cur.execute(sql)
    x = cur.fetchall()
    readids = []
    for ii in x:
        readid1 = ii[0]
        readid2 = ii[1]
        if readid1 not in readids:
            readids.append( readid1 )
        if readid2 not in readids:
            readids.append( readid2 )
        
    for readid in readids:
        sql = "select name, speciesid, fastqid from Reads where readid=" + readid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        library_name = x[0]
        fastqid = x[2]
        sql = "select filepath from FastqFiles where id=" + fastqid.__str__()
        cur.execute(sql)
        fastq = cur.fetchone()[0]
#         speciesid = x[1]
#         sql = "Select name from Species where id=" + speciesid.__str__()
#         cur.execute(sql)
#         species = cur.fetchone()[0]
        
        samoutpath = re.sub(".fastq", ".unique.sam", fastq)
        
        print "\n. I'm writing the perfect reads that are unique in", library_name, "to a new SAM file:"
        print "\t", samoutpath
        
        """Read the header from the original SAM file."""
        sql = "select sampath from BowtieOutput where readid=" + readid.__str__()
        cur.execute(sql)
        bowtie_sampath = cur.fetchone()[0]
        fin = open(bowtie_sampath, "r")
        header_lines = []
        for l in fin.xreadlines():
            if False == l.startswith("@"):
                break # stop parsing when the lines don't contain header hashes
            header_lines.append( l )
        fin.close()
         
        """readnames is a list of strings, each a unique name for a read."""
        readnames = Set([])
        sys.stdout.write(".")
        sys.stdout.flush()
        sql = "select readname from Reads" + readid.__str__() + " where Reads" + readid.__str__() + ".readid in (select readid from UniqueReads" + readid.__str__() + ") order by Reads" + readid.__str__() + ".order_seen"
        cur.execute(sql)
        x = cur.fetchall()
        count = 0
        for ii in x:
            readnames.add( ii[0] )
            count += 1
            if count%10000 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
        
        if get_setting("practice_mode", con) == "0":
            """Open the file and write the header lines."""
            fout = open(samoutpath, "w")
            for hl in header_lines:
                fout.write(hl)
             
            """Walk through the original SAM file and copy relevant lines into the new SAM file."""
            count = 0
            count_good_reads = 0
            fin = open(bowtie_sampath, "r")
            for erg_line in fin.xreadlines():
                count += 1
                if count%10000 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                
                if False == erg_line.startswith("@"):
                    this_readid = erg_line.split()[0]
                    
                    """WARNING: this next comparison is computationally expensive!"""
                    if this_readid in readnames:
                        """Copy the original line into the new SAM path, because this read satisfies what we're looking for."""
                        count_good_reads += 1
                        fout.write(erg_line)
             
            fin.close()
            fout.close()
        
        sql = "insert or replace into FilteredBowtieOutput (readid, sampath) VALUES("
        sql += readid.__str__() + ",'" + samoutpath + "')"
        cur.execute(sql)
        con.commit()
        
        print "\n. I found", count_good_reads, "reads for the new SAM file at", samoutpath
    