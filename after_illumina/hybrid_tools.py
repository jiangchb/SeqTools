from annotation_db import *
import re, os, sys
from sets import Set

def find_hybrid_unique_reads(con):
    """For hybrids X and Y, find the reads that uniquely
        mapped to X or Y, but not both. Before this method was called,
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
        annoid1 = pair[0]
        annoid2 = pair[1]

        sql = "CREATE TABLE IF NOT EXISTS UniqueReads" + annoid1.__str__() + "(readid INTEGER primary key)" # reads that are unique to an annotation
        cur.execute(sql)
        con.commit()
        sql = "CREATE TABLE IF NOT EXISTS UniqueReads" + annoid2.__str__() + "(readid INTEGER primary key)" # reads that are unique to an annotation
        cur.execute(sql)
        con.commit()

        sql = "delete from UniqueReads" + annoid1.__str__()
        cur.execute(sql)
        sql = "delete from UniqueReads" + annoid2.__str__()
        cur.execute(sql)
        con.commit()

        print "\n. Finding unique reads between", annoid1, annoid2
           
        sql = "select readname, readid from Reads" + annoid1.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        set1 = {}
        for ii in x:
            set1[ ii[0] ] = ii[1]
        
        sql = "select readname, readid from Reads" + annoid2.__str__()
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
        print "\n. N unique reads in", annoid1,"=", unique_read_names.__len__(), "of", names1.__len__(), "total reads."


#         """Sanity Check"""
        if len(unique_read_names) == 0:
            print "\n. Warning: There are no reads that map uniquely to annotation ", annoid1
            write_error(con, "There are no reads that map uniquely to annotation " + annoid1)       
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
                sql = "insert into UniqueReads" + annoid1.__str__() + "(readid) VALUES("
                sql += readid.__str__() + ")"
                cur.execute(sql)  
        
        sys.stdout.write("\r    --> %100.0\n")
        sys.stdout.flush()      
        
        """Repeat for the second species."""
        unique_read_names = names2.difference(names1)
        print "\n. N unique reads in", annoid2,"=", unique_read_names.__len__(), "of", names2.__len__(), "total reads."
        
        if len(unique_read_names) == 0:
            write_error(con, "There are no reads that map uniquely to the genome for annotation " + annoid2)
            print "\nWarning: There are no reads that map uniquely to annotation ", annoid2
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
                    sql = "insert into UniqueReads" + annoid2.__str__() + "(readid) VALUES("
                    sql += readid.__str__() + ")"
                cur.execute(sql)            
            except:
                con.rollback()
            con.commit()
        
        sys.stdout.write("\r    --> %100.0f\n")
        sys.stdout.flush()   
        
        """Now check what we inserted."""
        sql = "select count(*) from UniqueReads" + annoid1.__str__()
        cur.execute(sql)
        count1 = cur.fetchone()[0]
        
        sql = "select count(*) from UniqueReads" + annoid2.__str__()
        cur.execute(sql)
        count2 = cur.fetchone()[0]
        
        sql = "insert or replace into UniqueReadStats(annoid, nunique) VALUES("
        sql += annoid1.__str__() + "," + count1.__str__() + ")"
        cur.execute(sql)
        sql = "insert or replace into UniqueReadStats(annoid, nunique) VALUES("
        sql += annoid2.__str__() + "," + count2.__str__() + ")"
        cur.execute(sql) 
        con.commit()

def print_read_stats(con):
    """This method prints basic stats about the extent to which reads mapped to genomes."""    
    cur = con.cursor()
    sql = "select annoid from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    annoids = []
    for ii in x:
        annoids.append( ii[0] )

    fout = open("read_stats.xls", "w")        
    print "\nlibrary_name\tN_total\tN_perfect\tN_unique"
    fout.write("\nlibrary_name\tN_total\tN_perfect\tN_unique\n")
     
    for annoid in annoids:
                
        sql = "select nperfect, ntotal from ReadStats where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        nperfect = x[0]
        ntotal = x[1] # total reads mapped to this genome
        
        """Is this annotation part of a hybrid genome?"""
        sql = "select count(*) from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()[0]
        if x > 0:
            """Yes, it's part of a hybrid."""
            sql = "select nunique from UniqueReadStats where annoid=" + annoid.__str__()
            cur.execute(sql)
            x = cur.fetchone()
            nunique = x[0]
        else:
            """Nope, not a hybrid"""
            nunique = ntotal
        
        sql = "select library_name, species from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        rowname = x[0] + "-" + x[1]
        l = rowname + "\t" + ntotal.__str__() + "\t" + nperfect.__str__() + "\t" + nunique.__str__() + "\t" + "%.3f"%(nperfect/float(ntotal)) + "\t" + "%.3f"%(nunique/float(ntotal))
        print l
        fout.write(l + "\n")
    fout.close()

def print_read_histograms(con):
    print "\n. Plotting read matchcount histograms...."
    
    cur = con.cursor()

    #sql = "select * from Hybrids"
    sql = "select * from Annotations"
    cur.execute(sql)
    x = cur.fetchall()
    annoids = []
    for ii in x:
        annoids.append( ii[0] )
    
    #print "\n. Finding the maximum read mismatch..."
    #sql = "select max(mismatch) from Reads"
    #cur.execute(sql)
    #max_mismatch = cur.fetchone()[0] # the maximum X value
    max_mismatch = 18

#     """This next loop is the hard work -- get all the mismatch counts from the DB"""
#     annoid_mismatches = {}
#     for annoid in annoids:
#         sql = "select mismatch from Reads where annoid=" + annoid.__str__()
#         cur.execute(sql)
#         x = cur.fetchall()
#         """vals is a list of all mismatchd scores."""
#         vals = []
#         for ii in x:
#             vals.append( ii[0] )
#         annoid_mismatches[annoid] = vals
#         print "\n. Annotation", annoid

#     annoid_totalreads = {}
#     for annoid in annoids:
#         sql = "select count(*) from Reads where annoid=" + annoid.__str__()
#         cur.execute(sql)
#         totalreads = cur.fetchone()[0]
#         annoid_totalreads[annoid] = totalreads
#         print "\n. Annotation", annoid, "has", totalreads, "total reads."

    annoid_totalreads = {}
    annoid_bars = {} # key = annoid, value = list of bargraph Y values
    max_count = 0 # the maximum Y value
    for annoid in annoids:
        print "\n. Counting mismatches for annotation", annoid
        
        annoid_totalreads[annoid] = 0
        annoid_bars[annoid] = []

        for ii in range(0, max_mismatch+1):
            sql = "select count(*) from Reads" + annoid.__str__() + " where mismatch=" + ii.__str__()
            cur.execute(sql)
            this_count = cur.fetchone()[0]
            annoid_totalreads[annoid] += this_count
            annoid_bars[annoid].append( this_count )
            
            if max_count < this_count:
                max_count = this_count
    
        print "\n. Annotation", annoid, "has", annoid_totalreads[annoid], "total reads."
    

    """This next loop is for hybrids only."""
    annoid_countunique = {} # key = annoid, value count of unique reads
    annoid_unique_bars = {} # data for uniquely-mapped reads
    for annoid in annoids:
        print "\n. Counting unique reads for hybrid annotation", annoid
        sql = "select count(*) from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        if cur.fetchone()[0] > 0: # i.e., is this annotation a hybrid?
            annoid_countunique[annoid] = 0
            annoid_unique_bars[annoid] = []
            for ii in range(0, max_mismatch+1):
                sql = "select count(*) from Reads" + annoid.__str__() + " where mismatch=" + ii.__str__() + " and readid in (select readid from UniqueReads" + annoid.__str__() + ")"
                cur.execute(sql)
                this_count = cur.fetchone()[0]
                annoid_countunique[annoid] += this_count
                annoid_unique_bars[annoid].append( this_count )        

    bins = range(0, max_mismatch+1)

    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    pdfpath = get_setting("project_name", con) + ".read_histograms.pdf"
    pdf_pages = PdfPages(pdfpath)
    
    for ii in range(0, annoids.__len__()):
        annoid = annoids[ii]
        sql = "select library_name, species from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        this_library_name = x[0]
        this_species = x[1]
        
        fig = plt.figure(figsize=(8,4))
        p1 = plt.bar(bins, annoid_bars[annoid], 0.75, log=1, align="center", color="#99CCFF")
        if annoid in annoid_countunique: # i.e., is this annotation a hybrid?
            p2 = plt.bar(bins, annoid_unique_bars[annoid], 0.75, align="center", log=1, color="#0066CC")
        plt.xlim( -1,max_mismatch+1 )
        plt.ylim( 1,max_count    )
        plt.xlabel('Mismatched Sites')
        plt.ylabel('Reads')
        plt.title('Histogram of Read Mismatches -- ' + this_library_name)
        if annoid in annoid_countunique: # i.e., is this annotation a hybrid?
            plt.legend( (p1[0], p2[0]), ('Reads Mapped to ' + this_species + ' (' + annoid_totalreads[annoid].__str__() + ')', 'Reads Unique to ' + this_species + ' (' + annoid_countunique[annoid].__str__() + ')') )
        plt.tight_layout()
        pdf_pages.savefig(fig)
    pdf_pages.close()


def write_filtered_sam(con):
    """Writes a new SAM file containing only those reads whose mismatch level is
    equal or below the mismatch threshold,
    and (for hybrids) which are uniquely mapped to only one of the hybrid parent species."""
    cur = con.cursor()
    sql = "select * from Hybrids"
    cur.execute(sql)
    x = cur.fetchall()
    annoids = []
    for ii in x:
        annoids.append( ii[0] )
        
    for annoid in annoids:
        sql = "select library_name, species, fastqpath from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        library_name = x[0] + "-" + x[1]
        fastq = x[2]
        species = x[1]
        
        sql = "select count(*) from Hybrids where annoid=" + annoid.__str__()
        cur.execute(sql)
        count = cur.fetchone()[0]
        samoutpath = re.sub(".fastq", ".unique.sam", fastq)
        if count > 0:
            """Hybrids get a special SAM path"""
            samoutpath = get_setting("outdir", con) + "/" + re.sub(".fastq", "-" + species + ".unique.sam", fastq)
        
        print "\n. I'm writing the perfect reads that are unique in", library_name, "to a new SAM file:"
        print "\t", samoutpath
        
        """Read the header from the original SAM file."""
        sql = "select sampath from BowtieOutput where annoid=" + annoid.__str__()
        cur.execute(sql)
        bowtie_sampath = cur.fetchone()[0]
        fin = open(bowtie_sampath, "r")
        header_lines = []
        for l in fin.xreadlines():
            if False == l.startswith("@"):
                break # stop parsing when the lines don't contain header hashes
            header_lines.append( l )
        fin.close()
         
        """Finally, write the new SAM file."""
        sys.stdout.write(".")
        sys.stdout.flush()
        sql = "select readname from Reads" + annoid.__str__() + " where Reads" + annoid.__str__() + ".readid in (select readid from UniqueReads" + annoid.__str__() + ") order by Reads" + annoid.__str__() + ".order_seen"
        cur.execute(sql)
        x = cur.fetchall()
        
        """readnames is a list of strings, each a unique name for a read."""
        count = 0
        readnames = Set([])
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
        
        sql = "insert or replace into FilteredBowtieOutput (annoid, sampath) VALUES("
        sql += annoid.__str__() + ",'" + samoutpath + "')"
        cur.execute(sql)
        con.commit()
        
        print "\n. I found", count_good_reads, "reads for the new SAM file at", samoutpath
    