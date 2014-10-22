from annotation_db import *
import re, os, sys
from sets import Set

def extract_perfect_reads(annoid, con, chrom_filter = None):
    """Extracts reads that match perfectly, without mismatches."""
    """Writes a new (shorted) SAM file, and also returns a hash of read IDs."""
    cur = con.cursor()
    
    sql = "delete from Reads where annoid=" + annoid.__str__()
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
    
    print "\n. Extracting perfect reads from", sampath
    
    # writing is depricated here. the writing of the new sam file
    # now takes place AFTER we've determine which reads are unique
    # to each species.
    #print ". and writing new SAM output to", outsampath
    
    fin = open(sampath, "r")
    #fout = open(outsampath, "w")
    total = 0
    inserted = 0
    for line in fin.xreadlines():
        outline = ""
        if line.startswith("@"):
            continue
            #fout.write(line)
        elif line.__len__() > 5:
            total += 1
            if total%10000 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                con.commit()
            tokens = line.split()
            for t in tokens[10:]:
                if False == t.startswith("X"):
                    continue
                elif t == "XM:i:0":
                    readname = tokens[0]
                    
                    """The chrom_filter allows for reads to be considered
                    only if they match the filter, such "Chr1".
                    This is useful for generating toy-sized test databases
                    on a restricted number of reads."""
                    cflist = get_setting_list("chrom_filter",con)
                    if cflist.__len__() > 0:
                        chrom = tokens[2]
                        if chrom not in cflist:
                            continue
                    inserted += 1
                    sql = "insert into Reads(readname,annoid,mismatch,order_seen) VALUES("
                    sql += "'" + readname + "',"
                    sql += annoid.__str__() + ",0," + inserted.__str__() + ")"
                    cur.execute(sql)
                    
                    """We found all the stuff for this reads; let's skip remaining tokens for this
                    line, and just continue to the next read."""
                    continue # skip to the next line
    con.commit()
    sql = "select count(*) from Reads where annoid=" + annoid.__str__() + " and mismatch=0"
    cur.execute(sql)
    count = cur.fetchone()[0]
    ratio = float(count)/total
    print "\n\t--> I found", count, "perfect reads out of", total, "reads. (%.3f)"%ratio
    fin.close()
    #fout.close()
    
    sql = "insert or replace into ReadStats(annoid, nperfect, ntotal) VALUES("
    sql += annoid.__str__() + "," + count.__str__() + "," + total.__str__() + ")"
    cur.execute(sql)
    con.commit()


def find_hybrid_unique_reads(con):
    """For hybrids annoid1 and annoid2, find the reads that are unique to each,
    and write a new SAM file with only those reads."""
    cur = con.cursor()
    
    cur.execute("DELETE from UniqueReads")
    con.commit()
    
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

        print "\n. Finding unique reads between", annoid1, annoid2
           
        sql = "select readname, readid from Reads where annoid=" + annoid1.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        set1 = {}
        for ii in x:
            set1[ ii[0] ] = ii[1]
        
        sql = "select readname, readid from Reads where annoid=" + annoid2.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        set2 = {}
        for ii in x:
            set2[ ii[0] ] = ii[1]
        
        """Convert to Set, so that we can call the difference method"""
        names1 = set( set1.keys() )
        names2 = set( set2.keys() )
        
        """Process unique reads for the first species."""
        unique_read_names = names1.difference(names2)
        print "\n. N unique reads in", annoid1,"=", unique_read_names.__len__(), "of", names1.__len__(), "total reads."

        #
        # continue here
        # NOTE: the difference method could be eliminated, and sped-up, by taking advantage
        # of the new 'order' field in the Read table.
        # Basically, just walk two pointers through sorted order lists.
        #

        """Sanity Check"""
        """ The first read of unique_read_names should have just one entry in Reads"""
        if len(unique_read_names) == 0:
            print "\n. Warning: There are no reads that map uniquely to annotation ", annoid1
        else:
            rqqq = next(iter(unique_read_names))
            sql = "select count(*) from Reads where readname=\"" + rqqq + "\""
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count > 1:
                print "\n. We have a problem - python 153"
                print ".", rqqq, annoid1, annoid2
                exit()        
    
            """Write the Unique reads to the UniqueReads table."""
            print "\n. Updating the table UniqueReads"
            count = 0
            total_count = names1.__len__()
            for name in unique_read_names:
                count += 1
                if count%10000 == 0:
                    sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                    sys.stdout.flush()
                    con.commit()
                
                readid = set1[name]
                sql = "insert into UniqueReads(readid, annoid) VALUES("
                sql += readid.__str__() + ","
                sql += annoid1.__str__() + ")"
                cur.execute(sql)  
        
        sys.stdout.write("\r    --> %100.0f\n")
        sys.stdout.flush()      
        
        """Repeat for the second species."""
        unique_read_names = names2.difference(names1)
        print "\n. N unique reads in", annoid2,"=", unique_read_names.__len__(), "of", names2.__len__(), "total reads."
        
        if len(unique_read_names) == 0:
            print "\nWarning: There are no reads that map uniquely to annotation ", annoid2
        else:
         
            # see not above: this difference call could be eliminated.
            # continue here
            # to-do
        
            rqqq = next(iter(unique_read_names))
            sql = "select count(*) from Reads where readname=\"" + rqqq + "\""
            cur.execute(sql)
            count = cur.fetchone()[0]
            if count > 1:
                print "\n. We have a problem - python 153"
                print ".", rqqq, annoid1, annoid2
                exit()
        
            unique_readids = []
            count = 0
            total_count = names2.__len__()
            for name in unique_read_names:
                count += 1
                if count%10000 == 0:
                    sys.stdout.write("\r    --> %.1f%%" % (100*count/float(total_count)) )
                    sys.stdout.flush()
                    con.commit()
                    
                readid = set2[name]
                sql = "insert into UniqueReads(readid, annoid) VALUES("
                sql += readid.__str__() + ","
                sql += annoid2.__str__() + ")"
                cur.execute(sql)            
            con.commit()
        
        sys.stdout.write("\r    --> %100.0f\n")
        sys.stdout.flush()   
        
        """Now check what we inserted."""
        sql = "select count(*) from UniqueReads where annoid=" + annoid1.__str__()
        cur.execute(sql)
        count1 = cur.fetchone()[0]
        
        sql = "select count(*) from UniqueReads where annoid=" + annoid2.__str__()
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
    cur = con.cursor()
    sql = "select * from Hybrids"
    cur.execute(sql)
    x = cur.fetchall()
    annoids = []
    for ii in x:
        annoids.append( ii[0] )
        
    print "\nlibrary_name\tN_total\tN_perfect\tN_unique"
    fout = open("read_stats.xls", "w")
        
    for annoid in annoids:
        sql = "select nperfect, ntotal from ReadStats where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        nperfect = x[0]
        ntotal = x[1]
        sql = "select nunique from UniqueReadStats where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        nunique = x[0]
        sql = "select library_name, species from Annotations where annoid=" + annoid.__str__()
        cur.execute(sql)
        x = cur.fetchone()
        rowname = x[0] + "-" + x[1]
        l = rowname + "\t" + ntotal.__str__() + "\t" + nperfect.__str__() + "\t" + nunique.__str__() + "\t" + "%.3f"%(nperfect/float(ntotal)) + "\t" + "%.3f"%(nunique/float(ntotal))
        print l
        fout.write(l + "\n")
    fout.close()
        
def write_filtered_sam(con):
    """Writes a new SAM file containing only those reads that are non-mismatches
    and which are unique to a hybrid parent species."""
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
            samoutpath = get_setting("outdir") + "/" + re.sub(".fastq", "-" + species + ".unique.sam", fastq)
        
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
        sql = "select readname from Reads where Reads.readid in (select readid from UniqueReads where annoid=" + annoid.__str__() + ") order by Reads.order_seen"
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
    