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
                    if chrom_filter != None:
                        chrom = tokens[2]
                        if False == chrom.__contains__(chrom_filter):
                            continue
                    sql = "insert into Reads(readname,annoid,mismatch) VALUES("
                    sql += "'" + readname + "',"
                    sql += annoid.__str__() + ",0)"
                    cur.execute(sql)
                    
                    #fout.write(line)
                    continue # skip to the next line
    con.commit()
    sql = "select count(*) from Reads where annoid=" + annoid.__str__() + " and mismatch=0"
    cur.execute(sql)
    count = cur.fetchone()[0]
    ratio = float(count)/total
    print "\n. OK, I found", count, "perfect reads out of", total, "reads (",ratio,") in", sampath
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

        """Sanity Check"""
        """ The first read of unique_read_names should have just one entry in Reads"""
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
        for name in unique_read_names:
            count += 1
            if count%10000 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                con.commit()
            
            readid = set1[name]
            sql = "insert into UniqueReads(readid, annoid) VALUES("
            sql += readid.__str__() + ","
            sql += annoid1.__str__() + ")"
            cur.execute(sql)        
        
        """Repeat for the second species."""
        unique_read_names = names2.difference(names1)
        print "\n. N unique reads in", annoid2,"=", unique_read_names.__len__(), "of", names2.__len__(), "total reads."
         
        rqqq = next(iter(unique_read_names))
        sql = "select count(*) from Reads where readname=\"" + rqqq + "\""
        cur.execute(sql)
        count = cur.fetchone()[0]
        if count > 1:
            print "\n. We have a problem - python 153"
            print ".", rqqq, annoid1, annoid2
            exit()
    
        count = 0
        for name in unique_read_names:
            count += 1
            if count%10000 == 0:
                sys.stdout.write(".")
                sys.stdout.flush()
                con.commit()
                
            readid = set2[name]
            sql = "insert into UniqueReads(readid, annoid) VALUES("
            sql += readid.__str__() + ","
            sql += annoid2.__str__() + ")"
            cur.execute(sql)            
        con.commit()
        
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
            
    