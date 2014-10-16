import re, os, sys

def extract_perfect_reads(annoid, con):
    """Extracts reads that match perfectly, without mismatches."""
    """Writes a new (shorted) SAM file, and also returns a hash of read IDs."""
    cur = con.cursor()
    
    sql = "delete from AnnoReads where annoid=" + annoid.__str__()
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
                    sql = "insert or replace into Reads(readname) VALUES("
                    sql += "'" + readname + "')"
                    cur.execute(sql)
                    
                    sql = "select last_insert_rowid()"
                    cur.execute(sql)
                    x = cur.fetchone()
                    readid = x[0]
                    
                    sql = "insert into AnnoReads(readid,annoid,mismatch) VALUES("
                    sql += readid.__str__() + "," + annoid.__str__() + ",0)"
                    cur.execute(sql)
                    
                    #fout.write(line)
                    continue # skip to the next line
    con.commit()
    sql = "select count(*) from AnnoReads where annoid=" + annoid.__str__() + " and mismatch=0"
    cur.execute(sql)
    count = cur.fetchone()[0]
    print "\n. OK, I found", count, "perfect reads, out of", total, "reads."
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

        """Find reads unique to annoid1"""
        sql = "select readid from AnnoReads where annoid=" + annoid1.__str__() + " and readid not in (select readid from AnnoReads where annoid=" + annoid2.__str__() + ")"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            sql = "insert into UniqueReads(readid, annoid) VALUES("
            sql += ii[0].__str__() + ","
            sql += annoid1.__str__() + ")"
            #print sql
            cur.execute(sql)
        con.commit()

        sql = "select readid from AnnoReads where annoid=" + annoid2.__str__() + " and readid not in (select readid from AnnoReads where annoid=" + annoid1.__str__() + ")"
        cur.execute(sql)
        x = cur.fetchall()
        for ii in x:
            sql = "insert into UniqueReads(readid, annoid) VALUES("
            sql += ii[0].__str__() + ","
            sql += annoid2.__str__() + ")"
            #print sql
            cur.execute(sql)    
        con.commit()   
                
        sql = "select count(*) from UniqueReads where annoid=" + annoid1.__str__()
        cur.execute(sql)
        count1 = cur.fetchone()[0]
        
        sql = "select count(*) from UniqueReads where annoid=" + annoid2.__str__()
        cur.execute(sql)
        count2 = cur.fetchone()[0]
        
        print "\n. I found", count1, "unique reads for annotation", annoid1
        print "\n. I found", count2, "unique reads for annotation", annoid2

#     sql = "insert or replace into FilteredBowtieOutput(annoid, sampath) VALUES("
#     sql += annoid.__str__() + ",'" + outsampath + "')"
#     cur.execute(sql)
#     con.commit()
            
    