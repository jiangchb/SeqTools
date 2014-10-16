import re, os, sys

def extract_perfect_reads(annoid, con):
    """Extracts reads that match perfectly, without mismatches."""
    """Writes a new (shorted) SAM file, and also returns a hash of read IDs."""
    cur = con.cursor()
    
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
            tokens = line.split()
            for t in tokens[10:]:
                if False == t.startswith("X"):
                    continue
                elif t == "XM:i:0":
                    readname = tokens[0]
                    sql = "insert or replace into Reads(readname,annoid,mismatch) VALUES("
                    sql += "'" + readname + "'," + annoid.__str__() + ",0)"
                    cur.execute(sql)
                
                    #fout.write(line)
                    continue # skip to the next line
    con.commit()
    sql = "select count(*) from Reads where annoid=" + annoid.__str__() + " and mismatch=0"
    cur.execute(sql)
    count = cur.fetchone()[0]
    print "\n. OK, I found", count, "perfect reads, out of", total, "reads."
    fin.close()
    #fout.close()
    
    sql = "insert or replace into BowtieFilterStats(annoid, nperfect, ntotal) VALUES("
    sql += annoid.__str__() + "," + count.__str__() + "," + total.__str__() + ")"
    cur.execute(sql)
    con.commit()
    
#     sql = "insert or replace into FilteredBowtieOutput(annoid, sampath) VALUES("
#     sql += annoid.__str__() + ",'" + outsampath + "')"
#     cur.execute(sql)
#     con.commit()
            
    