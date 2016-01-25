from Bio import SeqIO

from chipseqdb import *
from chipseqdb_api import *

from argParser import ArgParser
ap = ArgParser(sys.argv)

def read_motifs(motif_path):
    """Reads PSWM data from a text file.
    Motifs are expected to look like this:
    
    /NAME=[SWISS] ABF1 TYPE=PSFM SPECIES=S.cer
    0.208782    0.155386    0.135415    0.500414
    0.745936    0.0    0.254063    0.0
    0.0    0.0    0.0    1.0
    0.0    1.0    0.0    0.0
    
    The motif file can contain multiple motifs
    """
    if False == os.path.exists(motif_path):
        print "ERROR: I can't find the motif path:", motif_path
        exit()
    
    gene_motif = {} # key = gene name, value = list of lists, each sublist if (A,C,G,T) weights.
    
    fin = open(motif_path, "r")
    curr_motif = []
    curr_gene = None
    for line in fin.xreadlines():
        if line.startswith("#"):
            continue
        
        if line.__len__() < 3:
            continue
        
        if line.startswith("/NAME"):
            tokens = line.split()
            gene = tokens[1]
            print "found gene:", gene
            if curr_gene != None:
                gene_motif[ curr_gene ] = curr_motif
            curr_gene = gene
            curr_motif = []
        else:
            tokens = line.split()
            if tokens.__len__() < 4:
                print "Hmmm, something is wrong with a line in the motif file. I was expecting four decimal numbers:"
                print line
                exit()
            wa = float(tokens[0])
            wc = float(tokens[1])
            wg = float(tokens[2])
            wt = float(tokens[3])
            curr_motif.append( (wa,wc,wg,wt) )
    
    gene_motif[ curr_gene ] = curr_motif
    
    return gene_motif

def build_motif_dbtables(con):
    cur = con.cursor()
    sql = "create table if not exists Motifs(id INTEGER primary key, name TEXT)"    
    cur.execute(sql)
    sql = "create table if not exists Summits2MotifScores(summitid INTEGER, motifid INTEGER, maxmotifscore FLOAT, maxmotifsite INT)"
    cur.execute(sql)
    con.commit()
    return con
    
def build_motif_table(con, motifid, motifdata):
    """motifdata is a list of lists. motifdata[site-1] = (wa,wc,wg,wt)"""
    cur = con.cursor()
    sql = "create table if not exists Motif" + motifid.__str__() + "(site INT, state TEXT, weight FLOAT)"
    cur.execute(sql)
    con.commit()

def read_genome(genome_path):
    from Bio import SeqIO
    handle = open(genome_path, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        print record.id
    handle.close()
    
    
##############################
#
# main
#

motifpath = ap.getArg("--motifpath")

readsdbpath = ap.getArg("--readdbpath")
rcon = lite.connect(readsdbpath, timeout=1)

vizdbpath = ap.getArg("--vizdbpath")
vcon = build_db(dbpath=vizdbpath)
vcon = build_motif_dbtables(vcon)

print "98:", vcon, vizdbpath

"""Import motifs"""
vcur = vcon.cursor()
gene_motif = read_motifs(motifpath)
for genename in gene_motif:
    sql = "select id from Motifs where name='" + genename + "'"
    vcur.execute(sql)
    xx = vcur.fetchall()
    if xx == None:
        sql = "insert into Motifs (name) values(" + genename + ")"
        vcur.execute(sql)
        vcon.commit()
        sql = "select id from Motifs where name='" + genename + "'"
        vcur.execute(sql)
    for ii in vcur.fetchall():
        motifid = ii[0]
        build_motif_table(con, motifid, gene_motif[genename])

rcur = rcon.cursor()
speciesid_name = {}
sql = "select id, name from Species"
rcur.execute(sql)
for ii in rcur.fetchall():
    id = ii[0]
    name = ii[1]
    speciesid_name[id] = name
    
speciesid_genomepath = {}
sql = "select speciesid, genomepath from SpeciesGenomepath"
rcur.execute(sql)
for ii in rcur.fetchall():
    id = ii[0]
    genomepath = ii[1]
    speciesid_genomepath[id] = genomepath

"""Parse each genome."""
for speciesid in speciesid_genomepath:
    genomepath = speciesid_genomepath[speciesid]
    handle = open(genomepath, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        """Check that this chromosom exists in the viz DB"""
        sql = "select id from Chromosomes where name='" + chromname + "' and species=" + speciesid.__str__()
        vcur.execute(sql)
        if vcur == None:
            continue
        else:
            chromid = vcur.fetchone()[0]
            
        maxchromsite = record.seq.__len__()
            
        """ Get a hold of peaks """
        sql = "select id, site from Summits where chrom=" + record.name.__str__() + " order by site ASC"
        vcur.execute(sql)
        for ii in vcur.fetchall():
            summitid = ii[0]
            summitsite = ii[1]
            lowersummitsite = summitsite - 50
            uppersummitsite = summitsite + 50
            if lowersummitsite < 50:
                lowersummitsite = 1
            if uppersummitsite > maxchromsite:
                uppersummitsite = maxchromsite
            summitseq = record.seq[lowersummitsite-1:uppersummitsite-1]
            
            for motifname in gene_motif:            
                score = score_motif_sequence(gene_motif[motifname], summitseq)
                print "161:", speciesid, chromid, summitid, motifname, score


def score_motif_sequence(motif, seq, startsite):
    """Motif is list of lists, motif[site-1] = (wa,wc,wg,wt)
        seq is a string of nts
        startsite is the site number of the first character in seq."""
        
    return 10.0

    


    