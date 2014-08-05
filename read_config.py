#
#
#
import os, sys, re

def read_config(path):
    """Returns a hashtable named params."""
    print "\n. Reading the configuration file", path
    
    params = {}
    if False == os.path.exists( path ):
        print "\n. Error, I can't find you configuration file at", path
        exit()
    fin = open(path, "r")
    lines = fin.readlines()
    fin.close()
    curr_species = None
    curr_rgroup = None
    curr_rep = None
    for l in lines:
        l = l.strip()
        
        if l.startswith("SPECIES") and False == l.startswith("SPECIES_UNION"):
            if "species" not in params:
                params["species"] = {}
            tokens = l.split()
            species = " ".join(tokens[1:])
            if species in params["species"]:
                print "\n. Error: your configuration file includes two species with the same name", species
                exit()
            params["species"][species] = {}
            curr_species = species
        
        if l.startswith("GFF"):
            gffpath = l.split("=")[1]
            gffpath = re.sub(" ", "", gffpath)
            #if False == os.path.exists(gffpath):
            #    print "\n. Error reading your configuration file", path
            #    print "--> The following GFF file doesn't exist:", gffpath
            #    exit()
            params["species"][species]["gff"] = gffpath
        
        if l.startswith("NAME"):
            name = l.split("=")[1]
            name = re.sub(" ", "", name)
            params["species"][species]["name"] = name
        
        if l.startswith("REPGROUP"):
            rgroup = l.split()[1]
            if "rgroups" not in params["species"][species]:
                params["species"][species]["rgroups"] = {}
            if rgroup not in params["species"][species]:   
                params["species"][species]["rgroups"][rgroup] = {}
            curr_rgroup = rgroup
        
        if l.startswith("REPLICATE"):
            repid = l.split("=")[0].split()[1]
            repid = re.sub(" ", "", repid)
            curr_rep = repid
            if "reps" not in params["species"][curr_species]["rgroups"][curr_rgroup]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"] = {}
        
        if l.startswith("SUMMITS"):
            summitpath = l.split("=")[1]
            summitpath = re.sub(" ", "", summitpath)
            if curr_rep not in params["species"][curr_species]["rgroups"][curr_rgroup]["reps"]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep] = {}
            params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep]["summitpath"] = summitpath
        
        if l.startswith("ENRICHMENTS"):
            bdgpath = l.split("=")[1]
            bdgpath = re.sub(" ", "", bdgpath)
            if curr_rep not in params["species"][curr_species]["rgroups"][curr_rgroup]["reps"]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep] = {}
            params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep]["bdgpath"] = bdgpath
        
        if l.startswith("UNION"):
            tokens = l.split("=")
            unionname = tokens[0].split()[1]
            unionname = re.sub(" ", "", unionname)                
            if "unions" not in params["species"][curr_species]:
                params["species"][curr_species]["unions"] = {}
            params["species"][curr_species]["unions"][unionname] = []
            x = tokens[1].split()
            for repgroup in x:
                params["species"][curr_species]["unions"][unionname].append( re.sub(" ", "", repgroup) )
            print "\n. Adding union", unionname
    
        if l.startswith("SPECIES_UNION"):
            tokens = l.split("=")
            spunionname = tokens[0].split()[1]
            spunionname = re.sub(" ", "", spunionname)                
            if "unionssp" not in params:
                params["unionssp"] = {}
            params["unionssp"][spunionname] = []
            x = tokens[1].split()
            for unionname in x:
                params["unionssp"][spunionname].append( re.sub(" ", "", unionname) )
            #print "read_config.py 96:", params["unionssp"]
    
    return params

