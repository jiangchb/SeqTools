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
    params["repgroup_memberorder"] = {}
    params["union_memberorder"] = {}
    params["spunion_memberorder"] = {}
    for l in lines:
        l = l.strip()
        
        j = l
        if re.sub(" ", "", j).startswith("#"):
            """Skip commented-out lines"""
            continue
        
        elif l.startswith("SPECIES") and False == l.startswith("SPECIES_UNION"):
            if "species" not in params:
                params["species"] = {}
            tokens = l.split()
            species = " ".join(tokens[1:])
            if species in params["species"]:
                print "\n. Error: your configuration file includes two species with the same name", species
                exit()
            params["species"][species] = {}
            curr_species = species
        
        elif l.startswith("GFF"):
            gffpath = l.split("=")[1]
            gffpath = re.sub(" ", "", gffpath)
            #if False == os.path.exists(gffpath):
            #    print "\n. Error reading your configuration file", path
            #    print "--> The following GFF file doesn't exist:", gffpath
            #    exit()
            params["species"][curr_species]["gff"] = gffpath
        
        elif l.startswith("NAME"):
            name = l.split("=")[1]
            name = re.sub(" ", "", name)
            params["species"][curr_species]["name"] = name
        
        elif l.startswith("REPGROUP"):
            rgroup = l.split()[1]
            if "rgroups" not in params["species"][curr_species]:
                params["species"][curr_species]["rgroups"] = {}
            if rgroup not in params["species"][curr_species]:   
                params["species"][curr_species]["rgroups"][rgroup] = {}
            curr_rgroup = rgroup
        
        elif l.startswith("REPLICATE"):
            repid = l.split("=")[0].split()[1]
            repid = re.sub(" ", "", repid)
            curr_rep = repid
            if "reps" not in params["species"][curr_species]["rgroups"][curr_rgroup]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"] = {}
        
        elif l.startswith("SUMMITS"):
            summitpath = l.split("=")[1]
            summitpath = re.sub(" ", "", summitpath)
            if curr_rep not in params["species"][curr_species]["rgroups"][curr_rgroup]["reps"]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep] = {}
            if False == os.path.exists( summitpath ):
                print "\n. Error: I can't find the SUMMITS file at", summitpath
                exit()
            params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep]["summitpath"] = summitpath
        
        elif l.startswith("ENRICHMENTS"):
            bdgpath = l.split("=")[1]
            bdgpath = re.sub(" ", "", bdgpath)
            if curr_rep not in params["species"][curr_species]["rgroups"][curr_rgroup]["reps"]:
                params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep] = {}
            if False == os.path.exists( bdgpath ):
                print "\n. Error: I can't find the ENRICHMENTS file file at", bdgpath
                exit()
            params["species"][curr_species]["rgroups"][curr_rgroup]["reps"][curr_rep]["bdgpath"] = bdgpath
        
        elif l.startswith("UNION"):
            tokens = l.split("=")
            if tokens.__len__() < 2:
                continue
            unionname = tokens[0].split()[1]
            unionname = re.sub(" ", "", unionname)                
            if "unions" not in params["species"][curr_species]:
                params["species"][curr_species]["unions"] = {}
            params["species"][curr_species]["unions"][unionname] = []
            x = tokens[1].split()
            for repgroup in x:
                params["species"][curr_species]["unions"][unionname].append( re.sub(" ", "", repgroup) )
            #print "\n. Adding union", unionname
    
        elif l.startswith("SPECIES_UNION"):
            tokens = l.split("=")
            spunionname = tokens[0].split()[1]
            spunionname = re.sub(" ", "", spunionname)                
            if "unionssp" not in params:
                params["unionssp"] = {}
            params["unionssp"][spunionname] = []
            x = tokens[1].split()
            for unionname in x:
                params["unionssp"][spunionname].append( re.sub(" ", "", unionname) )
    
    return params

def validate_config(params):
    """params is the hash table created by the function read_config."""
    """This method validates that the definitions within the configuration file
    are internally consistent. e.g., unions don't reference undefined replicate 
    groups, etc."""
    
#     for species in params["species"]:
#         for repgroup in params["species"][species]["rgroups"]:
#             for rep in params["species"][species]["rgroups"][repgroup]["reps"]:
#                 print "119:", species, repgroup, rep
#     exit()
    
    for species in params["species"]:
        if "rgroups" not in params["species"][species]:
            print "\n. Error: the species", species, "contains no replicate groups."
            exit()
        for repgroup in params["species"][species]["rgroups"]:
            """"Do all the replicates in this group exist?"""
            if "reps" not in params["species"][species]["rgroups"][repgroup]:
                print "\n. Error: the replicate group", repgroup, "contains no replicates."
                exit()
#             for rep in params["species"][species]["rgroups"][repgroup]["reps"]:
#                 if "summitpath" not in params["species"][species]["rgroups"][repgroup]["reps"][rep]:
#                     print "\n. Error: you didn't define a SUMMITS line for species", species, "repgroup", repgroup, "replicate", rep
#                     exit()
#                 if "bdgpath" not in params["species"][species]["rgroups"][repgroup]["reps"][rep]:
#                     print "\n. Error: you didn't define an ENRICHMENTS line for replicate", repgroup, "repgroup", repgroup, "replicate", rep
#                     exit()
        if "unions" in params["species"][species]:
            for unionname in params["species"][species]["unions"]:
                for repgroup in params["species"][species]["unions"][unionname]:
                    if repgroup not in params["species"][species]["rgroups"]:
                        print "\n. Error, repgroup", repgroup, "is references by union", unionname, "but it's undefined!"
                        exit()
    
    return True
