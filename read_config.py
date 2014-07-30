#
#
#
import os, sys, re

def read_config(path):
    """Returns a hashtable named params."""
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
        if l.startswith("SPECIES"):
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
            if False == os.path.exists(gffpath):
                print "\n. Error reading your configuration file", path
                print "--> The following GFF file doesn't exist:", gffpath
                exit()
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
            ids = l.split()[1:]
            params["species"][curr_species]["rgroups"][curr_rgroup]["union"] = ids
    return params

