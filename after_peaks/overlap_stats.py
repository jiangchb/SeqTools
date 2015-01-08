from chipseqdb import *


def test_species_overlap(con):
    cur = con.cursor()
    sql = "select id from Species"
    cur.execute(sql)
    x = cur.fetchall()
    speciesids = []
    for ii in x:
        speciesids.append( ii[0] )
      
    species_geneids = {}
    for spid in speciesids:
        genes = get_genes_for_species(con, spid)
        species_geneids[spid] = []
        for ii in genes:
            geneid = get_geneid_from_aliasid(ii[0], con)
            species_geneids[spid].append( geneid )
      
    species_species_overlap = {}
    for species in species_geneids:
        print species, get_species_name(species, con), species_geneids[species].__len__() 
         
    for ii in range(0, speciesids.__len__()):
        species_species_overlap[speciesids[ii]] = {}
        for jj in range(ii+1, speciesids.__len__()):
            species_species_overlap[speciesids[ii]][speciesids[jj]] = 0
             
            count = 0
            for gii in species_geneids[ speciesids[ii] ]:
                 if gii in species_geneids[ speciesids[jj] ]:
                     count += 1
                                  
            print speciesids[ii], speciesids[jj], "share", count, "genes."
         
            species_species_overlap[speciesids[ii]][speciesids[jj]] = count

    #print species_species_overlap
    
    sql = "select unionid from Speciesunions"
    cur.execute(sql)
    x = cur.fetchall()
    spunionids = []
    for ii in x:
        spunionids.append( ii[0] )
     
    unionid_speciesid = {}
    
     
    for spunionid in spunionids:
        sql = "select memunionid from SpeciesunionUnions where spunionid = " + spunionid.__str__()
        cur.execute(sql)
        x = cur.fetchall()
        unionids = []
        for ii in x:
            unionids.append( ii[0] )
        
        unionid_countpeaks = {}
        for uu in unionids:
            sql = "select count(*) from UnionGenes where unionid=" + uu.__str__()
            cur.execute(sql)
            unionid_countpeaks[uu] = cur.fetchone()[0]

            sql = "select species from Replicates where id in (select replicate from GroupReplicate where rgroup in (select id from ReplicateGroups where id in (select repgroupid from UnionRepgroups where unionid=" + uu.__str__() + ") ) )"
            cur.execute(sql)
            this_speciesid = cur.fetchone()[0]
            #print "86:", this_speciesid, get_species_name(this_speciesid, con)
            unionid_speciesid[uu] = this_speciesid
            
        sql = "select count(*) from SpeciesUnionGenes where unionid=" + spunionid.__str__()
        cur.execute(sql)
        intersection = cur.fetchone()[0]
        print intersection, unionid_countpeaks
        
        sql = "select count(*) from SpeciesunionGenes where unionid=" + spunionid.__str__()
        cur.execute(sql)
        count_intersect = cur.fetchone()[0] + 1
        
        print "80:", count_intersect
        
        prob_intersect = 1.0
        for uu in unionids:
            #p_x = unionid_countpeaks[uu] # the number of peaks in union uu
            print "85:", count_intersect, "/", species_geneids[unionid_speciesid[uu]].__len__()
            p_x = ( float(count_intersect) / species_geneids[unionid_speciesid[uu]].__len__() )
            prob_intersect *= p_x
            print uu, prob_intersect, p_x
        
        print "\n. P of this species intersection is about", (prob_intersect*count_intersect)

        
        
            
            
    