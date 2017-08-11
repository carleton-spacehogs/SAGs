#Author: Michael Hoffert
#Date: 6/15/2017
#Credit to Rika Anderson
#Python version 2.7.12

"""Constructs lists of contigs and SNVs from input files"""
from Contig import *
from SNV import *
from SAAV import *
import os

#example SAG fasta file name:
#AD-155-C09-many_assemblers_simple.CISA.ctg.fa
#example contig name .tsv file:
#C09_contig_names.tsv
#example names_map file:
#C09.names_map
#
#Note: The file names MUST be formatted in the above format

def make_contigs(SAG_fasta_file, contig_name_tsv, names_map):
    
    # uses SAG fasta file to load contig sequences
    # and mapped name/tsv files to update names
    # for contigs that were simplified using anvi-script
    
    SAG_fa_file = open(SAG_fasta_file)
    contig_names = open(contig_name_tsv)
    names_map_file = open(names_map)

    if SAG_fasta_file[7:10] == contig_name_tsv[:3]: #checkfilenames
        print 'filenames checked out'
        
        contig_list = [] #output list for contigs
        
        SAG_fa_lines = SAG_fa_file.readlines()
        contig_name_lines = contig_names.readlines()

        for i in range(0, len(SAG_fa_lines), 2): #make contigs from fasta file
            
            ID = SAG_fa_lines[i][1 : (len(SAG_fa_lines[i]) - 1)]
            sequence = SAG_fa_lines[i + 1]
            SAG = contig_name_tsv[:3]
            contig = Contig(ID, sequence, SAG)
            contig_list.append(contig)
            
        #update contig names from simplified versions
        #anvi-script -> SAG.fasta names
        for line in contig_name_lines:
            line = line.split()
            old_ID = line[0]
            new_ID = line[1]
            for i in range(len(contig_list)):
                if contig_list[i].update_ID(old_ID, new_ID):
                    break
        
        #SAG.fasta names -> JGI contig names
        names_map_lines = names_map_file.readlines()
        for line in names_map_lines:
            line = line.split()
            old_ID = line[0]
            new_ID = line[1]
            for i in range(len(contig_list)):
                if contig_list[i].update_ID(old_ID, new_ID):
                    break
    else:
        print "The SAG FASTA file and  contig name .tsv file do not match!"
    
    SAG_fa_file.close()
    contig_names.close()
    names_map_file.close()
    
    #print contig_list
    return contig_list


#makes list of SNVs from variability profile
def make_SNVs(variability_file):
    
    SNV_list = []
    
    variability_file = open(variability_file)
    variability_lines = variability_file.readlines()
    num_SNVs = len(variability_lines),'SNVs in file'
    print num_SNVs,'SNVs in file'
    
    count = 0
    for line in variability_lines:
        
        line = line.split()
        if not (line[0] == 'entry_id'):
            
            #print line
            snv = SNV(line)
            SNV_list.append(snv)
            
        count += 1
            
    #print SNV_list
    return SNV_list

def make_SAAVs(variability_file):
    
    SAAV_list = []
    
    variability_file = open(variability_file)
    variability_lines = variability_file.readlines()
    num_SAAVs = len(variability_lines),'SAAVs in file'
    print num_SAAVs,'SAAVs in file'
    
    count = 0
    for line in variability_lines:
        
        line = line.split()
        if not (line[0] == 'entry_id'):
            #print 'THIS IS THE LINE!',line
            saav = SAAV(line)
            SAAV_list.append(saav)
            
        count += 1
            
    return SAAV_list
    

def make_SNV_SAAV_contig_names_consistent(SNVs, Contigs, SAAVs):
    count = 0
    num_contigs = len(Contigs)
    
    for contig in Contigs:
        prev = count
        #os.system('cls')
        print float(count) / num_contigs * 100, 'percent complete!'
        for snv in SNVs:
            
            if snv.update_contig_name(contig.first_ID, contig.ID):
                pass         
        
        for saav in SAAVs:
            
            if saav.update_contig_name(contig.first_ID, contig.ID):
                pass
        
        count += 1
        
                    
    return [SNVs,Contigs,SAAVs]
                
    
        
if __name__ == '__main__':
    
    c = make_contigs('AD-155-C09-many_assemblers_simple.CISA.ctg.fa', 'C09_contig_names.tsv','C09.names_map')
    s = make_SNVs('C09_variability_new.txt')
    make_SNV_contig_names_consistent(s, c)
    
    
    
    
    
    
    
    
    
    
    